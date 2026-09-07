#!/usr/bin/env python3
"""
RAPID — scripts/run_primer3_parallel.py
Launch one primer3_core process per gene config file using GNU Parallel,
dispatched dynamically across --threads cores, with a live progress bar/ETA
and a per-gene wall-clock timeout.

Usage:
    python3 run_primer3_parallel.py <gene_configs_dir> <output_raw> <n_jobs>
                                     [timeout_seconds]

Input
-----
<gene_configs_dir> is the directory produced by fill_primer3_config.py: one
config file per gene (e.g. g10004.txt), each containing all of that gene's
junction blocks with their SEQUENCE_TEMPLATE already filled in.

timeout_seconds (optional, 4th arg): maximum wall-clock time allowed for a
single gene. A gene still running after this limit is killed and marked as
failed; the run continues with the other genes. Default: 7200 (2 hours).
Set to 0 to disable the timeout entirely.

How it works
------------
- GNU Parallel runs `primer3_core` on every gene file. It dispatches jobs
  dynamically: the moment a core finishes one gene, it immediately grabs the
  next gene from the queue. This keeps all n_jobs cores busy even though
  genes take very different amounts of time (a gene with many exons / a long
  transcript takes far longer than a short single-junction gene).
- Each primer3_core call is wrapped in `timeout <N>s` so that a single
  pathological gene (e.g. a very long fused transcript with thousands of
  junctions) can never stall the whole multi-day run: it is aborted after
  the limit, recorded as failed, and the pipeline moves on. `timeout`
  returns exit code 124 when it kills the job.
- `--bar` shows a live percentage + ETA, computed by Parallel itself from the
  number of gene files already finished, the number remaining, and the number
  of job slots (-j).
- `--will-cite` silences the one-time GNU Parallel citation prompt, so the run
  never blocks waiting for interactive input.
- `--joblog` records the exit status of every gene job, so failing / timed-out
  genes are reported in <output_raw>.failed.tsv instead of losing the run.
"""

import csv
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


DEFAULT_TIMEOUT = 7200   # seconds (2 hours)


def main():
    if len(sys.argv) not in (4, 5):
        print(f"Usage: {sys.argv[0]} <gene_configs_dir> <output_raw> <n_jobs> "
              f"[timeout_seconds]", file=sys.stderr)
        sys.exit(1)

    config_dir  = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    n_jobs      = max(1, int(sys.argv[3]))
    timeout_s   = int(sys.argv[4]) if len(sys.argv) == 5 else DEFAULT_TIMEOUT
    timeout_s   = max(0, timeout_s)   # negative → treated as disabled

    if shutil.which("parallel") is None:
        print("[primer3_parallel] ERROR: GNU Parallel ('parallel') not found in PATH.\n"
              "  It must be installed in the 'rapid' conda environment "
              "(see rapid.yaml) — run 'rapid install --update-only' to add it.",
              file=sys.stderr)
        sys.exit(1)

    if timeout_s > 0 and shutil.which("timeout") is None:
        print("[primer3_parallel] ERROR: 'timeout' (GNU coreutils) not found in PATH, "
              "but a per-gene timeout was requested.\n"
              "  Install coreutils in the 'rapid' environment, or pass a timeout of 0 "
              "to disable it.", file=sys.stderr)
        sys.exit(1)

    if not config_dir.is_dir():
        print(f"[primer3_parallel] ERROR: gene configs directory not found: {config_dir}",
              file=sys.stderr)
        sys.exit(1)

    gene_files = sorted(config_dir.glob("*.txt"))
    if not gene_files:
        print(f"[primer3_parallel] ERROR: no gene config (*.txt) files in {config_dir}",
              file=sys.stderr)
        sys.exit(1)

    n_genes = len(gene_files)
    print(f"[primer3_parallel] {n_genes} per-gene config file(s) in {config_dir}",
          flush=True)
    if timeout_s > 0:
        print(f"[primer3_parallel] Per-gene timeout: {timeout_s}s "
              f"({timeout_s / 3600:.2f} h) — genes exceeding it are excluded.",
              flush=True)
    else:
        print("[primer3_parallel] Per-gene timeout: DISABLED.", flush=True)
    print(f"[primer3_parallel] Launching GNU Parallel with {n_jobs} core(s) — "
          f"one primer3_core per gene, dynamically scheduled.", flush=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # List of gene files fed to Parallel (one path per line).
    joblist_path = output_path.parent / "primer3_gene_joblist.txt"
    with open(joblist_path, "w") as fh:
        for gf in gene_files:
            fh.write(str(gf) + "\n")

    joblog_path = output_path.parent / "primer3_gene_joblog.tsv"

    # Per-gene command run by Parallel.
    #   {}   = gene config file path       (e.g. .../g10004.txt)
    #   {.}  = same path without extension  (e.g. .../g10004)  → .out / .err
    # Wrap primer3_core in `timeout` when a limit is set. `timeout` kills the
    # job and returns 124 if the limit is reached.
    if timeout_s > 0:
        primer3_cmd = f"timeout {timeout_s}s primer3_core < {{}} > {{.}}.out 2> {{.}}.err"
    else:
        primer3_cmd = f"primer3_core < {{}} > {{.}}.out 2> {{.}}.err"

    # --will-cite : no interactive citation prompt
    # --bar       : live % + ETA on stderr
    # --joblog    : per-gene exit status
    # 2>&1 | tr '\r' '\n' : keep the progress bar readable both in a static
    #   log file and live under `tail -f`.
    cmd = (
        f"parallel --will-cite --bar -j {n_jobs} "
        f"--joblog {shlex.quote(str(joblog_path))} "
        f"{shlex.quote(primer3_cmd)} "
        f":::: {shlex.quote(str(joblist_path))} "
        f"2>&1 | tr '\\r' '\\n'"
    )
    subprocess.run(cmd, shell=True)

    # ── Read the joblog to know which genes succeeded ───────────────────
    if not joblog_path.is_file():
        print("[primer3_parallel] ERROR: GNU Parallel produced no joblog — "
              "it likely failed to start at all.", file=sys.stderr)
        sys.exit(1)

    exit_by_seq = {}
    with open(joblog_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                exit_by_seq[int(row["Seq"])] = int(row["Exitval"])
            except (KeyError, ValueError):
                pass

    # ── Merge successful outputs, in the same order as the joblist ──────
    failed_path = Path(str(output_path) + ".failed.tsv")
    failures      = []
    timed_out     = []

    with open(output_path, "w") as out_fh, open(failed_path, "w") as failed_fh:
        failed_fh.write("Gene_file\tExitval\tReason\tStderr\n")
        for seq, gf in enumerate(gene_files, start=1):   # joblist order == Seq order
            out_file = gf.with_suffix(".out")
            err_file = gf.with_suffix(".err")
            exitval  = exit_by_seq.get(seq)

            if exitval == 0 and out_file.is_file():
                out_fh.write(out_file.read_text())
            else:
                stderr_text = err_file.read_text().strip() if err_file.is_file() else ""
                # timeout returns 124 when it kills the job
                if exitval == 124:
                    reason = f"timeout (>{timeout_s}s)"
                    timed_out.append(gf.name)
                else:
                    reason = "error"
                failures.append(gf.name)
                failed_fh.write(f"{gf.name}\t{exitval}\t{reason}\t{stderr_text}\n")

    n_ok = n_genes - len(failures)
    print(f"[primer3_parallel] Done — {n_ok}/{n_genes} gene(s) succeeded.", flush=True)

    if timed_out:
        print(f"[primer3_parallel] {len(timed_out)} gene(s) excluded by the "
              f"{timeout_s}s timeout: {', '.join(timed_out[:10])}"
              + (" …" if len(timed_out) > 10 else ""), flush=True)

    if failures:
        print(f"[primer3_parallel] WARNING: {len(failures)} gene(s) failed or timed out. "
              f"Details in:\n  {failed_path}", flush=True)
        if n_ok == 0:
            print("[primer3_parallel] ERROR: every gene failed.", file=sys.stderr)
            sys.exit(1)

    sys.exit(0)


if __name__ == "__main__":
    main()
