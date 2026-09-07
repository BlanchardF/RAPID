#!/usr/bin/env python3
"""
RAPID — scripts/primer_report.py
Build an enriched, sortable report of the primer pairs kept by a `track` run,
to make picking the best candidates easy.

For every primer pair in best_primers.txt it collects:
  - amplicon (product) size, penalty, primer Tm and GC%
  - the gene's expression across ALL stages as CPM (depth-normalised)
  - the enrichment fold-change of the reference over the most-expressed other
    stage, and which stage is the limiting one

Outputs:
  <out_tsv>  : one row per primer pair with all metrics
  <out_html> : self-contained interactive table — sort/filter, hide (x) primers
               you don't want as you review, review/restore hidden ones, expand a
               row for a per-stage CPM bar chart, and download your selection.

Usage:
  python3 primer_report.py <results_dir> [--out-tsv FILE] [--out-html FILE]

<results_dir> is a `rapid track` output directory. The script reads:
  <results_dir>/rapid_track_config.yaml           (sample names / stage labels)
  <results_dir>/featurecounts/counts_<name>.txt   (per-sample counts)
  <results_dir>/primer3/results/best_primers.txt  (kept primer pairs)
"""

import argparse
import json
import math
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def load_samples(results_dir):
    """
    Return list of (name, role, label, counts_path). Uses rapid_track_config.yaml
    for nice stage labels (folder name of the reads); falls back to filenames.
    """
    cfg = Path(results_dir) / "rapid_track_config.yaml"
    fc_dir = Path(results_dir) / "featurecounts"
    samples = []
    if cfg.is_file():
        try:
            import yaml
            data = yaml.safe_load(cfg.read_text())
            for s in data.get("samples", []):
                name = s["name"]
                role = s.get("role", "other")
                r1 = s.get("r1", "")
                label = Path(r1).parent.name if r1 else name
                samples.append((name, role, label or name))
        except Exception as e:
            print(f"[primer_report] config unreadable ({e}); using filenames.",
                  file=sys.stderr)
            samples = []
    if not samples:
        if (fc_dir / "counts_ref.txt").is_file():
            samples.append(("ref", "reference", "ref"))
        for f in sorted(fc_dir.glob("counts_other_*.txt")):
            name = f.stem.replace("counts_", "")
            samples.append((name, "other", name))
    out = []
    for name, role, label in samples:
        p = fc_dir / f"counts_{name}.txt"
        if p.is_file():
            out.append((name, role, label, p))
        else:
            print(f"[primer_report] missing {p}, skipping sample {name}.", file=sys.stderr)
    return out


def read_counts(path):
    counts = {}
    header_seen = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True
                continue
            row = line.rstrip("\n")
            if not row:
                continue
            cols = row.split("\t")
            try:
                counts[cols[0]] = int(round(float(cols[-1])))
            except ValueError:
                counts[cols[0]] = 0
    return counts


def to_cpm(counts):
    lib = sum(counts.values())
    if lib == 0:
        return {g: 0.0 for g in counts}
    f = 1e6 / lib
    return {g: v * f for g, v in counts.items()}


def parse_best_primers(path):
    primers, block = [], {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "=":
                if block.get("SEQUENCE_ID"):
                    primers.append(block)
                block = {}
                continue
            if "=" in line:
                k, v = line.split("=", 1)
                block[k] = v
        if block.get("SEQUENCE_ID"):
            primers.append(block)
    return primers


def gene_of(seq_id):
    return seq_id.rsplit("_j", 1)[0]


def junction_of(seq_id):
    parts = seq_id.rsplit("_j", 1)
    return parts[1] if len(parts) == 2 else ""


def fnum(block, key, default=""):
    try:
        return float(block.get(key, default))
    except (ValueError, TypeError):
        return default


# ---------------------------------------------------------------------------
# Build records
# ---------------------------------------------------------------------------

def build_records(results_dir):
    samples = load_samples(results_dir)
    if not samples:
        print("[primer_report] ERROR: no featureCounts tables found.", file=sys.stderr)
        sys.exit(1)

    ref = next((s for s in samples if s[1] == "reference"), None)
    others = [s for s in samples if s[1] == "other"]
    if ref is None:
        print("[primer_report] ERROR: no reference sample found.", file=sys.stderr)
        sys.exit(1)

    ref_cpm = to_cpm(read_counts(ref[3]))
    others_cpm = [(lbl, to_cpm(read_counts(p))) for (_, _, lbl, p) in others]
    ref_label = ref[2]

    best = Path(results_dir) / "primer3" / "results" / "best_primers.txt"
    if not best.is_file():
        print(f"[primer_report] ERROR: {best} not found.", file=sys.stderr)
        sys.exit(1)

    records = []
    for b in parse_best_primers(best):
        sid = b["SEQUENCE_ID"]
        gene = gene_of(sid)
        rcpm = ref_cpm.get(gene, 0.0)
        stage_cpm = [(lbl, cpm.get(gene, 0.0)) for lbl, cpm in others_cpm]
        max_other = max((v for _, v in stage_cpm), default=0.0)
        limiting = max(stage_cpm, key=lambda x: x[1])[0] if stage_cpm else ""
        enrichment = float("inf") if max_other <= 0 else rcpm / max_other
        records.append({
            "gene": gene, "junction": junction_of(sid), "seq_id": sid,
            "amplicon_bp": int(fnum(b, "PRIMER_PAIR_0_PRODUCT_SIZE", 0) or 0),
            "penalty": round(fnum(b, "PRIMER_PAIR_0_PENALTY", 0) or 0, 3),
            "fwd": b.get("PRIMER_LEFT_0_SEQUENCE", ""),
            "rev": b.get("PRIMER_RIGHT_0_SEQUENCE", ""),
            "fwd_tm": round(fnum(b, "PRIMER_LEFT_0_TM", 0) or 0, 1),
            "rev_tm": round(fnum(b, "PRIMER_RIGHT_0_TM", 0) or 0, 1),
            "fwd_gc": round(fnum(b, "PRIMER_LEFT_0_GC_PERCENT", 0) or 0, 1),
            "rev_gc": round(fnum(b, "PRIMER_RIGHT_0_GC_PERCENT", 0) or 0, 1),
            "cpm_ref": round(rcpm, 1),
            "stage_cpm": [{"stage": lbl, "cpm": round(v, 1)} for lbl, v in stage_cpm],
            "enrichment": (None if math.isinf(enrichment) else round(enrichment, 1)),
            "limiting_stage": limiting,
        })
    return records, ref_label, [lbl for _, _, lbl, _ in others]


# ---------------------------------------------------------------------------
# TSV
# ---------------------------------------------------------------------------

def write_tsv(records, ref_label, other_labels, out_tsv):
    cols = (["gene", "junction", "amplicon_bp", "enrichment_fold", "limiting_stage",
             "penalty", "fwd", "rev", "fwd_tm", "rev_tm", "fwd_gc", "rev_gc",
             f"cpm_{ref_label}"] + [f"cpm_{l}" for l in other_labels])
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in records:
            stage_map = {d["stage"]: d["cpm"] for d in r["stage_cpm"]}
            enr = "inf" if r["enrichment"] is None else r["enrichment"]
            row = [r["gene"], r["junction"], r["amplicon_bp"], enr, r["limiting_stage"],
                   r["penalty"], r["fwd"], r["rev"], r["fwd_tm"], r["rev_tm"],
                   r["fwd_gc"], r["rev_gc"], r["cpm_ref"]]
            row += [stage_map.get(l, 0.0) for l in other_labels]
            fh.write("\t".join(str(x) for x in row) + "\n")


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RAPID — primer selection</title>
<style>
  :root{--bg:#0f1115;--panel:#171a21;--line:#262b36;--fg:#e6e9ef;--muted:#98a2b3;
        --accent:#5b9dff;--good:#2ea043;--warn:#d29922;--bad:#f85149;--bar:#33415a;--refbar:#2ea043;}
  *{box-sizing:border-box}
  body{margin:0;font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--fg)}
  header{padding:16px 22px;border-bottom:1px solid var(--line)}
  h1{margin:0 0 4px;font-size:18px}
  .sub{color:var(--muted);font-size:13px}
  .wrap{padding:16px 22px}
  .toolbar{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin-bottom:14px}
  label.ctl{font-size:13px;color:var(--muted);display:flex;gap:6px;align-items:center}
  select,input{background:#0d1017;color:var(--fg);border:1px solid var(--line);border-radius:8px;padding:6px 9px;font-size:13px}
  button{background:#1f2530;color:var(--fg);border:1px solid var(--line);border-radius:8px;padding:6px 11px;cursor:pointer;font-size:13px}
  button:hover{border-color:var(--accent)}
  table{border-collapse:collapse;width:100%}
  th,td{padding:8px 10px;text-align:left;border-bottom:1px solid var(--line);font-size:13px;white-space:nowrap}
  th{color:var(--muted);font-weight:600;position:sticky;top:0;background:var(--bg);cursor:pointer;user-select:none}
  th:hover{color:var(--fg)}
  tr.main:hover{background:#141821}
  .mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px}
  .enr{font-weight:700}.enr.hi{color:var(--good)}.enr.mid{color:var(--warn)}
  .pill{display:inline-block;padding:1px 8px;border-radius:999px;background:#1f2530;border:1px solid var(--line);font-size:12px}
  .expand{cursor:pointer;color:var(--accent)}
  .hidebtn{cursor:pointer;color:var(--muted);border:1px solid var(--line);border-radius:6px;padding:1px 7px;background:#1a1f28}
  .hidebtn:hover{color:var(--bad);border-color:var(--bad)}
  .hidebtn.restore:hover{color:var(--good);border-color:var(--good)}
  .detail td{background:#0d1017}
  .bars{display:flex;flex-direction:column;gap:4px;max-width:660px}
  .barrow{display:grid;grid-template-columns:150px 1fr 96px;align-items:center;gap:8px;font-size:12px}
  .barrow .name{color:var(--muted);text-align:right;font-style:italic;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .track{background:#0f1319;border:1px solid var(--line);border-radius:5px;height:14px;overflow:hidden}
  .fill{height:100%;background:var(--bar)}.fill.ref{background:var(--refbar)}
  .val{color:var(--muted);text-align:right}
  .seq{font-family:ui-monospace,monospace;font-size:12px}
</style>
</head>
<body>
<header>
  <h1>RAPID — primer selection</h1>
  <div class="sub">Sort and filter primer pairs. Click <b>x</b> to discard the ones you don't want; your selection is kept (even after reload).</div>
</header>
<div class="wrap">
  <div class="toolbar">
    <label class="ctl">Sort by:
      <select id="sortBy" onchange="render()">
        <option value="enr">enrichment &#8595;</option>
        <option value="amp">amplicon size &#8593;</option>
        <option value="pen">penalty &#8593;</option>
        <option value="gene">gene</option>
      </select>
    </label>
    <label class="ctl">Min enrichment: <input type="number" id="minEnr" value="0" style="width:80px" oninput="render()"></label>
    <label class="ctl">Amplicon: <input type="number" id="ampMin" value="0" style="width:70px" oninput="render()"> &#8211; <input type="number" id="ampMax" value="100000" style="width:70px" oninput="render()"> bp</label>
    <button id="toggleHidden" onclick="toggleHidden()"></button>
    <button onclick="restoreAll()">Restore all</button>
    <button onclick="downloadSelection()">&#8681; Download selection (TSV)</button>
    <span class="sub" id="count"></span>
  </div>
  <table>
    <thead><tr>
      <th></th>
      <th onclick="setSort('gene')">gene</th>
      <th onclick="setSort('amp')">amplicon</th>
      <th onclick="setSort('enr')">enrichment</th>
      <th>limiting stage</th>
      <th onclick="setSort('pen')">penalty</th>
      <th>Tm (F/R)</th>
      <th>GC% (F/R)</th>
      <th>primers</th>
      <th></th>
    </tr></thead>
    <tbody id="tbody"></tbody>
  </table>
</div>
<script>
const DATA = __DATA__;                 // {ref_label, other_labels, records:[...]}
const STORE_KEY = "rapid_primer_hidden_" + (DATA.ref_label||"ref");
let sortKey = "enr";
let showHidden = false;

function loadHidden(){
  try{ return new Set(JSON.parse(localStorage.getItem(STORE_KEY)||"[]")); }
  catch(e){ return new Set(); }
}
function saveHidden(s){
  try{ localStorage.setItem(STORE_KEY, JSON.stringify([...s])); }catch(e){}
}
let hidden = loadHidden();

function setSort(k){ document.getElementById('sortBy').value=k; sortKey=k; render(); }
function toggleHidden(){ showHidden=!showHidden; render(); }
function restoreAll(){ hidden.clear(); saveHidden(hidden); render(); }
function hide(id){ hidden.add(id); saveHidden(hidden); render(); }
function unhide(id){ hidden.delete(id); saveHidden(hidden); render(); }

function enrVal(r){ return r.enrichment===null ? Infinity : r.enrichment; }
function enrText(r){ return r.enrichment===null ? "\u221E" : ("x"+r.enrichment); }
function enrClass(r){ const v=enrVal(r); return v>=10?"enr hi":(v>=3?"enr mid":"enr"); }
function maxCpmFor(r){ let m=r.cpm_ref; for(const s of r.stage_cpm) m=Math.max(m,s.cpm); return m; }
function barWidth(v,max){ const lv=Math.log10(v+1), lm=Math.log10(max+1); return lm>0?Math.max(1,(lv/lm)*100):1; }

function visibleRows(){
  const minEnr=parseFloat(document.getElementById('minEnr').value)||0;
  const ampMin=parseFloat(document.getElementById('ampMin').value)||0;
  const ampMax=parseFloat(document.getElementById('ampMax').value)||1e9;
  let rows=DATA.records.filter(r=> enrVal(r)>=minEnr && r.amplicon_bp>=ampMin && r.amplicon_bp<=ampMax);
  rows=rows.filter(r=> showHidden ? hidden.has(r.seq_id) : !hidden.has(r.seq_id));
  rows.sort((a,b)=>{
    if(sortKey==="enr") return enrVal(b)-enrVal(a);
    if(sortKey==="amp") return a.amplicon_bp-b.amplicon_bp;
    if(sortKey==="pen") return a.penalty-b.penalty;
    if(sortKey==="gene") return a.gene.localeCompare(b.gene);
    return 0;
  });
  return rows;
}

function render(){
  sortKey = document.getElementById('sortBy').value;
  const rows=visibleRows();
  const kept = DATA.records.length - hidden.size;
  document.getElementById('toggleHidden').textContent =
    showHidden ? `\u21A9 Back to selection` : `View hidden (${hidden.size})`;
  document.getElementById('count').textContent =
    showHidden ? `${rows.length} hidden shown`
               : `${rows.length} shown — ${kept} kept / ${DATA.records.length}`;

  const tb=document.getElementById('tbody'); tb.innerHTML="";
  rows.forEach((r,idx)=>{
    const tr=document.createElement('tr'); tr.className="main";
    const btn = showHidden
      ? `<span class="hidebtn restore" title="Restore" data-un="${r.seq_id}">&#8635;</span>`
      : `<span class="hidebtn" title="Discard" data-hide="${r.seq_id}">&#10005;</span>`;
    tr.innerHTML=`
      <td>${btn}</td>
      <td class="mono">${r.gene}<span class="sub"> j${r.junction}</span></td>
      <td><b>${r.amplicon_bp}</b> bp</td>
      <td class="${enrClass(r)}">${enrText(r)}</td>
      <td><span class="pill">${r.limiting_stage||"\u2014"}</span></td>
      <td>${r.penalty}</td>
      <td>${r.fwd_tm} / ${r.rev_tm}</td>
      <td>${r.fwd_gc} / ${r.rev_gc}</td>
      <td class="seq">${r.fwd}<br>${r.rev}</td>
      <td class="expand" data-i="${idx}">details &#9662;</td>`;
    tb.appendChild(tr);

    const det=document.createElement('tr'); det.className="detail"; det.style.display="none";
    const max=maxCpmFor(r);
    let bars=`<div class="barrow"><div class="name">${DATA.ref_label} (ref)</div>`
           +`<div class="track"><div class="fill ref" style="width:${barWidth(r.cpm_ref,max)}%"></div></div>`
           +`<div class="val">${r.cpm_ref} CPM</div></div>`;
    for(const s of r.stage_cpm){
      bars+=`<div class="barrow"><div class="name">${s.stage}</div>`
          +`<div class="track"><div class="fill" style="width:${barWidth(s.cpm,max)}%"></div></div>`
          +`<div class="val">${s.cpm} CPM</div></div>`;
    }
    det.innerHTML=`<td colspan="10"><div class="bars">${bars}</div>
      <div class="sub" style="margin-top:8px">amplicon ${r.amplicon_bp} bp — enrichment ${enrText(r)} vs the most expressed stage (${r.limiting_stage||"none"})</div></td>`;
    tb.appendChild(det);

    tr.querySelector('.expand').addEventListener('click',()=>{ det.style.display = det.style.display==="none" ? "" : "none"; });
    const h=tr.querySelector('[data-hide]'); if(h) h.addEventListener('click',()=>hide(h.dataset.hide));
    const u=tr.querySelector('[data-un]');   if(u) u.addEventListener('click',()=>unhide(u.dataset.un));
  });
}

function downloadSelection(){
  // export currently kept (non-hidden) primers, ignoring the show-hidden view
  const wasShown=showHidden; showHidden=false;
  const rows=visibleRows(); showHidden=wasShown;
  const cols=["gene","junction","amplicon_bp","enrichment_fold","limiting_stage","penalty",
              "fwd","rev","fwd_tm","rev_tm","fwd_gc","rev_gc","cpm_"+DATA.ref_label]
              .concat(DATA.other_labels.map(l=>"cpm_"+l));
  let out=cols.join("\t")+"\n";
  for(const r of rows){
    const sm={}; r.stage_cpm.forEach(d=>sm[d.stage]=d.cpm);
    const enr = r.enrichment===null ? "inf" : r.enrichment;
    let row=[r.gene,r.junction,r.amplicon_bp,enr,r.limiting_stage,r.penalty,
             r.fwd,r.rev,r.fwd_tm,r.rev_tm,r.fwd_gc,r.rev_gc,r.cpm_ref]
             .concat(DATA.other_labels.map(l=>sm[l]!==undefined?sm[l]:0));
    out+=row.join("\t")+"\n";
  }
  const blob=new Blob([out],{type:"text/tab-separated-values"});
  const a=document.createElement("a");
  a.href=URL.createObjectURL(blob); a.download="primer_selection.tsv";
  document.body.appendChild(a); a.click(); a.remove();
}

render();
</script>
</body>
</html>
"""


def write_html(records, ref_label, other_labels, out_html):
    data = {"ref_label": ref_label, "other_labels": other_labels, "records": records}
    blob = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    Path(out_html).parent.mkdir(parents=True, exist_ok=True)
    Path(out_html).write_text(HTML_TEMPLATE.replace("__DATA__", blob), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Enriched primer report for RAPID track.")
    ap.add_argument("results_dir", help="rapid track output directory.")
    ap.add_argument("--out-tsv", default=None)
    ap.add_argument("--out-html", default=None)
    a = ap.parse_args()

    rd = Path(a.results_dir)
    out_tsv = a.out_tsv or str(rd / "primer3" / "results" / "primer_report.tsv")
    out_html = a.out_html or str(rd / "primer3" / "results" / "primer_report.html")

    records, ref_label, other_labels = build_records(rd)
    write_tsv(records, ref_label, other_labels, out_tsv)
    write_html(records, ref_label, other_labels, out_html)

    print(f"[primer_report] primers        : {len(records)}", file=sys.stderr)
    print(f"[primer_report] reference stage : {ref_label}", file=sys.stderr)
    print(f"[primer_report] other stages    : {', '.join(other_labels)}", file=sys.stderr)
    print(f"[primer_report] TSV  -> {out_tsv}", file=sys.stderr)
    print(f"[primer_report] HTML -> {out_html}", file=sys.stderr)


if __name__ == "__main__":
    main()
