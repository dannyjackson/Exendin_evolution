#!/usr/bin/env python3
"""
Curate one protein per exendin-family locus from the output of
2a_collect_exendin_projections.py.

Each query locus hit by an exendin or anchor (ANC) projection becomes one
sequence. Loci are named by where they sit in the conserved SLC4A8-ANC-AQP5
interval, measured from the ANC gene:

  <Query>_ANC       the ancestral exendin gene (LOC_00012323 in HelSus)
  <Query>_A038      38 kb from ANC, towards AQP5 (the side with the HelSus exendin-2/-4 cluster)
  <Query>_S118      118 kb from ANC, towards SLC4A8
  <Query>_U1, U2..  not on the same contig as ANC (cannot be placed)
  a letter suffix (a, b, ...) separates two loci that round to the same name;
  curation labels (e.g. _EX2) are appended at the end.

Direction comes from the AQP5 projection on the ANC contig. If AQP5 is not on
that contig, the ANC strand is used instead (ANC on the minus strand in HelSus
with AQP5 at higher coordinates); the table records which was used.

For each locus the projection with the best TOGA2 status (then longest CDS) is
kept. A locus enters the tree if its status is in --status, its CDS is a
multiple of 3, at least --min-cds bp, with no internal stop codon, and the
curation file does not drop it.

Curation file (tab-separated; '#' lines are comments). Loci are identified by
query, contig and locus start exactly as written in 2a's loci.tsv:
  query  q_chrom  q_start  action  value  note
  action = anchor : use this locus as ANC (when TOGA2 found no ANC projection)
           drop   : exclude from the tree (value ignored)
           keep   : override the sequence-quality filters
           replace: use the sequence in <value>.aa.fa / <value>.cds.fa instead of the
                    projection (the output of 2e_rescue_locus.py), for loci whose TOGA2
                    model is wrong because the reference isoform is mis-annotated
           label  : append _<value> to the name (e.g. EX2)

Outputs (--outdir)
  loci_curated.tsv   every locus: name, position, chosen projection, QC, decision
  curated_aa.fa      proteins that enter the tree
  curated_cds.fa     their CDS (no stop codon)

Usage
  python3 2b_curate_exendin_loci.py --collect-dir collect --curation 2_curation.tsv --outdir curated
"""
import argparse, csv, os, re, sys
from collections import defaultdict

RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
_B, _AA = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: _AA[16 * i + 4 * j + k]
        for i, a in enumerate(_B) for j, b in enumerate(_B) for k, c in enumerate(_B)}


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_collect_fasta(path):
    """2a headers: >QUERY|projection|status|q_gene"""
    out, key = {}, None
    if os.path.exists(path):
        for line in open(path):
            if line.startswith(">"):
                f = line[1:].strip().split("|")
                key = (f[0], f[1]); out[key] = []
            elif key:
                out[key].append(line.strip())
    return {k: re.sub(r"[^ACGTN]", "", "".join(v).upper()) for k, v in out.items()}


def read_one(path):
    """Return the single sequence in a FASTA file, or '' if it cannot be read."""
    try:
        return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()
    except OSError:
        return ""


def read_curation(path):
    cur = defaultdict(list)
    if path:
        for line in open(path):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] == "query":
                continue
            cur[(f[0], f[1], int(f[2]))].append((f[3], f[4] if len(f) > 4 else "", f[5] if len(f) > 5 else ""))
    return cur

def isoform(projection):
    """mRNA number, so the first (canonical) isoform wins ties."""
    m = re.search(r"-mRNA-(\d+)", projection)
    return int(m.group(1)) if m else 99
    
def mid(x):
    return (int(x["q_start"]) + int(x["q_end"])) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collect-dir", required=True, help="the 2a --outdir")
    ap.add_argument("--exendins", default="LOC_00012324,LOC_00012325",
                    help="reference exendin genes (HelSus exendin-2, exendin-4)")
    ap.add_argument("--anchor", default="LOC_00012323", help="reference ANC gene")
    ap.add_argument("--flank", default="LOC_00012326", help="reference AQP5 gene (orientation only)")
    ap.add_argument("--curation", help="curation TSV (see header)")
    ap.add_argument("--status", default="FI,I", help="TOGA2 statuses allowed into the tree")
    ap.add_argument("--min-cds", type=int, default=150)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    exn, allowed = set(a.exendins.split(",")), set(a.status.split(","))
    projs = read_tsv(os.path.join(a.collect_dir, "projections.tsv"))
    loci = read_tsv(os.path.join(a.collect_dir, "loci.tsv"))
    seqs = {}
    for g in exn | {a.anchor}:
        seqs.update(read_collect_fasta(os.path.join(a.collect_dir, f"{g}.nucleotide.fa")))
    if not seqs:
        sys.exit("no sequences: run 2a with --outdir and with the exendin and anchor genes in --genes")
    cur = read_curation(a.curation)
    used_cur = set()

    by_q = defaultdict(list)
    for p in projs:
        if p["q_chrom"] != "NA":
            by_q[p["query"]].append(p)

    # ---- 1. which loci are family loci, and which one is ANC in each query ----
    fam = []
    for L in loci:
        q, ch, st = L["query"], L["q_chrom"], L["q_strand"]
        s, e = int(L["q_start"]), int(L["q_end"])
        here = [p for p in by_q[q] if p["q_chrom"] == ch and p["q_strand"] == st
                and int(p["q_start"]) < e and int(p["q_end"]) > s]
        orth_anchor = any(p["ref_gene"] == a.anchor and p["q_gene"] not in ("NA", "None", "") for p in here)
        acts = [x[0] for x in cur.get((q, ch, s), [])]
        if "anchor" in acts or orth_anchor:
            kind = "anchor"
        elif any(p["ref_gene"] in exn for p in here):
            kind = "exendin"
        else:
            continue
        if acts:
            used_cur.add((q, ch, s))
        fam.append(dict(L=L, q=q, ch=ch, st=st, s=s, e=e, here=here, kind=kind,
                        forced="anchor" in acts))

    anchor = {}
    for f in fam:
        if f["kind"] == "anchor":
            old = anchor.get(f["q"])
            if old is None or (f["forced"] and not old["forced"]):
                anchor[f["q"]] = f
    for f in fam:                       # a second ANC-hit locus in one query is treated as an exendin copy
        if f["kind"] == "anchor" and anchor[f["q"]] is not f:
            f["kind"] = "exendin"

    # AQP5 (orientation): best-status projection on the ANC contig, per query
    flank = {}
    for p in projs:
        if p["ref_gene"] == a.flank and p["q_chrom"] != "NA" and p["query"] in anchor \
                and p["q_chrom"] == anchor[p["query"]]["ch"]:
            old = flank.get(p["query"])
            if old is None or RANK.get(p["status"], 99) < RANK.get(old["status"], 99):
                flank[p["query"]] = p

    # ---- 2. name, choose a projection, QC ----
    rows, used = [], set()
    for f in sorted(fam, key=lambda f: (f["q"], f["kind"] != "anchor", f["ch"], f["s"])):
        q, anc = f["q"], anchor.get(f["q"])
        orient, off = "NA", float("nan")
        if f["kind"] == "anchor":
            base, off, orient = f"{q}_ANC", 0.0, "self"
        elif anc and anc["ch"] == f["ch"]:
            amid = (anc["s"] + anc["e"]) / 2
            if q in flank:
                sign, orient = (1 if mid(flank[q]) > amid else -1), "AQP5"
            else:
                sign, orient = (1 if anc["st"] == "-" else -1), "ANC_strand"
            off = ((f["s"] + f["e"]) / 2 - amid) * sign / 1000
            base = f"{q}_{'A' if off >= 0 else 'S'}{abs(off):03.0f}"
        else:
            n = sum(1 for r in rows if r["query"] == q and r["orient"] == "unplaced") + 1
            base, orient = f"{q}_U{n}", "unplaced"
        name, k = base, 0
        while name in used:
            name = base + "abcdefghijklmnopqrstuvwxyz"[k]; k += 1
        used.add(name)

        want = {a.anchor} | exn if f["kind"] == "anchor" else exn
        cands = sorted((p for p in f["here"] if p["ref_gene"] in want),
                       key=lambda p: (p["ref_gene"] != a.anchor, (q, p["projection"]) not in seqs,
                                      RANK.get(p["status"], 99), isoform(p["projection"]),
                                      -int(p["cds_bp"] or 0)))
        best = cands[0] if cands else None
        nt = seqs.get((q, best["projection"]), "") if best else ""
        aa = translate(nt)
        if aa.endswith("*"):
            nt, aa = nt[:-3], aa[:-1]

        why = []
        if not best:
            why.append("no_projection")
        else:
            if not nt:
                why.append("no_sequence")
            if best["status"] not in allowed:
                why.append(f"status_{best['status']}")
            if len(nt) % 3:
                why.append("cds_not_mult3")
            if len(nt) < a.min_cds:
                why.append(f"cds<{a.min_cds}")
            if "*" in aa:
                why.append(f"internal_stops={aa.count('*')}")
        labels, notes, repl = [], [], None
        for act, val, note in cur.get((q, f["ch"], f["s"]), []):
            if act == "drop":
                why.append("curation_drop")
            elif act == "replace" and val:
                repl = val
            elif act == "keep":
                # override sequence-quality filters; a missing sequence cannot be overridden
                why = [w for w in why if w in ("no_projection", "no_sequence")]
            elif act == "label" and val:
                labels.append(val)
            if note:
                notes.append(note)
        if repl:
            raa, rnt = read_one(repl + ".aa.fa"), read_one(repl + ".cds.fa")
            if raa and rnt:
                aa, nt = raa.rstrip("*"), rnt
                why = [w for w in why if w == "curation_drop"]   # the rescued model replaces the projection
                notes.append(f"sequence from {repl}.aa.fa / .cds.fa")
            else:
                print(f"[warn] replace: cannot read {repl}.aa.fa or {repl}.cds.fa", file=sys.stderr)
                why.append("replace_failed")
        if labels:
            name = name + "_" + "_".join(labels)
        rows.append(dict(
            name=name, query=q, kind=f["kind"], q_chrom=f["ch"], q_start=f["s"], q_end=f["e"],
            q_strand=f["st"], kb_from_ANC=f"{off:.1f}", orient=orient,
            projection=best["projection"] if best else "NA", status=best["status"] if best else "NA",
            cds_bp=len(nt), included="yes" if not why else "no", reason=";".join(why),
            curation_note="; ".join(notes), protein=aa, _nt=nt))

    for k in cur:
        if k not in used_cur:
            print(f"[warn] curation entry matched no locus: {k}", file=sys.stderr)

    os.makedirs(a.outdir, exist_ok=True)
    cols = ["name", "query", "kind", "q_chrom", "q_start", "q_end", "q_strand", "kb_from_ANC", "orient",
            "projection", "status", "cds_bp", "included", "reason", "curation_note", "protein"]
    with open(os.path.join(a.outdir, "loci_curated.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    keep = [r for r in rows if r["included"] == "yes"]
    with open(os.path.join(a.outdir, "curated_aa.fa"), "w") as fa, \
         open(os.path.join(a.outdir, "curated_cds.fa"), "w") as fn:
        for r in keep:
            fa.write(f">{r['name']}\n{r['protein'].replace('*', 'X')}\n")
            fn.write(f">{r['name']}\n{r['_nt']}\n")
    nq = len({r["query"] for r in keep})
    print(f"[info] {len(rows)} loci; {len(keep)} in the tree ({sum(r['kind']=='anchor' for r in keep)} ANC, "
          f"{sum(r['kind']=='exendin' for r in keep)} exendin copies) from {nq} species", file=sys.stderr)


if __name__ == "__main__":
    main()
