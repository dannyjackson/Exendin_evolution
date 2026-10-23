#!/usr/bin/env python3
"""
Pick one representative coding sequence per query locus from the output of
2a_collect_exendin_projections.py, name each locus by its position relative to
the synteny-anchor gene, run basic CDS QC, and write alignment-ready FASTA.

Why one per locus: in a tandem array, the exendin-2 and exendin-4 reference
transcripts both project onto every copy, so one query gene can appear as
several projections. The gene tree needs each query copy once.

Tip names (safe for IQ-TREE / HyPhy / codeml: letters, digits, underscore):
  <Query>_R<kb><s|o>   locus on the anchor's contig, on the side where the
                       HelSus chr2 exendin cluster sits (R) or the other side
                       (L), <kb> kb from the anchor midpoint, on the same (s)
                       or opposite (o) strand to the anchor.
                       Orientation follows HelSus: anchor on the minus strand,
                       cluster at higher coordinates.
  <Query>_U<n>         locus not on the anchor's contig (or no anchor found)
  <Query>_ANC          the anchor gene itself (outgroup / rooting)
A suffix letter (a, b, ...) is added if two loci round to the same name.

Inputs  (--collect-dir = the 2a --outdir)
  projections.tsv, loci.tsv, <gene>.nucleotide.fa
Outputs (--outdir)
  loci_selected.tsv     every candidate locus, the projection chosen, QC, and
                        whether it was included
  exendin_cds.fa        included exendin loci, CDS without the final stop codon
  anchor_cds.fa         included anchor orthologs (same format)
  exendin_plus_anchor_cds.fa
  *_aa.fa               translations of the above (for checking peptides)

Usage
  python3 2b_select_exendin_loci.py \
      --collect-dir ${PROJDIR}/analyses/toga2/collect_HelSus \
      --genes LOC_00012324,LOC_00012325 --anchor LOC_00012323 \
      --outdir ${PROJDIR}/analyses/toga2/loci_HelSus [--status FI,I] [--min-cds 150]
"""
import argparse, csv, os, re, sys
from collections import defaultdict

RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODE[x + y + z] = _aa[16 * i + 4 * j + k]


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_collect_fasta(path):
    """2a headers: >QUERY|projection|status|q_gene"""
    out, key = {}, None
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].strip().split("|")
                key = (f[0], f[1])
                out[key] = []
            elif key:
                out[key].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def clean_cds(seq):
    s = re.sub(r"[^ACGTNacgtn]", "", seq).upper()   # drop gaps / exon separators if present
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--collect-dir", required=True)
    ap.add_argument("--genes", required=True, help="exendin reference genes, comma-separated")
    ap.add_argument("--anchor", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--status", default="FI,I", help="statuses allowed into the alignment")
    ap.add_argument("--min-cds", type=int, default=150, help="minimum CDS length (bp) to include")
    ap.add_argument("--max-internal-stops", type=int, default=0)
    a = ap.parse_args()

    genes = set(a.genes.split(","))
    allowed = set(a.status.split(","))
    projs = read_tsv(os.path.join(a.collect_dir, "projections.tsv"))
    loci = read_tsv(os.path.join(a.collect_dir, "loci.tsv"))
    seqs = {}
    for g in genes | {a.anchor}:
        seqs.update(read_collect_fasta(os.path.join(a.collect_dir, f"{g}.nucleotide.fa")))
    if not seqs:
        sys.exit("No sequences found: check that 2a was run with --outdir and wrote *.nucleotide.fa")

    # anchor position per query (best-status anchor projection)
    anchor = {}
    for p in projs:
        if p["ref_gene"] == a.anchor and p["q_chrom"] != "NA":
            cur = anchor.get(p["query"])
            if cur is None or RANK.get(p["status"], 99) < RANK.get(cur["status"], 99):
                anchor[p["query"]] = p

    by_q = defaultdict(list)
    for p in projs:
        if p["q_chrom"] != "NA":
            by_q[p["query"]].append(p)

    rows, used_names = [], set()
    for L in loci:
        hit = set(L["ref_genes_hitting"].split(","))
        is_anchor = hit == {a.anchor}
        if not (hit & genes or is_anchor):
            continue
        q, chrom, strand = L["query"], L["q_chrom"], L["q_strand"]
        s, e = int(L["q_start"]), int(L["q_end"])
        cands = [p for p in by_q[q] if p["q_chrom"] == chrom and p["q_strand"] == strand
                 and int(p["q_start"]) < e and int(p["q_end"]) > s
                 and (p["ref_gene"] in genes or (is_anchor and p["ref_gene"] == a.anchor))]

        def key(p):
            return (RANK.get(p["status"], 99), -int(p["cds_bp"] or 0), p["on_anchor_chain"] != "yes")
        cands.sort(key=key)
        best = next((p for p in cands if (q, p["projection"]) in seqs), cands[0] if cands else None)

        # position relative to anchor
        anc = anchor.get(q)
        if is_anchor:
            base, off_kb, rel = f"{q}_ANC", 0.0, "self"
        elif anc and anc["q_chrom"] == chrom:
            amid = (int(anc["q_start"]) + int(anc["q_end"])) / 2
            off_kb = ((s + e) / 2 - amid) / 1000 * (1 if anc["q_strand"] == "-" else -1)
            rel = "s" if strand == anc["q_strand"] else "o"
            base = f"{q}_{'R' if off_kb >= 0 else 'L'}{abs(off_kb):03.0f}{rel}"
        else:
            off_kb, rel = float("nan"), "NA"
            n = sum(1 for r in rows if r["query"] == q and r["tip"].startswith(f"{q}_U")) + 1
            base = f"{q}_U{n}"
        name, k = base, 0
        while name in used_names:
            name = base + "abcdefghijklmnopqrstuvwxyz"[k]; k += 1
        used_names.add(name)

        nt = clean_cds(seqs.get((q, best["projection"]), "")) if best else ""
        aa = translate(nt)
        final_stop = aa.endswith("*")
        if final_stop:
            nt, aa = nt[:-3], aa[:-1]
        istops = aa.count("*")
        reasons = []
        if not best:
            reasons.append("no_projection")
        elif best["status"] not in allowed:
            reasons.append(f"status_{best['status']}")
        if best and not nt:
            reasons.append("no_sequence")
        if nt and len(nt) % 3:
            reasons.append("len_not_mult3")
        if nt and len(nt) < a.min_cds:
            reasons.append(f"cds<{a.min_cds}")
        if istops > a.max_internal_stops:
            reasons.append(f"internal_stops={istops}")
        if nt and not aa.startswith("M"):
            reasons.append("no_start_M(warn)")
        include = not [r for r in reasons if not r.endswith("(warn)")]

        rows.append(dict(
            tip=name, query=q, q_chrom=chrom, q_start=s, q_end=e, q_strand=strand,
            offset_kb=f"{off_kb:.1f}", strand_vs_anchor=rel, ref_genes_hitting=L["ref_genes_hitting"],
            chosen_projection=best["projection"] if best else "NA",
            status=best["status"] if best else "NA", on_anchor_chain=best["on_anchor_chain"] if best else "NA",
            cds_bp=len(nt), final_stop=final_stop, internal_stops=istops,
            included="yes" if include else "no", notes=";".join(reasons), aa=aa, nt=nt, is_anchor=is_anchor))

    os.makedirs(a.outdir, exist_ok=True)
    cols = ["tip", "query", "q_chrom", "q_start", "q_end", "q_strand", "offset_kb", "strand_vs_anchor",
            "ref_genes_hitting", "chosen_projection", "status", "on_anchor_chain", "cds_bp", "final_stop",
            "internal_stops", "included", "notes", "aa"]
    rows.sort(key=lambda r: (r["query"], r["tip"]))
    with open(os.path.join(a.outdir, "loci_selected.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    def dump(fname, sel, field):
        with open(os.path.join(a.outdir, fname), "w") as fh:
            for r in sel:
                fh.write(f">{r['tip']}\n{r[field]}\n")
    ex = [r for r in rows if r["included"] == "yes" and not r["is_anchor"]]
    an = [r for r in rows if r["included"] == "yes" and r["is_anchor"]]
    dump("exendin_cds.fa", ex, "nt"); dump("exendin_aa.fa", ex, "aa")
    dump("anchor_cds.fa", an, "nt"); dump("anchor_aa.fa", an, "aa")
    dump("exendin_plus_anchor_cds.fa", ex + an, "nt"); dump("exendin_plus_anchor_aa.fa", ex + an, "aa")
    print(f"loci considered: {len(rows)}  exendin included: {len(ex)}  anchor included: {len(an)}", file=sys.stderr)


if __name__ == "__main__":
    main()
