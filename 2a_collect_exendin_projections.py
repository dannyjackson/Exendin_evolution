#!/usr/bin/env python3
"""
Collect TOGA2 projections of selected reference genes across all queries of one
reference run, and summarise them per projection and per query locus.

Read-only on the TOGA2 directory. Writes only to --outdir (default: stdout for
the projection table; locus table and FASTA only if --outdir is given).

Usage
  python3 collect_exendin_projections.py \
      --run ${TOGADIR}/target_HelSus \
      --genes LOC_00012323,LOC_00012324,LOC_00012325 \
      --anchor LOC_00012323 \
      [--outdir ${PROJDIR}/toga2_exendin] [--fasta-status FI,I,PI]

Outputs (with --outdir)
  projections.tsv  one row per projection (query x reference transcript x chain)
  loci.tsv         one row per distinct query locus (overlapping projections
                   on the same contig/strand merged), with the reference genes
                   that hit it and the best status
  <gene>.protein.fa / <gene>.nucleotide.fa  projections whose status is in
                   --fasta-status, headers  >QUERY|projection|status|q_gene
"""
import argparse, gzip, os, re, sys
from collections import defaultdict

STATUS_RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gene_of(proj_or_tx):
    return proj_or_tx.split("-mRNA")[0]


def read_loss(path, genes):
    status = {}
    with opener(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 3 and t[0] == "PROJECTION" and gene_of(t[1]) in genes:
                status[t[1]] = t[2]
    return status


def read_orth(path, genes):
    """q_transcript (projection) -> (q_gene, class); also ref genes with one2zero."""
    proj2 = {}
    zero = defaultdict(set)
    with opener(path) as fh:
        next(fh, None)
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 5 or t[0] not in genes:
                continue
            if t[3] == "None":
                zero[t[0]].add(t[4])
            else:
                proj2[t[3]] = (t[2], t[4])
    return proj2, zero


def read_bed(path, genes):
    out = {}
    with opener(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 12 or gene_of(t[3]) not in genes:
                continue
            sizes = [int(x) for x in t[10].strip(",").split(",") if x]
            out[t[3]] = (t[0], int(t[1]), int(t[2]), t[5], int(t[9]), sum(sizes))
    return out


def read_fasta(path, wanted):
    """Return {projection: seq} for headers whose first token equals or starts with a wanted projection id."""
    seqs, cur = {}, None
    if not os.path.exists(path):
        return seqs
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                key = re.split(r"[\s|]", line[1:].strip())[0]
                cur = key if key in wanted else None
                if cur is not None:
                    seqs.setdefault(cur, [])
            elif cur is not None:
                seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="e.g. ${TOGADIR}/target_HelSus")
    ap.add_argument("--genes", required=True, help="comma-separated reference gene IDs")
    ap.add_argument("--anchor", default=None, help="synteny-anchor gene; flags projections on its chain")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--fasta-status", default="FI,I,PI")
    a = ap.parse_args()

    genes = set(a.genes.split(","))
    if a.anchor:
        genes.add(a.anchor)
    fasta_ok = set(a.fasta_status.split(","))
    queries = sorted(d for d in os.listdir(a.run) if d.startswith("query_"))

    rows, zero_rows = [], []
    for qd in queries:
        q = qd.replace("query_", "")
        p = os.path.join(a.run, qd)
        status = read_loss(os.path.join(p, "loss_summary.tsv"), genes)
        orth, zero = read_orth(os.path.join(p, "orthology_classification.tsv"), genes)
        bed = read_bed(os.path.join(p, "query_annotation.bed"), genes)
        anchor_chains = {pr.split("#", 1)[1] for pr in status if a.anchor and gene_of(pr) == a.anchor}
        for g, txs in zero.items():
            for cls in txs:
                zero_rows.append((q, g, cls))
        for pr in sorted(set(status) | set(orth) | set(bed)):
            tx, chain = pr.split("#", 1)
            qg, cls = orth.get(pr, ("NA", "not_classified"))
            chrom, s, e, strand, nex, cds = bed.get(pr, ("NA", -1, -1, ".", 0, 0))
            shares = "yes" if set(chain.split(",")) & {c for ac in anchor_chains for c in ac.split(",")} else "no"
            rows.append([q, gene_of(tx), tx, pr, chain, status.get(pr, "NA"), cls, qg,
                         chrom, s, e, strand, nex, cds, shares])

    hdr = ["query", "ref_gene", "ref_transcript", "projection", "chain", "status", "orthology_class",
           "q_gene", "q_chrom", "q_start", "q_end", "q_strand", "n_exons", "cds_bp", "on_anchor_chain"]

    out = sys.stdout
    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)
        out = open(os.path.join(a.outdir, "projections.tsv"), "w")
    print("\t".join(hdr), file=out)
    for r in rows:
        print("\t".join(map(str, r)), file=out)
    for q, g, cls in zero_rows:
        print("\t".join([q, g, "NA", "NA", "NA", "NA", "one2zero", "None"] + ["NA"] * 7), file=out)
    if not a.outdir:
        return
    out.close()

    # ---- merge overlapping projections into query loci ----
    by = defaultdict(list)
    for r in rows:
        if r[8] != "NA":
            by[(r[0], r[8], r[11])].append(r)
    with open(os.path.join(a.outdir, "loci.tsv"), "w") as lo:
        print("\t".join(["query", "locus_id", "q_chrom", "q_start", "q_end", "q_strand", "ref_genes_hitting",
                         "best_status", "q_gene_labels", "chains", "on_anchor_chain", "n_projections"]), file=lo)
        for (q, chrom, strand), rs in sorted(by.items()):
            rs.sort(key=lambda r: r[9])
            clusters, cur, cur_end = [], [], -1
            for r in rs:
                if cur and r[9] > cur_end:
                    clusters.append(cur); cur = []
                cur.append(r); cur_end = max(cur_end, r[10]) if cur[:-1] else r[10]
            if cur:
                clusters.append(cur)
            for i, c in enumerate(clusters, 1):
                best = min((r[5] for r in c), key=lambda s: STATUS_RANK.get(s, 99))
                print("\t".join(map(str, [
                    q, f"{q}_{chrom}_{min(r[9] for r in c)}", chrom, min(r[9] for r in c), max(r[10] for r in c),
                    strand, ",".join(sorted({r[1] for r in c})), best,
                    ",".join(sorted({r[7] for r in c})), ",".join(sorted({r[4] for r in c})),
                    "yes" if any(r[14] == "yes" for r in c) else "no", len(c)])), file=lo)

    # ---- sequences ----
    wanted = defaultdict(dict)
    for r in rows:
        if r[5] in fasta_ok:
            wanted[r[0]][r[3]] = r
    for kind in ("protein", "nucleotide"):
        handles = {}
        for qd in queries:
            q = qd.replace("query_", "")
            if not wanted[q]:
                continue
            seqs = read_fasta(os.path.join(a.run, qd, f"{kind}.fa.gz"), set(wanted[q]))
            for pr, seq in seqs.items():
                r = wanted[q][pr]
                fh = handles.get(r[1])
                if fh is None:
                    fh = handles[r[1]] = open(os.path.join(a.outdir, f"{r[1]}.{kind}.fa"), "w")
                fh.write(f">{q}|{pr}|{r[5]}|{r[7]}\n{seq}\n")
        for fh in handles.values():
            fh.close()


if __name__ == "__main__":
    main()
