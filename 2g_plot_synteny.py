#!/usr/bin/env python3
"""
The synteny plot of 4d_plot_exendin_synteny.py, drawn from the step-2 curated loci
(2_exendin_origin.md): every exendin-family locus in curated/loci_curated.tsv from
2b_curate_exendin_loci.py, between the two conserved boundary genes, in species-tree
order.

Inputs
  - curated/loci_curated.tsv (2b). Every row is drawn: loci that entered the tree
    are filled, loci left out are hollow (reason in the hover text), and loci
    dropped by 2_curation.tsv are dashed.
  - the boundary genes (SLC4A8, AQP5) are taken from the TOGA2 run that step 2 used,
    so contig names and coordinates match 2b's: the TOGA2 ortholog projection of each
    boundary gene on each locus contig, nearest the loci (one locus per contig,
    overlapping isoforms/chains merged, best TOGA2 status). Give each as NAME=GENE_ID,
    or NAME alone with --gff to look it up in the HelSus annotation.

Colour = what 2b named the locus: ANC, a curation label (EX1, EX2, EX3, EX4 ...), or
'unnamed copy'. --groups TSV (name<TAB>group) overrides this, e.g. to colour by the
2d nomenclature verdict. Labels are the 2b name without the species prefix
(A076_EX2, S118, U1, ANC), so they read against the step-2 tree.

Drawing, orientation, clipping and the tree are shared with 4d (imported from the
same directory). Read-only on everything; writes the SVG, optional PNG and --layout.

Usage
  python3 2g_plot_synteny.py \
      --loci curated/loci_curated.tsv \
      --run ${TOGADIR}/target_HelSus \
      --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg exendin_origin_synteny.svg --png exendin_origin_synteny.png \
      --layout exendin_origin_synteny.tsv
"""
import argparse, csv, importlib.util, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("synteny4d", os.path.join(HERE, "4d_plot_exendin_synteny.py"))
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)


def group_of(name, kind):
    """2b name -> colour group: ANC, the last curation label (EX2 ...), or 'unnamed copy'."""
    if kind == "anchor" or name.endswith("_ANC"):
        return "ANC"
    m = re.search(r"_(EX\d+)$", name)
    return m.group(1) if m else "unnamed copy"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True, help="curated/loci_curated.tsv from 2b")
    ap.add_argument("--run", required=True, help="${TOGADIR}/target_HelSus (for the boundary genes)")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5=LOC_00012326",
                    help="two NAME[=REF_GENE_ID], left one first")
    ap.add_argument("--gff", help="HelSus reference GFF, to look up boundary names given without an ID")
    ap.add_argument("--groups", help="TSV name<TAB>group, overrides the colour group")
    ap.add_argument("--title", default="Exendin-family loci between the conserved boundary genes "
                                       "(step 2: curated TOGA2 loci)")
    S.add_common_args(ap)
    a = ap.parse_args()

    with open(a.loci) as fh:
        loci = [r for r in csv.DictReader(fh, delimiter="\t") if r["q_chrom"] not in ("NA", "")]
    if not loci:
        sys.exit("[error] no loci with coordinates in " + a.loci)
    override = {}
    if a.groups:
        for line in open(a.groups):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] not in ("name", "tip"):
                override[f[0]] = f[1]

    ref = [r for r in loci if r["query"] == a.reference_query]
    near = ((ref[0]["q_chrom"], sum(int(r["q_start"]) for r in ref) / len(ref)) if ref else ("", 0))
    bnames, bids = [], {}
    for item in a.boundaries.split(","):
        name, _, gid = item.partition("=")
        if not gid:
            if not a.gff:
                sys.exit(f"[error] boundary {name} has no gene ID; give {name}=GENE_ID or --gff")
            gid = S.gff_lookup(a.gff, name, near)
        bnames.append(name); bids[gid] = name
        print(f"[info] boundary {name} = {gid}", file=sys.stderr)
    if len(bnames) != 2:
        sys.exit("[error] --boundaries needs exactly two genes")

    contigs = defaultdict(lambda: defaultdict(list))
    for r in loci:
        nm = r["name"]
        lab = nm.split("_", 1)[1] if "_" in nm else nm
        reason = r.get("reason", "")
        contigs[r["query"]][r["q_chrom"]].append(dict(
            kind="exendin", name=nm, label=lab, group=override.get(nm, group_of(nm, r.get("kind", ""))),
            chrom=r["q_chrom"], start=int(r["q_start"]), end=int(r["q_end"]), strand=r["q_strand"],
            status=r.get("status", ""), filled=r.get("included") == "yes",
            dashed="curation_drop" in reason,
            hover=f'{r.get("kind", "")} | status {r.get("status", "")} | {r.get("cds_bp", "")} bp | '
                  f'{"in tree" if r.get("included") == "yes" else "not in tree: " + reason} | '
                  f'{r.get("kb_from_ANC", "")} kb from ANC | {r.get("curation_note", "")}'))

    queries = sorted(set(contigs) | {d[6:] for d in os.listdir(a.run) if d.startswith("query_")})
    for q in queries:
        bnd = S.boundary_projections(os.path.join(a.run, "query_" + q), set(bids))
        mid = {c: sum((f["start"] + f["end"]) / 2 for f in fs) / len(fs) for c, fs in contigs[q].items()}
        for gid, name in bids.items():
            ps = bnd.get(gid, [])
            on = defaultdict(list)
            for p in ps:
                if p[0] in mid:
                    on[p[0]].append(p)
            picks = ([S.choose_boundary(v, mid[c]) for c, v in on.items()] if on
                     else ([S.choose_boundary(ps, None)] if ps else []))
            for p in picks:
                contigs[q][p[0]].append(dict(kind="boundary", name=name, label=name, group=name, chrom=p[0],
                                             start=p[1], end=p[2], strand=p[3], status=p[4], filled=True,
                                             dashed=False, hover=f"status {p[4]} {p[5]}"))
    missing = [q for q in queries if not any(f["kind"] == "exendin" for fs in contigs[q].values() for f in fs)]
    if missing:
        print(f"[info] species with no step-2 locus: {', '.join(missing)}", file=sys.stderr)

    S.finish(a, contigs, queries, bnames, a.title,
             f"Tracks oriented with {bnames[0]} on the left; arrows give strand after orientation; "
             f"labels are 2b names; // = boundary gene clipped to {a.clip_kb:g} kb.",
             legend_extra=[("in the step-2 tree", True, False),
                           ("left out by 2b filters", False, False),
                           ("dropped in 2_curation.tsv", False, True)])


if __name__ == "__main__":
    main()
