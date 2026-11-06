#!/usr/bin/env python3
"""
The synteny plot of 4d_plot_exendin_synteny.py, drawn from the step-3 interval scans
(3_copy_number_survey.md) instead of TOGA2: every exendin-family locus that
2f_scan_interval.py found in each species, between the two boundary genes, in
species-tree order.

Inputs
  - scan_<label>_loci.tsv from 2f_scan_interval.py / 3a_scan_species.sh (--scan-dir).
    Coordinates in these files are already contig coordinates.
  - the species table used by 3a (label accession contig known_min known_max note),
    which says which contig each scan was run on. HelSus and HemCap, scanned in step 4
    of the survey, are in the built-in table below; pass your own with --species-table.
  - boundary genes (SLC4A8, AQP5) from each species' NCBI annotation: the GFF under
    ${GENOMEDIR}/<label>_*/ncbi_dataset/data/*/genomic.gff, or --gff-map. 3a downloads
    genomes without annotation, so fetch the GFFs once (see 3_copy_number_survey.md).
    The gene on the scanned contig whose gene/Name equals the boundary name, or whose
    gene/Name/description matches its pattern, nearest the scan loci, is used. A
    species with no GFF is drawn without boundary genes.

How loci are drawn (state from 2f):
  colour   lineage of the best-matching query (best_query, species prefix removed)
  filled   candidate_gene: both an exon-1 and a hormone-domain hit, no flag
  hollow   candidate_gene flagged premature_stop or domain_starts_X (disabled)
  dashed   N_term_only / domain_only: one part only (remnant, or a diverged exon)

Drawing, orientation, clipping and the tree are shared with 4d (imported from the same
directory). Read-only on everything; writes the SVG, optional PNG and --layout.

Usage
  python3 3b_plot_scan_synteny.py --scan-dir ${PROJDIR}/analyses/copy_number \
      --genomedir ${GENOMEDIR} \
      --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
      --svg scan_synteny.svg --png scan_synteny.png --layout scan_synteny.tsv
"""
import argparse, csv, glob, gzip, importlib.util, os, re, sys
from collections import defaultdict
from urllib.parse import unquote

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("synteny4d", os.path.join(HERE, "4d_plot_exendin_synteny.py"))
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)

# label  accession  contig  known_min  known_max  note   (3a's table, plus the step-4 species)
TABLE = """\
AnoCar  GCF_035594765.1  NC_085842.1  155087274  155138943  SLC4A8 not found; exendin + AQP5
CanAsp  GCF_035149785.1  NC_086154.1  135220948  135362723  AQP5 .. exendin .. SLC4A8
PogVit  GCF_051106095.1  NC_135784.1  163211717  163460527  AQP5 .. unknown LOC .. exendin .. SLC4A8
EubMac  GCF_028583425.1  NC_072808.1  3391923    3535608    AQP5 .. 3 annotated exendins .. SLC4A8
HelSus  GCF_057233095.1  NC_139459.1  188950000  189850000  step 4 positive control
HemCap  NA               NC_069658.1  429200000  429600000  step 4
"""

DEFAULT_BOUNDARIES = "SLC4A8=SLC4A8|solute carrier family 4 member 8,AQP5=AQP5|aquaporin[- ]?5"


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_table(text):
    rows = {}
    for line in text.splitlines():
        f = line.split()
        if not f or f[0].startswith("#") or len(f) < 3:
            continue
        rows[f[0]] = dict(label=f[0], accession=f[1], contig=f[2],
                          lo=int(f[3]) if len(f) > 3 and f[3].isdigit() else None,
                          hi=int(f[4]) if len(f) > 4 and f[4].isdigit() else None)
    return rows


def find_gff(genomedir, label, gff_map):
    if label in gff_map:
        return gff_map[label]
    if not genomedir:
        return None
    hits = sorted(glob.glob(os.path.join(genomedir, f"{label}_*", "ncbi_dataset", "data", "*", "genomic.gff*")))
    return hits[0] if hits else None


def gff_genes(gff, contig):
    """gene features on one contig: (start, end, strand, gene, Name, description, ID)."""
    out = []
    with opener(gff) as fh:
        for line in fh:
            if not line.startswith(contig + "\t"):
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 9 or t[2] != "gene":
                continue
            at = {k: unquote(v) for k, v in (kv.split("=", 1) for kv in t[8].strip(";").split(";") if "=" in kv)}
            out.append((int(t[3]), int(t[4]), t[6], at.get("gene", ""), at.get("Name", ""),
                        at.get("description", ""), at.get("ID", "")))
    return out


def pick_boundary(genes, name, pattern, near):
    pat = re.compile(rf"(?<![A-Za-z0-9])(?:{pattern})(?![A-Za-z0-9])", re.I)
    cands = []
    for g in genes:
        exact = name.upper() in (g[3].upper(), g[4].upper())
        if exact or any(pat.search(x) for x in g[3:6]):
            cands.append((not exact, abs((g[0] + g[1]) / 2 - near), g))
    return min(cands)[2] if cands else None


def short(q):
    """best_query -> lineage label: drop the species prefix, keep the last part (HelSus_A076_EX2 -> EX2)."""
    rest = q.split("_", 1)[1] if "_" in q else q
    return rest.split("_")[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", default=".", help="where scan_<label>_loci.tsv are (step 3 OUT)")
    ap.add_argument("--species-table", help="3a-format table (default: built-in, see source)")
    ap.add_argument("--genomedir", help="${GENOMEDIR}, to find each species' NCBI GFF")
    ap.add_argument("--gff-map", help="TSV label<TAB>gff path, overrides --genomedir")
    ap.add_argument("--boundaries", default=DEFAULT_BOUNDARIES,
                    help="two NAME=REGEX (matched to gene/Name/description), left one first")
    ap.add_argument("--title", default="Exendin-family loci between the conserved boundary genes "
                                       "(step-3 tblastn interval scans)")
    S.add_common_args(ap)
    a = ap.parse_args()

    table = read_table(open(a.species_table).read() if a.species_table else TABLE)
    gff_map = {}
    if a.gff_map:
        for line in open(a.gff_map):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and not line.startswith("#"):
                gff_map[f[0]] = f[1]
    bspec = []
    for item in a.boundaries.split(","):
        name, _, pat = item.partition("=")
        bspec.append((name, pat or re.escape(name)))
    if len(bspec) != 2:
        sys.exit("[error] --boundaries needs exactly two genes")
    bnames = [n for n, _ in bspec]

    found = {os.path.basename(p)[5:-len("_loci.tsv")]: p
             for p in glob.glob(os.path.join(a.scan_dir, "scan_*_loci.tsv"))}
    for lab in sorted(set(found) - set(table)):
        print(f"[warn] {found[lab]} has no row in the species table (contig unknown): skipped", file=sys.stderr)
    species = sorted(set(table) & set(found))
    for lab in sorted(set(table) - set(found)):
        print(f"[info] {lab}: no scan_{lab}_loci.tsv in {a.scan_dir}: skipped", file=sys.stderr)
    if not species:
        sys.exit(f"[error] no scan_<label>_loci.tsv in {a.scan_dir} matches the species table")

    contigs = defaultdict(lambda: defaultdict(list))
    for lab in species:
        ctg = table[lab]["contig"]
        with open(found[lab]) as fh:
            loci = list(csv.DictReader(fh, delimiter="\t"))
        for L in loci:
            state = L["state"]
            g = short(L["best_query"])
            contigs[lab][ctg].append(dict(
                kind="exendin", name=f"{lab}_locus{L['locus']}", label=g, group=g, chrom=ctg,
                start=int(L["g_start"]), end=int(L["g_end"]), strand=L["strand"], status=state,
                filled=state == "candidate_gene",
                dashed=state.split(";")[0] in ("N_term_only", "domain_only"),
                hover=f"{state} | parts {L['parts']} | best {L['best_query']} {L['best_pident']}% | "
                      f"stops {L['stops']} | domain start {L['domain_start']}"))
        print(f"[info] {lab}: {len(loci)} scan loci on {ctg}", file=sys.stderr)

        gff = find_gff(a.genomedir, lab, gff_map)
        if not gff:
            print(f"[warn] {lab}: no GFF found; drawn without boundary genes", file=sys.stderr)
            continue
        genes = gff_genes(gff, ctg)
        fs = contigs[lab][ctg]
        near = (sum((f["start"] + f["end"]) / 2 for f in fs) / len(fs) if fs else
                ((table[lab]["lo"] + table[lab]["hi"]) / 2 if table[lab]["lo"] else 0))
        for name, pat in bspec:
            g = pick_boundary(genes, name, pat, near)
            if not g:
                print(f"[info] {lab}: no {name} on {ctg} in {os.path.basename(gff)}", file=sys.stderr)
                continue
            print(f"[info] {lab}: {name} = {g[3] or g[4]} ({g[6]}) {ctg}:{g[0]}-{g[1]} {g[2]}"
                  f"{'  [description: ' + g[5] + ']' if g[5] else ''}", file=sys.stderr)
            fs.append(dict(kind="boundary", name=name, label=name, group=name, chrom=ctg, start=g[0], end=g[1],
                           strand=g[2], status="", filled=True, dashed=False,
                           hover=f"{g[3] or g[4]} {g[6]} {g[5]}"))

    S.finish(a, contigs, species, bnames, a.title,
             f"Tracks oriented with {bnames[0]} on the left; arrows give strand after orientation; "
             f"colour = best-matching lineage query; // = boundary gene clipped to {a.clip_kb:g} kb.",
             legend_extra=[("candidate gene, intact", True, False),
                           ("candidate gene, disabled (stop / domain start)", False, False),
                           ("remnant: one part only", False, True)])


if __name__ == "__main__":
    main()
