#!/usr/bin/env python3
"""
Plot, for every species, where each TOGA2-annotated exendin gene sits on the
orthologous chromosome, between the two conserved boundary genes (AQP5 and
SLC4A8). One track per contig, so tandem duplications show up as
several genes between the boundaries, and copies elsewhere as extra tracks.

Inputs are TOGA2 only:
  - the exendin genes: tips.tsv from 4b_collect_toga2_exendins.py (every row is
    drawn; rows with included = no are drawn hollow, so status filters stay visible);
  - the boundary genes: TOGA2 ortholog projections of AQP5 and SLC4A8 in each query
    (orthology_classification.tsv + query_annotation.bed + loss_summary.tsv).

Boundary genes are given as NAME=REFERENCE_GENE_ID. A NAME with no ID is looked up
in --gff: the gene whose attributes contain NAME as a whole word, nearest to the
HelSus exendin genes. The IDs actually used are printed, so check them.

Orientation: every track is drawn with the first --boundaries gene on the left
(default SLC4A8, which keeps HelSus chr2 in its own orientation, as in the
"SLC4A8 - ANC - AQP5" interval). If both boundaries are on a track, their order
decides it; otherwise a track is flipped when its boundary gene (or, failing that,
its exendin genes) lies on the opposite strand from HelSus. The flip is recorded in layout.tsv.

Read-only on the TOGA2 directory. Writes the SVG (stdout or --svg), optional PNG,
and --layout (one row per drawn feature).

Usage
  python3 4d_plot_exendin_synteny.py \
      --run ${TOGADIR}/target_HelSus \
      --tips toga2_exendins/tips.tsv \
      --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
      --svg exendin_synteny.svg --png exendin_synteny.png --layout exendin_synteny.tsv
"""
import argparse, csv, gzip, html, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict
from urllib.parse import unquote

RANK = {"FI": 0, "I": 1, "PI": 2, "UL": 3, "M": 4, "L": 5, "PG": 6, "N": 7}
PALETTE = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F3C8D", "#8C564B"]
BOUNDARY_COL = "#555555"


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def gene_of(pr):
    return pr.split("-mRNA")[0]


def gff_lookup(gff, name, near):
    """Reference gene ID whose gene/mRNA attributes contain NAME as a word, nearest to `near`."""
    pat = re.compile(rf"(?<![A-Za-z0-9]){re.escape(name)}(?![A-Za-z0-9])", re.I)
    genes, parent, hits = {}, {}, set()
    with opener(gff) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(t) < 9 or t[2] not in ("gene", "mRNA", "transcript"):
                continue
            at = dict(kv.split("=", 1) for kv in t[8].strip().strip(";").split(";") if "=" in kv)
            fid = at.get("ID", "")
            if t[2] == "gene":
                genes[fid] = (t[0], int(t[3]), int(t[4]))
            else:
                parent[fid] = at.get("Parent", "").split(",")[0]
            if any(pat.search(unquote(v)) for k, v in at.items() if k not in ("ID", "Parent")):
                hits.add(fid if t[2] == "gene" else parent[fid])
    hits = [h for h in hits if h in genes]
    if not hits:
        sys.exit(f"[error] no gene in {gff} matches '{name}'; give it as {name}=GENE_ID")
    chrom, pos = near

    def d(h):
        c, s, e = genes[h]
        return (c != chrom, abs((s + e) / 2 - pos))
    hits.sort(key=d)
    if len(hits) > 1:
        print(f"[info] {name}: {len(hits)} GFF matches; using the one nearest the HelSus exendins "
              f"({hits[0]} at {genes[hits[0]][0]}:{genes[hits[0]][1]}); others: {','.join(hits[1:6])}",
              file=sys.stderr)
    return hits[0]


def boundary_projections(qdir, ref_ids):
    """{ref_gene: [(chrom, start, end, strand, status, projection)]} for TOGA2 ortholog projections."""
    orth = set()
    p = os.path.join(qdir, "orthology_classification.tsv")
    if os.path.exists(p):
        with opener(p) as fh:
            next(fh, None)
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 5 and t[0] in ref_ids and t[3] != "None":
                    orth.add(t[3])
    status = {}
    p = os.path.join(qdir, "loss_summary.tsv")
    if os.path.exists(p):
        with opener(p) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 3 and t[0] == "PROJECTION" and t[1] in orth:
                    status[t[1]] = t[2]
    out = defaultdict(list)
    p = os.path.join(qdir, "query_annotation.bed")
    if os.path.exists(p):
        with opener(p) as fh:
            for line in fh:
                t = line.rstrip("\n").split("\t")
                if len(t) >= 6 and t[3] in orth:
                    out[gene_of(t[3])].append((t[0], int(t[1]), int(t[2]), t[5], status.get(t[3], "NA"), t[3]))
    return out


def write_png(svg, path, w):
    with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as fh:
        fh.write(svg); tmp = fh.name
    cmds = []
    if shutil.which("rsvg-convert"):
        cmds.append(["rsvg-convert", "-w", str(w), "-o", path, tmp])
    if shutil.which("inkscape"):
        cmds.append(["inkscape", tmp, "--export-type=png", f"--export-filename={path}", f"--export-width={w}"])
    if shutil.which("convert"):
        cmds.append(["convert", "-density", "150", "-background", "white", tmp, path])
    for c in cmds:
        if subprocess.run(c, capture_output=True).returncode == 0 and os.path.exists(path):
            print(f"[info] wrote {path} with {c[0]}", file=sys.stderr); os.unlink(tmp); return
    try:
        import cairosvg
        cairosvg.svg2png(bytestring=svg.encode(), write_to=path, output_width=w)
        print(f"[info] wrote {path} with cairosvg", file=sys.stderr); os.unlink(tmp); return
    except Exception:
        pass
    os.unlink(tmp)
    print("[warn] no SVG converter found; PNG not written (mamba install -c conda-forge cairosvg)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="${TOGADIR}/target_HelSus")
    ap.add_argument("--tips", required=True, help="tips.tsv from 4b")
    ap.add_argument("--boundaries", default="SLC4A8,AQP5=LOC_00012326",
                    help="two NAME[=REF_GENE_ID], left one first")
    ap.add_argument("--gff", help="reference GFF, to look up boundary names given without an ID")
    ap.add_argument("--reference-query", default="HelSus", help="query that sets the orientation")
    ap.add_argument("--order", help="comma-separated query order, or a file with one query per line")
    ap.add_argument("--pad-kb", type=float, default=20.0, help="margin around each track")
    ap.add_argument("--min-kb", type=float, default=0.0,
                    help="minimum axis length (kb), e.g. to match another figure")
    ap.add_argument("--title", default="Exendin genes between the conserved boundary genes (TOGA2, HelSus reference)")
    ap.add_argument("--width", type=int, default=1200)
    ap.add_argument("--svg", help="output SVG (default stdout)")
    ap.add_argument("--png")
    ap.add_argument("--layout", help="TSV of every drawn feature")
    a = ap.parse_args()

    with open(a.tips) as fh:
        tips = [r for r in csv.DictReader(fh, delimiter="\t") if r["q_chrom"] != "NA"]
    if not tips:
        sys.exit("[error] no tips with coordinates in " + a.tips)

    # ---- boundary gene IDs ----
    ref_tips = [t for t in tips if t["query"] == a.reference_query]
    near = ((ref_tips[0]["q_chrom"], sum(int(t["q_start"]) for t in ref_tips) / len(ref_tips))
            if ref_tips else ("", 0))
    bnames, bids = [], {}
    for item in a.boundaries.split(","):
        name, _, gid = item.partition("=")
        if not gid:
            if not a.gff:
                sys.exit(f"[error] boundary {name} has no gene ID; give {name}=GENE_ID or --gff")
            gid = gff_lookup(a.gff, name, near)
        bnames.append(name); bids[gid] = name
        print(f"[info] boundary {name} = {gid}", file=sys.stderr)
    if len(bnames) != 2:
        sys.exit("[error] --boundaries needs exactly two genes")
    left_name = bnames[0]

    # ---- per query: features grouped by contig ----
    by_q = defaultdict(list)
    for t in tips:
        by_q[t["query"]].append(dict(kind="exendin", name=t["tip"], group=t["group"], chrom=t["q_chrom"],
                                     start=int(t["q_start"]), end=int(t["q_end"]), strand=t["q_strand"],
                                     status=t["status"], included=t["included"] == "yes",
                                     orth=t.get("orthology_class", ""), cls=t.get("kind", "ortholog")))
    queries = sorted(d[6:] for d in os.listdir(a.run) if d.startswith("query_"))
    bnd = {q: boundary_projections(os.path.join(a.run, "query_" + q), set(bids)) for q in queries}

    # reference strands, for orientation
    ref_strand = {}
    for gid, name in bids.items():
        hs = bnd.get(a.reference_query, {}).get(gid, [])
        if hs:
            ref_strand[name] = sorted(hs, key=lambda x: RANK.get(x[4], 99))[0][3]
    ex_ref_strand = {t["group"]: t["q_strand"] for t in ref_tips}

    tracks = []
    for q in queries:
        ex = by_q.get(q, [])
        contigs = defaultdict(list)
        for f in ex:
            contigs[f["chrom"]].append(f)
        # boundary projections: keep those on an exendin contig; if none, the best one anywhere
        for gid, name in bids.items():
            ps = bnd[q].get(gid, [])
            on = [p for p in ps if p[0] in contigs]
            for p in (on or sorted(ps, key=lambda x: RANK.get(x[4], 99))[:1]):
                contigs[p[0]].append(dict(kind="boundary", name=name, group=name, chrom=p[0], start=p[1],
                                          end=p[2], strand=p[3], status=p[4], included=True, orth="", cls=""))
        if not contigs:
            tracks.append(dict(query=q, chrom="", feats=[], flip=False, why="no exendin or boundary ortholog"))
            continue
        main = max(contigs, key=lambda c: (sum(f["kind"] == "boundary" for f in contigs[c]),
                                            sum(f["kind"] == "exendin" for f in contigs[c])))
        for c in sorted(contigs, key=lambda c: (c != main, c)):
            fs = contigs[c]
            b = {f["name"]: f for f in fs if f["kind"] == "boundary"}
            flip, why = False, "as is"
            if len(b) == 2:
                flip = (b[left_name]["start"] > b[bnames[1]]["start"]); why = "boundary order"
            elif b:
                n, f = next(iter(b.items()))
                if n in ref_strand:
                    flip = f["strand"] != ref_strand[n]; why = f"{n} strand"
            else:
                votes = [f["strand"] != ex_ref_strand[f["group"]] for f in fs if f["group"] in ex_ref_strand]
                if votes:
                    flip = sum(votes) > len(votes) / 2; why = "exendin strand"
            tracks.append(dict(query=q, chrom=c, feats=fs, flip=flip, why=why, main=c == main))

    # ---- order ----
    if a.order:
        order = ([l.strip() for l in open(a.order) if l.strip()] if os.path.exists(a.order)
                 else a.order.split(","))
    else:
        order = [a.reference_query] + [q for q in queries if q != a.reference_query]
    rank = {q: i for i, q in enumerate(order)}
    tracks.sort(key=lambda t: (rank.get(t["query"], len(rank)), not t.get("main", True), t["chrom"]))

    # ---- coordinates: kb from the left edge of each track (after orientation) ----
    pad = a.pad_kb * 1000
    for t in tracks:
        if not t["feats"]:
            t["span"] = 0; continue
        lo = min(f["start"] for f in t["feats"]) - pad
        hi = max(f["end"] for f in t["feats"]) + pad
        t["lo"], t["hi"], t["span"] = lo, hi, hi - lo
        for f in t["feats"]:
            if t["flip"]:
                f["x0"], f["x1"] = hi - f["end"], hi - f["start"]
                f["dir"] = {"+": "-", "-": "+"}.get(f["strand"], ".")
            else:
                f["x0"], f["x1"] = f["start"] - lo, f["end"] - lo
                f["dir"] = f["strand"]
    maxspan = max([t["span"] for t in tracks] + [a.min_kb * 1000, 1])

    # ---- draw ----
    groups = sorted({f["group"] for t in tracks for f in t["feats"] if f["kind"] == "exendin"})
    col = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(groups)}
    lab_w, left, top, row = 230, 15, 60, 46
    plot_w = a.width - lab_w - left - 30
    X = lambda bp: lab_w + left + bp / maxspan * plot_w
    H = top + row * len(tracks) + 90
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{a.width}" height="{H}" '
           f'font-family="Helvetica, Arial, sans-serif" font-size="11">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="{left}" y="24" font-size="15" font-weight="bold">{html.escape(a.title)}</text>',
           f'<text x="{left}" y="42" font-size="10" fill="#555">Each track drawn with {html.escape(left_name)} '
           f'on the left; arrows give strand after orientation. Hollow = not in the tree '
           f'(status/length filter); dashed outline = paralogous projection.</text>']
    lay = []
    prev_q = None
    for i, t in enumerate(tracks):
        y = top + i * row + 26
        if t["query"] != prev_q and i:
            out.append(f'<line x1="{left}" x2="{a.width - 15}" y1="{y - 30}" y2="{y - 30}" stroke="#ddd"/>')
        qlab = t["query"] if t["query"] != prev_q else ""
        prev_q = t["query"]
        bold = ' font-weight="bold"' if t["query"].startswith("Hel") else ""
        out.append(f'<text x="{left}" y="{y + 4}"{bold}>{html.escape(qlab)}</text>')
        if not t["feats"]:
            out.append(f'<text x="{lab_w + left}" y="{y + 4}" fill="#999" font-style="italic">{t["why"]}</text>')
            continue
        flip_txt = " (rev)" if t["flip"] else ""
        out.append(f'<text x="{left + 62}" y="{y + 4}" font-size="9" fill="#666">'
                   f'{html.escape(t["chrom"][:24])}{flip_txt}</text>')
        out.append(f'<line x1="{X(0):.1f}" x2="{X(t["span"]):.1f}" y1="{y}" y2="{y}" stroke="#999" stroke-width="2"/>')
        # coordinate of the track ends, in Mb on the contig
        lo_mb, hi_mb = (t["hi"], t["lo"]) if t["flip"] else (t["lo"], t["hi"])
        out.append(f'<text x="{X(0):.1f}" y="{y + 16}" font-size="8" fill="#888">{max(lo_mb, 0) / 1e6:.3f} Mb</text>')
        out.append(f'<text x="{X(t["span"]):.1f}" y="{y + 16}" font-size="8" fill="#888" text-anchor="end">'
                   f'{max(hi_mb, 0) / 1e6:.3f} Mb</text>')
        exs = sorted((f for f in t["feats"]), key=lambda f: f["x0"])
        for k, f in enumerate(exs):
            x0, x1 = X(f["x0"]), X(f["x1"])
            if x1 - x0 < 4:
                m = (x0 + x1) / 2; x0, x1 = m - 2, m + 2
            h, tip = 12, 5
            c = BOUNDARY_COL if f["kind"] == "boundary" else col[f["group"]]
            if f["dir"] == "-":
                pts = [(x0, y), (x0 + tip, y - h / 2), (x1, y - h / 2), (x1, y + h / 2), (x0 + tip, y + h / 2)]
            elif f["dir"] == "+":
                pts = [(x0, y - h / 2), (x1 - tip, y - h / 2), (x1, y), (x1 - tip, y + h / 2), (x0, y + h / 2)]
            else:
                pts = [(x0, y - h / 2), (x1, y - h / 2), (x1, y + h / 2), (x0, y + h / 2)]
            fill = c if f["included"] else "white"
            dash = ' stroke-dasharray="2,1.5"' if f["cls"] == "paralog" else ""
            title = (f'{f["name"]} {f["chrom"]}:{f["start"]}-{f["end"]} ({f["strand"]}) status {f["status"]} '
                     f'{f["orth"]}')
            out.append(f'<polygon points="{" ".join(f"{px:.1f},{py:.1f}" for px, py in pts)}" fill="{fill}" '
                       f'stroke="{c}" stroke-width="1.2"{dash}><title>{html.escape(title)}</title></polygon>')
            short = f["name"].split("_", 1)[1] if f["kind"] == "exendin" and "_" in f["name"] else f["name"]
            ly = y - 9 - (k % 2) * 9
            out.append(f'<text x="{(x0 + x1) / 2:.1f}" y="{ly:.1f}" font-size="8" text-anchor="middle" '
                       f'fill="{c}">{html.escape(short)}</text>')
            lay.append([t["query"], t["chrom"], "rev" if t["flip"] else "fwd", t["why"], f["kind"], f["name"],
                        f["group"], f["start"], f["end"], f["strand"], f["status"],
                        "yes" if f["included"] else "no", f"{f['x0'] / 1000:.1f}", f"{f['x1'] / 1000:.1f}"])

    # scale bar and legend
    yb = top + row * len(tracks) + 20
    step = next(s for s in (1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7)
                if s >= maxspan / 6)
    out.append(f'<line x1="{X(0):.1f}" x2="{X(step):.1f}" y1="{yb}" y2="{yb}" stroke="#333" stroke-width="2"/>')
    out.append(f'<text x="{X(0):.1f}" y="{yb + 13}" font-size="9">{step / 1000:g} kb</text>')
    lx = X(step) + 40
    for g, c in list(col.items()) + [("boundary gene", BOUNDARY_COL)]:
        out.append(f'<rect x="{lx:.1f}" y="{yb - 6}" width="12" height="10" fill="{c}"/>')
        out.append(f'<text x="{lx + 16:.1f}" y="{yb + 3}" font-size="10">{html.escape(g)}</text>')
        lx += 30 + 6.5 * len(g)
    out.append("</svg>")
    svg = "\n".join(out)

    if a.svg:
        open(a.svg, "w").write(svg)
    else:
        sys.stdout.write(svg)
    if a.png:
        write_png(svg, a.png, a.width * 2)
    if a.layout:
        with open(a.layout, "w") as fh:
            fh.write("\t".join(["query", "contig", "orientation", "oriented_by", "kind", "name", "group",
                                "start", "end", "strand", "status", "in_tree", "x0_kb", "x1_kb"]) + "\n")
            for r in lay:
                fh.write("\t".join(map(str, r)) + "\n")
    print(f"[info] {len(tracks)} tracks from {len({t['query'] for t in tracks})} queries; "
          f"axis {maxspan / 1000:.0f} kb", file=sys.stderr)


if __name__ == "__main__":
    main()
