#!/usr/bin/env python3
"""
Find and extract glucagon/secretin-superfamily hormone domains.

The family's mature peptides all begin with a recognisable N-terminus
(HSQGTFT glucagon, HAEGTFT GLP-1, HADGSFS GLP-2, YAEGTFI GIP, HGEGTFT exendin-4,
HSDAIFT helodermin, HSDAVFT VIP, HSDGIFT PACAP, HSDGTFT secretin, YADAIFT GHRH).
A precursor can carry one domain (exendins, GIP, GHRH) or several (proglucagon:
glucagon + GLP-1 + GLP-2), so the family tree is built from domains, not precursors.

A domain is called where the regex below matches and the next --width residues
reach --min-id percent identity (ungapped, from residue 1) to at least one
peptide in the reference panel. The panel is only used to find and tentatively
label domains; the tree decides the relationships.

Modes
  scan     list domains in any protein FASTA (plain or .gz), e.g. the HelSus
           self-projection, to discover which reference genes are family members
  extract  write each domain as protein (and codon, if --cds is given) FASTA,
           named <seq_id>_<k><best_ref>, plus a table

Usage
  python3 2d_glucagon_family_domains.py scan --aa protein.fa.gz > domains_scan.tsv
  python3 2d_glucagon_family_domains.py extract --aa family_aa.fa --cds family_cds.fa \
          --out-prefix family_domains
"""
import argparse, gzip, re, sys

MOTIF = re.compile(r"(?=([HY][AGS][DEQ][AG][ATVIGS]F[TSI]))")

# Human mature peptides except where noted (Heloderma for exendin-4 / helodermin).
PANEL = {
    "GCG":   "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
    "GLP1":  "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGRG",
    "GLP2":  "HADGSFSDEMNTILDNLAARDFINWLIQTKITD",
    "GIP":   "YAEGTFISDYSIAMDKIHQQDFVNWLLAQKGKKNDWKHNITQ",
    "EX4":   "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
    "EX2":   "HSDAIFTEEYSKLLAKLALQKYLASILGSRTSPPP",
    "VIP":   "HSDAVFTDNYTRLRKQMAVKKYLNSILN",
    "PHI":   "HADGVFTSDFSKLLGQLSAKKYLESLM",
    "PACAP": "HSDGIFTDSYSRYRKQMAVKKYLAAVLGKRYKQRVKNK",
    "SCT":   "HSDGTFTSELSRLREGARLQRLLQGLV",
    "GHRH":  "YADAIFTNSYRKVLGQLSARKLLQDIMSRQQGESNQERGARARL",
}


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_fasta(p):
    name, buf, out = None, [], []
    with opener(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].strip().split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def ident(a, b):
    n = min(len(a), len(b))
    return 100.0 * sum(1 for x, y in zip(a[:n], b[:n]) if x == y) / n if n else 0.0


def find_domains(seq, width, min_id):
    hits, last = [], -10 ** 9
    for m in MOTIF.finditer(seq):
        s = m.start()
        if s < last + width:
            continue
        win = seq[s:s + width]
        if len(win) < width or "*" in win:
            continue
        scores = sorted(((ident(win, p), k) for k, p in PANEL.items()), reverse=True)
        if scores[0][0] >= min_id:
            hits.append((s, win, scores[0], scores[1]))
            last = s
    return hits


def fallback_domain(seq, width, min_id):
    """No motif match: slide a window along the whole sequence and keep the best
    panel match, if it reaches min_id. For diverged copies whose N-terminus has
    drifted from the canonical motif. Only used in extract mode (cheap on short
    precursors, too slow for a whole proteome)."""
    best = None
    for s in range(0, len(seq) - width + 1):
        win = seq[s:s + width]
        if "*" in win:
            continue
        sc = sorted(((ident(win, p), k) for k, p in PANEL.items()), reverse=True)
        if best is None or sc[0][0] > best[2][0]:
            best = (s, win, sc[0], sc[1])
    return [best] if best and best[2][0] >= min_id else []


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    for m in ("scan", "extract"):
        p = sub.add_parser(m)
        p.add_argument("--aa", required=True)
        p.add_argument("--width", type=int, default=28,
                       help="domain length in aa from the first residue (default 28: fits inside glucagon's 29)")
        p.add_argument("--min-id", type=float, default=40.0)
        if m == "extract":
            p.add_argument("--cds", help="CDS FASTA with the same IDs (no stop codon); adds codon output")
            p.add_argument("--out-prefix", required=True)
            p.add_argument("--fallback-min-id", type=float, default=45.0,
                           help="identity needed for the sliding-window fallback (0 disables)")
    a = ap.parse_args()

    seqs = read_fasta(a.aa)
    if a.mode == "scan":
        print("seq_id\tdomain\taa_start\tmotif\tdomain_seq\tbest_ref\tbest_id\tsecond_ref\tsecond_id\tn_domains_in_seq")
        for name, seq in seqs:
            hits = find_domains(seq.upper(), a.width, a.min_id)
            for k, (s, win, b1, b2) in enumerate(hits, 1):
                print(f"{name}\t{k}\t{s + 1}\t{win[:7]}\t{win}\t{b1[1]}\t{b1[0]:.0f}\t{b2[1]}\t{b2[0]:.0f}\t{len(hits)}")
        return

    cds = dict(read_fasta(a.cds)) if a.cds else {}
    with open(a.out_prefix + "_aa.fa", "w") as fa, open(a.out_prefix + ".tsv", "w") as tab:
        fc = open(a.out_prefix + "_codon.fa", "w") if a.cds else None
        tab.write("domain_id\tsource\tdomain\taa_start\tbest_ref\tbest_id\tsecond_ref\tsecond_id\tdomain_seq\tmethod\n")
        n = 0
        for name, seq in seqs:
            hits, method = find_domains(seq.upper(), a.width, a.min_id), "motif"
            if not hits and a.fallback_min_id > 0:
                hits, method = fallback_domain(seq.upper(), a.width, a.fallback_min_id), "fallback"
            if not hits:
                print(f"[warn] no domain found in {name}", file=sys.stderr)
            for k, (s, win, b1, b2) in enumerate(hits, 1):
                did = f"{name}_{k}{b1[1]}"
                fa.write(f">{did}\n{win}\n")
                tab.write(f"{did}\t{name}\t{k}\t{s + 1}\t{b1[1]}\t{b1[0]:.0f}\t{b2[1]}\t{b2[0]:.0f}\t{win}\t{method}\n")
                if fc is not None:
                    nt = cds.get(name)
                    if nt is None:
                        print(f"[warn] no CDS for {name}", file=sys.stderr)
                    else:
                        fc.write(f">{did}\n{nt[3 * s:3 * (s + a.width)]}\n")
                n += 1
        if fc:
            fc.close()
    print(f"[info] {n} domains from {len(seqs)} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
