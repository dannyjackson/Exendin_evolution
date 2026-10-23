#!/usr/bin/env python3
"""
Prepare a gene tree for HyPhy / codeml:
  1. prune tips not in the alignment (e.g. the anchor outgroup used only for rooting),
     collapsing the unary nodes this leaves and summing their branch lengths
  2. strip internal support values (IQ-TREE writes "SH-aLRT/UFBoot" labels, which
     HyPhy would otherwise read as node names)
  3. label branches: a tip whose name matches --fg is {Foreground}; an internal
     branch is {Foreground} only if every tip below it matches (so the branches
     leading into and within an all-Heloderma clade are foreground, but the stem
     of a mixed clade is not). With --tips-only, only terminal branches are labelled.
     Optionally label every other branch {Reference} (for HyPhy RELAX).

Pruning keeps the rooting implied by the input tree, so root the tree first
(IQ-TREE -o <anchor tips>) and prune the anchor afterwards.

Usage
  python3 2c_label_tree.py --tree exendin.treefile --keep exendin_cds.aln.fa \
      --fg '^Hel(Sus|Cha)_' [--reference] [--tips-only] [--codeml] > exendin.fg.nwk

--codeml writes PAML-style marks (#1) instead of {Foreground} and omits branch
lengths (codeml re-estimates them).
"""
import argparse, re, sys


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self):
        self.name, self.length, self.children, self.parent = "", None, [], None


def parse(s):
    s = s.strip().rstrip(";")
    i, root = 0, Node()
    stack, cur = [], root
    while i < len(s):
        c = s[i]
        if c == "(":
            n = Node(); n.parent = cur; cur.children.append(n); stack.append(cur); cur = n; i += 1
        elif c == ",":
            p = stack[-1]; n = Node(); n.parent = p; p.children.append(n); cur = n; i += 1
        elif c == ")":
            cur = stack.pop(); i += 1
        elif c == ":":
            m = re.match(r":([^,():;\[]+)", s[i:]); cur.length = float(m.group(1)); i += len(m.group(0))
        elif c == "[":
            i = s.index("]", i) + 1           # drop comments
        else:
            m = re.match(r"('[^']*'|[^,():;\[]+)", s[i:]); cur.name = m.group(1).strip("'"); i += len(m.group(0))
    # first '(' created a child of the dummy root; unwrap
    return root.children[0] if len(root.children) == 1 and not root.name else root


def tips(n):
    return [n] if not n.children else [t for c in n.children for t in tips(c)]


def prune(root, keep):
    def rec(n):
        if not n.children:
            return n if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            k = kids[0]
            if n.length is not None or k.length is not None:
                k.length = (n.length or 0.0) + (k.length or 0.0)
            return k
        n.children = kids
        for k in kids:
            k.parent = n
        return n
    r = rec(root)
    if r is None:
        sys.exit("no tips left after pruning")
    r.length = None
    return r


def write(n, label, is_root=True, lengths=True):
    if n.children:
        s = "(" + ",".join(write(c, label, False, lengths) for c in n.children) + ")"
        name = ""                              # drop support values
    else:
        s, name = "", n.name
    s += name + (label(n) if not is_root else "")
    if lengths and n.length is not None and not is_root:
        s += f":{n.length:.10g}"
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--keep", help="FASTA (aligned or not); prune tips not in it")
    ap.add_argument("--fg", required=True, help="regex for foreground tip names")
    ap.add_argument("--tips-only", action="store_true")
    ap.add_argument("--reference", action="store_true", help="label all other branches {Reference}")
    ap.add_argument("--codeml", action="store_true", help="PAML marks (#1) instead of HyPhy {Foreground}")
    ap.add_argument("--keep-rooted", action="store_true", help="do not unroot a bifurcating root")
    a = ap.parse_args()

    root = parse(open(a.tree).read())
    if a.keep:
        keep = {l[1:].split()[0] for l in open(a.keep) if l.startswith(">")}
        missing = keep - {t.name for t in tips(root)}
        if missing:
            print(f"[warn] in alignment but not in tree: {sorted(missing)}", file=sys.stderr)
        root = prune(root, keep)
    fg = re.compile(a.fg)
    fgmark = "#1" if a.codeml else "{Foreground}"
    refmark = "" if a.codeml else "{Reference}"

    # HyPhy and codeml work on unrooted trees and merge the two root branches;
    # if those carry different labels the merged branch is ambiguous. Unroot here
    # by dissolving the root child that is NOT an all-foreground clade, so the
    # branch separating an all-foreground clade from the rest stays foreground.
    if not a.keep_rooted and len(root.children) == 2:
        def allfg(n):
            return all(fg.search(t.name) for t in tips(n))
        internal = [c for c in root.children if c.children]
        if internal:
            dis = next((c for c in internal if not allfg(c)), internal[0])
            other = root.children[0] if root.children[1] is dis else root.children[1]
            other.length = (other.length or 0.0) + (dis.length or 0.0)
            root.children = [other] + dis.children
            for c in root.children:
                c.parent = root

    def label(n):
        ts = tips(n)
        is_fg = all(fg.search(t.name) for t in ts) and (not a.tips_only or not n.children)
        if is_fg:
            return fgmark
        return refmark if a.reference else ""

    nfg = sum(1 for t in tips(root) if fg.search(t.name))
    print(f"[info] tips: {len(tips(root))}  foreground tips: {nfg}", file=sys.stderr)
    print(write(root, label, lengths=not a.codeml) + ";")


if __name__ == "__main__":
    main()
