# Exendin gene tree, codon alignment and selection (TOGA2, HelSus reference)

Builds on `2a_collect_exendin_projections.py`. Steps:
1. collect projections → 2. one CDS per locus → 3. check the loci by hand →
4. codon alignment → 5. gene tree and duplication-timing test →
6. recombination / gene conversion → 7. selection → 8. ancestral sequences

Reference genes (collaborator's curated HelSus annotation):

| ID | What it is |
|---|---|
| `LOC_00012324` | exendin-2 (helodermin-type), HelSus chr2:189,564,066–189,569,010 (+) |
| `LOC_00012325` | exendin-4, HelSus chr2:189,613,335–189,621,417 (+) |
| `LOC_00012323` | glucagon-family neighbour, chr2:189,489,482–189,492,340 (−); intact 1:1 in 13/14 queries. Used as the **synteny anchor** and **outgroup** |
| `LOC_00001483` | proglucagon, chr1 ~180.5 Mb |
| `LOC_00000122` | glucagon-family, chr1 ~13.3 Mb, 1:1 everywhere; identity still to confirm |

## Set up environment
```
source myconda
mamba activate exendin_env

PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution

cd ${SCRIPTDIR}
git pull
cd ${PROJDIR}

sinteractive --cpus-per-task=8 --mem=16g --time=8:00:00

OUT=${PROJDIR}/analyses/toga2
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
```

## 1. Collect projections (one row per projection, one per query locus, FASTA)
```
python3 ${SCRIPTDIR}/2a_collect_exendin_projections.py \
    --run ${TOGADIR}/target_HelSus \
    --genes LOC_00012324,LOC_00012325,LOC_00001483,LOC_00000122 \
    --anchor LOC_00012323 \
    --outdir ${OUT}/collect_HelSus 2> ${OUT}/collect_HelSus.log

# Any "[warn] ... cannot tell which is the query" means TOGA2 stores reference and
# query records under the same ID; look at the headers and tell me the format.
cat ${OUT}/collect_HelSus.log
zcat ${TOGADIR}/target_HelSus/query_ShiCro/nucleotide.fa.gz | grep '>' | head -4
```

## 2. One representative CDS per locus
Each locus is named by its position relative to the anchor, oriented like HelSus
(anchor on the minus strand, exendin cluster at higher coordinates):
`ShiCro_R028o` = 28 kb from the anchor on the cluster side (R), opposite strand (o).
`L` = the other side; `s` = same strand; `U<n>` = not on the anchor's contig; `ANC` = the anchor.
```
python3 ${SCRIPTDIR}/2b_select_exendin_loci.py \
    --collect-dir ${OUT}/collect_HelSus \
    --genes LOC_00012324,LOC_00012325 --anchor LOC_00012323 \
    --status FI,I --min-cds 150 \
    --outdir ${OUT}/loci_HelSus

cut -f1,7,8,11,13,15-17 ${OUT}/loci_HelSus/loci_selected.tsv
```
Use `--status FI,I` for selection analyses. A second run with `--status FI,I,PI` can be
used for the gene tree only (more copies, but partial sequences).

## 3. Check the loci by hand before aligning
```
cd ${OUT}/loci_HelSus

# All-vs-all protein identity: catches haplotig duplicates (~100% identity between two
# loci of one species on different ptg contigs), mislabelled anchors, and paralogs
blastp -query exendin_plus_anchor_aa.fa -subject exendin_plus_anchor_aa.fa \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" -evalue 1e-3 \
  | awk '$1!=$2' | sort -k1,1 -k6,6gr > allvsall_blastp.tsv

# Mature-peptide region of each copy (first ~45 aa after the signal peptide)
awk '/^>/{n=$0; next} {print n"\t"substr($0,1,90)}' exendin_plus_anchor_aa.fa
```
Things to decide here (edit the FASTA by hand, or rerun 2b with other filters):
- **HelCha anchor.** TOGA2 found no anchor projection in HelCha, but HelCha's exendin
  copies are on one contig with the same spacing as in HelSus. The `LOC_00012324_b` locus
  (JANEZZ010000386__1:4,521,249, + strand) sits exactly where the anchor should be
  (~35 kb past the Frag1-equivalent copy, on the anchor's strand). If it matches
  `HelSus_ANC` better than any exendin in `allvsall_blastp.tsv`, it is the anchor:
  rename it `HelCha_ANC` and move it to `anchor_cds.fa`.
- **Off-cluster HelSus copies** `HelSus_L118s`, `HelSus_L363s`, `HelSus_L379s` (same strand
  as the anchor, 118–379 kb on the far side). They aren't in the curated annotation. Keep
  them only if the translation looks like a real precursor (signal peptide, an intact
  peptide region, a C-terminal Gly before a basic residue for amidation).
- **ptg-contig assemblies** (BarImb, ElgCoe, GerLio) can carry haplotigs. Two copies of one
  species with ~100% identity on different contigs are probably one gene; keep one.
- **Pogona** (`#retro` projection, 479 bp single block): probably a processed copy; exclude it from selection.

## 4. Codon alignment (MACSE)
MACSE aligns at the codon level and tolerates frameshifts in projected sequences.
Two alignments: with the anchor (to root the tree) and exendin-only (for selection).
```
cd ${OUT}
mkdir -p aln && cd aln

for SET in exendin_plus_anchor exendin; do
  macse -prog alignSequences -seq ../loci_HelSus/${SET}_cds.fa \
        -out_NT ${SET}.macse_NT.fa -out_AA ${SET}.macse_AA.fa
  # frameshifts (!) and internal stops become NNN; the final stop is removed
  macse -prog exportAlignment -align ${SET}.macse_NT.fa \
        -codonForFinalStop --- -codonForInternalStop NNN \
        -codonForExternalFS --- -codonForInternalFS NNN \
        -out_NT ${SET}.codon.fa -out_AA ${SET}.aa.fa
done

# Look at it (signal peptide, peptide region, and spacer should line up across copies)
seqkit fx2tab -n -l exendin.codon.fa
```
No column trimming at this stage: the precursor is ~80–110 codons and trimming would
remove the peptide region, which is the part of interest. If one sequence badly
disrupts the alignment, remove that sequence instead of trimming columns.

## 5. Gene tree, and when did the duplications happen?
```
cd ${OUT} && mkdir -p tree && cd tree

ANC=$(grep '>' ../loci_HelSus/anchor_cds.fa | sed 's/>//' | paste -sd,)

# codon model, rooted on the anchor orthologs
$IQTREE -s ../aln/exendin_plus_anchor.codon.fa -st CODON -m MFP \
        -B 1000 -alrt 1000 -o ${ANC} -T AUTO --prefix exendin_anchor.codon

# protein tree as a cross-check
$IQTREE -s ../aln/exendin_plus_anchor.aa.fa -m MFP \
        -B 1000 -alrt 1000 -o ${ANC} -T AUTO --prefix exendin_anchor.aa
```
**Duplication-timing test.** H1: the *Heloderma* copies arose by lineage-specific
duplication, so all *Heloderma* exendins form one clade. The alternative: copies predate
*Heloderma*, so each *Heloderma* copy groups with the copy at the same position in
*Shinisaurus* / anguids (compare the R/L + kb tip names). Fit the constrained tree and
compare it with the unconstrained one (AU test):
```
# Heloderma exendin copies (R/L/U tips) vs everything else, including all anchors
HEL=$(grep '>' ../aln/exendin_plus_anchor.codon.fa | sed 's/>//' | grep -E  '^Hel(Sus|Cha)_[RLU]' | paste -sd,)
OTH=$(grep '>' ../aln/exendin_plus_anchor.codon.fa | sed 's/>//' | grep -vE '^Hel(Sus|Cha)_[RLU]' | paste -sd,)
echo "((${HEL}),${OTH});" > H1_heloderma_monophyletic.constraint.nwk

MODEL=$(grep "Best-fit model according to" exendin_anchor.codon.iqtree | awk '{print $NF}')
echo ${MODEL}
$IQTREE -s ../aln/exendin_plus_anchor.codon.fa -st CODON -m ${MODEL} \
        -g H1_heloderma_monophyletic.constraint.nwk -o ${ANC} -T AUTO --prefix H1

cat exendin_anchor.codon.treefile H1.treefile > candidate_trees.nwk
$IQTREE -s ../aln/exendin_plus_anchor.codon.fa -st CODON -m ${MODEL} \
        -z candidate_trees.nwk -n 0 -zb 10000 -au -T AUTO --prefix AUtest
grep -A 6 "USER TREES" AUtest.iqtree
```
The alignment is short, so expect weak support. Give synteny (tip position names) equal
weight with the tree when deciding positional orthology.

Optional: reconcile with the species tree (`Anguimorpha_species_tree_iteration-01_5-23-26.nwk`)
to place duplications and losses on branches (GeneRax or NOTUNG; neither is in `exendin_env`).

## 6. Recombination / gene conversion (GARD)
Tandem copies often exchange sequence by gene conversion, which distorts both the tree
and dN/dS. Check before running selection tests:
```
cd ${OUT} && mkdir -p selection && cd selection
hyphy gard --alignment ../aln/exendin.codon.fa --type codon --rate-classes 3 \
     --output gard.json > gard.log
grep -i "breakpoint" gard.log | tail
```
If GARD finds a supported breakpoint, run step 7 on each segment with its own tree (or
use the GARD partitioned output with HyPhy), and note it as possible gene conversion.

## 7. Selection (HyPhy)
Tree for selection = the anchor-rooted codon tree, with the anchor pruned and
*Heloderma* branches labelled:
```
python3 ${SCRIPTDIR}/2c_label_tree.py \
    --tree ../tree/exendin_anchor.codon.treefile \
    --keep ../aln/exendin.codon.fa \
    --fg '^Hel(Sus|Cha)_' --reference > exendin.HelFG.nwk

# gene-wide episodic diversifying selection: anywhere / on Heloderma branches
hyphy busted --alignment ../aln/exendin.codon.fa --tree exendin.HelFG.nwk \
      --output busted_all.json > busted_all.log
hyphy busted --alignment ../aln/exendin.codon.fa --tree exendin.HelFG.nwk \
      --branches Foreground --output busted_HelFG.json > busted_HelFG.log

# which branches (exploratory, every branch tested)
hyphy absrel --alignment ../aln/exendin.codon.fa --tree exendin.HelFG.nwk \
      --output absrel.json > absrel.log

# intensified vs relaxed selection on Heloderma (venom) vs other branches
hyphy relax --alignment ../aln/exendin.codon.fa --tree exendin.HelFG.nwk \
      --test Foreground --reference Reference --output relax_HelFG.json > relax_HelFG.log

# which codons (pervasive: FEL; episodic: MEME)
hyphy fel  --alignment ../aln/exendin.codon.fa --tree exendin.HelFG.nwk --output fel.json  > fel.log
hyphy meme --alignment ../aln/exendin.codon.fa --tree exendin.HelFG.nwk --output meme.json > meme.log
```
Other foregrounds are worth testing once the tree is settled, e.g. only the branches
leading to the exendin-4 copies (relabel with a different `--fg` regex), or only
*Heloderma* terminal branches (`--tips-only`). For a codeml branch-site cross-check,
write the tree with `--codeml` (PAML `#1` marks, no branch lengths).

## 8. Ancestral sequences
```
cd ${OUT}/tree
$IQTREE -s ../aln/exendin_plus_anchor.codon.fa -st CODON -m ${MODEL} \
        -te exendin_anchor.codon.treefile -asr -T AUTO --prefix ASR
# ASR.state: per-node, per-codon posterior probabilities. Next script maps the
# changes in the mature-peptide region onto branches.
```

## Where the Cactus alignment still matters
- Varanids (no projections) and *Lanthanotus* (10-bp remnant): use `halLiftover` on the
  HelSus cluster interval (chr2:189,450,000–189,700,000) to tell a deletion from a gap.
- *Pogona* / CanAsp: check the anchor's neighbourhood for diverged copies TOGA2 can't project.
- HelCha: confirm that the anchor sits at JANEZZ010000386__1:~4.52 Mb.
