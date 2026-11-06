# Exendin gene tree from TOGA2 annotations only

**Scope.** This is a clean rebuild that replaces the exploratory analyses in
`2_exendin_origin.md` and `3_copy_number_survey.md`. It uses **only TOGA2**:

- **Gene set:** genes that the HelSus reference annotation (the one TOGA2 projected) calls an
  exendin. No other genes are added: no ANC anchor, no AQP5, and no proglucagon unless the
  annotation itself calls them exendin.
- **Tips:** TOGA2 query genes (`q_gene` in `orthology_classification.tsv`) that TOGA2 assigns
  as orthologs of those reference genes, in every query of `target_HelSus`. HelSus enters
  through its own self-projection (`query_HelSus`).
- **Sequences:** TOGA2's projected CDS (`nucleotide.fa.gz`). There is one per `q_gene`,
  chosen by best TOGA2 status, then the lowest mRNA number, then the longest CDS.

**Excluded:** BLAST/tBLASTn scans, rescued models (`2e`), the curation file
(`2_curation.tsv`), interval scans (`2f`, `3a`), NCBI or other outgroups, and nomenclature
assignment. A TOGA2 model that is wrong stays in the tree as TOGA2 built it. Its row in
`tips.tsv` lets you identify it.

**Tip names:** `<Query>_<label>`. The label is that of the reference exendin(s) the gene is
orthologous to, so `ShiCro_EX2_EX4` is a many2many gene. A letter suffix (`b`, `c` …)
separates several genes with the same orthology in one species. `tips.tsv` maps every tip
back to its TOGA2 `q_gene`, projection and coordinates.

## Set up environment
```
sinteractive --cpus-per-task=8 --mem=16g --time=4:00:00

source myconda
mamba activate exendin_env

PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution
HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3

cd ${SCRIPTDIR}
git pull

OUT=${PROJDIR}/analyses/toga2_exendin_tree
mkdir -p ${OUT} && cd ${OUT}
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
```

## 1. Which reference genes are annotated as exendin
```
python3 ${SCRIPTDIR}/4a_list_reference_exendins.py \
    --gff ${HELSUS_EVIANN_GFF} \
    --run ${TOGADIR}/target_HelSus \
    --out reference_exendins.tsv 2> reference_exendins.log

cat reference_exendins.log
cut -f1-6,9,10 reference_exendins.tsv
cut -f1,8 reference_exendins.tsv        # the annotation text that put each gene in the set
```
Check:
- **Every gene has `toga2_queries` > 0.** A warning means the GFF gene ID is not a reference
  gene in this TOGA2 run, i.e. the GFF is not the annotation TOGA2 used.
- **Which genes came back.** From earlier work, expect `LOC_00012324` (exendin-2) and
  `LOC_00012325` (exendin-4). Whether `LOC_00012323` (ANC) is in the set depends only on
  whether its annotation says "exendin". Step 4 depends on this.
- **Labels.** `label` is taken from the number in the annotation text (`EX2`, `EX4`; `EXN`
  if there is none). To rename tips, edit that column in `reference_exendins.tsv` before
  step 2.

## 2. One TOGA2 sequence per exendin gene per query
```
python3 ${SCRIPTDIR}/4b_collect_toga2_exendins.py \
    --run ${TOGADIR}/target_HelSus \
    --gene-table reference_exendins.tsv \
    --status FI,I,PI --min-cds 150 \
    --outdir toga2_exendins 2> toga2_exendins.log

cat toga2_exendins.log
cut -f1,3,4,6,7,9,15-19 toga2_exendins/tips.tsv | column -t
grep -c '>' toga2_exendins/exendin_cds.fa
```
The only filters are the TOGA2 status and the minimum CDS. Things to read in `tips.tsv`:
- **`orthology_class`.** one2one, one2many and many2many come straight from TOGA2.
- **`aa_len`.** The precursors are about 84–113 aa. A tip far outside that range is a TOGA2
  model problem, not biology. Under this protocol it stays in (or you drop it by changing
  `--status`/`--min-cds`), but note it. The known case is HelCha exendin-3: TOGA2 had only
  reference isoform `LOC_00012325-mRNA-2`, so expect a ~198 aa model there.
- **`internal_stops`.** Can be > 0 for PI projections. They become `X`/`NNN` in the alignment.
- **Species with no row.** TOGA2 assigned no ortholog there. `toga2_exendins/projections.tsv`
  lists the `one2zero` calls and any projections TOGA2 did not classify.

Two optional variants (separate `--outdir`, so they don't overwrite the main run):
```
# stricter: intact projections only (use this set for any selection test later)
python3 ${SCRIPTDIR}/4b_collect_toga2_exendins.py --run ${TOGADIR}/target_HelSus \
    --gene-table reference_exendins.tsv --status FI,I --outdir toga2_exendins_intact

# also keep TOGA2's paralogous projections (tips get a 'p' suffix). These are TOGA2 output,
# but TOGA2 itself does not report them as genes
python3 ${SCRIPTDIR}/4b_collect_toga2_exendins.py --run ${TOGADIR}/target_HelSus \
    --gene-table reference_exendins.tsv --include-paralogs --outdir toga2_exendins_paralogs
```

## 3. Codon alignment (MACSE)
```
mkdir -p aln && cd aln
macse -prog alignSequences -seq ../toga2_exendins/exendin_cds.fa \
      -out_NT exendin.macse_NT.fa -out_AA exendin.macse_AA.fa
macse -prog exportAlignment -align exendin.macse_NT.fa \
      -codonForFinalStop --- -codonForInternalStop NNN \
      -codonForExternalFS --- -codonForInternalFS NNN \
      -out_NT exendin.codon.fa -out_AA exendin.aa.fa
seqkit fx2tab -n -l exendin.codon.fa
cd ..
```
There is no column trimming, for the same reason as before: the precursor is short and the
peptide region is the part that matters.

## 4. Gene trees
```
mkdir -p tree && cd tree
$IQTREE -s ../aln/exendin.codon.fa -st CODON -m MFP -B 1000 -alrt 1000 -T AUTO --prefix exendin.codon
$IQTREE -s ../aln/exendin.aa.fa              -m MFP -B 1000 -alrt 1000 -T AUTO --prefix exendin.aa
grep "Best-fit model according to" exendin.codon.iqtree exendin.aa.iqtree
```
**Rooting.** The gene set has no outgroup unless the annotation calls ANC an exendin (step 1).
- **ANC is not in the set (expected):** midpoint-root the tree. This is a display choice,
  not a test of where the root is.
  ```
  for T in exendin.codon exendin.aa; do
    python3 ${SCRIPTDIR}/4c_midpoint_root.py ${T}.treefile > ${T}.rooted.nwk
  done
  ```
- **ANC is in the set:** its orthologs are the natural root. Re-run IQ-TREE with `-o` set to
  them, then use the treefile directly.
  ```
  ANCLAB=$(awk -F'\t' '$1=="LOC_00012323"{print $2}' ../reference_exendins.tsv)
  OG=$(grep '>' ../aln/exendin.codon.fa | sed 's/>//' | grep -E "_${ANCLAB}[b-z]?$" | paste -sd,)
  $IQTREE -s ../aln/exendin.codon.fa -st CODON -m MFP -B 1000 -alrt 1000 -o ${OG} -T AUTO \
          --prefix exendin.codon.ANCroot --redo
  cp exendin.codon.ANCroot.treefile exendin.codon.rooted.nwk
  ```

## 5. Draw
Tips are coloured by the reference exendin they are TOGA2 orthologs of.
```
for T in exendin.codon exendin.aa; do
  python3 ${SCRIPTDIR}/2c_draw_tree.py --tree ${T}.rooted.nwk \
      --groups ../toga2_exendins/groups.tsv \
      --title "TOGA2 exendin orthologs (HelSus reference), ${T#exendin.}" \
      --png ${T}.png > ${T}.svg
done
```
Node labels are SH-aLRT/UFBoot; trust a node at ≥80 / ≥95. Bold = *Heloderma*.

## Reading the result
- **Within each reference exendin's orthologs**, the tree should broadly follow the species
  tree if TOGA2's orthology calls are right. A species placed far from where the species tree
  puts it points to a TOGA2 call or model to look at in `tips.tsv`.
- **many2many / one2many tips** are where TOGA2 could not split the copies. Where they sit
  relative to the one2one tips is the main thing this tree adds over the orthology table.
- **Compare with the exploratory tree** (`2_exendin_origin.md`) only after this one is read on
  its own terms. Differences come from what this protocol deliberately leaves out: ANC, the
  zebra finch root, rescued models, and paralogous projections.
