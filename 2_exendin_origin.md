# What are the "exendins"? Curated proteins, published names, and a rooted gene tree

**Question:** the name *exendin* was applied to proteins prior to understanding their genomic 
and phylogenetic context. Some proteins were collected from Heloderma horridium samples that
predated the split of this species into four distinct lineages: **Mexican Beaded Lizard** 
(*Heloderma horridum*), **Chiapan Beaded Lizard** (*Heloderma alvarezi*), **Guatemalan Beaded Lizard**,
(*Heloderma charlesbogerti*), and **Río Fuerte Beaded Lizard** (*Heloderma exasperatum*). 
These proteins were all purified and named based on distinct protein coding sequences, but
could easily be orthologs from distinct populations and species, or could be paralogous genes
within the same species, or (most likely) could be a mix of both options.

Additionally, reference genomes have applied diverse nomenclatures to a gene that appears to
be orthologous to the ancestral exendin gene in preliminary analyses. It is sometimes named
exendin-2 like, exendin-3 like, pro-glucagon (zebra finch), and glucagon-like.

Before embarking on an analysis of the exendin gene family, I will first work out a naming system
to keep track of genes across this complex and messy nomenclature system

| Use of "exendin" | Where it comes from | What it actually denotes |
|---|---|---|
| exendin-1, -2, -3, -4 | peptides purified from *Heloderma* venom, 1984–1992 | **mature peptide sequences**, named before any genome existed. Nothing in the names refers to a locus |
| "exendin-like", "exendin-3-like" | NCBI annotations | **sequence similarity only**, not an orthology claim |
| the ancestral exendin gene | Irwin 2012; = **GCGL** of Wang et al. 2012 | a **single-copy gene** found in birds, reptiles, amphibians and fish, with its own receptor (GCGLR), apparently absent from mammals |

**Working hypothesis to test here:** the third sense is gene (`LOC_00012323`, tentatively 
referred to as "ANC" until we properly identify it), has 1:1 orthologs in all 14 genomes, 
all in the conserved **SLC4A8 – ANC – AQP5** interval, and is the same gene as the zebra 
finch locus annotated "pro-glucagon". The *Heloderma* venom exendins are tandem duplicates 
beside it. (Supporting circumstantial evidence from the first-pass GFF notes: *Elgaria* 
carries a `GCGR-like` receptor between CLUAP1 and ANKS3, which is where GCGLR would sit 
— a receptor for this ligand lineage, not for glucagon.)

**Approach:** collect every copy in the conserved **SLC4A8 – ANC – AQP5** interval 
cross the 14 squamate genomes annotated with TOGA2, curate one protein per copy, 
**assign each one to the published nomenclature by mature-peptide sequence**, then 
build a protein tree rooted on the zebra finch ortholog.

**Copy names** give the position relative to ANC: `HelSus_A076` = 76 kb from ANC on the
AQP5 side; `S` = SLC4A8 side; `U` = not on ANC's contig. Curation labels are appended
(`HelSus_A076_EX2`). These are provisional handles, not nomenclature — step 3 decides
what each gene should be called.

## Set up environment
```
sinteractive --cpus-per-task=8 --mem=16g --time=4:00:00

source myconda
mamba activate exendin_env

PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution

cd ${SCRIPTDIR}
git pull

OUT=${PROJDIR}/analyses/exendin_origin
mkdir -p ${OUT} && cd ${OUT}
IQTREE=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree)
```

## 1. Collect TOGA2 projections
Projections of ANC, the two annotated exendins and AQP5 in every query, merged into loci.
Because TOGA2 projects each annotated exendin through every paralogous chain, this finds
copies that are in no annotation.
```
python3 ${SCRIPTDIR}/2a_collect_exendin_projections.py \
    --run ${TOGADIR}/target_HelSus \
    --genes LOC_00012324,LOC_00012325,LOC_00012326 \
    --anchor LOC_00012323 \
    --outdir collect 2> collect.log
cat collect.log
```

## 2. Curate one protein per locus
Judgment calls live in `2_curation.tsv` (HelCha's ANC, which TOGA2 missed; provisional
labels; copies left out). Each line carries a note saying why.
```
python3 ${SCRIPTDIR}/2b_curate_exendin_loci.py \
    --collect-dir collect \
    --curation ${SCRIPTDIR}/2_curation.tsv \
    --outdir curated

cut -f1,3,8,9,11,12,13,14 curated/loci_curated.tsv
```
Filters: TOGA2 status intact (FI or I), CDS ≥150 bp, no internal stop. Check that every
species has an `_ANC` row and that no curation warnings were printed.

## 3. Assign each gene to the published nomenclature
Compares the **mature peptide** of every curated protein with the published peptides in
`2_reference_peptides.fa` (exendin-1 to -4, and the glucagon, GLP-1, GLP-2 and PACAP
domains as controls), by global alignment.
```
python3 ${SCRIPTDIR}/2d_assign_nomenclature.py \
    --query curated/curated_aa.fa \
    --reference ${SCRIPTDIR}/2_reference_peptides.fa \
    --out-prefix nomenclature

cut -f1,2,4,5,6,7,10 nomenclature_assignments.tsv     # verdict per gene
cut -f1-5 nomenclature_identity.tsv                    # identity to each published peptide
```
Read it like this:
- **"encodes published EXENDIN<i> exactly"** — this locus is the gene behind that
  published peptide. Where two species each have such a gene and they are each other's
  closest relatives, the two published names are one gene in two species.
- **"EXENDIN<i>-type, N difference(s)"** — same peptide class, not the published isolate:
  another species, another individual, or another copy.
- **"no published match"** — a copy that has never been described as a peptide. Expect
  several; the exendin-1/-2 class and the exendin-3/-4 class are only ~30–40% identical
  to each other, so anything in that range is genuinely something else.
- **`how_found`** says how the mature peptide was located. Check any row that is not
  `motif` by eye; those precursors have a diverged N-terminus.

Note that the panel is *Heloderma* venom peptides, so ANC orthologs from other species
will come back as "no published match". That is expected: ANC's published name is **GCGL**,
not an exendin number, and no GCGL peptide sequence is in the panel. Add one to
`2_reference_peptides.fa` if you obtain it.

## 4. Outgroup: zebra finch GCGL
The zebra finch gene between AQP5 and SLC4A8, annotated by NCBI as "pro-glucagon".
Confirmed as GCGL: a single hormone domain of exactly 29 aa before its dibasic cleavage
site, 66% identical to glucagon — matching the "cGCGL1–29 ... 62–66% identity with GCG"
of Wang et al. 2012. Re-run the check below if you change the accession.
```
ZEFI=XP_004175947.5      # its RefSeq protein accession

datasets download gene accession ${ZEFI} --include protein --filename zefi.zip
unzip -o -q zefi.zip -d zefi
# keep ONE record: NCBI often returns several identical isoforms, and duplicate
# tips break the tree step
seqkit seq -w 0 zefi/ncbi_dataset/data/protein.faa | seqkit head -n 1 \
  | sed '1s/^>.*/>TaeGut_GCGL/' > outgroup.faa
grep -c '>' outgroup.faa     # must be 1
cat outgroup.faa

# how many hormone domains? one = ANC/GCGL-type; three = real proglucagon
python3 ${SCRIPTDIR}/2d_assign_nomenclature.py --query outgroup.faa \
    --reference ${SCRIPTDIR}/2_reference_peptides.fa --out-prefix zefi_check
cut -f1,2,4,5,6,10 zefi_check_assignments.tsv
```
If it has more than one isoform, keep the longest.

## 5. Align
```
cat curated/curated_aa.fa outgroup.faa > exendin_origin_aa.fa
mafft --localpair --maxiterate 1000 exendin_origin_aa.fa > exendin_origin.aln.fa
```

## 6. Tree, rooted on zebra finch
```
$IQTREE -s exendin_origin.aln.fa -m MFP -B 1000 -alrt 1000 -o TaeGut_ANC -T AUTO \
        --prefix exendin_origin
grep "Best-fit model according to" exendin_origin.iqtree
```

## 7. Draw
```
python3 ${SCRIPTDIR}/2c_draw_tree_svg.py --tree exendin_origin.treefile --outgroup TaeGut_ANC \
    --title "Exendin-family proteins, SLC4A8–AQP5 interval, rooted on zebra finch ANC/GCGL" \
    > exendin_origin.svg
```
Node labels are SH-aLRT/UFBoot; trust a node at ≥80 / ≥95. Bold = *Heloderma*.

## Reading the result
Put `nomenclature_assignments.tsv` beside the tree; together they answer the question.

1. **ANC orthologs** should form a clade following the species tree, with zebra finch
   outside it. That identifies the ancestral gene and dates it to before the
   bird–squamate split.
2. **Where the exendin copies attach.** Sister to squamate ANC = the array arose by
   duplication of ANC within squamates. Outside all ANC copies = the duplication is older.
3. **Which published name goes with which gene.** Each "exactly" row names a gene; genes
   that pair in the tree and carry two different published names are one gene named twice.
4. **Copies with no published name** are genes that exist but have never been purified
   from venom.

## Renaming
Once steps 3 and 6 agree, the genes can be renamed on evidence rather than on which
peptide was purified first. A scheme that would survive new genomes:
- keep **GCGL** (or *exendin*, per Irwin) for the ANC lineage, since that name is already
  in the literature for this gene in birds and fish;
- name the venom copies by lineage and position rather than by discovery order, and
  record the historical peptide names as synonyms of the genes step 3 identifies.

Do not rename anything before step 3 is run on the full curated set: the published names
are peptide names, and only the mature-peptide comparison can say which gene each belongs to.
