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
(`HelSus_A076_EX2`). These are provisional handles, not nomenclature — step 4 decides
what each gene should be called.

**A caution carried through the whole pipeline.** TOGA2 projects the *reference* transcript's
exon structure onto each query, so a mis-annotated reference isoform produces a wrong model in
every query that uses it. That happened here to the HelCha exendin-3 gene (step 3), and the
same failure mode could affect any locus. Step 2's QC table is where it shows up: a CDS length
far from its orthologs' is the symptom to chase.

## Set up environment
```
sinteractive --cpus-per-task=8 --mem=16g --time=4:00:00

source myconda
mamba activate exendin_env

PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution
GENOMEDIR=/data/Wilson_Lab/data/Exendin_Evolution/genomes

# genomes used for the rescue and annotation-evidence checks (see 1_download_files.md)
HELSUS_FASTA=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/GCF_057233095.1_16-S_ref_genomic.fna
HELCHA_FASTA=${GENOMEDIR}/HelCha_GCA_026122225.1/ncbi_dataset/data/GCA_026122225.1/GCA_026122225.1_Hcha_1_genomic.fna
HELSUS_NCBI_GFF=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/genomic.gff
HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3

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
labels; copies left out; rescued models from step 3). Each line carries a note saying why.
```
python3 ${SCRIPTDIR}/2b_curate_exendin_loci.py \
    --collect-dir collect \
    --curation ${SCRIPTDIR}/2_curation.tsv \
    --outdir curated

cut -f1,3,8,9,11,12,13,14,15 curated/loci_curated.tsv
```
Filters: TOGA2 status intact (FI or I), CDS ≥150 bp, no internal stop. Check that every
species has an `_ANC` row and that no curation warnings were printed.

Then look at CDS lengths. Anything far from its orthologs is a broken model, not biology:
```
awk -F'\t' 'NR==1||$13=="yes"' curated/loci_curated.tsv | cut -f1,12 | sort -k2,2n
```
These precursors are 84–113 aa (252–339 bp). A locus well outside that range goes to step 3.

## 3. Rescue loci whose TOGA2 model is wrong
`2e_rescue_locus.py` realigns a trusted query (a published cDNA, or a protein) to the
genome, chains the BLAST hits into exons, splices them and reports the longest ORF.

**The case this was written for.** `HelCha_A129_EX3` came out at 594 bp / 198 aa against
87 aa for its HelSus ortholog. Three lines of evidence showed the model was wrong:

- the two proteins only become homologous at `LSSEDSASSESFASKIKR`; HelCha's first 42
  residues are Ala/Pro/Arg-rich with no hydrophobic core, so not a signal peptide;
- the published cDNA AJ580309 aligns to the HelCha model only from that point on;
- TOGA2 had only `LOC_00012325-mRNA-2` for HelCha, whose first exon is 419 bp against
  86 bp in `mRNA-1` (see the annotation-evidence check at the end of this step).

```
# the published H. horridum exendin-3 cDNA

cat << EOF > AJ580309.fa
>AJ580309.1 Heloderma horridum mRNA for exendin 3
GGTGCCTGTGCAAAGGAGAAGATGAAAATCATCCTGTGGCTGTGTGTTTTCGGGCTGTTCCTTGCAACTT
TATTCCCTGTCAGCTGGCAAATGCCTGTTGAATCTGGGTTGTCTTCTGAGGATTCTGCAAGCTCAGAAAG
CTTTGCTTCGAAGATTAAGCGACATAGTGATGGAACATTTACCAGTGACTTGTCAAAACAGATGGAAGAG
GAGGCAGTGCGGTTATTTATTGAGTGGCTTAAGAACGGAGGACCAAGTAGCGGGGCACCTCCGCCATCGG
GTTAAGGTCTTTCAATTGTGGAACAAGACACACACCTGATGTTTGACGACCATTTTAAAGAAATGTTTCC
AGCAATAAGTCACATGTCTTTGTTTCCAATTAGTGAGCGACACAGCC
EOF

python3 ${SCRIPTDIR}/2e_rescue_locus.py \
    --genome ${HELCHA_FASTA} --contig JANEZZ010000386.1 \
    --query AJ580309.fa --name HelCha_A129_EX3 --outdir rescued
```
Check the printed exon table: hits collinear, one strand, query covered end to end. For this
locus it should report 3 exons on the minus strand, exon 1 at ~4,399,218–4,399,326 — **5.9 kb
upstream of where TOGA2 placed it** — and an 87 aa protein identical to the published
exendin-3 precursor (P20394).

Then record it in `2_curation.tsv` and re-run step 2:
```
# HelCha  JANEZZ010000386__1  4391184  replace  rescued/HelCha_A129_EX3  <note>
python3 ${SCRIPTDIR}/2b_curate_exendin_loci.py \
    --collect-dir collect --curation ${SCRIPTDIR}/2_curation.tsv --outdir curated
```
The `replace` action substitutes `<value>.aa.fa` and `<value>.cds.fa` for the projection, so
the rescue is reproducible from the repo rather than by hand-editing FASTA files.

### Why the reference isoform is wrong (annotation evidence)
Worth recording, and worth passing to whoever maintains the annotation:
```
grep "LOC_00012325" ${HELSUS_EVIANN_GFF} | awk -F'\t' '$3=="mRNA"'

awk -F'\t' '$1=="NC_139459.1" && $4<189625000 && $5>189610000 && ($3=="mRNA"||$3=="exon")' \
    ${HELSUS_NCBI_GFF} | cut -f3,4,5,7,9 | sed 's/;Dbxref.*product=/ | /; s/;.*//'
```
What this showed: `mRNA-1` carries `EvidenceTranscriptID=MSTRG_00103018:1:4423.74` and
`Class==` (exact intron-chain match) with the full exendin-3 protein as evidence; `mRNA-2`
carries `MSTRG_00103020:1:106.45` and `Class=k`, with only a PDB chain of the **mature
peptide** as protein evidence — so its 111-residue N-terminal extension has no protein
support, at ~2% of the locus's expression. NCBI's independent annotation of the same assembly
(`XM_084836889.1`, `XR_014788128.1`) uses the canonical first exon at 189,613,255–189,613,421
and never the alternative one. Conclusion: the transcript is probably a real minor variant,
but its CDS call is an artifact, and `mRNA-1` is the correct model.

## 4. Assign each gene to the published nomenclature
Compares the **mature peptide** of every curated protein with the published peptides in
`2_reference_peptides.fa` (exendin-1 to -4, GCGL, and the glucagon, GLP-1, GLP-2 and PACAP
domains as controls), by global alignment. `n_hormone_domains` counts the domains in the whole
precursor: **1** for an exendin or GCGL, **3** for a real proglucagon.
```
python3 ${SCRIPTDIR}/2d_assign_nomenclature.py \
    --query curated/curated_aa.fa \
    --reference ${SCRIPTDIR}/2_reference_peptides.fa \
    --out-prefix nomenclature

# select columns by name, so this does not drift if a column is added
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
  {print $h["gene"]"\t"$h["mature_peptide"]"\t"$h["n_hormone_domains"]"\t"$h["how_found"]"\t"\
         $h["best_match"]"\t"$h["identity"]"\t"$h["verdict"]}' nomenclature_assignments.tsv

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

The panel's GCGL entry is the zebra finch peptide confirmed in step 5, so squamate ANC
orthologs should come back as **GCGL-type**. Note that the *Heloderma* ANC peptides do not
trim at a dibasic site: they carry `ERYS..` where every other squamate has `KR`, so the
K→E substitution has removed the processing signal immediately after the hormone domain.
Confirm against the full precursors before making anything of it:
```
cut -f1,16 curated/loci_curated.tsv | grep ANC     # look at what follows WLINT in each
```

## 5. Outgroup: zebra finch GCGL
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
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i}
  {print $h["gene"]"\t"$h["peptide_len"]"\t"$h["n_hormone_domains"]"\t"$h["verdict"]}' \
  zefi_check_assignments.tsv
```
Expect `n_hormone_domains` = 1 and "encodes published GCGL exactly". Three domains would mean
the accession is the real proglucagon and the wrong gene was picked.

**This outgroup only orients the tree; it does not test where the duplicates came from.**
With one outgroup sequence the root is placed on that branch by construction. To test whether
the duplicates descend from GCGL, add GCGL orthologs from several vertebrates (chicken
**EU718628**; also *Xenopus*, coelacanth and a teleost, per Wang et al. 2012, each confirmed to
sit between AQP5 and SLC4A8) plus a single-domain non-GCGL root such as GIP. Then:
*(bird/frog/fish GCGL, (squamate GCGL, duplicates))* = duplication from GCGL within reptiles;
*((all GCGL), duplicates)* = older duplication with losses elsewhere; duplicates attaching
outside GCGL = they did not come from GCGL.

## 6. Align
```
cat curated/curated_aa.fa outgroup.faa > exendin_origin_aa.fa
seqkit fx2tab -n -l exendin_origin_aa.fa | sort -k2,2n | tail -3     # no outliers
mafft --localpair --maxiterate 1000 exendin_origin_aa.fa > exendin_origin.aln.fa
```

## 7. Tree, rooted on zebra finch
```
$IQTREE -s exendin_origin.aln.fa -m MFP -B 1000 -alrt 1000 -o TaeGut_GCGL -T AUTO \
        --prefix exendin_origin

grep "Best-fit model according to" exendin_origin.iqtree
```
Sequences IQ-TREE reports as identical are kept and are expected here: several anguid ANC
peptides are identical, as are the two *Shinisaurus* copies.

## 8. Draw
```
python3 ${SCRIPTDIR}/2c_draw_tree.py --tree exendin_origin.treefile --outgroup TaeGut_GCGL \
    --title "Exendin-family proteins, SLC4A8–AQP5 interval, rooted on zebra finch GCGL" \
    --png exendin_origin.png > exendin_origin.svg
```
Node labels are SH-aLRT/UFBoot; trust a node at ≥80 / ≥95. Bold = *Heloderma*.
## 8b
```
ANC=$(grep '>' curated/curated_aa.fa | sed 's/>//' | grep _ANC | paste -sd,)
echo "((${ANC}));" > ANC_mono.constraint.nwk
MODEL=$(grep "Best-fit model according to" exendin_origin.iqtree | awk '{print $NF}')

$IQTREE -s exendin_origin.aln.fa -m ${MODEL} -g ANC_mono.constraint.nwk \
        -o TaeGut_GCGL -T AUTO --prefix ANC_mono

cat exendin_origin.treefile ANC_mono.treefile > cand.nwk
$IQTREE -s exendin_origin.aln.fa -m ${MODEL} -z cand.nwk -n 0 -zb 10000 -au -T AUTO --prefix AUtest
grep -A6 "USER TREES" AUtest.iqtree
```
## Reading the result
Put `nomenclature_assignments.tsv` beside the tree; together they answer the question.

1. **ANC orthologs** should form a clade following the species tree, with zebra finch
   outside it. That identifies the ancestral gene and dates it to before the
   bird–squamate split.
2. **Where the exendin copies attach.** Sister to squamate ANC = the array arose by
   duplication of ANC within squamates. Outside all ANC copies = the duplication is older.
   Remember the caveat in step 5: with a single outgroup this is oriented, not tested.
3. **Which published name goes with which gene.** Each "exactly" row names a gene; genes
   that pair in the tree and carry two different published names are one gene named twice.
4. **Copies with no published name** are genes that exist but have never been purified
   from venom.

## Renaming
Once steps 4 and 7 agree, the genes can be renamed on evidence rather than on which
peptide was purified first. A scheme that would survive new genomes:
- keep **GCGL** (or *exendin*, per Irwin) for the ANC lineage, since that name is already
  in the literature for this gene in birds and fish;
- name the venom copies by lineage and position rather than by discovery order, and
  record the historical peptide names as synonyms of the genes step 4 identifies.

Do not rename anything before step 4 is run on the full curated set: the published names
are peptide names, and only the mature-peptide comparison can say which gene each belongs to.



# Copy number across squamates: how many exendin-family genes, annotated or not

**Why this is needed.** The array is systematically under-annotated. In HelSus four of
six copies had no gene model; in HemCap, which the gene tables list as having *one*
exendin, a direct sequence scan found three more loci. So presence/absence cannot be
taken from a GFF — it has to come from the sequence.

**Why it matters.** The gene tree's backbone is unresolved (116 alignable columns, ~79
informative sites), so it cannot date the duplication. Copy number mapped onto the species
tree can. The HemCap result already suggests the array predates the Scinciformata /
Anguimorpha split.

**Method.** tblastn full precursors — one per lineage — against the AQP5–SLC4A8 interval,
then classify every hit by *where in the query* it falls. A hit in the signal/pro region is
exon 1; a hit after the hormone-domain motif is exon 2. A locus with both is a candidate
gene; a locus with one is a remnant. Premature stops and a hormone domain not starting
His/Tyr mark disabled copies.

## Set up
```
source myconda && mamba activate exendin_env
PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution
GENOMEDIR=/data/Wilson_Lab/data/Exendin_Evolution/genomes
OUT=${PROJDIR}/analyses/copy_number
mkdir -p ${OUT} ${GENOMEDIR} && cd ${OUT}
```

## 1. Query set
One full precursor per lineage identified in stage 2. The scan classifies hits by position
within these, so they must be complete precursors, not mature peptides.
```
seqkit grep -p HemCap_ANC,HelSus_A076_EX2,HelSus_A126_EX4,HelSus_A038,HelSus_S379 \
    ${PROJDIR}/analyses/exendin_origin/exendin_origin_aa.fa > lineage_queries.faa
grep -c '>' lineage_queries.faa      # 5
```

## 2. Species table
Coordinates are the known flanks and annotated candidates; the scan searches the whole
interval regardless, so these are supporting evidence rather than targets.
```
cat > species.tsv << 'EOF'
#label	accession	contig	known_min	known_max	note
AnoCar	GCF_035594765.1	NC_085842.1	155087274	155138943	SLC4A8 not found; exendin + AQP5 only
CanAsp	GCF_035149785.1	NC_086154.1	135220948	135362723	AQP5 .. exendin .. SLC4A8
PogVit	GCF_051106095.1	NC_135784.1	163211717	163460527	AQP5 .. unknown LOC .. exendin .. SLC4A8
EubMac	GCF_028583425.1	NC_072808.1	3391923	3535608	AQP5 .. 3 annotated exendins .. SLC4A8
EOF
```

## 3. Download, scan
`PAD` extends past the known flanks, because HelSus has copies 118–379 kb from ANC —
well outside the AQP5–SLC4A8 interval as drawn from the annotation.
```
PAD=150000

while IFS=$'\t' read -r LABEL ACC CTG LO HI NOTE; do
  [[ "$LABEL" == \#* ]] && continue
  echo "=== $LABEL ($ACC) $NOTE"

  if [[ ! -d ${GENOMEDIR}/${LABEL}_${ACC} ]]; then
    datasets download genome accession ${ACC} --include genome \
        --filename ${GENOMEDIR}/${LABEL}.zip
    unzip -o -q ${GENOMEDIR}/${LABEL}.zip -d ${GENOMEDIR}/${LABEL}_${ACC}
  fi
  FA=$(ls ${GENOMEDIR}/${LABEL}_${ACC}/ncbi_dataset/data/*/*_genomic.fna)
  [[ -f ${FA}.fai ]] || samtools faidx ${FA}

  # clamp the window to the contig
  CLEN=$(awk -v C=$CTG '$1==C{print $2}' ${FA}.fai)
  S=$(( LO - PAD )); [[ $S -lt 1 ]] && S=1
  E=$(( HI + PAD )); [[ $E -gt $CLEN ]] && E=$CLEN
  echo "    window ${CTG}:${S}-${E}  (contig ${CLEN} bp)"

  python3 ${SCRIPTDIR}/2f_scan_interval.py \
      --genome ${FA} --contig ${CTG} --region ${S}-${E} \
      --query lineage_queries.faa --out-prefix scan_${LABEL}
done < species.tsv
```

## 4. Same scan on the species already analysed
So the table is comparable end to end. HelSus is the positive control: it must recover all
six known loci.
```
HELSUS=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/GCF_057233095.1_16-S_ref_genomic.fna
python3 ${SCRIPTDIR}/2f_scan_interval.py --genome ${HELSUS} --contig NC_139459.1 \
    --region 188950000-189850000 --query lineage_queries.faa --out-prefix scan_HelSus

HEMCAP=$(ls ${GENOMEDIR}/HemCap_*/ncbi_dataset/data/*/*_genomic.fna)
python3 ${SCRIPTDIR}/2f_scan_interval.py --genome ${HEMCAP} --contig NC_069658.1 \
    --region 429200000-429600000 --query lineage_queries.faa --out-prefix scan_HemCap
```
Check `scan_HelSus_loci.tsv` recovers ANC (189,489,482), the Frag1 copy (189,527,255),
exendin-2 (189,564,066), exendin-4 (189,613,335) and the two distal copies (189,371,236
and 189,110,176). If the scan misses known loci, lower `--min-identity` before trusting
any absence elsewhere.

## 5. Summarise
```
{ echo -e "species\tloci\tcandidate_genes\tremnants\tdisabled"
  for F in scan_*_loci.tsv; do
    S=${F#scan_}; S=${S%_loci.tsv}
    awk -F'\t' -v S=$S 'NR>1{
      n++
      if ($7 ~ /^candidate_gene/) c++
      else r++
      if ($7 ~ /premature_stop|domain_starts_/) d++
    } END{printf "%s\t%d\t%d\t%d\t%d\n", S, n, c, r, d}' $F
  done
} | column -t

# every locus, side by side
for F in scan_*_loci.tsv; do S=${F#scan_}; S=${S%_loci.tsv}
  awk -F'\t' -v S=$S 'NR>1{print S"\t"$0}' $F; done \
  | sort -k1,1 -k3,3n | column -t
```

## Reading it
- **`candidate_gene`** — both an exon-1 and a hormone-domain hit at one locus. A real gene
  unless flagged disabled.
- **`N_term_only` / `domain_only`** — a remnant, or a gene whose other exon is too diverged
  to hit. Check the `subject_translation` column in `_hits.tsv` before calling it either way.
- **`premature_stop`** — a stop codon inside the aligned region: pseudogene.
- **`domain_starts_X`** — the hormone domain does not begin His or Tyr, so the peptide
  cannot function even if the ORF is open.

Then map counts onto the species tree. The question is whether non-anguimorph squamates
carry remnants (array ancient, repeatedly lost) or genuinely nothing (array young, arising
within Anguimorpha). *Eublepharis* with three annotated copies and HemCap with three
unannotated remnants both point to the former.

Loci found here are **not** added to the stage 2 tree: remnants and pseudogenes sit on long
branches and destabilise an already weakly resolved backbone. They are presence/absence
evidence. Any locus that comes out a clean, intact candidate gene can be rebuilt with
`2e_rescue_locus.py` and considered for the tree on its own merits.
