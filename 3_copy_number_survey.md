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
