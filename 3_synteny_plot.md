# Synteny plot of the step-3 scan loci
This draws the same figure as step 6 of `4_toga2_exendin_tree.md`, but from the tblastn
interval scans of `3_copy_number_survey.md` instead of TOGA2. Compare the two figures side by
side: the TOGA2 one shows what TOGA2 annotates, and this one shows what the sequence scan finds.
Append it to `3_copy_number_survey.md` as step 6, or keep it as its own file.

`3b_plot_scan_synteny.py` imports its drawing code from `4d_plot_exendin_synteny.py`, so keep
the two in the same directory.

## Set up
```
source myconda && mamba activate exendin_env
PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution
GENOMEDIR=/data/Wilson_Lab/data/Exendin_Evolution/genomes
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026
OUT=${PROJDIR}/analyses/copy_number
cd ${OUT}
```

## 1. Boundary-gene annotations
`3a_scan_species.sh` downloads genomes without annotation, so fetch each species' NCBI GFF
once, into the genome directory 3a already made. HelSus already has its GFF
(`1_download_files.md`).
```
for ROW in AnoCar:GCF_035594765.1 CanAsp:GCF_035149785.1 PogVit:GCF_051106095.1 EubMac:GCF_028583425.1; do
  L=${ROW%%:*}; ACC=${ROW#*:}
  datasets download genome accession ${ACC} --include gff3 --filename ${GENOMEDIR}/${L}_gff.zip \
    && unzip -o -q ${GENOMEDIR}/${L}_gff.zip -d ${GENOMEDIR}/${L}_${ACC}
done
ls ${GENOMEDIR}/*/ncbi_dataset/data/*/genomic.gff
```
HemCap is in the script's built-in table (contig `NC_069658.1`). Its GFF is found under
`${GENOMEDIR}/HemCap_*`, wherever you put it in step 4.

## 2. Plot
```
python3 ${SCRIPTDIR}/3b_plot_scan_synteny.py \
    --scan-dir ${OUT} --genomedir ${GENOMEDIR} \
    --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
    --svg scan_synteny.svg --png scan_synteny.png --layout scan_synteny.tsv 2> scan_synteny.log
cat scan_synteny.log
```
Check the log:
- **Each species reports its loci, plus the gene it used for SLC4A8 and for AQP5.** The ID
  and description are printed. If no gene is called SLC4A8 exactly, the script uses one whose
  description matches, e.g. a "solute carrier family 4 member 8-like" LOC. If that is the
  wrong gene, pass your own pattern, e.g.
  `--boundaries "SLC4A8=SLC4A8|solute carrier family 4 member 8,AQP5=AQP5|aquaporin[- ]?5"`.
  A species with no GFF, or no boundary gene on the scanned contig, is drawn without it.
  This is expected for AnoCar, where SLC4A8 was not found.
- **Some `scan_*_loci.tsv` files have no row in the table** (their contig is unknown). Add
  them with your own `--species-table`, in the same format as 3a
  (`label accession contig known_min known_max note`).
- **Tree matching:** tree tips are matched as `Genus_species` → `GenSpe`. Species not in the
  Anguimorpha tree (e.g. AnoCar, EubMac or HemCap, if absent) are drawn below it. Place them
  with `--tree-map`, or use a broader tree.

## Reading it
- **Colour** is the lineage query that hit best (`best_query` without its species prefix:
  ANC, EX2, EX4, A038, S379).
- **Filled** = `candidate_gene` (both an exon-1 and a hormone-domain hit, no flag).
- **Hollow** = a candidate gene flagged `premature_stop` or `domain_starts_X`.
- **Dashed** = `N_term_only` or `domain_only`: a remnant, or a gene whose other exon is too
  diverged to hit.
- Hover over a locus in the SVG to see its state, parts, best identity and domain start.

The scan only covers the window 3a searched (`known_min`–`known_max` ± `PAD`). An empty
stretch beyond that window is unsearched, not empty.
