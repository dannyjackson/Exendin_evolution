# Synteny plot of the step-2 curated loci
This shows where every locus from `2_exendin_origin.md` sits on its chromosome, between
SLC4A8 and AQP5, in species-tree order. The labels are the 2b names (`A076_EX2`, `S118`, `U1`,
`ANC`), so the plot reads directly against the step-2 tree. Append it to
`2_exendin_origin.md` after step 8, or keep it as its own file.

`2g_plot_synteny.py` imports its drawing code from `4d_plot_exendin_synteny.py`, so keep the
two in the same directory.

```
PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution
HELSUS_EVIANN_GFF=${PROJDIR}/datasets/genomes/Heloderma_suspectum/Heloderma_suspectum.final.reference.gff3
OUT=${PROJDIR}/analyses/exendin_origin
cd ${OUT}

python3 ${SCRIPTDIR}/2g_plot_synteny.py \
    --loci curated/loci_curated.tsv \
    --run ${TOGADIR}/target_HelSus \
    --boundaries SLC4A8,AQP5=LOC_00012326 --gff ${HELSUS_EVIANN_GFF} \
    --species-tree ${CACTUSDIR}/Anguimorpha_species_tree_iteration-01_5-23-26.nwk \
    --svg exendin_origin_synteny.svg --png exendin_origin_synteny.png \
    --layout exendin_origin_synteny.tsv 2> exendin_origin_synteny.log
cat exendin_origin_synteny.log
```
- **Boundary genes:** these are TOGA2's ortholog projections of SLC4A8 and AQP5 from the same
  `target_HelSus` run 2a used, so contig names and coordinates match `loci_curated.tsv`.
  Check the SLC4A8 ID the log prints, then pass it as `SLC4A8=LOC_...`.
- **How loci are drawn:**
  - *Filled:* the locus entered the step-2 tree.
  - *Hollow:* left out by 2b's filters (TOGA2 status, CDS length, internal stops).
  - *Dashed:* dropped in `2_curation.tsv`.
  - The reason and curation note are in the SVG hover text.
- **Colour** is the 2b name: ANC, a curation label (EX1–EX4), or unnamed copy. To colour by
  the 2d verdict instead, build a groups file from the nomenclature table:
  ```
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{print $h["gene"]"\t"$h["best_match"]}' \
      nomenclature_assignments.tsv > nomenclature_groups.tsv
  # then add: --groups nomenclature_groups.tsv
  ```
- **Coordinates of rescued loci:** these come from the TOGA2 locus in 2a, not from the
  `2e` model. HelCha EX3's rebuilt exon 1 lies about 5.9 kb from where TOGA2 placed it. At the
  plot's scale this doesn't show.
- **Species with no step-2 locus** are still drawn, with their boundary genes only, so an
  absence is visible rather than silent.
