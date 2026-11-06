# Download files for analysis of exendins

## TOGA2 
***downloaded on October 23rd, 2026***
```
ssh -L 53682:localhost:53682 helix.nih.gov

module load rclone

rclone config


rclone lsd gdrive: --drive-shared-with-me 

# Top-level contents (should show integrated_annotations, query_*, target_*)
rclone lsd gdrive: --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug

# Full recursive listing with sizes
rclone lsf -R --format "ps" --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug gdrive: > toga2_listing.txt
wc -l toga2_listing.txt
head toga2_listing.txt

# pull all files, excluding the genome browser files and the Nextflow working files
module load rclone

cat > ~/toga2_filters.txt << 'EOF'
- **/ucsc_browser_files/**
- **/nextflow/**
- .command.*
- .exitcode
+ **
EOF

# Check the size first; it should be about 26-27 GB
rclone size gdrive: --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug \
  --filter-from ~/toga2_filters.txt

# Then copy, inside tmux so it survives a disconnect
tmux new -s toga2
rclone copy gdrive: /data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/ \
  --drive-root-folder-id 1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug \
  --filter-from ~/toga2_filters.txt \
  --progress --transfers 8 --checkers 16
```
Readme:
```
# TOGA2 annotations (Toxicofera / Anguimorpha)

- **Downloaded:** 2026-10-23, on helix.nih.gov
- **Source:** Google Drive folder shared by Elise Lauterbur, who ran the analysis
  https://drive.google.com/drive/folders/1PqZeWAav5iu-Baj7kkXPoPnYT30n29ug
- **Destination:** `/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/`
- **Size:** 25G
- **Excluded:** UCSC browser files (`ucsc_browser_files/`, ~11.5 GB) and Nextflow
  working files (`nextflow/`, `.command.*`, `.exitcode`): not needed for analysis

Script for downloading:
https://github.com/dannyjackson/Exendin_evolution/blob/main/1_download_files.md
```

## Reference genomes (for model rescue and annotation-evidence checks)
***downloaded on October 24th, 2026***

Needed by `2_exendin_origin.md` step 3: the HelCha genome to rebuild the exendin-3 model,
the HelSus genome and both annotations (EviAnn and NCBI RefSeq) to show why the reference
isoform used by TOGA2 is mis-annotated. Both NCBI assemblies are of the same material as
the TOGA2 runs, so coordinates transfer: `GCF_057233095.1` is the Wilson-lab HelSus
(isolate Gila_16-S) and `GCA_026122225.1` is the HelCha JANEZZ assembly TOGA2 used.

```
GENOMEDIR=/data/Wilson_Lab/data/Exendin_Evolution/genomes
mkdir -p ${GENOMEDIR} && cd ${GENOMEDIR}

# HelCha — JANEZZ contigs, same names as the TOGA2 query
datasets download genome accession GCA_026122225.1 --include genome,gff3,seq-report \
    --filename HelCha_GCA_026122225.1.zip
unzip -o -q HelCha_GCA_026122225.1.zip -d HelCha_GCA_026122225.1

# HelSus — RefSeq chromosomes (NC_139459.1 = chr2) plus NCBI's own annotation
datasets download genome accession GCF_057233095.1 --include genome,gff3,seq-report \
    --filename HelSus_GCF_057233095.1.zip
unzip -o -q HelSus_GCF_057233095.1.zip -d HelSus_GCF_057233095.1

for A in HelCha_GCA_026122225.1 HelSus_GCF_057233095.1; do
  samtools faidx $(ls ${A}/ncbi_dataset/data/*/*_genomic.fna)
done

# confirm NCBI and TOGA2 coordinates agree: this must be exendin-4 exon 1
HELSUS_FASTA=${GENOMEDIR}/HelSus_GCF_057233095.1/ncbi_dataset/data/GCF_057233095.1/GCF_057233095.1_16-S_ref_genomic.fna
samtools faidx ${HELSUS_FASTA} NC_139459.1:189613335-189613421   # ATGAAAATCATCCTGTGG...
```

## Published exendin cDNAs
Used as trusted queries for `2e_rescue_locus.py` and as the reference peptides in
`2_reference_peptides.fa`.

| Accession | What it is |
|---|---|
| AJ580309.1 | *H. horridum* (sensu lato, 1990) exendin-3 mRNA |
| U77613.1 | *H. suspectum* exendin-4 mRNA |
| EU790960.1 | *H. suspectum cinctum* exendin-2 mRNA |

```
cd ${PROJDIR}/datasets/query_fastas/nucleotide
for ACC in AJ580309.1 U77613.1 EU790960.1; do
  efetch -db nuccore -id ${ACC} -format fasta > ${ACC%%.*}.fa
done
```
If EDirect is unavailable, download them from NCBI Nucleotide by hand; `AJ580309.fa` is the
one `2e_rescue_locus.py` needs.
