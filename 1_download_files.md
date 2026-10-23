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
https://github.com/dannyjackson/Exendin_evolution/1_download_files.md
```
