#!/usr/bin/env bash
# Scan the AQP5-SLC4A8 interval of each species for exendin-family copies.
#
#   bash 3a_scan_species.sh [species_table]
#
# Expects in the environment (or edit the defaults below):
#   SCRIPTDIR  repo with 2f_scan_interval.py
#   GENOMEDIR  where genomes are downloaded
# and in the working directory:
#   lineage_queries.faa   one full precursor per lineage
#
# Species table: whitespace-separated, '#' comments, columns
#   label  accession  contig  known_min  known_max  [note...]
# Defaults to the built-in table if no file is given.

set -u

SCRIPTDIR=${SCRIPTDIR:-/data/Wilson_Lab/projects/GLP1_evolution/scripts/Exendin_evolution}
GENOMEDIR=${GENOMEDIR:-/data/Wilson_Lab/data/Exendin_Evolution/genomes}
QUERY=${QUERY:-lineage_queries.faa}
PAD=${PAD:-150000}
MINID=${MINID:-35}

SCAN=${SCRIPTDIR}/2f_scan_interval.py
[[ -f $SCAN  ]] || { echo "!! not found: $SCAN"; exit 1; }
[[ -f $QUERY ]] || { echo "!! not found: $QUERY (see step 1)"; exit 1; }
mkdir -p "$GENOMEDIR"

# label  accession  contig  known_min  known_max  note
read -r -d '' TABLE << 'EOF' || true
AnoCar  GCF_035594765.1  NC_085842.1  155087274  155138943  SLC4A8 not found; exendin + AQP5
CanAsp  GCF_035149785.1  NC_086154.1  135220948  135362723  AQP5 .. exendin .. SLC4A8
PogVit  GCF_051106095.1  NC_135784.1  163211717  163460527  AQP5 .. unknown LOC .. exendin .. SLC4A8
EubMac  GCF_028583425.1  NC_072808.1  3391923    3535608    AQP5 .. 3 annotated exendins .. SLC4A8
EOF
[[ $# -ge 1 ]] && TABLE=$(cat "$1")

printf '%s\n' "$TABLE" | while read -r LABEL ACC CTG LO HI REST; do
  [[ -z "${LABEL:-}" || "$LABEL" == \#* ]] && continue
  echo "=== $LABEL ($ACC) ${REST:-}"

  DEST=${GENOMEDIR}/${LABEL}_${ACC}
  if [[ ! -d $DEST ]]; then
    echo "    downloading $ACC ..."
    datasets download genome accession "$ACC" --include genome \
        --filename "${GENOMEDIR}/${LABEL}.zip" || { echo "    !! download failed"; continue; }
    unzip -o -q "${GENOMEDIR}/${LABEL}.zip" -d "$DEST" || { echo "    !! unzip failed"; continue; }
  fi

  FA=$(ls "$DEST"/ncbi_dataset/data/*/*_genomic.fna 2>/dev/null | head -1)
  if [[ -z "$FA" || ! -f "$FA" ]]; then echo "    !! no FASTA under $DEST"; continue; fi
  [[ -f ${FA}.fai ]] || samtools faidx "$FA"

  CLEN=$(awk -v C="$CTG" '$1==C{print $2}' "${FA}.fai")
  if [[ -z "$CLEN" ]]; then
    echo "    !! contig $CTG not in ${FA}.fai; first few names:"
    cut -f1 "${FA}.fai" | head -3 | sed 's/^/       /'
    continue
  fi

  S=$(( LO - PAD )); (( S < 1 )) && S=1
  E=$(( HI + PAD )); (( E > CLEN )) && E=$CLEN
  echo "    window ${CTG}:${S}-${E}  (contig ${CLEN} bp)"

  python3 "$SCAN" --genome "$FA" --contig "$CTG" --region "${S}-${E}" \
      --query "$QUERY" --min-identity "$MINID" --out-prefix "scan_${LABEL}" \
      || echo "    !! scan failed for $LABEL"
done

echo
echo "=== summary"
{ echo -e "species\tloci\tcandidate_genes\tremnants\tflagged\tintact"
  for F in scan_*_loci.tsv; do
    [[ -f $F ]] || continue
    S=${F#scan_}; S=${S%_loci.tsv}
    awk -F'\t' -v S="$S" 'NR>1{
        n++
        if ($7 ~ /^candidate_gene/) c++; else r++
        if ($7 ~ /premature_stop|domain_starts_/) d++
        if ($7 ~ /^candidate_gene$/) i++
      } END{printf "%s\t%d\t%d\t%d\t%d\t%d\n", S, n+0, c+0, r+0, d+0, i+0}' "$F"
  done
} | { command -v column >/dev/null && column -t || cat; }
