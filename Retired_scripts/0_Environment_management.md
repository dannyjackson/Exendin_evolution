# Environment management

## Make working directory
```
PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026

mkdir -p ${PROJDIR}/scripts
cd ${PROJDIR}/scripts
git clone https://github.com/dannyjackson/Exendin_evolution

SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution

cd /data/Wilson_Lab/projects/GLP1_evolution/analyses
mkdir toga2 cactus

```
## Software

### Cactus / HAL tools (halLiftover, halStats, hal2fasta, cactus-hal2maf)
Use the Biowulf module (or the Cactus container), not conda.
```bash
module avail cactus
module load cactus/3.3.0
```

### mamba environment: TOGA2 output, NCBI genomes, sequence/evolution analyses
```bash
source myconda
mamba create -n exendin_env -c conda-forge -c bioconda --strict-channel-priority \
  python=3.11 pandas biopython \
  ncbi-datasets-cli \
  samtools bedtools seqkit gffread \
  ucsc-twobittofa ucsc-fatotwobit ucsc-bigbedtobed ucsc-liftover \
  blast miniprot exonerate \
  mafft macse trimal pal2nal \
  iqtree paml hyphy

mamba activate exendin_env

# record exact versions for the methods section / reproducibility
mamba env export --no-builds > exendin_env.yml
```

| Purpose | Tools |
|---|---|
| Download NCBI genomes/annotations | `datasets` (ncbi-datasets-cli) |
| Coordinates, intervals, FASTA handling | samtools (faidx), bedtools, seqkit, gffread |
| TOGA2 files (`exon_seqs.2bit`, UCSC formats, chains) | twoBitToFa, faToTwoBit, bigBedToBed, liftOver |
| Find exendin/GCG homologs in raw genomes | blast (tblastn), miniprot, exonerate |
| Codon alignments | mafft, macse, pal2nal, trimal |
| Gene trees and selection (dN/dS) | iqtree, paml (codeml), hyphy |
```

## Set up environment
Run this at the top of every markdown file
```
soure myconda
mamba activate hal_toga_env

PROJDIR=/data/Wilson_Lab/projects/GLP1_evolution/
TOGADIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_TOGA2_annotations/
CACTUSDIR=/data/Wilson_Lab/projects/squamates/Arizona_squamates/misc/Toxicofera_cactus_alignment_2026
SCRIPTDIR=${PROJDIR}/scripts/Exendin_evolution

cd ${SCRIPTDIR}
git pull
cd ${PROJDIR}
```