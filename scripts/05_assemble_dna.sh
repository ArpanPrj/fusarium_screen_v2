#!/usr/bin/env bash
# Exploratory DNA assembly using all reads, including the host.
# MEGAHIT is useful here for non-target contigs, not a polished fungal genome.
set -Eeuo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
WD="/scratch/aubaxp004/fusarium_screen_v2"
THREADS=8
MEM_GB=48
RUN=SRR27384773
OUT="$WD/assembly/dna_megahit"
mkdir -p "$WD/assembly" "$WD/logs"
if [[ -d "$OUT" ]]; then
    echo "Output directory exists. Inspect it; do not overwrite a prior assembly."
    echo "To resume an interrupted MEGAHIT run: megahit --continue -o $OUT"
    exit 1
fi
megahit -1 "$WD/raw/${RUN}_1.fastq.gz" -2 "$WD/raw/${RUN}_2.fastq.gz" \
    -t "$THREADS" -m "$((MEM_GB * 1024 * 1024 * 1024))" \
    --min-contig-len 500 -o "$OUT" > "$WD/logs/$RUN.megahit.log" 2>&1
samtools faidx "$OUT/final.contigs.fa"
echo "DNA contigs: $OUT/final.contigs.fa"
