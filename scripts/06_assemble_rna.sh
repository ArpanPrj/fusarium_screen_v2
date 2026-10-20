#!/usr/bin/env bash
# rnaSPAdes: the RNA assembler suggested in the assignment.
# Uses all paired RNA reads. Do not apply DNA assembly coverage assumptions.
set -Eeuo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
WD="/scratch/aubaxp004/fusarium_screen_v2"
THREADS=8
MEM_GB=48
RUN=SRR1185280
OUT="$WD/assembly/rna_rnaspades"
mkdir -p "$WD/assembly" "$WD/logs"
if [[ -d "$OUT" ]]; then
    echo "Output exists. Inspect it before restarting."
    echo "To resume: rnaspades.py --continue -o $OUT"
    exit 1
fi
# No strand-specific flag because library orientation has not been verified.
rnaspades.py -1 "$WD/raw/${RUN}_1.fastq.gz" -2 "$WD/raw/${RUN}_2.fastq.gz" \
    -t "$THREADS" -m "$MEM_GB" -o "$OUT" > "$WD/logs/$RUN.rnaspades.log" 2>&1
samtools faidx "$OUT/transcripts.fasta"
echo "RNA assembly: $OUT/transcripts.fasta"
