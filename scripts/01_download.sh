#!/usr/bin/env bash
# Independent downloads into the new project; no symbolic links.
set -Eeuo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
WD="/scratch/aubaxp004/fusarium_screen_v2"
THREADS=4
mkdir -p "$WD"/{raw/sra,tmp,logs,metadata}
trap 'echo "ERROR at line $LINENO; inspect $WD/logs" >&2' ERR
for RUN in SRR27384773 SRR1185280; do
    R1="$WD/raw/${RUN}_1.fastq.gz"
    R2="$WD/raw/${RUN}_2.fastq.gz"
    [[ ! -L "$R1" && ! -L "$R2" ]] || { echo "Refusing symbolic links: $RUN"; exit 1; }
    if [[ -s "$R1" && -s "$R2" ]]; then
        gzip -t "$R1" "$R2"
        echo "Reusing complete new-project FASTQs: $RUN"
        continue
    fi
    [[ ! -e "$R1" && ! -e "$R2" ]] || { echo "Partial pair exists: inspect $WD/raw"; exit 1; }
    prefetch "$RUN" --output-directory "$WD/raw/sra" > "$WD/logs/$RUN.prefetch.log" 2>&1
    JOBTMP=$(mktemp -d "$WD/tmp/$RUN.XXXXXX")
    mkdir -p "$JOBTMP/fastq" "$JOBTMP/work"
    fasterq-dump "$WD/raw/sra/$RUN/$RUN.sra" --split-files --threads "$THREADS" --outdir "$JOBTMP/fastq" --temp "$JOBTMP/work" > "$WD/logs/$RUN.fasterq.log" 2>&1
    for MATE in 1 2; do
        test -s "$JOBTMP/fastq/${RUN}_${MATE}.fastq"
        gzip "$JOBTMP/fastq/${RUN}_${MATE}.fastq"
        gzip -t "$JOBTMP/fastq/${RUN}_${MATE}.fastq.gz"
    done
    mv "$JOBTMP/fastq/${RUN}_1.fastq.gz" "$R1"
    mv "$JOBTMP/fastq/${RUN}_2.fastq.gz" "$R2"
done
echo "Both complete readsets are ready."
