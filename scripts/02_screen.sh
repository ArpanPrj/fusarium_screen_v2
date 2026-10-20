#!/usr/bin/env bash
# ONLY SendSketch and Kraken2; NO screening database downloads.
set -Eeuo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
WD="/scratch/aubaxp004/fusarium_screen_v2"
KDB="/datasets/unzipped/kraken2_full"
THREADS=8
mkdir -p "$WD"/{screen/sendsketch,screen/kraken2,logs,metadata}
trap 'echo "ERROR at line $LINENO; inspect $WD/logs" >&2' ERR
for TOOL in sendsketch.sh kraken2; do command -v "$TOOL" >/dev/null; done
for FILE in hash.k2d opts.k2d taxo.k2d; do
    [[ -r "$KDB/$FILE" && -s "$KDB/$FILE" ]] || { echo "Missing shared index: $KDB/$FILE. No download will be attempted."; exit 1; }
done
printf 'path\t%s\ncoverage\tunverified; inspect cluster database documentation\n' "$KDB" > "$WD/metadata/kraken_database.tsv"
stat -c '%n %s %Y' "$KDB/hash.k2d" "$KDB/opts.k2d" "$KDB/taxo.k2d" > "$WD/metadata/kraken_index_files.txt"
STATUS="$WD/screen/screen_status.tsv"
printf 'run\tdatabase\tcommand_status\tutc\n' > "$STATUS"
FAILED=0
for RUN in SRR27384773 SRR1185280; do
    R1="$WD/raw/${RUN}_1.fastq.gz"; R2="$WD/raw/${RUN}_2.fastq.gz"
    test -s "$R1"; test -s "$R2"
    for DATABASE in nt refseq silva; do
        echo "SendSketch: $RUN $DATABASE"
        # All input reads, but comparisons still use sketches.
        if gzip -dc "$R1" "$R2" | sendsketch.sh -Xmx4g in=stdin.fq interleaved=f threads="$THREADS" "$DATABASE" samplerate=1 minkeycount=2 records=100 color=f > "$WD/screen/sendsketch/$RUN.$DATABASE.txt" 2> "$WD/logs/$RUN.$DATABASE.sendsketch.log"; then
            STATE=returned_success_inspect_report
        else
            STATE=failed; FAILED=1
        fi
        printf '%s\t%s\t%s\t%s\n' "$RUN" "$DATABASE" "$STATE" "$(date -u +%FT%TZ)" >> "$STATUS"
    done
    # Shared database is read-only. Ensure RAM allocation accommodates it.
    # Paired mode reports fragments, not individual-read fractions.
    if kraken2 --db "$KDB" --threads "$THREADS" --paired --gzip-compressed --confidence 0.1 --report "$WD/screen/kraken2/$RUN.report.tsv" --output "$WD/screen/kraken2/$RUN.classifications.tsv" "$R1" "$R2" > "$WD/logs/$RUN.kraken.stdout.log" 2> "$WD/logs/$RUN.kraken.stderr.log"; then
        STATE=completed
    else
        STATE=failed; FAILED=1
    fi
    printf '%s\tkraken2_shared\t%s\t%s\n' "$RUN" "$STATE" "$(date -u +%FT%TZ)" >> "$STATUS"
done
echo "Review screen reports and SRA STAT. Select candidates before step 03."
exit "$FAILED"
