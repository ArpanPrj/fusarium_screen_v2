#!/usr/bin/env bash
# Competitive DNA and splice-aware RNA alignment; candidate measurements.
# Reads = individual mates. Denominator = every raw input read.
# MAPQ >=20; retain primary mapped alignments only. No duplicate removal.
# Depth counts overlapping mates separately and uses base quality >=0.
set -Eeuo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
WD="/scratch/aubaxp004/fusarium_screen_v2"
THREADS=8
mkdir -p "$WD"/{mapping,results,logs}
test -s "$WD/refs/combined.fa.fai"
MAP="$WD/metadata/reference_sequences.tsv"
TABLE="$WD/results/measurements.tsv"
printf 'run\tcandidate\tinput_reads\tmapped_reads\tread_fraction\tread_percent\treference_bp\tcovered_bp\tbreadth_percent\tmean_depth\n' > "$TABLE"
for RUN in SRR27384773 SRR1185280; do
    R1="$WD/raw/${RUN}_1.fastq.gz"; R2="$WD/raw/${RUN}_2.fastq.gz"
    PREFIX="$WD/mapping/$RUN"
    # Count four-line FASTQ records. These are SRA-extracted FASTQs.
    N1=$(gzip -dc "$R1" | awk 'END {if (NR==0 || NR%4) exit 1; printf "%.0f",NR/4}')
    N2=$(gzip -dc "$R2" | awk 'END {if (NR==0 || NR%4) exit 1; printf "%.0f",NR/4}')
    [[ "$N1" == "$N2" ]] || { echo "Unequal mate counts: $RUN"; exit 1; }
    TOTAL=$((N1 + N2))
    if [[ "$RUN" == SRR27384773 ]]; then
        bowtie2 --very-sensitive -p "$THREADS" -x "$WD/refs/index/dna" -1 "$R1" -2 "$R2" 2> "$WD/logs/$RUN.mapping.log" |
            samtools sort -@ 1 -m 1G -o "$PREFIX.all.bam" -
    else
        # Unstranded default: do not assume RNA library strand orientation.
        hisat2 -p "$THREADS" -x "$WD/refs/index/rna" -1 "$R1" -2 "$R2" 2> "$WD/logs/$RUN.mapping.log" |
            samtools sort -@ 1 -m 1G -o "$PREFIX.all.bam" -
    fi
    samtools index "$PREFIX.all.bam"
    # 2308 = unmapped + secondary + supplementary flags.
    samtools view -b -q 20 -F 2308 "$PREFIX.all.bam" > "$PREFIX.filtered.bam"
    samtools index "$PREFIX.filtered.bam"
    samtools coverage --ff 0 -q 0 -Q 0 "$PREFIX.filtered.bam" > "$PREFIX.coverage.tsv"
    samtools idxstats "$PREFIX.filtered.bam" > "$PREFIX.idxstats.tsv"
    # Exact coverage sums; zero-covered bases are included in the denominator
    # using full reference lengths. D/N CIGAR operations do not count as bases.
    samtools depth -q 0 -Q 0 -g 1792 "$PREFIX.filtered.bam" |
        awk 'BEGIN{OFS="\t"} {if($3>0) covered[$1]++; depth[$1]+=$3}
             END{for(s in depth) print s,covered[s]+0,depth[s]}' > "$PREFIX.depth_sums.tsv"
    awk -v run="$RUN" -v total="$TOTAL" -v map="$MAP" -v dep="$PREFIX.depth_sums.tsv" '
    BEGIN {FS=OFS="\t";
        while ((getline < map)>0) if($1!="sequence_id") {group[$1]=$2; role[$1]=$3}
        close(map);
        while ((getline < dep)>0) {cov[$1]=$2; ds[$1]=$3}
        close(dep)
    }
    $1!="*" && role[$1]=="candidate" {
        g=group[$1]; len[g]+=$2; reads[g]+=$3; covered[g]+=cov[$1]; sums[g]+=ds[$1]
    }
    END {for(g in len) printf "%s\t%s\t%.0f\t%.0f\t%.9g\t%.6f\t%.0f\t%.0f\t%.6f\t%.6f\n",run,g,total,reads[g],reads[g]/total,100*reads[g]/total,len[g],covered[g],100*covered[g]/len[g],sums[g]/len[g]}
    ' "$PREFIX.idxstats.tsv" >> "$TABLE"
done
echo "Saved: $TABLE. RNA genomic breadth is expression-dependent; interpret accordingly."
