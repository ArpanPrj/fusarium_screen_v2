#!/usr/bin/env bash
# Follow-up assembly validation for EVERY candidate in the reference panel.
# Does not rerun read screening/mapping or download any reference database.
# makeblastdb indexes the already downloaded combined FASTA locally.
# Fusarium stays in the search so host matches can be compared with candidates.
set -Eeuo pipefail
WD="/scratch/aubaxp004/fusarium_screen_v2"
trap 'echo "ERROR at line $LINENO; inspect $WD/logs" >&2' ERR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
export LC_ALL=C
THREADS=8
MAX_PER_CANDIDATE=2
for TOOL in blastn makeblastdb samtools awk sort sha256sum flock; do
    command -v "$TOOL" >/dev/null || { echo "Missing tool: $TOOL" >&2; exit 1; }
done
MAP="$WD/metadata/reference_sequences.tsv"
REF="$WD/refs/combined.fa"
DNA="$WD/assembly/dna_megahit/final.contigs.fa"
RNA="$WD/assembly/rna_rnaspades/transcripts.fasta"
for FILE in "$MAP" "$REF" "$DNA" "$RNA"; do
    [[ -s "$FILE" ]] || { echo "Missing/empty input: $FILE" >&2; exit 1; }
done
mkdir -p "$WD/validation/all_candidates" "$WD/logs"
exec 9> "$WD/validation/all_candidates/selection.lock"
flock -n 9 || { echo "Another selection job is active" >&2; exit 1; }
OUT=$(mktemp -d "$WD/validation/all_candidates/run.XXXXXX")
mkdir -p "$OUT/queries" "$OUT/remote"
cp "$MAP" "$OUT/reference_sequences.tsv"
awk -F '\t' 'NR>1 && $3=="candidate" {print $2}' "$MAP" | sort -u > "$OUT/candidates.txt"
[[ -s "$OUT/candidates.txt" ]] || { echo "No candidates in reference map" >&2; exit 1; }

# Reuse a local BLAST index only for this exact reference content.
# Version 4 avoids the large virtual-address reservation of LMDB databases.
HASH=$(sha256sum "$REF" | awk '{print $1}')
DBDIR="$WD/validation/all_candidates/db_$HASH"
if [[ ! -s "$DBDIR/complete" ]]; then
    BUILD=$(mktemp -d "$WD/validation/all_candidates/dbbuild.XXXXXX")
    makeblastdb -in "$REF" -dbtype nucl -parse_seqids -blastdb_version 4 \
        -out "$BUILD/panel" > "$BUILD/build.log" 2>&1
    printf '%s\n' "$HASH" > "$BUILD/complete"
    [[ ! -e "$DBDIR" ]] || { echo "Incomplete DB directory: $DBDIR" >&2; exit 1; }
    mv "$BUILD" "$DBDIR"
fi
printf 'run\tcandidate\tcontigs_passing\tcontigs_selected\n' > "$OUT/selection_summary.tsv"
printf 'run\tcandidate\tquery_fasta\n' > "$OUT/query_batches.tsv"
printf 'query_id\trun\tcandidate\toriginal_contig\tassembly\n' > "$OUT/selected_contigs.tsv"

for RUN in SRR27384773 SRR1185280; do
    if [[ "$RUN" == SRR27384773 ]]; then ASSEMBLY="$DNA"; else ASSEMBLY="$RNA"; fi
    samtools faidx "$ASSEMBLY"
    echo "Searching all assembled sequences from $RUN against host + candidates."
    # A nucleotide search, not a splice-aware transcript-to-genome alignment.
    # RNA hits may cover individual exons only. No hits != organism absent.
    # Word size 20 balances runtime and sensitivity for these reference matches.
    blastn -task blastn -word_size 20 -num_threads "$THREADS" \
        -query "$ASSEMBLY" -db "$DBDIR/panel" -evalue 1e-10 \
        -max_target_seqs 100 -max_hsps 5 \
        -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp' \
        -out "$OUT/$RUN.raw.tsv" 2> "$OUT/$RUN.blast.log"

    # Record candidate hits alongside the strongest reported host HSP.
    # A higher host score is a warning, not an automatic discard: chimeric
    # contigs and different aligned regions require inspection of coordinates.
    awk -v map="$MAP" 'BEGIN{FS=OFS="\t";
        while((getline < map)>0) if($1!="sequence_id"){group[$1]=$2;role[$1]=$3}
        close(map)
    }
    FNR==NR {if(role[$2]=="host" && $12>host[$1]) host[$1]=$12; next}
    role[$2]=="candidate" {
        flag=(host[$1]>=$12 ? "HOST_HSP_GE_CANDIDATE" : "CANDIDATE_HSP_GT_HOST");
        print group[$2],$0,host[$1]+0,flag
    }' "$OUT/$RUN.raw.tsv" "$OUT/$RUN.raw.tsv" > "$OUT/$RUN.candidate_hits.raw.tsv"
    {
        printf 'candidate\tcontig\treference\tidentity_percent\talignment_bp\tcontig_bp\treference_bp\tquery_start\tquery_end\tref_start\tref_end\tevalue\tbitscore\tquery_coverage_hsp_percent\tbest_host_hsp_bitscore\thost_comparison\n'
        cat "$OUT/$RUN.candidate_hits.raw.tsv"
    } > "$OUT/$RUN.candidate_hits.tsv"

    while IFS= read -r CANDIDATE; do
        # Selection thresholds, NOT thresholds for biological significance.
        # Include candidates even when step 04 reported zero mapped reads.
        awk -F '\t' -v c="$CANDIDATE" '$1==c && $4>=80 && $5>=300 && $12<=1e-20' \
            "$OUT/$RUN.candidate_hits.raw.tsv" |
            sort -t $'\t' -k13,13gr -k5,5nr |
            awk -F '\t' '!seen[$2]++ {print $2}' > "$OUT/$RUN.$CANDIDATE.passing_ids.txt"
        awk -v n="$MAX_PER_CANDIDATE" 'NR<=n' "$OUT/$RUN.$CANDIDATE.passing_ids.txt" \
            > "$OUT/$RUN.$CANDIDATE.selected_ids.txt"
        PASS=$(awk 'END{print NR+0}' "$OUT/$RUN.$CANDIDATE.passing_ids.txt")
        N=$(awk 'END{print NR+0}' "$OUT/$RUN.$CANDIDATE.selected_ids.txt")
        printf '%s\t%s\t%s\t%s\n' "$RUN" "$CANDIDATE" "$PASS" "$N" >> "$OUT/selection_summary.tsv"
        (( N > 0 )) || continue
        QUERY="$OUT/queries/$RUN.$CANDIDATE.fa"
        samtools faidx -r "$OUT/$RUN.$CANDIDATE.selected_ids.txt" "$ASSEMBLY" |
            awk -v run="$RUN" '/^>/{sub(/^>/, ">" run "_")} {print}' > "$QUERY"
        printf '%s\t%s\t%s\n' "$RUN" "$CANDIDATE" "$QUERY" >> "$OUT/query_batches.tsv"
        awk -v r="$RUN" -v c="$CANDIDATE" -v a="$ASSEMBLY" 'BEGIN{OFS="\t"}
            {print r "_" $1,r,c,$1,a}' "$OUT/$RUN.$CANDIDATE.selected_ids.txt" >> "$OUT/selected_contigs.tsv"
    done < "$OUT/candidates.txt"
done
printf '%s\n' "$OUT" > "$WD/validation/all_candidates/latest_selection.txt.part"
mv "$WD/validation/all_candidates/latest_selection.txt.part" "$WD/validation/all_candidates/latest_selection.txt"
echo "Completed: $OUT"
echo "Review selection_summary.tsv, then submit 07_validate_all_candidate_contigs.sh."
echo "No selected contigs for a candidate is a valid result, not an error."
