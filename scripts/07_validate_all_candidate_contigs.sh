#!/usr/bin/env bash
# Validate selected contigs against online nt, one small candidate batch at a time.
# Uses the successful output of 06b_select_all_candidate_contigs.sh.
# No database download. Local matches select queries; remote matches test identity.
set -Eeuo pipefail
WD="/scratch/aubaxp004/fusarium_screen_v2"
trap 'echo "ERROR at line $LINENO; inspect the remote/*.log files" >&2' ERR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
for TOOL in blastn sha256sum awk flock; do
    command -v "$TOOL" >/dev/null || { echo "Missing tool: $TOOL" >&2; exit 1; }
done
POINTER="$WD/validation/all_candidates/latest_selection.txt"
[[ -s "$POINTER" ]] || { echo "Complete the new step 06b first" >&2; exit 1; }
IFS= read -r OUT < "$POINTER"
[[ "$OUT" == "$WD/validation/all_candidates/run."* && -s "$OUT/query_batches.tsv" ]] || exit 1
exec 9> "$OUT/remote.lock"
flock -n 9 || { echo "Validation already running for $OUT" >&2; exit 1; }
printf 'run\tcandidate\tstatus\tresult_file\n' > "$OUT/remote_status.tsv"
while IFS=$'\t' read -r RUN CANDIDATE QUERY; do
    [[ "$RUN" == run ]] && continue
    [[ -s "$QUERY" ]] || { echo "Missing query: $QUERY" >&2; exit 1; }
    N=$(awk '/^>/{n++} END{print n+0}' "$QUERY")
    (( N>=1 && N<=5 )) || { echo "Expected 1-5 contigs per batch" >&2; exit 1; }
    PREFIX="$OUT/remote/$RUN.$CANDIDATE"
    KEY=$(sha256sum "$QUERY" | awk '{print $1}')
    # A completed search with no hits is still a valid result (outfmt 7 has comments).
    if [[ -s "$PREFIX.complete" && -s "$PREFIX.blast.tsv" ]] && [[ "$(cat "$PREFIX.complete")" == "$KEY" ]]; then
        STATUS=REUSED
    else
        echo "Remote nt validation: $RUN / $CANDIDATE ($N contigs)"
        blastn -remote -db nt -task blastn -query "$QUERY" -evalue 1e-10 \
            -max_target_seqs 25 \
            -outfmt '7 qseqid saccver pident length qlen qstart qend sstart send evalue bitscore staxids stitle' \
            -out "$PREFIX.blast.tsv.part" 2> "$PREFIX.log"
        mv "$PREFIX.blast.tsv.part" "$PREFIX.blast.tsv"
        date -u +%FT%TZ > "$PREFIX.date.txt"
        printf '%s\n' "$KEY" > "$PREFIX.complete"
        STATUS=COMPLETED
        # Run serially and leave a short gap between remote requests.
        sleep 10
    fi
    printf '%s\t%s\t%s\t%s\n' "$RUN" "$CANDIDATE" "$STATUS" "$PREFIX.blast.tsv" >> "$OUT/remote_status.tsv"
done < "$OUT/query_batches.tsv"

# Small review bundle: measurements, selections, local candidate alignments,
# remote outputs and queries. Excludes genomes, indexes and full assemblies.
BUNDLE="$OUT/review_bundle"
mkdir -p "$BUNDLE"
cp "$OUT/selection_summary.tsv" "$OUT/selected_contigs.tsv" "$OUT/remote_status.tsv" "$OUT/reference_sequences.tsv" "$BUNDLE/"
cp "$OUT/"*.candidate_hits.tsv "$BUNDLE/"
cp -r "$OUT/queries" "$OUT/remote" "$BUNDLE/"
if [[ -s "$WD/results/measurements.tsv" ]]; then cp "$WD/results/measurements.tsv" "$BUNDLE/"; fi
tar -czf "$OUT/candidate_validation_results.tar.gz" -C "$OUT" review_bundle
echo "Finished. Upload: $OUT/candidate_validation_results.tar.gz"
echo "Do not call a conserved/host-like match contamination without reviewing alternatives."
