#!/usr/bin/env bash
# Prepare reviewable submission templates; never fabricate interpretations.
# Does not publish, push to GitHub, or submit coursework.
set -Eeuo pipefail
WD="/scratch/aubaxp004/fusarium_screen_v2"
mkdir -p "$WD/results" "$WD/submission"
test -s "$WD/results/measurements.tsv"
TARGET="$WD/results/detection_table.tsv"
if [[ ! -e "$TARGET" ]]; then
    awk 'BEGIN{FS=OFS="\t"} NR==1{print $0,"status","bin","likely_source","supporting_evidence","alternatives_and_limits";next}
    {print $0,"REVIEW_REQUIRED","UNASSIGNED","UNRESOLVED","ADD_EVIDENCE","ADD_LIMITATIONS"}' "$WD/results/measurements.tsv" > "$TARGET"
else
    echo "Preserving existing detection_table.tsv; update measurements manually if remapped."
fi
if [[ ! -e "$WD/results/writeup.md" ]]; then
cat > "$WD/results/writeup.md" <<'EOF'
# Most interesting detection — one-page writeup

## Sample and question
State the run, host, library type, and candidate accession.version.

## Methods
State screening references/date, candidate selection, competitive mapping,
MAPQ filter, and that read fractions count individual reads from raw FASTQs.

## Results
Insert measured read fraction, mean depth, breadth, and coverage distribution.
Explain contig/BLAST support if available. Do not treat screen hits as proof.

## Interpretation
Choose noise, trace, consequential, or biological and explain why.
Give the likely source, evidence over alternatives, and remaining uncertainty.
RNA depth reflects expression; genomic breadth is not expected to be uniform.
Do not claim reagent contamination or index hopping without supporting evidence.
EOF
fi
if [[ ! -e "$WD/results/agent_transcript.md" ]]; then
cat > "$WD/results/agent_transcript.md" <<'EOF'
# Agent record
The assistant drafted Bash scripts. It did not execute the HPC analyses.
The student executes the workflow and makes the final bin/source decisions.
Paste/export the actual conversation here; this note is not the full transcript.
EOF
fi
cp "$TARGET" "$WD/results/writeup.md" "$WD/results/agent_transcript.md" "$WD/submission/"
cp -r "$WD/scripts" "$WD/metadata" "$WD/submission/"
tar -czf "$WD/submission_bundle.tar.gz" -C "$WD" submission
echo "Edit interpretations and writeup, then rerun to refresh the submission bundle."
echo "Include a repository URL separately if your instructor requires one."
