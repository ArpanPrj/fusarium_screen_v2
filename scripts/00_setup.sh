#!/usr/bin/env bash
# Run once. Installs software only, never screening databases.
set -Eeuo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
WD="/scratch/aubaxp004/fusarium_screen_v2"
mkdir -p "$WD"/{scripts,raw/sra,tmp,logs,metadata,screen/sendsketch,screen/kraken2,refs/candidates,mapping,assembly,results,validation}
if conda env list | awk '{print $1}' | grep -qx contamination; then
    conda install -y -n contamination -c conda-forge -c bioconda --strict-channel-priority sra-tools bbmap kraken2 bowtie2 hisat2 samtools megahit spades blast seqtk curl
else
    conda create -y -n contamination -c conda-forge -c bioconda --strict-channel-priority sra-tools bbmap kraken2 bowtie2 hisat2 samtools megahit spades blast seqtk curl
fi
conda activate contamination
conda list --explicit > "$WD/metadata/environment.txt"
printf 'run\treadset\nSRR27384773\tDNA_WGS\nSRR1185280\tRNA_seq\n' > "$WD/metadata/samples.tsv"
# This manifest is deliberately not populated with guessed detections.
# One row per reference, TAB separated. role = host or candidate.
# id: letters/digits/single underscores; no double underscores.
# url: direct nucleotide FASTA or FASTA.gz download URL.
# accession: exact accession.version; include all segments of segmented viruses.
if [[ ! -e "$WD/metadata/references.tsv" ]]; then
    printf 'role\tid\taccession\turl\n' > "$WD/metadata/references.tsv"
fi
cat > "$WD/metadata/tool_rationale.tsv" <<'EOF'
tool	reason_for_selection
SRA_Toolkit	Download original reads by accession.
SendSketch	Remote candidate screening without database downloads.
Kraken2	Read-pair classification using the existing cluster database.
Bowtie2	Competitive alignment of short DNA reads.
HISAT2	Splice-aware competitive alignment of RNA reads.
SAMtools	Alignment filtering and read/depth/breadth measurements.
MEGAHIT	Exploratory assembly of the DNA readset.
rnaSPAdes	Transcriptome assembly of the RNA readset as suggested in the assignment.
BLASTn	Remote validation of a few selected assembled contigs.
EOF
echo "Setup complete. Scripts use $WD and never the previous project."
