#!/usr/bin/env bash
# Step 03: selected reference panel from the two screening reports.
# Downloads individual mapping references, NOT a screening database.
# No subsampling. Outputs are compatible with steps 04-08.
# Plan: validate software; write panel; fetch/validate missing FASTAs;
# combine host + candidates; build DNA and splice-aware RNA indexes.
# curl: HTTPS downloads, avoiding the blocked rsync port.
# samtools: reference lengths and FASTA indexing.
# Bowtie2: sensitive short-read DNA mapping reference.
# HISAT2: splice-aware short-read RNA mapping reference.
set -Eeuo pipefail
trap 'echo "ERROR at line $LINENO. Check $WD/logs and the URL printed above." >&2' ERR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
WD="/scratch/aubaxp004/fusarium_screen_v2"
THREADS="${THREADS:-8}"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid THREADS" >&2; exit 1; }
mkdir -p "$WD"/{refs/downloads,metadata,logs}
for TOOL in curl gzip awk sha256sum samtools bowtie2-build hisat2-build flock; do
    command -v "$TOOL" >/dev/null || { echo "Missing software: $TOOL (environment contamination)" >&2; exit 1; }
done
# Prevent simultaneous jobs from altering the same reference panel.
exec 9> "$WD/refs/prepare.lock"
flock -n 9 || { echo "Another reference preparation job is active." >&2; exit 1; }
MANIFEST="$WD/metadata/references.tsv"
cat > "$WD/metadata/references.selected.tsv" <<'PANEL'
role	id	accession	url
host	Fgraminearum_PH1	GCF_000240135.3	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/135/GCF_000240135.3_ASM24013v3/GCF_000240135.3_ASM24013v3_genomic.fna.gz
host	Fgraminearum_mito	NC_009493.1	https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_009493.1&rettype=fasta&retmode=text
candidate	FgV1	NC_006937.1	https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_006937.1&rettype=fasta&retmode=text
candidate	Bacillus_cereus	GCF_000007825.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/825/GCF_000007825.1_ASM782v1/GCF_000007825.1_ASM782v1_genomic.fna.gz
candidate	Homo_sapiens	GCF_000001405.40	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
candidate	Peribacillus_frigoritolerans	GCA_024169475.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/169/475/GCA_024169475.1_ASM2416947v1/GCA_024169475.1_ASM2416947v1_genomic.fna.gz
candidate	Stenotrophomonas_maltophilia	GCF_000072485.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/072/485/GCF_000072485.1_ASM7248v1/GCF_000072485.1_ASM7248v1_genomic.fna.gz
candidate	Enterobacter_T2	GCA_010692925.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/010/692/925/GCA_010692925.1_ASM1069292v1/GCA_010692925.1_ASM1069292v1_genomic.fna.gz
candidate	Bacillus_velezensis	GCF_000015785.2	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/015/785/GCF_000015785.2_ASM1578v2/GCF_000015785.2_ASM1578v2_genomic.fna.gz
candidate	Spirosoma_pollinicola	GCA_002831565.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/831/565/GCA_002831565.1_ASM283156v1/GCA_002831565.1_ASM283156v1_genomic.fna.gz
PANEL
# Preserve a previous custom manifest before installing this selected panel.
if [[ -f "$MANIFEST" ]] && ! cmp -s "$MANIFEST" "$WD/metadata/references.selected.tsv"; then
    cp "$MANIFEST" "$MANIFEST.backup.$(date -u +%Y%m%dT%H%M%S).$$"
fi
cp "$WD/metadata/references.selected.tsv" "$MANIFEST"
BUILD=$(mktemp -d "$WD/refs/build.XXXXXX")
COMBINED="$BUILD/combined.fa"
MAP="$BUILD/reference_sequences.tsv"
: > "$COMBINED"
printf 'sequence_id\tgroup\trole\taccession\toriginal_id\n' > "$MAP"
valid_fasta() {
    awk '
    /^>/ {if(n && !bases) bad=1; n++; bases=0; next}
    {gsub(/[[:space:]]/, ""); if(length($0)) {
        if(!n || toupper($0) ~ /[^ACGTRYSWKMBDHVNU]/) bad=1;
        bases+=length($0)
    }}
    END {exit (bad || n==0 || bases==0)}' "$1"
}
while IFS=$'\t' read -r ROLE ID ACCESSION URL; do
    [[ "$ROLE" == role ]] && continue
    FILE="$WD/refs/downloads/${ID}.${ACCESSION}.fa"
    KEY=$(printf '%s\n' "$ACCESSION" "$URL" | sha256sum | awk '{print $1}')
    EXPECTED="$KEY"
    if [[ -s "$FILE" && -s "$FILE.source" && -s "$FILE.sha256" ]] &&
       [[ "$(cat "$FILE.source")" == "$EXPECTED" ]] &&
       sha256sum --status -c "$FILE.sha256" && valid_fasta "$FILE"; then
        echo "Reusing validated reference: $ID ($ACCESSION)"
    else
        echo "Downloading $ID ($ACCESSION): $URL"
        curl --fail --location --retry 5 --connect-timeout 60 "$URL" -o "$FILE.download"
        if [[ "$URL" == *.gz ]]; then
            gzip -t "$FILE.download"
            gzip -dc "$FILE.download" > "$FILE.part"
        else
            cp "$FILE.download" "$FILE.part"
        fi
        valid_fasta "$FILE.part" || { echo "Invalid/incomplete FASTA: $ID" >&2; exit 1; }
        # Fixed nucleotide accessions must match the returned FASTA header.
        if [[ "$ACCESSION" == NC_* ]]; then
            awk -v acc="$ACCESSION" 'NR==1 {exit $1 != ">"acc}' "$FILE.part"
        fi
        mv "$FILE.part" "$FILE"
        printf '%s\n' "$EXPECTED" > "$FILE.source"
        sha256sum "$FILE" > "$FILE.sha256"
        rm -f "$FILE.download"
        printf '%s\t%s\t%s\t%s\n' "$ID" "$ACCESSION" "$URL" "$(date -u +%FT%TZ)" >> "$WD/metadata/reference_downloads.tsv"
    fi
    # Mitochondria are supplied separately; remove any mitochondrial records
    # from the PH-1 assembly to avoid duplicate competing host sequences.
    awk -v group="$ID" -v role="$ROLE" -v acc="$ACCESSION" -v map="$MAP" '
    /^>/ {
        skip=(group=="Fgraminearum_PH1" && (tolower($0) ~ /mitochond/ || $1 ~ /^>NC_009493[.]/));
        if(skip) next;
        n++; original=substr($0,2); sub(/[ \t].*/,"",original);
        id=group "__" n; print ">" id;
        print id "\t" group "\t" role "\t" acc "\t" original >> map; next
    }
    !skip {gsub(/\r/,""); print}
    END {if(n==0) exit 1}
    ' "$FILE" >> "$COMBINED"
done < "$MANIFEST"
valid_fasta "$COMBINED"
samtools faidx "$COMBINED"
mkdir -p "$BUILD/index"
echo "Building DNA index; human reference dominates time and memory."
bowtie2-build --threads "$THREADS" "$COMBINED" "$BUILD/index/dna" > "$WD/logs/dna.index.log" 2>&1
echo "Building splice-aware RNA index."
hisat2-build -p "$THREADS" "$COMBINED" "$BUILD/index/rna" > "$WD/logs/rna.index.log" 2>&1
# Publish only after both builders succeed. Keep an older index if present.
if [[ -e "$WD/refs/index" ]]; then
    mv "$WD/refs/index" "$WD/refs/index.previous.$(date -u +%Y%m%dT%H%M%S).$$"
fi
mv "$BUILD/index" "$WD/refs/index"
mv "$COMBINED" "$WD/refs/combined.fa"
mv "$COMBINED.fai" "$WD/refs/combined.fa.fai"
mv "$MAP" "$WD/metadata/reference_sequences.tsv"
sha256sum "$WD/refs/combined.fa" > "$WD/metadata/combined_reference.sha256"
sha256sum "$MANIFEST" > "$WD/metadata/references_manifest.sha256"
rmdir "$BUILD"
echo "Step 03 complete. Run 04_map_measure.sh after this job finishes."
echo "If previously mapped, rerun step 04 against this panel."
