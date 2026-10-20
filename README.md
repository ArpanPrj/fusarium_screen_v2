# Fusarium non-target sequence screening

This project compares a DNA readset and an RNA readset, checks screening candidates by competitive reference mapping, and examines assembled sequences with local BLAST and GC-versus-depth plots. It uses the complete readsets, without subsampling.

## Samples and project location

| Run | Library used in this project | Individual reads in the uploaded measurement table |
|---|---|---:|
| SRR27384773 | DNA/WGS | 12,446,790 |
| SRR1185280 | RNA-seq | 27,608,126 |

Working directory: `/scratch/aubaxp004/fusarium_screen_v2`.
Conda environment: `contamination`.
Shared Kraken2 database: `/datasets/unzipped/kraken2_full`.

The accessions identify two separate datasets. They should not be treated as paired DNA/RNA measurements from the same biological specimen. The figures and results here are based on the supplied files; library strandedness, treatment history, and controls have not been independently verified.

## What is complete, and what is not

The supplied results support completion of read screening, selected-reference mapping, DNA/RNA assembly, local panel BLAST, and original blobplots. Corrected tables and plots were generated from those existing outputs, and `10_correct_blobplots.sh` was provided to reproduce the correction on ASC. That Bash/R script was syntax-checked but was not executed by the assistant on ASC; the student's corrected job completion has not yet been reported.

Remote nt BLAST was attempted but skipped because of its long runtime. The older Kraken-based assembly plotting job was replaced by the local-BLAST plotting workflow. Neither abandoned job is treated as completed validation.

No Bracken, sourmash, UNITE, geNomad, DIAMOND, CheckV, FCS-GX, or BlobToolKit analysis is included in these final results. The plots are custom GC-versus-depth plots, not BlobToolKit datasets. No spike-in experiment was performed, so no empirical detection limit is claimed.

## Files to keep in the repository

Place `.gitignore`, `environment.yml`, this `README.md`, and `summary.md` in the repository root. Keep the actual Bash scripts in `scripts/`. Preserve small screening reports in `screen/`, reference manifests and software records in `metadata/`, and selected final tables and plots in `results/`.

The `.gitignore` excludes raw reads, downloaded references, alignment files, assemblies, BLAST working directories, plotting workspaces, and large per-read classifications. It deliberately keeps small text reports and curated results. Copy desired final plots from the ignored working directories into `results/figures/` before committing them. It does not remove files that Git already tracks.

This documentation package does not contain a newly published Git repository or the raw data. Add it to the existing project with the real scripts and results. The PDF is a one-page version of `summary.md`.

## Software environment

For a NEW installation, from the directory containing `environment.yml`:

```bash
source ~/miniconda3/etc/profile.d/conda.sh
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
conda activate contamination
```

The command applies strict channel priority for this installation. If the `contamination` environment already works, keep it for finishing the project. To intentionally update it later:

```bash
CONDA_CHANNEL_PRIORITY=strict conda env update -n contamination -f environment.yml
```

The YAML lists required packages; it is not a lockfile and does not claim to reproduce the exact original versions. No environment installation or dependency solve was run while preparing these files. Before changing a working environment, capture its actual versions:

```bash
WD="/scratch/aubaxp004/fusarium_screen_v2"
mkdir -p "$WD/metadata"
conda list -n contamination --explicit > "$WD/metadata/conda-linux-explicit.txt"
conda env export -n contamination --no-builds > "$WD/metadata/environment-resolved.yml"
```

The explicit specification is platform-specific. Review exported paths/URLs before publishing. Database versions and reference FASTAs are separate from the software environment; keep their provenance too.

| Tool | Why it was included |
|---|---|
| SRA Toolkit | Downloads runs and extracts paired FASTQ reads. |
| BBTools/SendSketch | Quickly compares read sketches with online reference collections. |
| Kraken2 | Screens reads against the existing shared database without downloading another database. |
| Bowtie2 | Maps DNA reads to the combined panel and both libraries to their own assemblies. |
| HISAT2 | Allows splice-aware mapping of RNA reads to genomic references. |
| SAMtools | Sorts, filters, indexes, counts alignments, and measures coverage. |
| MEGAHIT | Provides an exploratory assembly of the DNA readset. |
| SPAdes/rnaSPAdes | Provides an RNA-specific transcript assembly. |
| BLAST+ | Compares contigs/transcripts to the selected local reference panel. |
| Base R | Draws plots and recalculates query-interval support without extra R packages. |
| curl / Java | Provides HTTPS reference downloads and the BBTools runtime. |

ASC also supplies standard shell tools including Bash, awk, gzip, tar, sha256sum, and flock. SPAdes and other packages may bring Python as an internal dependency; the user-facing pipeline scripts are Bash with embedded R.

## Running jobs and choosing resources

Each computational script activates Conda explicitly and uses an absolute project path, because ASC may run it in a temporary working directory. Submit one script with:

```bash
cd /scratch/aubaxp004/fusarium_screen_v2/scripts
/scripts/run_script 04_map_measure.sh
```

The submission wrapper is specific to ASC. Submission order alone does not create dependencies: wait for prerequisite jobs to finish successfully, or use scheduler dependencies appropriate to the cluster. A job disappearing from `qstat` does not prove success; check its logs and expected outputs.

Most existing analysis scripts use eight threads. The assembly scripts specify 48 GB of memory; allocate more than that limit to allow overhead. The human-inclusive index build and full Kraken database can require substantial memory. Allocating 64 cores does not make a script with `THREADS=8` use 64 cores. The correction script only reads tables and plots and does not need the large assembly/database allocation.

Use `environment.yml` for dependency setup; the older `00_setup.sh` installs packages directly and does not include every plotting dependency in the new YAML. Do not run it merely to recreate an already working environment. Create the project folders if starting fresh:

```bash
WD="/scratch/aubaxp004/fusarium_screen_v2"
mkdir -p "$WD"/{scripts,raw/sra,tmp,logs,metadata,screen/sendsketch,screen/kraken2,refs,mapping,assembly,results,validation}
```

## Workflow and reasons for each step

### 01 - Download both complete readsets

Script: `01_download.sh`.

`prefetch` downloads SRA archives; `fasterq-dump --split-files` extracts mates. Each FASTQ is compressed and checked with `gzip -t`. Existing complete pairs are reused; partial pairs stop the script for inspection. Inputs are the two run IDs. Outputs are `raw/SRR27384773_1.fastq.gz`, its mate, and the corresponding RNA pair.

This step does not subsample or perform a separate adapter-trimming or read-quality filtering stage. Do not describe an unperformed QC stage in the methods. The mapping denominator is the raw extracted read count.

### 02 - Screen with Kraken2 and SendSketch

Script: `02_screen.sh`.

Kraken2 reads the shared `/datasets/unzipped/kraken2_full` index. It checks for `hash.k2d`, `opts.k2d`, and `taxo.k2d`; it does not download a replacement if they are absent. The run uses paired mode and confidence 0.1. The database name alone does not establish its content or build date; these require cluster documentation.

SendSketch processes both mates against the online `nt`, `refseq`, and `silva` services, with `samplerate=1`, `minkeycount=2`, and up to 100 returned records. It examines all input reads but represents their sequence content as sketches. Its scores are not direct read percentages. Save each returned report and the command-status table; an exit status alone is not a biological finding.

Outputs: `screen/kraken2/*.report.tsv`, per-fragment classification files, `screen/sendsketch/*.txt`, and `screen/screen_status.tsv`. Kraken's paired-mode counts are fragments/read pairs, whereas step 04 counts individual mates. Do not compare their raw counts without accounting for this difference.

### 03 - Prepare a host-plus-candidate reference panel

Script: `03_prepare_references.sh`.

Selected FASTAs are downloaded over HTTPS, validated, and reused when checksums/source records match. This is a small selected mapping panel relative to a broad screening database, although the human genome makes it large in absolute terms. The panel contains:

| Role | Group | Accession.version |
|---|---|---|
| Host | Fusarium graminearum PH-1 | GCF_000240135.3 |
| Host | F. graminearum mitochondrion | NC_009493.1 |
| Candidate | FgV1 | NC_006937.1 |
| Candidate | Bacillus cereus | GCF_000007825.1 |
| Candidate | Homo sapiens | GCF_000001405.40 |
| Candidate | Peribacillus frigoritolerans | GCA_024169475.1 |
| Candidate | Stenotrophomonas maltophilia | GCF_000072485.1 |
| Candidate | Enterobacter T2 | GCA_010692925.1 |
| Candidate | Bacillus velezensis | GCF_000015785.2 |
| Candidate | Spirosoma pollinicola | GCA_002831565.1 |

Reference headers are renamed to unique group-specific IDs. The correspondence to original accession IDs is saved in `metadata/reference_sequences.tsv`. Host mitochondrial duplication is avoided. The script builds Bowtie2 and HISAT2 indexes for `refs/combined.fa`.

Keeping host and candidates together makes mapping competitive: a read can match the host instead of being forced onto a candidate-only reference. This reduces some false assignments but cannot compensate for organisms or strains missing from the panel. The human FASTA includes alternate/patch regions, which can introduce ambiguous matches and affect the total reference-length denominator.

### 04 - Measure read fraction, breadth, and depth

Script: `04_map_measure.sh`; output: `results/measurements.tsv`.

DNA uses Bowtie2 `--very-sensitive`. RNA uses HISAT2 with no assumed strand orientation. Alignments are sorted, then filtered with `samtools view -q 20 -F 2308` to retain primary mapped alignments. Supplementary and secondary records are removed. The filter does not require proper pairs and does not remove duplicate or QC-fail flags. Overlapping mates count separately.

SAMtools generates coverage reports and index statistics. The final measurement table uses depth sums from `samtools depth -q 0 -Q 0 -g 1792` on the filtered BAM, plus full reference lengths; it is not simply an average of per-sequence coverage rows. With modern SAMtools, depth is uncapped. Deletions and skipped RNA introns are not counted as sequenced bases. Reference positions with no coverage remain in the mean-depth denominator.

Definitions, for a candidate with reference length L:

- Read fraction = retained primary mapped individual reads / all raw individual input reads.
- Read percent = read fraction times 100.
- Breadth percent = reference positions with depth greater than zero / L times 100.
- Mean depth = sum of per-position read depths / L, including zero-covered positions.

These quantities answer different questions. A few reads can match a repeat without supporting an entire organism. RNA genome breadth depends on expression and should not be judged by DNA uniformity expectations. RNA read percentages are not percentages of organisms or cells.

### 05 and 06 - Assemble DNA and RNA separately

Scripts: `05_assemble_dna.sh` and `06_assemble_rna.sh`.

MEGAHIT assembles all DNA reads with a minimum contig length of 500 bp, producing `assembly/dna_megahit/final.contigs.fa`. This is exploratory assembly, not a claim of a polished fungal genome.

rnaSPAdes assembles all RNA reads and produces `assembly/rna_rnaspades/transcripts.fasta`. No strand-specific setting was selected without verified library orientation. Its transcript labels and coverage fields should not be confused with independent organisms or mapped base depth.

The assemblies can run independently after downloads, provided each has its own resource allocation. Their scripts stop if the output directory already exists; inspect a prior run before using the documented resume command.

### 06b - Compare assembled sequences with the local panel

Use `06b_select_all_candidate_contigs.sh`, which supersedes the earlier FgV1-only selection script.

A local nucleotide BLAST index is built from the already downloaded combined panel and keyed by its checksum. The search uses word size 20, E-value 1e-10, at most 100 target sequences, and at most five HSPs per query/subject pair. An HSP is one local alignment segment. Raw output preserves query/reference coordinates, identity, bit score, and lengths.

The script records candidate matches and the strongest host HSP for comparison. Selection for optional follow-up requires identity >=80%, alignment length >=300 bp, and E-value <=1e-20. It selects up to two distinct contigs per candidate per sample. These are selection rules, not biological-significance thresholds. Candidate matches are retained for inspection even if the host has a stronger match.

Outputs are in `validation/all_candidates/run.*/`. The pointer `validation/all_candidates/latest_selection.txt` identifies a completed selection run. Preserve its `reference_sequences.tsv` and `SRR*.raw.tsv`; the plotting steps depend on them.

### 07 - Optional broader remote validation: skipped

`07_validate_all_candidate_contigs.sh` sends selected contigs to remote nt BLAST. This was too slow and was skipped. Empty remote logs or submitted requests are not completed validation. This step is not required to draw the local-panel plots. Do not claim that the local panel search confirms a globally best species match.

### 09 - Build plots using existing local BLAST results

Use `09_blobplots_local.sh`, not the superseded `09_blobplots.sh` Kraken-based version.

Both readsets are mapped to their own assemblies with Bowtie2, because the RNA targets here are assembled transcripts rather than genomic references. Primary alignments with MAPQ >=20 are retained. SAMtools coverage provides per-sequence depth and breadth; GC is calculated using A/C/G/T bases. The script combines these measurements with local BLAST evidence and writes `blobplots_local/SRR*_blobplot.tsv` and figures.

The original script classified full/partial support using only the strongest HSP's query coverage. Multiple exon-sized or fragmented matches could therefore make a well-supported sequence look partial. Step 10 corrects this display assessment.

### 10 - Correct combined alignment support and plot layout

Script: `10_correct_blobplots.sh`.

This step reads the original plot tables and the recorded step-06b BLAST source. It validates query IDs/lengths against the tables and subject IDs against the reference key. It does not rerun mapping, BLAST, screening, or assembly.

For the original best-reference group, it retains the same qualifying HSP filters, merges overlapping query-coordinate intervals, and counts each query position once. Combined query support is this union length divided by full query length. Groups within 95% of the best individual bit score are flagged as competing. Broad support means at least 80% query union; lower values are partial. Host nuclear and mitochondrial matches share a host group.

New fields are `query_union_bp`, `query_union_percent`, `support_status`, `reported_hsp_cap_reached`, and `best_group_hsp_cap_reached`. Original fields remain unchanged for audit: use the new support fields, not the old `reference_group`/`plot_group`, to interpret the corrected figures.

The figure includes sequences >=300 bp, with readable legends, visible small points, and a depth axis that is linear from zero to one and logarithmic above one. Gray means no qualifying match; it does not mean contaminant. Colored filled circles mean broad support, open circles partial support, and triangles competing groups. Marker area is deliberately not proportional to length because a minimum size is applied. Tables retain all sequences, including shorter ones.

Outputs go to a new `blobplots_corrected/run.*/` folder. `latest_corrected.txt` records the completed folder. PDF is always produced; PNG/SVG require Cairo support in R. Original depth/read counts remain unchanged. Base-R rendering differs slightly from the assistant's earlier figure previews.

### 08 - Package the final submission, after interpretation

`08_prepare_submission.sh` creates a detection-table template, a writeup template, and an agent-record placeholder, then packages scripts and metadata. Run or refresh it after the final interpretations are ready, despite its earlier number. It preserves existing files, so it does not replace old values after remapping. Check for `REVIEW_REQUIRED` and other placeholders before submission.

The package does not automatically include the newer corrected plots or these documentation files. Add selected final tables, figures, and the real conversation transcript explicitly. A short assistant-use note is not a complete agent transcript. No GitHub push or coursework submission has been performed by these scripts.

## Findings and proposed significance assessments

The table below reports reference-mapping measurements, not blobplot depths. Percentages use individual-read denominators. Bins and sources are provisional judgments based on the supplied results.

| Dataset / candidate | Retained reads | Read % | Reference breadth % | Mean depth (x) | Assessment and likely source |
|---|---:|---:|---:|---:|---|
| RNA / FgV1 | 14,497,528 | 52.511815 | 100.000000 | 220,615.475996 | Biological; also consequential. Likely virus associated with the sampled culture; natural versus experimental infection is unresolved. |
| DNA / human | 1,034 | 0.008307 | 0.003989 | 0.000047 | Provisional trace. Handling or library carryover is plausible; shared sequences/reference ambiguity remain alternatives. Harmlessness is not demonstrated for every downstream analysis. |
| RNA / human | 645 | 0.002336 | 0.001463 | 0.000019 | Provisional noise/unresolved. Sparse matching is not enough to establish human contamination or its source. |
| RNA / S. maltophilia | 3 | 0.000011 | 0.004617 | 0.000062 | Noise. Too few reads and too little breadth for a credible biological detection; incidental matching is more plausible than a substantial bacterial component. |

DNA had zero retained mapped reads for FgV1 and all six tested bacterial groups. RNA had zero for B. cereus, P. frigoritolerans, Enterobacter T2, B. velezensis, and S. pollinicola. These are non-detections under this panel/filter, not evidence that the samples are sterile or that related strains are absent. Any initial screening suggestions for these candidates lack mapping support here.

The four assignment bins are: noise (not credible), trace (credible but negligible for the analysis), consequential (changes an analysis), and biological (worth reporting in its own right). A finding can have both biological interest and analytical consequences; if one bin is required, use biological for FgV1 and state its consequences separately.

### Why FgV1 is the strongest finding

Over half the RNA reads mapped to the 6,624-base viral reference, and every reference base had retained coverage. This is much stronger than a few reads matching a conserved region. Assembly comparison provides additional sequence support: 54 transcripts had FgV1 as their best panel group, and 53 reached at least 80% combined query support. They are redundant assembled sequences, not 54 independent viruses.

The likely source is a viral component of the sampled fungal culture. These results do not establish infectivity, symptoms, natural infection, or exclude deliberate infection or substantial sample carryover. Those distinctions require sample metadata and controls. Do not use failure to detect FgV1 in the separate DNA run as a matched-sample test for integration.

### What changed in the corrected plots

| Measure | DNA | RNA |
|---|---:|---:|
| Assembly sequences in the tables | 260 | 20,510 |
| Broad reference support after correction | 67 | 12,403 |
| Partial support after correction | 65 | 1,050 |
| No qualifying panel hit | 128 | 7,057 |
| Partial labels that gained >=80% combined support | 38 | 2,164 |
| Query/subject pairs reaching the five-HSP limit | 402 | 706 |

There were 123 DNA sequences with a host best group and nine with a human best group; all nine human-best sequences reached >=80% query union. This strengthens a human-like sequence signal within the panel, but is not independent broad taxonomic validation. The RNA table had 13,399 host-best sequences. Of its 7,057 unmatched sequences, 6,614 were shorter than 300 bp and could not meet the minimum alignment-length rule.

## Remaining limits

- **Limited reference coverage:** a best panel match can be misleading if alternatives are absent. The exact shared Kraken database content/build remains unverified.
- **BLAST reporting caps:** five HSPs and 100 targets limit reported evidence. Reaching a cap suggests possible truncation, not proof of omitted hits. Union support may be incomplete; low support alone does not establish chimerism.
- **RNA mapping ambiguity:** similar assembled transcripts compete for reads. MAPQ filtering can make individual viral transcripts appear weakly covered even though the separate viral reference is deeply covered. Do not substitute the reference-genome depth into transcript points.
- **RNA interpretation:** expression, transcript redundancy, and library selection prevent direct conversion of read fraction to cell abundance. Assembly header coverage is not the measured mapped-base depth.
- **Source attribution:** there are no evaluated negative controls, batch comparisons, or verified treatment histories here. Reagent contamination, index hopping, handling, and genuine biology cannot always be separated.
- **Detection limit:** no spike-ins were performed. A zero count is conditional on this panel and these filters.

## Final repository and submission checks

Keep `results/measurements.tsv`, a completed detection table, `summary.md`, selected final plots, software/database records, and the actual script versions. Add the true conversation/agent transcript if required; no agent execution log should be invented. The assistant drafted scripts, reviewed uploaded output tables, and generated corrected figures locally; the student ran the HPC workflow.

Review what Git would add before committing:

```bash
git status --short --ignored
git add .gitignore environment.yml README.md summary.md scripts metadata results
git diff --cached --stat
```

Inspect the staged file list before committing large or machine-specific outputs. `.gitignore` does not erase previously committed data. Publishing a repository and submitting the assignment are separate actions not performed here.

## Documentation sources

Methods and results above were checked against the provided scripts, `measurements.tsv`, original blobplot tables, raw local BLAST files, and the corrected review tables. Software references:

- Kraken2 manual: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
- rnaSPAdes documentation: https://ablab.github.io/spades/rna.html
- SAMtools depth manual: https://www.htslib.org/doc/samtools-depth.html

The current documentation may describe newer software than the original analysis. Preserve the actual environment export to identify the versions used.
