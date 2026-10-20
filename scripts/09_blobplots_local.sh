#!/usr/bin/env bash
# Local-reference blobplots: no Kraken, remote BLAST or database download.
# Reuse completed step 06b alignments; compute GC/length from assemblies;
# map each readset to its OWN assembly and plot per-contig coverage.
# Bowtie2 is appropriate for RNA here because the targets are transcripts.
# Base R plotting is embedded below; no Python or extra R packages are needed.
# Colors are provisional best-reference groups, not verified species identities.
set -Eeuo pipefail
WD="/scratch/aubaxp004/fusarium_screen_v2"
trap 'echo "ERROR at line $LINENO. See $WD/blobplots_local for logs." >&2' ERR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
export LC_ALL=C
THREADS=8
MAPQ=20
ROOT="$WD/blobplots_local"
mkdir -p "$ROOT"
for TOOL in bowtie2 bowtie2-build samtools Rscript awk sha256sum stat flock; do
    command -v "$TOOL" >/dev/null || {
        echo "Missing software: $TOOL in contamination environment." >&2
        echo "If Rscript is missing, install r-base before submitting again." >&2
        exit 1
    }
done
# PDF always works with base R. PNG/SVG are also saved when Cairo is available.
Rscript -e 'cat("R plotting ready; Cairo:", capabilities("cairo"), "\n")'
# This pointer is written only after step 06b completes both local searches.
POINTER="$WD/validation/all_candidates/latest_selection.txt"
[[ -s "$POINTER" ]] || { echo "Missing completed local-BLAST pointer: $POINTER" >&2; exit 1; }
IFS= read -r BLASTDIR < "$POINTER"
[[ -s "$BLASTDIR/reference_sequences.tsv" ]] || { echo "Missing reference key in $BLASTDIR" >&2; exit 1; }
for RUN in SRR27384773 SRR1185280; do
    # Empty local BLAST output is valid (no hits); absent output is not.
    [[ -f "$BLASTDIR/$RUN.raw.tsv" ]] || { echo "Missing local BLAST output for $RUN" >&2; exit 1; }
done
printf '%s\n' "$BLASTDIR" > "$ROOT/local_blast_source.txt"
exec 9> "$ROOT/job.lock"
flock -n 9 || { echo "Another blobplot job is active" >&2; exit 1; }
MANIFEST="$ROOT/current_samples.tsv.part"
printf 'run\tlibrary\tdirectory\n' > "$MANIFEST"
for RUN in SRR27384773 SRR1185280; do
    if [[ "$RUN" == SRR27384773 ]]; then
        TYPE=DNA; ASSEMBLY="$WD/assembly/dna_megahit/final.contigs.fa"
    else
        TYPE=RNA; ASSEMBLY="$WD/assembly/rna_rnaspades/transcripts.fasta"
    fi
    R1="$WD/raw/${RUN}_1.fastq.gz"
    R2="$WD/raw/${RUN}_2.fastq.gz"
    for FILE in "$ASSEMBLY" "$R1" "$R2"; do
        [[ -s "$FILE" ]] || { echo "Missing/empty: $FILE" >&2; exit 1; }
    done
    # Resume completed stages for identical assemblies, settings and input file
    # metadata. Raw-read identity is checked by path/size/mtime, not full hashing.
    KEY=$({
        sha256sum "$ASSEMBLY"
        stat -c '%n %s %Y' "$R1" "$R2"
        printf 'local-blob-v1 MAPQ=%s\n' "$MAPQ"
        bowtie2 --version; samtools --version
    } | sha256sum | awk '{print $1}')
    OUT="$ROOT/$RUN.$KEY"
    mkdir -p "$OUT/index"
    printf '%s\t%s\t%s\n' "$RUN" "$TYPE" "$OUT" >> "$MANIFEST"
    echo "Working on $RUN ($TYPE): $OUT"

    # GC denominator is A+C+G+T only; retain ambiguous-base count separately.
    awk 'BEGIN{OFS="\t";print "contig","length_bp","gc_percent","ambiguous_bp"}
        function emit(){if(id!="") print id,len,(acgt ? 100*gc/acgt : "NA"),len-acgt}
        /^>/{emit();id=substr($0,2);sub(/[ \t].*/,"",id);len=gc=acgt=0;next}
        {s=toupper($0);gsub(/[[:space:]]/,"",s);len+=length(s);
         t=s;gc+=gsub(/[GC]/,"",t);t=s;acgt+=gsub(/[ACGT]/,"",t)}
        END{emit()}' "$ASSEMBLY" > "$OUT/gc.tsv"

    if [[ ! -f "$OUT/index.complete" ]]; then
        bowtie2-build --threads "$THREADS" "$ASSEMBLY" "$OUT/index/assembly" \
            > "$OUT/index.log" 2>&1
        touch "$OUT/index.complete"
    fi
    if [[ ! -f "$OUT/mapping.complete" ]]; then
        echo "Mapping $RUN reads back to its own assembly..."
        # Align transcripts directly: a splice-aware mapper is not needed here.
        # Primary mapped reads only, MAPQ>=20. Retain duplicates and overlapping
        # mates. RNA isoforms/repeats may lose ambiguous reads under this filter.
        bowtie2 --very-sensitive -p "$THREADS" -x "$OUT/index/assembly" \
            -1 "$R1" -2 "$R2" 2> "$OUT/mapping.log" |
            samtools view -u -q "$MAPQ" -F 2308 - |
            samtools sort -@ 2 -m 1G -o "$OUT/assembly.filtered.bam.part" -
        mv "$OUT/assembly.filtered.bam.part" "$OUT/assembly.filtered.bam"
        samtools quickcheck "$OUT/assembly.filtered.bam"
        samtools index "$OUT/assembly.filtered.bam"
        touch "$OUT/mapping.complete"
    fi
    # No further flag/base-quality filtering. Use full sequence lengths as the
    # coverage denominator. The table join below retains zero-coverage contigs.
    samtools coverage --ff 0 -q 0 -Q 0 -d 0 "$OUT/assembly.filtered.bam" > "$OUT/coverage.tsv"
done
mv "$MANIFEST" "$ROOT/current_samples.tsv"

# A 95%-of-best bit-score rule flags near ties between reference groups.
# Ranking uses the best individual HSP per group, not summed overlap counts.
# Full sequence labels require >=80% query coverage by that HSP; shorter
# matches are labeled partial/unresolved. These are display rules, not bins.
Rscript --vanilla - "$ROOT" "$MAPQ" "$BLASTDIR" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]; mapq <- args[2]; blastdir <- args[3]
readtab <- function(path, header=TRUE) read.delim(path, header=header,
    sep="\t", quote="", comment.char="", stringsAsFactors=FALSE,
    check.names=FALSE, strip.white=FALSE)
manifest <- readtab(file.path(root,"current_samples.tsv"))
tables <- list()
for(i in seq_len(nrow(manifest))) {
    out <- manifest$directory[i]
    g <- readtab(file.path(out,"gc.tsv"))
    if(anyDuplicated(g$contig)) stop("Duplicate assembly IDs")
    refmap <- readtab(file.path(blastdir,"reference_sequences.tsv"))
    g$reference_group <- "No qualifying reference hit"
    g$best_reference_group <- NA_character_
    g$near_best_groups <- NA_character_
    g$best_hsp_bitscore <- NA_real_
    g$best_hsp_identity <- NA_real_
    g$best_hsp_query_coverage_percent <- NA_real_
    g$best_hsp_alignment_bp <- NA_real_
    hitfile <- file.path(blastdir,paste0(manifest$run[i],".raw.tsv"))
    if(file.info(hitfile)$size>0) {
        h <- readtab(hitfile,FALSE)
        if(ncol(h)!=13) stop("Unexpected local BLAST format: ",hitfile)
        names(h) <- c("query","subject","identity","aligned_bp","qlen","slen",
                       "qstart","qend","sstart","send","evalue","bitscore","qcovhsp")
        # Verify results belong to the current assembly, not an unrelated file.
        qi <- match(h$query,g$contig)
        if(anyNA(qi) || any(h$qlen != g$length_bp[qi]))
            stop("Local BLAST query IDs/lengths do not match current assembly")
        si <- match(h$subject,refmap$sequence_id)
        if(anyNA(si)) stop("Some BLAST subject IDs are missing from reference key")
        h$group <- refmap$group[si]
        h$group[refmap$role[si]=="host"] <- "Fusarium host reference"
        h <- h[h$identity>=80 & h$aligned_bp>=300 & h$evalue<=1e-20,]
        if(nrow(h)) {
            h <- h[order(-h$bitscore,-h$aligned_bp),]
            h <- h[!duplicated(h[c("query","group")]),]
            for(rows in split(seq_len(nrow(h)),h$query)) {
                x <- h[rows,]; x <- x[order(-x$bitscore,-x$aligned_bp),]
                ix <- match(x$query[1],g$contig)
                near <- unique(x$group[x$bitscore>=0.95*x$bitscore[1]])
                g$best_reference_group[ix] <- x$group[1]
                g$near_best_groups[ix] <- paste(near,collapse="; ")
                g$best_hsp_bitscore[ix] <- x$bitscore[1]
                g$best_hsp_identity[ix] <- x$identity[1]
                g$best_hsp_query_coverage_percent[ix] <- x$qcovhsp[1]
                g$best_hsp_alignment_bp[ix] <- x$aligned_bp[1]
                g$reference_group[ix] <- if(length(near)>1) "Ambiguous reference match" else
                    if(x$qcovhsp[1]<80) "Partial reference match" else x$group[1]
            }
        }
    }
    cov <- readtab(file.path(out,"coverage.tsv"))
    ci <- match(g$contig,cov[[1]])
    g$mapped_reads <- cov$numreads[ci]
    g$covered_bp <- cov$covbases[ci]
    g$breadth_percent <- cov$coverage[ci]
    g$mean_depth <- cov$meandepth[ci]
    for(nm in c("mapped_reads","covered_bp","breadth_percent","mean_depth"))
        g[[nm]][is.na(g[[nm]])] <- 0
    g$run <- manifest$run[i]; g$library <- manifest$library[i]
    tables[[i]] <- g
}
all <- do.call(rbind,tables)
special <- c("Partial reference match","Ambiguous reference match","No qualifying reference hit")
named <- sort(setdiff(unique(all$reference_group),special))
groups <- c(named,special)
cols <- setNames(c(grDevices::hcl.colors(max(length(named),1),"Dark 3")[seq_along(named)],
                  "#B7A9CF","#555555","#BBBBBB"),groups)
all$plot_group <- all$reference_group
write.table(all,file.path(root,"all_contigs_blobplot.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
ymax <- max(1,log10(max(all$mean_depth)+1)*1.06)
draw <- function(d,title) {
    good <- is.finite(d$gc_percent) & is.finite(d$mean_depth)
    z <- d[good,]; z <- z[order(z$length_bp,decreasing=TRUE),]
    par(mar=c(5,5,4.5,14),xpd=FALSE)
    plot(NA,xlim=c(0,100),ylim=c(0,ymax),xaxt="n",yaxt="n",
         xlab="GC content (%)",ylab="Mean read depth (log10(depth + 1) axis)",main=title)
    axis(1); ticks <- 0:ceiling(ymax)
    axis(2,at=ticks,labels=format(10^ticks-1,scientific=FALSE,trim=TRUE),las=1,cex.axis=.8)
    abline(h=ticks,col="#ECECEC",lwd=.6)
    if(nrow(z)) symbols(z$gc_percent,log10(z$mean_depth+1),
        circles=sqrt(z$length_bp),inches=.30,add=TRUE,
        bg=adjustcolor(cols[z$plot_group],alpha.f=.55),fg=NA)
    present <- groups[groups %in% d$plot_group]
    legend("topright",inset=c(-.48,0),legend=present,pch=21,
           pt.bg=cols[present],col=NA,pt.cex=1.4,bty="n",cex=.74,xpd=NA,
           title="Provisional reference match")
    mtext(sprintf("%s sequences; %s with zero retained coverage; %s GC-undefined omitted",
          nrow(d),sum(d$mean_depth==0),sum(!good)),side=3,line=.5,cex=.72)
    mtext(paste0("Point area proportional to length within panel; MAPQ >= ",mapq,
          "; gray = no qualifying reference hit"),side=1,line=3.6,cex=.67)
}
for(i in seq_len(nrow(manifest))) {
    d <- all[all$run==manifest$run[i],]
    prefix <- file.path(root,paste0(manifest$run[i],"_",manifest$library[i],"_blobplot"))
    title <- paste(manifest$run[i],if(manifest$library[i]=="RNA") "RNA transcripts" else "DNA contigs")
    pdf(paste0(prefix,".pdf"),width=12,height=7.5); draw(d,title); dev.off()
    if(capabilities("cairo")) {
        png(paste0(prefix,".png"),width=2400,height=1500,res=200,type="cairo")
        draw(d,title); dev.off()
        svg(paste0(prefix,".svg"),width=12,height=7.5); draw(d,title); dev.off()
    }
    write.table(d,paste0(prefix,".tsv"),sep="\t",row.names=FALSE,quote=FALSE)
}
capture.output(sessionInfo(),file=file.path(root,"R_session.txt"))
RSCRIPT

cat > "$ROOT/INTERPRETATION.txt" <<'EOF'
Custom GC-versus-coverage plots, not a BlobToolKit dataset.
Colors reuse existing step 06b LOCAL BLAST alignments; no new search is run.
Reference labels are restricted to the selected panel, not broad taxonomy.
Best group uses its highest individual HSP (>=80% identity, >=300 aligned bp,
E<=1e-20). Groups within 95% of best bit score are ambiguous. Below 80% query
coverage by the best HSP, color is partial/unresolved, not a full-contig call.
The TSV retains best-group evidence for partial and ambiguous matches too.
These are plotting heuristics, not validated significance thresholds. A best
hit does not prove species, exclude organisms absent from the panel, or rule
out a chimera. RNA exon matches can fall into the partial category.
Fusarium nuclear and mitochondrial references share a host color.
Gray is no qualifying panel hit, NOT evidence of contamination.
GC uses A/C/G/T only. Depth includes zero-covered positions, primary mapped
reads with MAPQ>=20, duplicate reads, and overlapping mates counted separately.
RNA targets are transcripts, so Bowtie2 is used without splice-aware alignment.
RNA coverage reflects expression and is reduced by ambiguity among isoforms.
The axis is log10(depth+1), labeled in original depth units. Zero-depth points
are retained. GC-undefined sequences appear in tables but not plots.
Point area is proportional to sequence length within each panel; maximum
bubble size is rescaled separately for DNA and RNA.
Old step 04 reference-genome coverage is NOT substituted for assembly coverage.
Script resumes its own completed mapping/index stages using assembly checksum,
settings/tool versions and raw-read path/size/mtime. It does not automatically
reuse older Kraken-workflow BAMs, whose provenance is not recorded sufficiently.
EOF
(
    cd "$ROOT"
    tar -czf blobplot_results.tar.gz \
        SRR27384773_DNA_blobplot.* SRR1185280_RNA_blobplot.* \
        all_contigs_blobplot.tsv INTERPRETATION.txt R_session.txt current_samples.tsv local_blast_source.txt
)
echo "Finished. Plots/tables: $ROOT"
echo "Upload $ROOT/blobplot_results.tar.gz"
