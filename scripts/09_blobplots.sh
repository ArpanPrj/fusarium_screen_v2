#!/usr/bin/env bash
# Custom blobplots for DNA contigs and RNA transcripts (not a BlobToolKit dataset).
# Plan: preflight -> GC/length -> classify assemblies -> map own reads back
# -> per-sequence coverage -> plot and export tables. No remote services.
# Tool choices: Kraken2 reuses the shared taxonomy DB; Bowtie2 aligns short reads
# to assembled sequences (RNA transcripts have no genomic introns);
# samtools measures alignment coverage; base R draws plots without packages.
# This Bash entry point contains its plotting R code; no Python is used.
set -Eeuo pipefail
WD="/scratch/aubaxp004/fusarium_screen_v2"
trap 'echo "ERROR at line $LINENO. See $WD/blobplots for logs." >&2' ERR
source ~/miniconda3/etc/profile.d/conda.sh
conda activate contamination
export LC_ALL=C
THREADS=8
MAPQ=20
CONFIDENCE=0.1
DB="/datasets/unzipped/kraken2_full"
ROOT="$WD/blobplots"
mkdir -p "$ROOT"
for TOOL in bowtie2 bowtie2-build samtools kraken2 Rscript awk sha256sum stat flock; do
    command -v "$TOOL" >/dev/null || {
        echo "Missing software: $TOOL in contamination environment." >&2
        echo "If Rscript is missing, install r-base before submitting again." >&2
        exit 1
    }
done
# PDF always works with base R. PNG/SVG are also saved when Cairo is available.
Rscript -e 'cat("R plotting ready; Cairo:", capabilities("cairo"), "\n")'
for FILE in hash.k2d opts.k2d taxo.k2d; do
    [[ -r "$DB/$FILE" ]] || { echo "Missing Kraken DB file: $DB/$FILE" >&2; exit 1; }
done
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
        stat -c '%n %s %Y' "$R1" "$R2" "$DB/hash.k2d" "$DB/opts.k2d" "$DB/taxo.k2d"
        printf 'blob-v1 MAPQ=%s confidence=%s\n' "$MAPQ" "$CONFIDENCE"
        bowtie2 --version; samtools --version; kraken2 --version
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

    if [[ ! -f "$OUT/kraken.complete" ]]; then
        echo "Classifying assembled sequences against shared Kraken DB..."
        # --memory-mapping avoids eagerly loading the entire shared DB into RAM;
        # it may still need substantial memory/I/O. Use a compute-node job.
        kraken2 --db "$DB" --memory-mapping --threads "$THREADS" \
            --confidence "$CONFIDENCE" --report "$OUT/kraken.report.tsv.part" \
            --output "$OUT/kraken.classifications.tsv.part" "$ASSEMBLY" \
            > "$OUT/kraken.stdout.log" 2> "$OUT/kraken.log"
        mv "$OUT/kraken.report.tsv.part" "$OUT/kraken.report.tsv"
        mv "$OUT/kraken.classifications.tsv.part" "$OUT/kraken.classifications.tsv"
        touch "$OUT/kraken.complete"
    fi
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
    samtools coverage --ff 0 -q 0 -Q 0 "$OUT/assembly.filtered.bam" > "$OUT/coverage.tsv"
done
mv "$MANIFEST" "$ROOT/current_samples.tsv"

# Base R only. Full per-contig taxonomy is saved even when plot colors combine
# rare phyla. Colors represent Kraken evidence, not confirmed contamination.
Rscript --vanilla - "$ROOT" "$MAPQ" "$CONFIDENCE" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]; mapq <- args[2]; confidence <- args[3]
readtab <- function(path, header=TRUE) read.delim(path, header=header,
    sep="\t", quote="", comment.char="", stringsAsFactors=FALSE,
    check.names=FALSE, strip.white=FALSE)
manifest <- readtab(file.path(root,"current_samples.tsv"))
tables <- list()
for(i in seq_len(nrow(manifest))) {
    out <- manifest$directory[i]
    g <- readtab(file.path(out,"gc.tsv"))
    if(anyDuplicated(g$contig)) stop("Duplicate assembly IDs")
    k <- readtab(file.path(out,"kraken.classifications.tsv"),FALSE)
    if(ncol(k)<5) stop("Unexpected Kraken output")
    names(k)[1:5] <- c("status","contig","taxid","query_length","matches")
    if(anyDuplicated(k$contig)) stop("Duplicate Kraken IDs")
    report <- readtab(file.path(out,"kraken.report.tsv"),FALSE)
    if(ncol(report)!=6) stop("Expected standard six-column Kraken report")
    names(report) <- c("percent","clade","direct","rank","taxid","name")
    # Kraken standard reports are depth-first taxonomy trees. Infer ancestry
    # from leading spaces, preserving unclassified and above-phylum calls.
    stack <- integer(); broad <- character(nrow(report))
    for(j in seq_len(nrow(report))) {
        nm <- report$name[j]
        depth <- nchar(nm)-nchar(sub("^ +","",nm))
        olddepth <- if(length(stack)) nchar(report$name[stack])-nchar(sub("^ +","",report$name[stack])) else numeric()
        stack <- c(stack[olddepth < depth],j)
        phylum <- stack[report$rank[stack]=="P"]
        domains <- stack[trimws(report$name[stack]) %in% c("Bacteria","Archaea","Eukaryota","Viruses")]
        broad[j] <- if(report$taxid[j]==0) "Unclassified" else if(length(phylum))
            trimws(report$name[tail(phylum,1)]) else if(length(domains))
            paste0(trimws(report$name[tail(domains,1)])," (no phylum)") else "Above phylum / unresolved"
    }
    ki <- match(g$contig,k$contig)
    if(anyNA(ki)) stop("Some assembly IDs missing from Kraken output")
    ri <- match(k$taxid[ki],report$taxid)
    g$kraken_status <- k$status[ki]
    g$taxid <- k$taxid[ki]
    g$assigned_taxon <- trimws(report$name[ri])
    g$taxonomy_group <- broad[ri]
    g$taxonomy_group[is.na(g$taxonomy_group)] <- "Above phylum / unresolved"
    g$taxonomy_group[g$kraken_status=="U"] <- "Unclassified"
    g$assigned_taxon[g$kraken_status=="U"] <- "Unclassified"
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
bp <- sort(tapply(all$length_bp,all$taxonomy_group,sum),decreasing=TRUE)
named <- setdiff(names(bp),c("Unclassified","Above phylum / unresolved"))
keep <- head(named,10)
groups <- c(keep,"Other assigned groups","Above phylum / unresolved","Unclassified")
cols <- setNames(c(grDevices::hcl.colors(max(length(keep),1),"Dark 3")[seq_along(keep)],
                  "#A6761D","#555555","#BBBBBB"),groups)
all$plot_group <- ifelse(all$taxonomy_group %in% groups,all$taxonomy_group,"Other assigned groups")
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
           title="Kraken taxonomy (mostly phylum)")
    mtext(sprintf("%s sequences; %s with zero retained coverage; %s GC-undefined omitted",
          nrow(d),sum(d$mean_depth==0),sum(!good)),side=3,line=.5,cex=.72)
    mtext(paste0("Point area proportional to length within panel; MAPQ >= ",mapq,
          "; Kraken confidence ",confidence,"; gray = unclassified"),side=1,line=3.6,cex=.67)
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
Each point is an assembled DNA contig or RNA transcript, not an organism.
GC uses A/C/G/T bases only. Depth includes zero-covered positions and counts
overlapping mates separately; primary mapped reads with MAPQ>=20 are retained.
The depth axis uses log10(depth+1), with ticks labeled in original depth units.
Bubble area is proportional to sequence length within each panel; maximum
bubble size is rescaled separately for DNA/RNA. All sequences are in the TSV;
GC-undefined sequences are omitted only from the plots.
RNA coverage reflects expression and is affected by isoforms and multimapping.
Low/zero coverage can reflect MAPQ filtering rather than biological absence.
Kraken colors are preliminary taxonomic assignments, not significance bins.
The shared DB may lack adequate fungal representation: gray is not contamination.
Rare taxonomy groups are combined only for plot colors; TSV preserves labels.
This custom plot is not a BlobToolKit dataset. No Bracken or remote BLAST used.
No plots from reference-genome measurements were substituted for contig coverage.
EOF
(
    cd "$ROOT"
    tar -czf blobplot_results.tar.gz \
        SRR27384773_DNA_blobplot.* SRR1185280_RNA_blobplot.* \
        all_contigs_blobplot.tsv INTERPRETATION.txt R_session.txt current_samples.tsv
)
echo "Finished. Plots/tables: $ROOT"
echo "Upload $ROOT/blobplot_results.tar.gz"
