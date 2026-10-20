#!/usr/bin/env bash
# Correct existing blobplots using existing step 06b LOCAL BLAST output.
# No downloads, BLAST, Kraken, assembly, or read mapping are run.
# Bash + embedded base R. No Python or additional R packages required.
# Submit: /scripts/run_script 10_correct_blobplots.sh
set -Eeuo pipefail
trap 'printf "ERROR at line %s: %s\n" "$LINENO" "$BASH_COMMAND" >&2' ERR

# ---- Project and existing results ----
WD="${WD:-/scratch/aubaxp004/fusarium_screen_v2}"
INPUT_DIR="${INPUT_DIR:-$WD/blobplots_local}"
# Usually leave BLAST_DIR empty: use the source recorded by step 09.
BLAST_DIR="${BLAST_DIR:-}"
MIN_LENGTH=300
MAPQ_LABEL=20                    # Describes EXISTING depth; does not refilter.

# Activate the same software environment inside the batch job.
[[ -r "$HOME/miniconda3/etc/profile.d/conda.sh" ]] || {
    echo 'Missing ~/miniconda3/etc/profile.d/conda.sh' >&2; exit 1;
}
set +u
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate contamination
set -u
export LC_ALL=C
command -v Rscript >/dev/null || {
    echo 'Missing Rscript in contamination. Install r-base before resubmitting.' >&2
    exit 1
}

# Reuse exactly the BLAST directory recorded when the original plots were made.
if [[ -z "$BLAST_DIR" ]]; then
    POINTER="$INPUT_DIR/local_blast_source.txt"
    if [[ ! -s "$POINTER" ]]; then
        POINTER="$WD/validation/all_candidates/latest_selection.txt"
    fi
    [[ -s "$POINTER" ]] || {
        echo 'Missing BLAST source pointer. Set BLAST_DIR to your step 06b result directory.' >&2
        exit 1
    }
    BLAST_DIR=$(head -n 1 "$POINTER")
    BLAST_DIR=${BLAST_DIR%$'\r'}
fi
[[ -s "$BLAST_DIR/reference_sequences.tsv" ]] || {
    echo "Missing: $BLAST_DIR/reference_sequences.tsv" >&2; exit 1;
}
for SAMPLE in SRR27384773_DNA SRR1185280_RNA; do
    [[ -s "$INPUT_DIR/${SAMPLE}_blobplot.tsv" ]] || {
        echo "Missing: $INPUT_DIR/${SAMPLE}_blobplot.tsv" >&2; exit 1;
    }
    RUN=${SAMPLE%_*}
    [[ -f "$BLAST_DIR/$RUN.raw.tsv" ]] || {
        echo "Missing: $BLAST_DIR/$RUN.raw.tsv" >&2; exit 1;
    }
done

# Every invocation gets its own directory. Original results are untouched.
ROOT="$WD/blobplots_corrected"
mkdir -p "$ROOT"
OUT=$(mktemp -d "$ROOT/run.$(date -u +%Y%m%dT%H%M%SZ).XXXXXX")
exec > >(tee "$OUT/correction.log") 2>&1
printf 'Input tables: %s\nLocal BLAST: %s\nOutput: %s\n' "$INPUT_DIR" "$BLAST_DIR" "$OUT"
printf '%s\n' "$INPUT_DIR" > "$OUT/input_table_source.txt"
printf '%s\n' "$BLAST_DIR" > "$OUT/local_blast_source.txt"
sha256sum "$INPUT_DIR/SRR27384773_DNA_blobplot.tsv" \
    "$INPUT_DIR/SRR1185280_RNA_blobplot.tsv" \
    "$BLAST_DIR/SRR27384773.raw.tsv" "$BLAST_DIR/SRR1185280.raw.tsv" \
    "$BLAST_DIR/reference_sequences.tsv" > "$OUT/input_sha256.txt"

Rscript --vanilla - "$INPUT_DIR" "$BLAST_DIR" "$OUT" "$MIN_LENGTH" "$MAPQ_LABEL" <<'RSCRIPT'
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]; blastdir <- args[2]; out <- args[3]
min_length <- as.integer(args[4]); mapq <- args[5]
read_tsv <- function(p, header=TRUE) read.delim(p, header=header,
    sep="\t", quote="", comment.char="", stringsAsFactors=FALSE,
    check.names=FALSE)
write_tsv <- function(x,p) write.table(x,p,sep="\t",row.names=FALSE,
    quote=FALSE,na="NA")
need <- function(x, cols, label) {
    absent <- setdiff(cols,names(x))
    if(length(absent)) stop(label,": missing ",paste(absent,collapse=", "))
}
ref <- read_tsv(file.path(blastdir,"reference_sequences.tsv"))
need(ref,c("sequence_id","group","role"),"Reference key")
if(anyDuplicated(ref$sequence_id)) stop("Duplicate reference sequence IDs")

# Inclusive query intervals; overlap is counted only once.
union_bp <- function(a,b) {
    starts <- pmin(a,b); ends <- pmax(a,b)
    ix <- order(starts,ends); starts <- starts[ix]; ends <- ends[ix]
    total <- 0; last_end <- 0
    for(j in seq_along(starts)) {
        total <- total + max(0, ends[j] - max(last_end,starts[j]-1))
        last_end <- max(last_end,ends[j])
    }
    total
}

# A continuous axis: linear for 0..1, log10 above 1; zero is retained.
y_transform <- function(x) ifelse(x<=1,x,1+log10(pmax(x,1)))

# Separate plotting and legend panels keep the legend within the canvas.
draw_plot <- function(d,run,lib) {
    valid <- is.finite(d$gc_percent) & is.finite(d$mean_depth)
    z <- d[valid & d$length_bp>=min_length, ,drop=FALSE]
    palette <- c("Fusarium host reference"="#26836b", "FgV1"="#c73377",
                 "Homo_sapiens"="#347bb5", "No qualifying hit"="#999999")
    extra <- setdiff(unique(na.omit(z$best_reference_group)),names(palette))
    if(length(extra)) palette <- c(palette,setNames(hcl.colors(length(extra),"Dark 3"),extra))
    group <- z$best_reference_group
    group[is.na(group)] <- "No qualifying hit"
    ymax <- if(nrow(z)) max(1,max(y_transform(z$mean_depth))*1.05) else 1
    maxdepth <- if(nrow(z)) max(z$mean_depth) else 1
    ticks <- c(0,.1,1)
    if(maxdepth>=10) ticks <- c(ticks,10^seq_len(floor(log10(maxdepth))))
    layout(matrix(c(1,2),1,2),widths=c(3.1,1.35))
    par(oma=c(5,0,0,0),mar=c(4.5,5,4,1),xpd=FALSE)
    plot(NA,xlim=c(0,100),ylim=c(0,ymax),yaxt="n",
         xlab="GC (%)",ylab="Mean mapped depth (x; linear 0-1, log above 1)",
         main=paste0(run," | ",lib,"\nAssembly sequences >= ",min_length," bp"))
    abline(h=y_transform(ticks),v=seq(0,100,20),col="#eeeeee")
    axis(2,at=y_transform(ticks),labels=format(ticks,scientific=FALSE,trim=TRUE),las=1,cex.axis=.8)
    # Draw no-hit first, followed by host, then other matches for visibility.
    plot_order <- c("No qualifying hit","Fusarium host reference",
                    setdiff(unique(group),c("No qualifying hit","Fusarium host reference")))
    for(g in plot_order) {
        ix <- which(group==g)
        if(!length(ix)) next
        ix <- ix[order(z$length_bp[ix],decreasing=TRUE)]
        color <- adjustcolor(palette[g],alpha.f=.55)
        # Diameter scaling; minimum size makes short contigs visible.
        # Point area scales approximately with 12+100*sqrt(length/10000).
        cex <- sqrt((12+100*sqrt(z$length_bp[ix]/10000))/36)
        filled <- z$support_status[ix] %in% c("Broad support","No qualifying hit")
        bg <- rep(NA_character_,length(ix)); bg[filled] <- color
        symbol <- ifelse(z$support_status[ix]=="Competing groups",24,21)
        points(z$gc_percent[ix],y_transform(z$mean_depth[ix]),pch=symbol,
               cex=cex,col=color,bg=bg,lwd=.7)
    }
    par(mar=c(1,0,4,1)); plot.new()
    present <- names(palette)[names(palette) %in% group]
    legend("topleft",legend=c(present,"Broad: >=80% union","Partial: <80% union","Competing groups"),
        pch=c(rep(21,length(present)),21,21,24),
        pt.bg=c(unname(palette[present]),"#333333",NA,NA),
        col=c(unname(palette[present]),rep("#333333",3)),
        pt.cex=1.2,cex=.8,bty="n",title="Provisional reference support")
    caption <- paste0("Existing primary-mapping depth; MAPQ >= ",mapq,". No remapping.\n",
        "Minimum point size applied; area is not proportional to sequence length.\n",
        "Gray = no qualifying hit. Reference colors do not establish biological origin.")
    if(lib=="RNA") caption <- paste0(caption,"\nRNA ambiguity among redundant transcripts can reduce retained depth.")
    mtext(caption,side=1,outer=TRUE,line=1,cex=.78,adj=.05)
}

summaries <- list()
for(k in 1:2) {
    run <- c("SRR27384773","SRR1185280")[k]
    lib <- c("DNA","RNA")[k]
    d <- read_tsv(file.path(input,paste0(run,"_",lib,"_blobplot.tsv")))
    need(d,c("contig","length_bp","gc_percent","mean_depth","best_reference_group",
             "reference_group"),run)
    if(anyDuplicated(d$contig) || anyNA(d$contig)) stop("Invalid contig IDs: ",run)
    if(any(!is.finite(d$length_bp) | d$length_bp<=0)) stop("Invalid lengths: ",run)
    if(any(d$mean_depth<0,na.rm=TRUE)) stop("Negative depth: ",run)
    # All old columns stay unchanged for an audit trail.
    d$query_union_bp <- 0
    d$query_union_percent <- 0
    d$support_status <- "No qualifying hit"
    d$reported_hsp_cap_reached <- FALSE
    d$best_group_hsp_cap_reached <- FALSE
    raw <- file.path(blastdir,paste0(run,".raw.tsv"))
    cap_pairs <- 0L; cap_queries <- character()
    if(file.info(raw)$size>0) {
        h <- read_tsv(raw,FALSE)
        if(ncol(h)!=13) stop("Expected 13 BLAST columns: ",raw)
        names(h) <- c("query","subject","identity","aligned_bp","qlen","slen",
                      "qstart","qend","sstart","send","evalue","bitscore","qcovhsp")
        qi <- match(h$query,d$contig); ri <- match(h$subject,ref$sequence_id)
        if(anyNA(qi) || any(h$qlen!=d$length_bp[qi]))
            stop("BLAST query IDs/lengths do not match table: ",run)
        if(anyNA(ri)) stop("BLAST subjects missing from reference key: ",run)
        if(anyNA(h[,3:13]) || any(pmin(h$qstart,h$qend)<1) ||
           any(pmax(h$qstart,h$qend)>h$qlen)) stop("Invalid BLAST values: ",run)
        h$group <- ref$group[ri]
        h$group[ref$role[ri]=="host"] <- "Fusarium host reference"
        # Cap assessment uses raw HSP counts, BEFORE filtering.
        pair <- paste(h$query,h$subject,sep="\034")
        counts <- table(pair)
        h$at_cap <- as.integer(counts[pair])>=5
        cap_pairs <- sum(counts>=5)
        cap_queries <- unique(h$query[h$at_cap])
        d$reported_hsp_cap_reached <- d$contig %in% cap_queries
        h <- h[h$identity>=80 & h$aligned_bp>=300 & h$evalue<=1e-20,,drop=FALSE]
        chunks <- split(seq_len(nrow(h)),h$query)
        for(q in names(chunks)) {
            ix <- match(q,d$contig); x <- h[chunks[[q]],,drop=FALSE]
            scores <- tapply(x$bitscore,x$group,max)
            near <- names(scores)[scores>=.95*max(scores)]
            # Preserve the prior best group, including its original tie handling.
            best <- d$best_reference_group[ix]
            if(is.na(best) || !(best %in% names(scores)) || scores[best]!=max(scores))
                stop("Best group conflicts with original table for ",q,
                     ". Check that the BLAST source matches step 09.")
            bx <- x[x$group==best,,drop=FALSE]
            bp <- union_bp(bx$qstart,bx$qend)
            d$query_union_bp[ix] <- bp
            d$query_union_percent[ix] <- 100*bp/d$length_bp[ix]
            d$best_group_hsp_cap_reached[ix] <- any(bx$at_cap)
            d$support_status[ix] <- if(length(near)>1) "Competing groups" else
                if(d$query_union_percent[ix]>=80) "Broad support" else "Partial support"
        }
    }
    if(any(!is.na(d$best_reference_group) & d$support_status=="No qualifying hit"))
        stop("Original best groups have no supporting BLAST rows: ",run)
    prefix <- file.path(out,paste0(run,"_",lib))
    write_tsv(d,paste0(prefix,"_reviewed.tsv"))
    pdf(paste0(prefix,"_revised_blobplot.pdf"),width=12,height=8)
    draw_plot(d,run,lib); dev.off()
    if(capabilities("cairo")) {
        png(paste0(prefix,"_revised_blobplot.png"),width=2400,height=1600,res=200,type="cairo")
        draw_plot(d,run,lib); dev.off()
        svg(paste0(prefix,"_revised_blobplot.svg"),width=12,height=8)
        draw_plot(d,run,lib); dev.off()
    } else warning("Cairo unavailable: saved PDF but skipped PNG/SVG for ",run)
    summaries[[k]] <- data.frame(run=run,library=lib,sequences=nrow(d),
        broad=sum(d$support_status=="Broad support"),
        partial=sum(d$support_status=="Partial support"),
        competing=sum(d$support_status=="Competing groups"),
        no_hit=sum(d$support_status=="No qualifying hit"),
        partial_to_broad=sum(d$reference_group=="Partial reference match" &
                            d$query_union_percent>=80,na.rm=TRUE),
        pairs_at_hsp_cap=cap_pairs,queries_at_hsp_cap=length(cap_queries),
        shorter_than_plot_cutoff=sum(d$length_bp<min_length))
    print(summaries[[k]],row.names=FALSE)
}
write_tsv(do.call(rbind,summaries),file.path(out,"review_summary.tsv"))
capture.output(sessionInfo(),file=file.path(out,"R_session.txt"))
RSCRIPT

cat > "$OUT/INTERPRETATION.txt" <<'NOTES'
Custom GC-versus-depth plots, not a BlobToolKit dataset.
The new support fields use the union of qualifying query alignment intervals.
Filters: identity >=80%, alignment length >=300 bp, E-value <=1e-20.
Overlapping query bases count once. Both host reference groups are combined.
Broad support means >=80% query union; it is not a validated taxonomic cutoff.
Intervals may span multiple references within a group; they do not prove a
continuous alignment, one organism of origin, or absence of a chimera.
All old table columns remain unchanged. Use best_reference_group plus the
NEW support_status for these plots. Old reference_group/plot_group are audit
columns and still describe the original single-HSP labels.
At least five raw HSPs for a query/subject pair flags possible reporting-cap
truncation; it does NOT prove additional alignments were omitted. A targeted
uncapped local search is needed to investigate affected partial matches.
No alignment search or read mapping was rerun. Depth, breadth, and read counts
are unchanged. RNA transcript ambiguity under MAPQ filtering is unresolved.
Colors identify provisional best matches in a restricted reference panel.
Gray means no qualifying hit, not contamination. Sequences shorter than the
300-bp alignment cutoff cannot meet it; they remain in the tables but are
excluded from the main figures. Human-like matches require specificity checks.
No new significance bins, species identifications, or source claims are made.
This Bash/base-R implementation reproduces the corrected interval calculation
and plot conventions. Rendering differs slightly from the supplied previews.
NOTES
# Include this script for reproducibility where the batch wrapper exposes it.
if [[ -f "${BASH_SOURCE[0]}" ]]; then
    cp -- "${BASH_SOURCE[0]}" "$OUT/10_correct_blobplots.sh"
fi
# Package completed deliverables; omit the live log from the archive.
(
    cd "$OUT"
    tar -czf corrected_blobplots.tar.gz \
        SRR27384773_DNA_* SRR1185280_RNA_* review_summary.tsv \
        INTERPRETATION.txt R_session.txt input_sha256.txt \
        input_table_source.txt local_blast_source.txt
)
printf '%s\n' "$OUT" > "$ROOT/latest_corrected.txt.part"
mv "$ROOT/latest_corrected.txt.part" "$ROOT/latest_corrected.txt"
printf '\nFinished: %s\nArchive: %s/corrected_blobplots.tar.gz\n' "$OUT" "$OUT"
