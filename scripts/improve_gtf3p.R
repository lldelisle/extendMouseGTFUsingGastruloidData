if (!all(c("rtracklayer", "plyr") %in% installed.packages())) {
  if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
  }
  devtools::install_github("lldelisle/usefulLDfunctions", upgrade = "never")
  library(usefulLDfunctions)
  safelyLoadAPackageInCRANorBioconductor("rtracklayer")
  safelyLoadAPackageInCRANorBioconductor("plyr")
} else {
  library(rtracklayer)
  library(plyr)
}

rm(list = ls())

if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help") {
  cat("Usage: Rscript improve_gtf3p.R pathForInputGtf pathForOutputGtf pathForMergedStringtie\n")
  stop()
}
input.gtf <- commandArgs(TRUE)[1]
output.gtf <- commandArgs(TRUE)[2]
stringtie.gtf <- commandArgs(TRUE)[3]

cat("Loading input gtf file...")
input.gr <- readGFFAsGRanges(input.gtf)
cat("Done.\n")
# I store the CDS:
input.gr.cds <- subset(input.gr, type == "CDS")
# I only work with exons
input.gr <- subset(input.gr, type == "exon")
# Build a dataframe with all exon info:
input.gr.tr.exon <- data.frame(i = 1:length(input.gr), transcript_id = input.gr$transcript_id,
                               strand = strand(input.gr), start = start(input.gr), end = end(input.gr),
                               exon_id = input.gr$exon_id)
# Summarize by transcript_id to keep only the left and right exon
input.gr.tr.exon.first.last <- ddply(input.gr.tr.exon, .(transcript_id), summarize,
                                     left.exon = i[start == min(start)],
                                     right.exon = i[end == max(end)],
                                     strand = unique(strand))
# Adjust the last to corresponds to 3'
input.gr.tr.exon.first.last$last <- input.gr.tr.exon.first.last$right.exon
input.gr.tr.exon.first.last$last[input.gr.tr.exon.first.last$strand == "-"] <- input.gr.tr.exon.first.last$left.exon[input.gr.tr.exon.first.last$strand == "-"]

# Create a GRange with only last exons
input.gr.last.exon <- input.gr[input.gr.tr.exon.first.last$last]

# Load stringtie output
cat("Loading input merged Stringtie gtf file...")
stringtie.gr <- subset(readGFFAsGRanges(stringtie.gtf), type == "exon")
cat("Done.\n")

# Get overlap between stringtie output and last exon
# Get start/end/strand for both
overlap.last <- as.data.frame(findOverlaps(input.gr.last.exon, stringtie.gr))
overlap.last$exon_id <- input.gr.last.exon$exon_id[overlap.last$queryHits]
overlap.last$transcript_id <- input.gr.last.exon$transcript_id[overlap.last$queryHits]
overlap.last$start <- start(input.gr.last.exon)[overlap.last$queryHits]
overlap.last$end <- end(input.gr.last.exon)[overlap.last$queryHits]
overlap.last$strand <- as.character(strand(input.gr.last.exon)[overlap.last$queryHits])
overlap.last$start.stringtie <- start(stringtie.gr)[overlap.last$subjectHits]
overlap.last$end.stringtie <- end(stringtie.gr)[overlap.last$subjectHits]
overlap.last$strand.stringtie <- as.character(strand(stringtie.gr)[overlap.last$subjectHits])

# If we have incoherent strand we remove the overlap:
if (!all(overlap.last$strand.stringtie == "*" | overlap.last$strand == overlap.last$strand.stringtie)) {
  overlaps.to.remove <- which(overlap.last$strand.stringtie != "*" & overlap.last$strand != overlap.last$strand.stringtie)
  cat("Remove", length(overlaps.to.remove), "overlaps between last exon and stringtie found with 2 orientations.\n")
  overlap.last <- overlap.last[-overlaps.to.remove, ]
}

# Select overlap which would extend the 3'
# First if the exon directly extend
overlap.last.extend <- subset(overlap.last, (strand == "+" & end < end.stringtie) | (strand == "-" & start > start.stringtie))
overlap.last.extend$new.start <- overlap.last.extend$start
overlap.last.extend$new.start[overlap.last.extend$strand == "-"] <- overlap.last.extend$start.stringtie[overlap.last.extend$strand == "-"]
overlap.last.extend$new.end <- overlap.last.extend$end
overlap.last.extend$new.end[overlap.last.extend$strand == "+"] <- overlap.last.extend$end.stringtie[overlap.last.extend$strand == "+"]
# Order to only keep one overlap per 3'
overlap.last.extend <- overlap.last.extend[order(overlap.last.extend$new.start, -overlap.last.extend$new.end), ]
overlap.last.extend.uniq <- overlap.last.extend[!duplicated(overlap.last.extend[, c("exon_id", "transcript_id")]), ]

overlap.last.extend.uniq$id <- paste0(overlap.last.extend.uniq$exon_id, "_", overlap.last.extend.uniq$transcript_id) 
rownames(overlap.last.extend.uniq) <- overlap.last.extend.uniq$id
cat("Found", length(unique(overlap.last.extend.uniq$exon_id)), "potential exon extensions. Will now check for collision.\n")

# Try to do the change
input.gr.nochange <- subset(input.gr, !paste0(exon_id, "_", transcript_id) %in% overlap.last.extend.uniq$id)
input.gr.changed <- subset(input.gr, paste0(exon_id, "_", transcript_id) %in% overlap.last.extend.uniq$id)
start(input.gr.changed) <- overlap.last.extend.uniq[paste0(input.gr.changed$exon_id, "_", input.gr.changed$transcript_id), "new.start"]
end(input.gr.changed) <- overlap.last.extend.uniq[paste0(input.gr.changed$exon_id, "_", input.gr.changed$transcript_id), "new.end"]
# Change the exon id
input.gr.changed$exon_id <- paste0(input.gr.changed$exon_id, ".ext")
# Check that identical exon id has identical start and end:
uniq.start.end.exon.id <- unique(data.frame(seqid = as.character(seqnames(input.gr.changed)), start = start(input.gr.changed), end = end(input.gr.changed), exon_id = input.gr.changed$exon_id))
if (anyDuplicated(uniq.start.end.exon.id$exon_id) != 0) {
  stop("I did not expect that")
}

# plot(width(input.gr.changed), width(subset(input.gr, exon_id %in% overlap.last.extend.uniq$exon_id)))
input.gr.combined <- c(input.gr.nochange, input.gr.changed)

# Test the overlap and check the gene_id
overlap.new <- as.data.frame(findOverlaps(input.gr.changed, input.gr.combined))
overlap.new$gene_id <- input.gr.changed$gene_id[overlap.new$queryHits]
overlap.new$gene_id.overlapped <- input.gr.combined$gene_id[overlap.new$subjectHits]

# Get the one with overlap between different gene_id
overlap.new.problem <- subset(overlap.new, gene_id != gene_id.overlapped)
overlap.new.problem$exon_id <- gsub(".ext$", "", input.gr.changed$exon_id[overlap.new.problem$queryHits])
overlap.new.problem$new.width <- width(pintersect(input.gr.changed[overlap.new.problem$queryHits], input.gr.combined[overlap.new.problem$subjectHits]))
overlap.new.problem$subject.start <- start(input.gr.combined)[overlap.new.problem$subjectHits]
overlap.new.problem$subject.end <- end(input.gr.combined)[overlap.new.problem$subjectHits]
overlap.new.problem <- overlap.new.problem[order(overlap.new.problem$new.width, decreasing = T), ]
overlap.new.problem.uniq <- overlap.new.problem[!duplicated(overlap.new.problem[, c("exon_id", "gene_id", "gene_id.overlapped")]), ]
cat("Found", length(unique(overlap.new.problem.uniq$exon_id)), "exon extensions which collide.\n")

# Check the one which were already overlapping in original version
old.exons.problem <- subset(input.gr, exon_id %in% overlap.new.problem$exon_id)
overlap.old <- as.data.frame(findOverlaps(old.exons.problem, input.gr))
overlap.old$gene_id <- old.exons.problem$gene_id[overlap.old$queryHits]
overlap.old$gene_id.overlapped <- input.gr$gene_id[overlap.old$subjectHits]
overlap.old$exon_id <- old.exons.problem$exon_id[overlap.old$queryHits]
overlap.old.problem <- subset(overlap.old, gene_id != gene_id.overlapped)
overlap.old.problem$old.width <- width(pintersect(old.exons.problem[overlap.old.problem$queryHits], input.gr[overlap.old.problem$subjectHits]))
overlap.old.problem <- overlap.old.problem[order(overlap.old.problem$old.width, decreasing = T), ]
overlap.old.problem.uniq <- overlap.old.problem[!duplicated(overlap.old.problem[, c("exon_id", "gene_id", "gene_id.overlapped")]), ]

# Identify the one which increase the size of overlap
compare.prob <- merge(overlap.new.problem.uniq[, c("exon_id", "gene_id", "gene_id.overlapped", "new.width")],
                      overlap.old.problem.uniq[, c("exon_id", "gene_id", "gene_id.overlapped", "old.width")],
                      all = TRUE)
# summary(compare.prob)
compare.prob[is.na(compare.prob)] <- 0
# plot(compare.prob$new.width, compare.prob$old.width)
# abline(a = 0, b = 1)
compare.prob.increased <- subset(compare.prob, new.width > old.width)
exons.increased.prob <- unique(compare.prob.increased$exon_id)
cat("Among them", length(exons.increased.prob), "exon extensions overlap a larger region.\n")

# Adapt extension for these new exons or remove extension
overlap.last.extend.uniq$gene_id <- input.gr.last.exon$gene_id[overlap.last.extend.uniq$queryHits]
overlap.last.extend.uniq.filtered <- subset(overlap.last.extend.uniq, ! exon_id %in% exons.increased.prob)
collisions <- ddply(overlap.new.problem, .(exon_id), summarise,
                    subject.start = min(subject.start),
                    subject.end = max(subject.end))

to.change.fwd <- subset(overlap.last.extend.uniq, exon_id %in% exons.increased.prob & strand == "+")
to.change.fwd <- merge(to.change.fwd, collisions, all.x = T)
# For this case:
# gene1:      ------->
# gene2:                    -------->
# stringtie:  ---------------------->
# or
# gene2:                    -------------->
to.change.fwd.stop.at.collision <- subset(to.change.fwd, subject.start > (end + 1))
to.change.fwd.stop.at.collision$new.end <- to.change.fwd.stop.at.collision$subject.start - 1
# For this case:
# gene1:      ------->                (will be extended by previous step)
# gene2:                    ---->
# stringtie:  ---------------------->
# Or
# gene1:      ------->
# gene2:      ------------->
# But not this case:
# gene1:     ----------------/\--------
# gene2:        ----->
# stringtie: ----------------/\--------
fwd.exons <- subset(input.gr, strand == "+")
overlap.fwd <- as.data.frame(findOverlaps(fwd.exons, stringtie.gr))
overlap.fwd$end <- end(fwd.exons)[overlap.fwd$queryHits]
overlap.fwd$gene_id <- fwd.exons$gene_id[overlap.fwd$queryHits]
prefered.gene.id <- unique(ddply(subset(overlap.fwd, subjectHits %in% to.change.fwd$subjectHits),
                                 .(subjectHits), summarise,
                                 prefered.gene.id = gene_id[end == max(end)]))
dup.sub <- prefered.gene.id$subjectHits[duplicated(prefered.gene.id$subjectHits)]
prefered.gene.id <- subset(prefered.gene.id, !subjectHits %in% dup.sub)
to.change.fwd <- merge(to.change.fwd, prefered.gene.id, all.x = TRUE)
to.change.fwd.favour.last <- subset(to.change.fwd, subject.start <= (end + 1) & gene_id == prefered.gene.id)

to.change.rev <- subset(overlap.last.extend.uniq, exon_id %in% exons.increased.prob & strand == "-")
to.change.rev <- merge(to.change.rev, collisions, all.x = T)
to.change.rev.stop.at.collision <- subset(to.change.rev, subject.end < (start - 1))
to.change.rev.stop.at.collision$new.start <- to.change.rev.stop.at.collision$subject.end + 1
rev.exons <- subset(input.gr, strand == "-")
overlap.rev <- as.data.frame(findOverlaps(rev.exons, stringtie.gr))
overlap.rev$start <- start(rev.exons)[overlap.rev$queryHits]
overlap.rev$gene_id <- rev.exons$gene_id[overlap.rev$queryHits]
prefered.gene.id <- unique(ddply(subset(overlap.rev, subjectHits %in% to.change.rev$subjectHits),
                                 .(subjectHits), summarise,
                                 prefered.gene.id = gene_id[start == min(start)]))
dup.sub <- prefered.gene.id$subjectHits[duplicated(prefered.gene.id$subjectHits)]
prefered.gene.id <- subset(prefered.gene.id, !subjectHits %in% dup.sub)
to.change.rev <- merge(to.change.rev, prefered.gene.id, all.x = TRUE)
to.change.rev.favour.last <- subset(to.change.rev, subject.end >= (start - 1) & gene_id == prefered.gene.id)

overlap.last.extend.uniq.filtered <- rbind(overlap.last.extend.uniq.filtered,
                                           to.change.fwd.stop.at.collision[, colnames(overlap.last.extend.uniq.filtered)],
                                           to.change.fwd.favour.last[, colnames(overlap.last.extend.uniq.filtered)],
                                           to.change.rev.stop.at.collision[, colnames(overlap.last.extend.uniq.filtered)],
                                           to.change.rev.favour.last[, colnames(overlap.last.extend.uniq.filtered)])
rownames(overlap.last.extend.uniq.filtered) <- overlap.last.extend.uniq.filtered$id
cat("Therefore only", length(unique(overlap.last.extend.uniq.filtered$exon_id)), "exons will be extended.\n")

cat(setdiff(exons.increased.prob, overlap.last.extend.uniq.filtered$exon_id), sep = "\n", file = gsub(".gtf$", "_rejected.txt", output.gtf))

# Combine both gtf (unchanged and changed)
input.gr.nochange <- subset(input.gr, !paste0(exon_id, "_", transcript_id) %in% overlap.last.extend.uniq.filtered$id)
input.gr.changed <- subset(input.gr, paste0(exon_id, "_", transcript_id) %in% overlap.last.extend.uniq.filtered$id)
start(input.gr.changed) <- overlap.last.extend.uniq.filtered[paste0(input.gr.changed$exon_id, "_", input.gr.changed$transcript_id), "new.start"]
end(input.gr.changed) <- overlap.last.extend.uniq.filtered[paste0(input.gr.changed$exon_id, "_", input.gr.changed$transcript_id), "new.end"]
input.gr.changed$exon_id <- paste0(input.gr.changed$exon_id, ".ext")
# Check that identical exon id has identical start and end:
uniq.start.end.exon.id <- unique(data.frame(seqid = as.character(seqnames(input.gr.changed)), start = start(input.gr.changed), end = end(input.gr.changed), exon_id = input.gr.changed$exon_id))
if (anyDuplicated(uniq.start.end.exon.id$exon_id) != 0) {
  stop("I did not expect that")
}

input.gr.combined <- c(input.gr.nochange, input.gr.changed)

## Add new exons from stringtie if belong to same transcript but not overlap existing genes
cat("Trying to add new exons in 3'\n")
# Build a dataframe with all exon info:
input.gr.tr.exon <- data.frame(i = 1:length(input.gr.combined), transcript_id = input.gr.combined$transcript_id,
                               strand = strand(input.gr.combined), start = start(input.gr.combined), end = end(input.gr.combined),
                               exon_id = input.gr.combined$exon_id)
# Summarize by transcript_id to keep only the left and right exon
input.gr.tr.exon.first.last <- ddply(input.gr.tr.exon, .(transcript_id), summarize,
                                     left.exon = i[start == min(start)],
                                     right.exon = i[end == max(end)],
                                     strand = unique(strand))
# Adjust the last to corresponds to 3'
input.gr.tr.exon.first.last$last <- input.gr.tr.exon.first.last$right.exon
input.gr.tr.exon.first.last$last[input.gr.tr.exon.first.last$strand == "-"] <- input.gr.tr.exon.first.last$left.exon[input.gr.tr.exon.first.last$strand == "-"]

# Create a GRange with only last exons
input.gr.last.exon <- input.gr.combined[input.gr.tr.exon.first.last$last]

# Get overlap between stringtie output and last exon
# Get start/end/strand for both
overlap.last <- as.data.frame(findOverlaps(input.gr.last.exon, stringtie.gr))
overlap.last$exon_id <- input.gr.last.exon$exon_id[overlap.last$queryHits]
overlap.last$transcript_id <- input.gr.last.exon$transcript_id[overlap.last$queryHits]
overlap.last$start <- start(input.gr.last.exon)[overlap.last$queryHits]
overlap.last$end <- end(input.gr.last.exon)[overlap.last$queryHits]
overlap.last$strand <- as.character(strand(input.gr.last.exon)[overlap.last$queryHits])
overlap.last$start.stringtie <- start(stringtie.gr)[overlap.last$subjectHits]
overlap.last$end.stringtie <- end(stringtie.gr)[overlap.last$subjectHits]
overlap.last$strand.stringtie <- as.character(strand(stringtie.gr)[overlap.last$subjectHits])
overlap.last$transcript.stringtie <- stringtie.gr$transcript_id[overlap.last$subjectHits]

# Get extremities of transcript.stringtie
transcript.stringtie.char <- ddply(subset(as.data.frame(stringtie.gr), transcript_id %in% overlap.last$transcript.stringtie),
                                   .(transcript_id), summarize,
                                   transcript.stringtie.start = min(start),
                                   transcript.stringtie.end = max(end))
colnames(transcript.stringtie.char)[1] <- "transcript.stringtie"

overlap.last <- merge(overlap.last, transcript.stringtie.char)

# Sort to get only the last overlapping exon:
overlap.last <- overlap.last[order(-overlap.last$start.stringtie), ]
overlap.last.uniq.fwd <- overlap.last[!duplicated(overlap.last[, c("transcript_id", "transcript.stringtie")]), ]

# Select last exons which overlap non-last stringtie exons but have same extremity (splicing reason)
non.last.fwd <- subset(overlap.last.uniq.fwd, strand == "+" & end.stringtie != transcript.stringtie.end & end.stringtie == end)
non.last.fwd$exon.nb.stringtie <- as.numeric(stringtie.gr$exon_number[non.last.fwd$subjectHits])
if (nrow(non.last.fwd) > 0) {
  rownames(non.last.fwd) <- with(non.last.fwd,  paste0(transcript_id, "__", transcript.stringtie))
  
  # Will recursively check if the following exon of stringtie does not overlap any exon.
  # If it is the case, the exon will be added and the next exon checked
  transcript.ids.to.check <- unique(non.last.fwd[,c("transcript_id", "transcript.stringtie", "exon.nb.stringtie")])
  transcript.ids.to.check$exon.nb.stringtie <- transcript.ids.to.check$exon.nb.stringtie + 1
  transcript.ids.to.check$id <- with(transcript.ids.to.check,
                                     paste0(transcript.stringtie, "__", exon.nb.stringtie))
  selected.stringtie <- subset(stringtie.gr, paste0(transcript_id, "__", exon_number) %in% transcript.ids.to.check$id)
  exons.to.add <- NULL
  
  while (length(selected.stringtie) > 0) {
    names(selected.stringtie) <- with(selected.stringtie, paste0(transcript_id, "__", exon_number))
    selected.stringtie.overlap <- as.data.frame(findOverlaps(selected.stringtie, input.gr.combined))
    selected.stringtie.overlap$id <- names(selected.stringtie)[selected.stringtie.overlap$queryHits]
    
    transcript.ids.to.check <- subset(transcript.ids.to.check, !id %in% selected.stringtie.overlap$id &
                                        id %in% names(selected.stringtie))
    if (nrow(transcript.ids.to.check) > 0) {
      exons.to.add.this.step <- transcript.ids.to.check
      exons.to.add.this.step$start <- start(selected.stringtie[exons.to.add.this.step$id])
      exons.to.add.this.step$end <- end(selected.stringtie[exons.to.add.this.step$id])
      exons.to.add.this.step$chr <- as.character(seqnames(selected.stringtie[exons.to.add.this.step$id]))
      exons.to.add.this.step$strand <- "+"
      exons.to.add <- rbind(exons.to.add,
                            exons.to.add.this.step)
      transcript.ids.to.check$exon.nb.stringtie <- transcript.ids.to.check$exon.nb.stringtie + 1
      transcript.ids.to.check$id <- with(transcript.ids.to.check,
                                         paste0(transcript.stringtie, "__", exon.nb.stringtie))
    }
    selected.stringtie <- subset(stringtie.gr, paste0(transcript_id, "__", exon_number) %in% transcript.ids.to.check$id)
  }
  
  # Do the same for reverse
  
  # Sort to get only the last overlapping exon:
  overlap.last <- overlap.last[order(overlap.last$start.stringtie), ]
  overlap.last.uniq.rev <- overlap.last[!duplicated(overlap.last[, c("transcript_id", "transcript.stringtie")]), ]
  
  # Select last exons which overlap non-last stringtie exons but have same extremity (splicing reason)
  non.last.rev <- subset(overlap.last.uniq.rev, strand == "-" & start.stringtie != transcript.stringtie.start & start.stringtie == start)
  non.last.rev$exon.nb.stringtie <- as.numeric(stringtie.gr$exon_number[non.last.rev$subjectHits])
  rownames(non.last.rev) <- with(non.last.rev,  paste0(transcript_id, "__", transcript.stringtie))
  
  # Determine exon ordering.
  test.transcript.id <- subset(stringtie.gr, exon_number == 2 & strand == "-")$transcript_id[1]
  start.exon1 <- start(subset(stringtie.gr, exon_number == 1 & transcript_id == test.transcript.id))
  start.exon2 <- start(subset(stringtie.gr, exon_number == 2 & transcript_id == test.transcript.id))
  if (start.exon1 < start.exon2) {
    increment <- -1
  } else {
    increment <- 1
  }
  # Will recursively check if the following exon of stringtie does not overlap any exon.
  # If it is the case, the exon will be added and the next exon checked
  transcript.ids.to.check <- unique(non.last.rev[,c("transcript_id", "transcript.stringtie", "exon.nb.stringtie")])
  transcript.ids.to.check$exon.nb.stringtie <- transcript.ids.to.check$exon.nb.stringtie + increment
  transcript.ids.to.check$id <- with(transcript.ids.to.check,
                                     paste0(transcript.stringtie, "__", exon.nb.stringtie))
  selected.stringtie <- subset(stringtie.gr, paste0(transcript_id, "__", exon_number) %in% transcript.ids.to.check$id)
  
  while (length(selected.stringtie) > 0) {
    names(selected.stringtie) <- with(selected.stringtie, paste0(transcript_id, "__", exon_number))
    selected.stringtie.overlap <- as.data.frame(findOverlaps(selected.stringtie, input.gr.combined))
    selected.stringtie.overlap$id <- names(selected.stringtie)[selected.stringtie.overlap$queryHits]
    
    transcript.ids.to.check <- subset(transcript.ids.to.check, !id %in% selected.stringtie.overlap$id &
                                        id %in% names(selected.stringtie))
    if (nrow(transcript.ids.to.check) > 0) {
      exons.to.add.this.step <- transcript.ids.to.check
      exons.to.add.this.step$start <- start(selected.stringtie[exons.to.add.this.step$id])
      exons.to.add.this.step$end <- end(selected.stringtie[exons.to.add.this.step$id])
      exons.to.add.this.step$chr <- as.character(seqnames(selected.stringtie[exons.to.add.this.step$id]))
      exons.to.add.this.step$strand <- "-"
      exons.to.add <- rbind(exons.to.add,
                            exons.to.add.this.step)
      transcript.ids.to.check$exon.nb.stringtie <- transcript.ids.to.check$exon.nb.stringtie + increment
      transcript.ids.to.check$id <- with(transcript.ids.to.check,
                                         paste0(transcript.stringtie, "__", exon.nb.stringtie))
    }
    selected.stringtie <- subset(stringtie.gr, paste0(transcript_id, "__", exon_number) %in% transcript.ids.to.check$id)
  }
  
  cat("Found", length(unique(exons.to.add$transcript_id)), "transcripts to extend.\n")
  
  # Now we add the exons to transcript.
  # Extract the exon keyword:
  exon.id.no.ext <- grep("\\.ext", input.gr$exon_id, invert = T, value = T)[1]
  exon.prefix <- regmatches(exon.id.no.ext, regexpr("^\\D*", exon.id.no.ext))
  exon.ids <- regmatches(input.gr.combined$exon_id, regexpr("\\d+", input.gr.combined$exon_id))
  max.exon.id <- max(as.numeric(exon.ids))
  exon.id.90pc <- round(unname(quantile(as.numeric(exon.ids), probs = 0.9)))
  if (2 * exon.id.90pc > max.exon.id) {
    current.exon.id <- 2 * exon.id.90pc
  } else {
    current.exon.id <- 2 * max.exon.id
  }
  nchar.exon.id <- median(sapply(exon.ids, nchar))
  exons.to.add.tr <- split(exons.to.add[, c("chr", "start", "end", "strand")], f = exons.to.add$transcript_id)
  exons.to.add.tr.gr <- lapply(exons.to.add.tr, makeGRangesFromDataFrame)
  exons.to.add.tr.gr.red <- lapply(exons.to.add.tr.gr, reduce)
  new.exons.annotated <- NULL
  for (tr.id in names(exons.to.add.tr.gr.red)) {
    current.last.exon <- subset(input.gr.last.exon, transcript_id == tr.id)
    new.exons <- exons.to.add.tr.gr.red[[tr.id]]
    metacols <- mcols(current.last.exon)
    metacols <- metacols[, grep("exon", colnames(metacols), invert = T)]
    mcols(new.exons) <- metacols
    new.exons$source <- "improve_gtf3p"
    if (as.character(strand(current.last.exon)) == "+") {
      new.exons$exon_number <- as.character(as.numeric(current.last.exon$exon_number) + 1:length(new.exons))
    } else {
      new.exons$exon_number <- as.character(as.numeric(current.last.exon$exon_number) + length(new.exons):1)
    }
    new.exons$exon_id <- sprintf(paste0(exon.prefix, "%0", nchar.exon.id, "d.new"), current.exon.id + 0:(length(new.exons) - 1))
    current.exon.id <- current.exon.id + length(new.exons)
    new.exons.annotated <- c(new.exons.annotated, new.exons)
  }
  new.exons.annotated.unlist <- suppressWarnings(unlist(as(new.exons.annotated, "GRangesList")))
  cat("Added", length(new.exons.annotated.unlist), "exons.\n")
  input.gr.combined <- c(input.gr.combined, new.exons.annotated.unlist)
} else {
  cat("Found no transcript to extend.\n")
}

# Add the CDS which we do not extend
input.gr.combined <- c(input.gr.combined, input.gr.cds)

# Export to file
export.gff(sort(input.gr.combined), output.gtf)
