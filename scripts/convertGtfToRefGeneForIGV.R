options(stringsAsFactors = F, scipen = 999)

if (!"rtracklayer" %in% installed.packages()) {
  if (!"devtools" %in% installed.packages()) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
  }
  devtools::install_github("lldelisle/usefulLDfunctions", upgrade = "never")
  library(usefulLDfunctions)
  safelyLoadAPackageInCRANorBioconductor("rtracklayer")
} else {
  library(rtracklayer)
}

rm(list = ls())

if (length(commandArgs(TRUE)) == 0) {
  cat("Choose the gtf file.\n")
  gtfFile <- file.choose()
} else {
  if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help") {
    cat("Usage: Rscript convertGtfToRefGeneForIGV.R pathForGtf \n")
    stop()
  }
  gtfFile <- commandArgs(TRUE)[1]
}

##################################################
cat("Loading gtf file...")
gtfInput <- readGFF(gtfFile)
cat("Done.\n")
###################################################
#The format of RefGene.txt is
#Col1 number (?)
#Col2 ID
#Col3 chr
#Col4 strand
#Col5 start 0-based included
#Col6 end 0-based excluded
#Col7 start of coding (or end) 0-based included
#Col8 end of coding (or end) 0-based excluded
#Col9 nb of exons
#Col10 begin of exons separated by comma, finished by comma (not necessary)
#Col11 end of exons idem
#Col12 0
#Col13 Name to display
#Col14 cmpl, incmpl or unk
#Col15 cmpl, incmpl or unk
#Col16 phase? for each exon 0,1,2 or -1

# I store the CDS:
gtfInput.cds <- subset(gtfInput, type == "CDS")
# I only work with exons
gtfInput <- subset(gtfInput, type == "exon")
gtfInput <- gtfInput[order(gtfInput$start, gtfInput$end), ]
cat("Building the transcripts from the exons.\n")
transcripts <- unique(gtfInput[, c("seqid", "strand", "gene_name", "gene_id", "transcript_id")])
if (anyDuplicated(transcripts$transcript_id) != 0) {
  print(transcripts[duplicated(transcripts$transcript_id), ])
  stop("The gtf file contains different information for the same transcript.")
}
# Get the start in 0-based
start <- aggregate(gtfInput$start - 1, by = list(transcript_id = gtfInput$transcript_id), FUN = min)
# Add it to the transcript data frame
transcripts$start <- start[match(transcripts$transcript_id, start$transcript_id), "x"]
# Get the end
end <- aggregate(gtfInput$end, by = list(transcript_id = gtfInput$transcript_id), FUN = max)
# Add it to the transcript data frame
transcripts$end <- end[match(transcripts$transcript_id, end$transcript_id), "x"]
# Get nb of exons
nbExons <- aggregate(gtfInput$exon_id, by = list(transcript_id = gtfInput$transcript_id), FUN = length)
# Add it to the transcript data frame
transcripts$nbExons <- nbExons[match(transcripts$transcript_id, nbExons$transcript_id), "x"]
# Get the position of exon starts
startsExons <- aggregate(gtfInput$start - 1, by = list(transcript_id = gtfInput$transcript_id),
                         FUN = paste, collapse = ",")
transcripts$startsExons <- startsExons[match(transcripts$transcript_id, startsExons$transcript_id), "x"]
# Get the position of exon ends
endsExons <- aggregate(gtfInput$end, by = list(transcript_id = gtfInput$transcript_id), FUN = paste, collapse = ",")
transcripts$endsExons <- endsExons[match(transcripts$transcript_id, endsExons$transcript_id), "x"]
# Get the 3' of the transcript
transcripts$endOriented <- transcripts$end
transcripts$endOriented[transcripts$strand == "-"] <- transcripts$start[transcripts$strand == "-"]
# I don't fill the phase of exons
transcripts$phaseExons <- sapply(transcripts$nbExons, function(n) {
  paste(rep("-1", n), collapse = ",")
})
# Try to get the coding part:
if (nrow(gtfInput.cds) > 0) {
  # Get the position of start coding
  leftCoding <- aggregate(gtfInput.cds$start, by = list(transcript_id = gtfInput.cds$transcript_id), FUN = min)
  transcripts$leftCoding <- leftCoding[match(transcripts$transcript_id, leftCoding$transcript_id), "x"]
  transcripts$leftCoding[is.na(transcripts$leftCoding)] <- transcripts$endOriented[is.na(transcripts$leftCoding)]
  # Get the position of end coding
  rightCoding <- aggregate(gtfInput.cds$end, by = list(transcript_id = gtfInput.cds$transcript_id), FUN = max)
  transcripts$rightCoding <- rightCoding[match(transcripts$transcript_id, leftCoding$transcript_id), "x"]
  transcripts$rightCoding[is.na(transcripts$rightCoding)] <- transcripts$endOriented[is.na(transcripts$rightCoding)]
  # I put all phase to 0
  transcripts$phaseExons[transcripts$transcript_id %in% gtfInput.cds$transcript_id] <-
    sapply(transcripts$nbExons[transcripts$transcript_id %in% gtfInput.cds$transcript_id], function(n) {
    paste(rep("0", n), collapse = ",")
  })
} else {
  transcripts$leftCoding <- transcripts$endOriented
  transcripts$rightCoding <- transcripts$endOriented
}

# Put id if name does not exists:
transcripts$gene_name[is.na(transcripts$gene_name)] <- transcripts$gene_id[is.na(transcripts$gene_name)]

finalDataFrame <- data.frame(nb = 1:nrow(transcripts), id = transcripts$transcript_id,
                             chr = transcripts$seqid, strand = transcripts$strand,
                             start = transcripts$start, end = transcripts$end,
                             scod = transcripts$leftCoding, ecod = transcripts$rightCoding,
                             nbE = transcripts$nbExons, bE = transcripts$startsExons,
                             eE = transcripts$endsExons, v12 = rep(0, nrow(transcripts)),
                             name = transcripts$gene_name, v14 = rep("unk", nrow(transcripts)),
                             v15 = rep("unk", nrow(transcripts)), pE = transcripts$phaseExons)
cat("Writting the RefGene.txt file...")
write.table(finalDataFrame, paste0(dirname(gtfFile), "/RefGene.txt"), sep = "\t", quote = F,
            row.names = F, col.names = F)
cat("Done.\n")
