library("dada2")
n_cores <- 17

pathF="/scratch/workspace/caroline_desouza_uml_edu-microbe-Run3/Run3_trim/16S/R1"
pathR="/scratch/workspace/caroline_desouza_uml_edu-microbe-Run3/Run3_trim/16S/R2"

list.files(pathF)
list.files(pathR)

fnFs = sort(list.files(pathF, pattern="_R1_16S.fastq", full.names = TRUE))
fnRs = sort(list.files(pathR, pattern="_R2_16S.fastq", full.names = TRUE))

sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#look at quality profiles to decide filtration parameters
#have already looked at multiqc

#set up where the filtered outputs will be stored
filtFs = file.path(pathF, "filtered", paste0(sample.names, "_F_16S_filt.fastq.gz"))
filtRs = file.path(pathR, "filtered", paste0(sample.names, "_R_16S_filt.fastq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

#filtering and truncation
#truncLen must be large enough to maintain 20 + biological.length.variation nucleotides of overlap between them.
#since mine were run 2x250 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 10, truncLen=c(210,200),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, verbose=TRUE) # On Windows set multithread=FALSE
head(out)

#learn non-random error rates made during sequencing in R1 and R2
errF = learnErrors(filtFs, nbases=1e8, multithread=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

#removing pcr and sequencing noise
#dada2 doing its actual job
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#manual loop for merging paired ends
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}

#make sequence tables
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

write.csv(seqtab, file="/work/pi_sarah_gignouxwolfsohn_uml_edu/caroline/ASVs/Run3_16S_ASV.csv")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#frequency of chimeras
#most reads should remain after chimera removal
sum(seqtab.nochim)/sum(seqtab)

#track number of reads that made it through each step ofthe pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, file="/work/pi_sarah_gignouxwolfsohn_uml_edu/caroline/ASVs/Run3_Seqnums_16S.csv")

#taxonomy
taxa = assignTaxonomy(seqtab.nochim, "/scratch/workspace/caroline_desouza_uml_edu-microbe-Run3/Run3_trim/16S/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)

dim(taxa)

write.csv(seqtab.nochim, file="/work/pi_sarah_gignouxwolfsohn_uml_edu/caroline/ASVs/Run3_16S_ASV_nochim.csv")

write.csv(taxa, file="/work/pi_sarah_gignouxwolfsohn_uml_edu/caroline/ASVs/Run3_16S_taxa.csv")












