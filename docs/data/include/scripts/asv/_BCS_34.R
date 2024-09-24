################################################################
#####           TRANSISTHMIAN SHRIMP MICROBIOME            #####
#####                 MISEQ RUN: BCS_34                    #####
#####               PLATES ISTHMO S01, S02                 #####
################################################################
#
# SCRIPT prepared by Matthieu Leray
# last modified January 28th 2023
#
# READ PROCESSING
#
# Modified by Jarrod Scott
# last modified September 11 2024
################################################################

library(dada2)
library(tidyverse)
library(ff)
library(phyloseq)
library(gridExtra)
library(dplyr)
library(decontam)

##Creating filepaths to data 
path <- "RAW_DATA/BCS_34"
head(list.files(path)) #eventually to check if the path works

##File preparation
#extracting Forward (fnFs) and Reverse (fnRs) reads from files
fnFs <- sort(list.files(path, pattern = "_R1_001.trimmed.fastq"))
fnRs <- sort(list.files(path, pattern = "_R2_001.trimmed.fastq"))
sample.names <- sapply(strsplit(fnFs, "_R1_"), `[`, 1)
fnFs <-file.path(path, fnFs)
fnRs <-file.path(path, fnRs)

x <- length(list.files(path, pattern = "_R1_001.trimmed.fastq"))

qprofile_fwd <- plotQualityProfile(fnFs[1:x], aggregate = TRUE)
qprofile_rev <- plotQualityProfile(fnRs[1:x], aggregate = TRUE)

qprofile <- grid.arrange(qprofile_fwd, qprofile_rev, nrow = 1)
ggsave("figures/BCS_34_filt_plot_qscores.png", qprofile, width = 7, height = 3)

#placing filtered files in a new filtered subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#filtering and trimming, here truncation at 220 (Fwd) and 180 (Rev) bp, 
#2expected errors max (N discarded automatically)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,180),
                        maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=20)

head(out)

#learning error rates
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

## ----plot_errF------------------------------
p3 <- plotErrors(errF, nominalQ = TRUE)
ggsave("figures/BCS_34_plot_errorF_1.png", p3, width = 7, height = 5)
ggsave("figures/BCS_34_plot_errorF_2.png", p3)
## ----plot_errR------------------------------
p4 <- plotErrors(errR, nominalQ = TRUE)
ggsave("figures/BCS_34_plot_errorR_1.png", p4, width = 7, height = 5)
ggsave("figures/BCS_34_plot_errorR_2.png", p4)

##Dereplicating reads
sam.names <- sapply(strsplit(basename(filtFs), "_F_"), `[`, 1)
derepFs <- derepFastq(filtFs)
names(derepFs) <- sam.names
derepRs <- derepFastq(filtRs)
names(derepRs) <- sam.names

##Infering Sequence Variants
dadaFs <- dada(derepFs, err = errF, pool = "pseudo", multithread = TRUE)
dadaFs[[1]]
dadaRs <- dada(derepRs, err = errR, pool = "pseudo", multithread = TRUE)
dadaRs[[1]]

##Merging paired ends
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
BCS_34 <- makeSequenceTable(mergers)
dim(BCS_34)
# [1]   384 29481
table(nchar(getSequences(BCS_34)))

#exporting files to use in the next part of the workflow
saveRDS(BCS_34, "BCS_34/BCS_34.rds")

#tracking changes through each step
getN <- function(x) sum(getUniques(x))
track <-    cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
                  sapply(mergers, getN))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR",
                        "merged")
rownames(track) <- sam.names
write.table(track, "BCS_34/BCS_34_read_changes.txt", sep = "\t", quote = FALSE,
            col.names=NA)

read_length <-  data.frame(nchar(getSequences(BCS_34)))

colnames(read_length) <- "length"

plot_BCS_34 <- qplot(length, data = read_length, 
                      geom = "histogram", binwidth = 1, 
                      xlab = "read length", 
                      ylab = "total variants", 
                      xlim = c(225,275)) 
ggsave("figures/read_length_before_pseudo_BCS_34.png", plot_BCS_34, width = 7, height = 3)

save.image("BCS_34.rdata")

