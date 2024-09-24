#!/usr/bin/env Rscript
set.seed(919191)
library(dada2); packageVersion("dada2")
library(ggplot2)
library(ff)
library(phyloseq)
library(gridExtra)
library(dplyr)
library(decontam)
library(grid)
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")
library(DECIPHER); packageVersion("DECIPHER")

########################################
#
# 2. MERGE ALL SEQ TABS
#
########################################

############################################
# PROBLEM: Duplicated sample names 
# detected in the sequence table row names: 
# 7512-G, 7512-H, 7512-M, 7512-S
# BCS_30 and BCS_35
# FOR BCS_30 changed name to 7512A...
# FOR BCS_35 changed name to 7512B...
############################################

BCS_26 <- readRDS("BCS_26/BCS_26.rds")
BCS_28 <- readRDS("BCS_28/BCS_28.rds")
BCS_29 <- readRDS("BCS_29/BCS_29.rds")
BCS_30 <- readRDS("BCS_30/BCS_30.rds")
BCS_34 <- readRDS("BCS_34/BCS_34.rds")
BCS_35 <- readRDS("BCS_35/BCS_35.rds")

seqtab.merge <- mergeSequenceTables(BCS_26, BCS_28, BCS_29, 
                                    BCS_30, BCS_34, BCS_35
                                    )
dim(seqtab.merge)
table(nchar(getSequences(seqtab.merge)))

read_length_all <-  data.frame(nchar(getSequences(seqtab.merge)))
colnames(read_length_all) <- "length"
plot_all <- qplot(length, data = read_length_all, 
                  geom = "histogram", binwidth = 1, 
                  xlab = "read length", 
                  ylab = "total variants", 
                  xlim = c(200,400)) 
ggsave("figures/read_length_before_collapse.png", 
        plot_all, width = 7, height = 3)
saveRDS(seqtab.merge, "2.seqtab.merge.rds")

save.image("rdata/2.merge.seqtabs.rdata")

########################################
#
# collapseNoMismatch
# TESTED, only minor differences in
# 13 samples. Takes long time to run
#
########################################

# seqtab_to_collapse <- collapseNoMismatch(st_all, minOverlap = 20, orderBy = "abundance",
#   identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)

# dim(seqtab_to_collapse)
# table(nchar(getSequences(seqtab_to_collapse)))

# read_length_all_collapse <-  data.frame(nchar(getSequences(seqtab_to_collapse)))
# colnames(read_length_all_collapse) <- "length"
# plot_all_collapse <- qplot(length, data = read_length_all_collapse, geom = "histogram", binwidth = 1, xlab = "read length", ylab = "total variants", xlim = c(200,400)) 
# ggsave("figures/read_length_after_collapse.png", plot_all_collapse, width = 7, height = 3)
# saveRDS(seqtab_to_collapse, "seqtab_after_collapse.rds")

# save.image("rdata/2_merge_seqtabs_collapsed.rdata")

########################################
#
# 3. REMOVING CHIMERAS
#
########################################

############################################
# PROBLEM: Duplicated sample names 
# detected in the sequence table row names: 
# 7512-G, 7512-H, 7512-M, 7512-S
# BCS_30 and BCS_35
# FOR BCS_30 changed name to 7512A...
# FOR BCS_35 changed name to 7512B...
############################################

## REMVOE OUTLIER READ LENGTHS
seqtab <- seqtab.merge
#seqtab.merge <- readRDS("seqtab_before_collapse.rds")
table(nchar(getSequences(seqtab)))

################################################################################# 
## 
## 220   221   222   223   224   225   226   227   228   229   230   231   232
##   125    67    14    36    20    13    10    25     9     6     4    27     2
##   234   235   236   237   238   239   240   241   242   243   244   245   246
##     9     8  1373   151    46     6    99   407   298   452    31    14    13
##   247   248   249   250   251   252   253   254   255   256   257   258   259
##    26    23    19    49   159  3587 84485  3772   319   123    96    20    10
##   260   261   262   263   264   265   266   267   268   269   270   271   272
##     8    16     9     4     2     1     1     1     4     2     9     8     4
##   273   274   275   276   277   278   279   280   281   282   284   285   286
##     7     3     2     5     1     7     4     1     1     2     1     4     4
##   288   289   290   291   292   293   294   295   296   297   298   300   303
##     1     3     1     2     4     8     7     2     3     2     3     2     3
##   304   305   307   308   309   310   311   312   313   315   316   317   318
##     1     5     2     3     2     1     3     1     3     1     4     1     3
##   319   320   321   322   323   324   325   326   328   329   330   332   333
##     3     1     2     3     2     1     3     1     3     3     2     1     3
##   334   335   336   337   338   339   340   341   342   343   344   345   346
##    13     6     7    18     5    25    16    70    52     8     7     8     4
##   347   348   349   350   351   352   353   354   355   356   357   358   359
##    25    17    21    10     2    11     1     1     7     6    31     6    15
##   360   361   362   363   364   365   366   367   368   369   370   371   372
##    21   161   188    43   141   108    19     9    26     5     3     3     8
##   373   374   376   377   378   379   380   384   385   386   387   388
##    11     2     3     5     2     3     1     1     1     1     2     1
##
#################################################################################

#######################################################
## ----REMOVE OUTLIER READ LENGTHS------------------ ##
#######################################################

seqtab.trim <- seqtab[,nchar(colnames(seqtab)) %in% seq(252, 254)]
dim(seqtab.trim)
table(nchar(getSequences(seqtab.trim)))

#####################
##   252   253   254
##  3587 84485  3772
#####################

#######################################################
## ----chimera pooled------------------------------- ##
#######################################################

seqtab.trim.nochim.pool <- 
          removeBimeraDenovo(seqtab.trim, 
                             method = "pooled", 
                             multithread = 20, 
                             verbose = TRUE)
dim(seqtab.trim.nochim.pool)
sum(seqtab.trim.nochim.pool)/sum(seqtab.trim)

saveRDS(seqtab.trim.nochim.pool, "3.seqtab.trim.nochim.pool.rds")

table(nchar(getSequences(seqtab.trim.nochim.pool)))

table(colSums(seqtab.trim.nochim.pool > 0))
table(rowSums(seqtab.trim.nochim.pool > 0))

##########################################################
## ----chimera consensus------------------------------- ##
##########################################################

seqtab.trim.nochim.consensus <- 
           removeBimeraDenovo(seqtab.trim, 
                              method = "consensus", 
                              multithread = 20, 
                              verbose = TRUE)
dim(seqtab.trim.nochim.consensus)
sum(seqtab.trim.nochim.consensus)/sum(seqtab.trim)

saveRDS(seqtab.trim.nochim.consensus, "3.seqtab.trim.nochim.consensus.rds")

table(nchar(getSequences(seqtab.trim.nochim.consensus)))

table(colSums(seqtab.trim.nochim.consensus > 0))
table(rowSums(seqtab.trim.nochim.consensus > 0))

##########################################################
## ----tracking changes-------------------------------- ##
##########################################################

getN <- function(x) sum(getUniques(x))
track <- cbind(rowSums(seqtab), 
               rowSums(seqtab.trim), 
               rowSums(seqtab.trim.nochim.pool), 
               rowSums(seqtab.trim.nochim.consensus))

colnames(track) <- c("merged", "trim", "chimera_pool", "chimera_concensus")
write.table(track, "3.chimera_read_changes_pipeline.txt", 
            sep = "\t", quote = FALSE, col.names=NA)

save.image("rdata/3.trim.chimera.rdata")

########################################
#
# 4. ASSIGNING TAXONOMY
#
########################################

###########################################################
# reference datasets formatted for DADA2 can be found here: 
# https://benjjneb.github.io/dada2/training.html
###########################################################

########################################
#
# TAXONOMY chimera = pooled
#
########################################

# seqtab <- readRDS("3.seqtab.trim.nochim.pool.rds")

########################################
# TAXONOMY = silva
########################################
seqtab.pool <- seqtab.trim.nochim.pool

tax_silva_v138.pool <- 
              assignTaxonomy(seqtab.pool, 
              "TAXONOMY_FILES/silva_nr99_v138.1_train_set.fa.gz", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_silva_v138.pool, "4.tax_silva_v138.pool.rds")

tax_silva_v132.pool <- 
              assignTaxonomy(seqtab.pool, 
              "TAXONOMY_FILES/silva_nr_v132_train_set.fa.gz", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_silva_v132.pool, "4.tax_silva_v132.pool.rds")

########################################
# TAXONOMY = RDP
########################################

tax_rdp_v138.pool <- 
             assignTaxonomy(seqtab.pool, 
             "TAXONOMY_FILES/rdp_train_set_18.fa.gz", 
             multithread = TRUE, verbose = TRUE)
saveRDS(tax_rdp_v138.pool, "4.tax_rdp_v138.pool.rds")

########################################
#
# TAXONOMY chimera = consensus
#
########################################

#remove(list = ls())
#seqtab <- readRDS("3.seqtab.trim.nochim.consensus.rds")
#objects()

########################################
# TAXONOMY = silva
########################################
seqtab.consensus <- seqtab.trim.nochim.consensus

tax_silva_v138.consensus <- 
              assignTaxonomy(seqtab.consensus, 
              "TAXONOMY_FILES/silva_nr99_v138.1_train_set.fa.gz", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_silva_v138.consensus, "4.tax_silva_v138.consensus.rds")

tax_silva_v132.consensus <- 
              assignTaxonomy(seqtab.consensus, 
              "TAXONOMY_FILES/silva_nr_v132_train_set.fa.gz", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_silva_v132.consensus, "4.tax_silva_v132.consensus.rds")

########################################
# TAXONOMY = RDP
########################################

tax_rdp_v138.consensus <- 
              assignTaxonomy(seqtab.consensus, 
              "TAXONOMY_FILES/rdp_train_set_18.fa.gz", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_rdp_v138.consensus, "4.tax_rdp_v138.consensus.rds")

########################################
# TAXONOMY = ITGDB
########################################

tax_itgdb.consensus <- 
              assignTaxonomy(seqtab.consensus, 
              "TAXONOMY_FILES/itgdb_dada2.fa", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_itgdb.consensus, "4.tax_itgdb.consensus.rds")

########################################
# TAXONOMY = GSRDB
########################################

tax_gsrdb.consensus <- 
              assignTaxonomy(seqtab.consensus, 
              "TAXONOMY_FILES/gsrdb_dada2.fa", 
              multithread = TRUE, verbose = TRUE)
saveRDS(tax_gsrdb.consensus, "4.tax_gsrdb.consensus.rds")

save.image("rdata/4.dada2.pipeline.rdata")

sessionInfo()
devtools::session_info()

quit()