#title: DESeq2 analysis 
#description: Differential expression analysis (using DESeq2) of RNA-seq data that was quantified with kallisto.
#author: Maridel Fredericksen
#date: 03 1 2021


#follow this vignette https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#and this one which has more info http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#preparing-quantification-input-to-deseq2


##############################################################
#SETUP WORKING ENVIRONMENT
##############################################################

setwd("C:/PathToFiles")
getwd()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")#version 1.28.1
BiocManager::install("tximport")#version 1.16.1
BiocManager::install("readr")#version 1.3.1
BiocManager::install("tximportData")#version 1.16.0
BiocManager::install("rhdf5")#version 2.32.2


#If you have performed transcript quantification (with Salmon, kallisto, RSEM, etc.) you could import the data with tximport, which produces a list, and then you can use DESeqDataSetFromTximport()


library("tximport")
library("readr")
library("tximportData")
library("rhdf5")


##############################################################
#IMPORT AND FILTER DATA
##############################################################

#specify path to kallisto results
dir <- "C:/PathToFiles"
list.files(dir)
list.files(file.path(dir, "iinb1_ref_all_reads"))
list.files(file.path(dir, "iinb1_ref_all_reads", "HS_BN_I_r1_21"))

#load samples file and call it coldata (same table for X and I, but I will modify in next step)
coldata_X <- read.table(file.path(dir, "samples.tsv"), header = TRUE, stringsAsFactors = FALSE)
coldata_I <- read.table(file.path(dir, "samples.tsv"), header = TRUE, stringsAsFactors = FALSE)


#make "names" column for sample names, and "files" column for count data files
coldata_X$names <- coldata_X$sample
coldata_X$files <- file.path(dir, "xinb3_ref_all_reads", coldata_X$names, "abundance.h5")
coldata_I$names <- coldata_I$sample
coldata_I$files <- file.path(dir, "iinb1_ref_all_reads", coldata_I$names, "abundance.h5")
#check that all files exist
file.exists(coldata_X$files)
file.exists(coldata_I$files)


#Load files using this tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#kallisto

#create named vector pointing to quantification files 
samples_X <- read.table(file.path(dir, "samples.tsv"), header = TRUE, stringsAsFactors = FALSE)
samples_I <- read.table(file.path(dir, "samples.tsv"), header = TRUE, stringsAsFactors = FALSE)


files_X <- file.path(dir, "xinb3_ref_all_reads", samples_X$sample, "abundance.h5")
files_I <- file.path(dir, "iinb1_ref_all_reads", samples_I$sample, "abundance.h5")

#give names to the files, corresponding to the sample names
names(files_X) <- paste0(samples_X$sample)
names(files_I) <- paste0(samples_I$sample)


file.exists(files_X)
file.exists(files_I)

#read in transcript-level info from abundance.h5 files 
txi.kallisto_X <- tximport(files_X, type = "kallisto", txOut = TRUE)
txi.kallisto_I <- tximport(files_I, type = "kallisto", txOut = TRUE)
head(txi.kallisto_X$counts)
head(txi.kallisto_I$counts)

#create DESeqDataSet object (custom class for storing data) from txi object and sample information
library("DESeq2")
ddsTxi_X <- DESeqDataSetFromTximport(txi.kallisto_X,
                                     colData = samples_X,
                                     design = ~ treatment + clone)

ddsTxi_I <- DESeqDataSetFromTximport(txi.kallisto_I,
                                     colData = samples_I,
                                     design = ~ treatment + clone)

#pre-filter to keep only rows that have at least 10 reads total
keep_X <- rowSums(counts(ddsTxi_X)) >= 10
keep_I <- rowSums(counts(ddsTxi_I)) >= 10
ddsTxi_X <- ddsTxi_X[keep_X,]
ddsTxi_I <- ddsTxi_I[keep_I,]

#set reference level (controls) for both "treatment" and "clone"
ddsTxi_X$treatment <- relevel(ddsTxi_X$treatment, ref = "CO")
ddsTxi_I$treatment <- relevel(ddsTxi_I$treatment, ref = "CO")

ddsTxi_X$clone <- relevel(ddsTxi_X$clone, ref = "Xinb3")
ddsTxi_I$clone <- relevel(ddsTxi_I$clone, ref = "Iinb1")



##############################################################
#DIFFERENTIAL EXPRESSION ANALYSIS
##############################################################
ddsTxi_X <- DESeq(ddsTxi_X)
ddsTxi_I <- DESeq(ddsTxi_I)

#output:
#estimating size factors
#using 'avgTxLength' from assays(dds), correcting for library size
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

#generate results tables (don't need to use contrast argument to explicitly tell results which comparisons to make, because I did this previously with relevel function)
res_X <- results(ddsTxi_X)
res_I <- results(ddsTxi_I)


#alpha is set to 0.1 by default. Set adjusted p value cutoff to 0.05 instead 
res_X_05 <- results(ddsTxi_X, alpha = 0.05)
res_I_05 <- results(ddsTxi_I, alpha = 0.05)
summary(res_X_05)
  #out of 17587 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 4880, 28%
  #LFC < 0 (down)     : 5471, 31%
  #outliers [1]       : 21, 0.12%
  #low counts [2]     : 682, 3.9%
  #(mean count < 0)
  #[1] see 'cooksCutoff' argument of ?results
  #[2] see 'independentFiltering' argument of ?results

summary(res_I_05)
  #out of 17890 with nonzero total read count
  #adjusted p-value < 0.05
  #LFC > 0 (up)       : 5662, 32%
  #LFC < 0 (down)     : 5495, 31%
  #outliers [1]       : 1, 0.0056%
  #low counts [2]     : 0, 0%
  #(mean count < 0)
  #[1] see 'cooksCutoff' argument of ?results
  #[2] see 'independentFiltering' argument of ?results


sum(res_X_05$padj < 0.05, na.rm = TRUE)
#[1] 10351
sum(res_I_05$padj < 0.05, na.rm = TRUE)
#[1] 11157


##############################################################
#PLOT RESULTS 
##############################################################
#compare gene expression for genes occurring in F-locus region of both Xinb3 and Iinb1

#PlotCounts in ggplot for genes in region of interest
#normalizes counts by sequencing depth and adds pseudocount of 1/2 to allow for log scale plotting
library("ggplot2")#version 3.3.2



#MAPPED TO Xinb3

xT_2392 <- plotCounts(ddsTxi_X, gene = "maker-000011F-augustus-gene-23.92-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_2392_all <- xT_2392
xT_2392_all$treatment = "all"

xT_2392_comb <- rbind(xT_2392, xT_2392_all)
xT_2392_comb$separation <- "treatments"
xT_2392_comb[xT_2392_comb$treatment == "all",]$separation <- "all" 
xT_2392_comb$separation <- factor(xT_2392_comb$separation, levels = c("treatments", "all"))
xT_2392_comb$separation


png(file = "DESeq2_curated_rplots/2392_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_2392_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.92")+
  ylab("normalized count")

dev.off()





xT_238 <- plotCounts(ddsTxi_X, gene = "augustus_masked-000011F-processed-gene-23.8-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_238_all <- xT_238
xT_238_all$treatment = "all"

xT_238_comb <- rbind(xT_238, xT_238_all)
xT_238_comb$separation <- "treatments"
xT_238_comb[xT_238_comb$treatment == "all",]$separation <- "all" 
xT_238_comb$separation <- factor(xT_238_comb$separation, levels = c("treatments", "all"))
xT_238_comb$separation


png(file = "DESeq2_curated_rplots/238_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_238_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.8")+
  ylab("normalized count")

dev.off()





xT_2393 <- plotCounts(ddsTxi_X, gene = "maker-000011F-augustus-gene-23.93-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_2393_all <- xT_2393
xT_2393_all$treatment = "all"

xT_2393_comb <- rbind(xT_2393, xT_2393_all)
xT_2393_comb$separation <- "treatments"
xT_2393_comb[xT_2393_comb$treatment == "all",]$separation <- "all" 
xT_2393_comb$separation <- factor(xT_2393_comb$separation, levels = c("treatments", "all"))
xT_2393_comb$separation


png(file = "DESeq2_curated_rplots/2393_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_2393_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.93")+
  ylab("normalized count")

dev.off()





xT_2373 <- plotCounts(ddsTxi_X, gene = "maker-000011F-snap-gene-23.73-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_2373_all <- xT_2373
xT_2373_all$treatment = "all"

xT_2373_comb <- rbind(xT_2373, xT_2373_all)
xT_2373_comb$separation <- "treatments"
xT_2373_comb[xT_2373_comb$treatment == "all",]$separation <- "all" 
xT_2373_comb$separation <- factor(xT_2373_comb$separation, levels = c("treatments", "all"))
xT_2373_comb$separation


png(file = "DESeq2_curated_rplots/2373_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_2373_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.73")+
  ylab("normalized count")

dev.off()





xT_239 <- plotCounts(ddsTxi_X, gene = "augustus_masked-000011F-processed-gene-23.9-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_239_all <- xT_239
xT_239_all$treatment = "all"

xT_239_comb <- rbind(xT_239, xT_239_all)
xT_239_comb$separation <- "treatments"
xT_239_comb[xT_239_comb$treatment == "all",]$separation <- "all" 
xT_239_comb$separation <- factor(xT_239_comb$separation, levels = c("treatments", "all"))
xT_239_comb$separation


png(file = "DESeq2_curated_rplots/239_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_239_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.9")+
  ylab("normalized count")

dev.off()





xT_2323 <- plotCounts(ddsTxi_X, gene = "augustus_masked-000011F-processed-gene-23.23-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_2323_all <- xT_2323
xT_2323_all$treatment = "all"

xT_2323_comb <- rbind(xT_2323, xT_2323_all)
xT_2323_comb$separation <- "treatments"
xT_2323_comb[xT_2323_comb$treatment == "all",]$separation <- "all" 
xT_2323_comb$separation <- factor(xT_2323_comb$separation, levels = c("treatments", "all"))
xT_2323_comb$separation


png(file = "DESeq2_curated_rplots/2323_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_2323_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.23")+
  ylab("normalized count")

dev.off()





xT_2310 <- plotCounts(ddsTxi_X, gene = "augustus_masked-000011F-processed-gene-23.10-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_2310_all <- xT_2310
xT_2310_all$treatment = "all"

xT_2310_comb <- rbind(xT_2310, xT_2310_all)
xT_2310_comb$separation <- "treatments"
xT_2310_comb[xT_2310_comb$treatment == "all",]$separation <- "all" 
xT_2310_comb$separation <- factor(xT_2310_comb$separation, levels = c("treatments", "all"))
xT_2310_comb$separation


png(file = "DESeq2_curated_rplots/2310_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_2310_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.10")+
  ylab("normalized count")

dev.off()





xT_2398 <- plotCounts(ddsTxi_X, gene = "maker-000011F-augustus-gene-23.98-mRNA-1", intgroup = c("clone", "treatment"), returnData=TRUE)
xT_2398_all <- xT_2398
xT_2398_all$treatment = "all"

xT_2398_comb <- rbind(xT_2398, xT_2398_all)
xT_2398_comb$separation <- "treatments"
xT_2398_comb[xT_2398_comb$treatment == "all",]$separation <- "all" 
xT_2398_comb$separation <- factor(xT_2398_comb$separation, levels = c("treatments", "all"))
xT_2398_comb$separation


png(file = "DESeq2_curated_rplots/2398_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(xT_2398_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("red", "blue")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("23.98")+
  ylab("normalized count")

dev.off()




#################################################################################



#MAPPED TO IINB1

iT_2589 <- plotCounts(ddsTxi_I, gene = "FUN_002589-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2589_all <- iT_2589
iT_2589_all$treatment = "all"

iT_2589_comb <- rbind(iT_2589, iT_2589_all)
iT_2589_comb$separation <- "treatments"
iT_2589_comb[iT_2589_comb$treatment == "all",]$separation <- "all" 
iT_2589_comb$separation <- factor(iT_2589_comb$separation, levels = c("treatments", "all"))
iT_2589_comb$separation


png(file = "DESeq2_curated_rplots/2589_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2589_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2589")+
  ylab("normalized count")

dev.off()





iT_2590 <- plotCounts(ddsTxi_I, gene = "FUN_002590-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2590_all <- iT_2590
iT_2590_all$treatment = "all"

iT_2590_comb <- rbind(iT_2590, iT_2590_all)
iT_2590_comb$separation <- "treatments"
iT_2590_comb[iT_2590_comb$treatment == "all",]$separation <- "all" 
iT_2590_comb$separation <- factor(iT_2590_comb$separation, levels = c("treatments", "all"))
iT_2590_comb$separation


png(file = "DESeq2_curated_rplots/2590_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2590_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2590")+
  ylab("normalized count")

dev.off()





iT_2591 <- plotCounts(ddsTxi_I, gene = "FUN_002591-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2591_all <- iT_2591
iT_2591_all$treatment = "all"

iT_2591_comb <- rbind(iT_2591, iT_2591_all)
iT_2591_comb$separation <- "treatments"
iT_2591_comb[iT_2591_comb$treatment == "all",]$separation <- "all" 
iT_2591_comb$separation <- factor(iT_2591_comb$separation, levels = c("treatments", "all"))
iT_2591_comb$separation


png(file = "DESeq2_curated_rplots/2591_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2591_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2591")+
  ylab("normalized count")

dev.off()





iT_2592 <- plotCounts(ddsTxi_I, gene = "FUN_002592-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2592_all <- iT_2592
iT_2592_all$treatment = "all"

iT_2592_comb <- rbind(iT_2592, iT_2592_all)
iT_2592_comb$separation <- "treatments"
iT_2592_comb[iT_2592_comb$treatment == "all",]$separation <- "all" 
iT_2592_comb$separation <- factor(iT_2592_comb$separation, levels = c("treatments", "all"))
iT_2592_comb$separation


png(file = "DESeq2_curated_rplots/2592_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2592_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2592")+
  ylab("normalized count")

dev.off()





iT_2593 <- plotCounts(ddsTxi_I, gene = "FUN_002593-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2593_all <- iT_2593
iT_2593_all$treatment = "all"

iT_2593_comb <- rbind(iT_2593, iT_2593_all)
iT_2593_comb$separation <- "treatments"
iT_2593_comb[iT_2593_comb$treatment == "all",]$separation <- "all" 
iT_2593_comb$separation <- factor(iT_2593_comb$separation, levels = c("treatments", "all"))
iT_2593_comb$separation


png(file = "DESeq2_curated_rplots/2593_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2593_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2593")+
  ylab("normalized count")

dev.off()





iT_2595 <- plotCounts(ddsTxi_I, gene = "FUN_002595-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2595_all <- iT_2595
iT_2595_all$treatment = "all"

iT_2595_comb <- rbind(iT_2595, iT_2595_all)
iT_2595_comb$separation <- "treatments"
iT_2595_comb[iT_2595_comb$treatment == "all",]$separation <- "all" 
iT_2595_comb$separation <- factor(iT_2595_comb$separation, levels = c("treatments", "all"))
iT_2595_comb$separation


png(file = "DESeq2_curated_rplots/2595_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2595_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2595")+
  ylab("normalized count")

dev.off()





iT_2596 <- plotCounts(ddsTxi_I, gene = "FUN_002596-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2602_all <- iT_2602
iT_2602_all$treatment = "all"

iT_2602_comb <- rbind(iT_2602, iT_2602_all)
iT_2602_comb$separation <- "treatments"
iT_2602_comb[iT_2602_comb$treatment == "all",]$separation <- "all" 
iT_2602_comb$separation <- factor(iT_2602_comb$separation, levels = c("treatments", "all"))
iT_2602_comb$separation


png(file = "DESeq2_curated_rplots/2602_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2602_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2602")+
  ylab("normalized count")

dev.off()





iT_2604 <- plotCounts(ddsTxi_I, gene = "FUN_002604-T1", intgroup = c("clone", "treatment"), returnData=TRUE)
iT_2604_all <- iT_2604
iT_2604_all$treatment = "all"

iT_2604_comb <- rbind(iT_2604, iT_2604_all)
iT_2604_comb$separation <- "treatments"
iT_2604_comb[iT_2604_comb$treatment == "all",]$separation <- "all" 
iT_2604_comb$separation <- factor(iT_2604_comb$separation, levels = c("treatments", "all"))
iT_2604_comb$separation


png(file = "DESeq2_curated_rplots/2604_all.png", 
    width = 12, height = 9, units = "cm", res = 1000)

ggplot(iT_2604_comb, aes(x=treatment, y=count, fill=clone)) +
  geom_boxplot(outlier.shape=NA) +
  scale_fill_manual(values = c("blue", "red")) +
  geom_jitter(size=0.5)+
  facet_grid(~clone + separation, space = "free", scale = "free") +
  scale_y_log10(limits = c(0.5,400), breaks = c(1, 10, 100, 300), labels = c(1, 10, 100, 300)) +
  theme_bw()+
  theme(legend.position = "none") +
  ggtitle("2604")+
  ylab("normalized count")

dev.off()




#results info: get more info on which variables and tests were used
mcols(res_X)$description
#[1] "mean of normalized counts for all samples"   
#[2] "log2 fold change (MLE): clone Iinb1 vs Xinb3"
#[3] "standard error: clone Iinb1 vs Xinb3"        
#[4] "Wald statistic: clone Iinb1 vs Xinb3"        
#[5] "Wald test p-value: clone Iinb1 vs Xinb3"     
#[6] "BH adjusted p-values"
mcols(res_I)$description
#[1] "mean of normalized counts for all samples"   
#[2] "log2 fold change (MLE): clone Xinb3 vs Iinb1"
#[3] "standard error: clone Xinb3 vs Iinb1"        
#[4] "Wald statistic: clone Xinb3 vs Iinb1"        
#[5] "Wald test p-value: clone Xinb3 vs Iinb1"     
#[6] "BH adjusted p-values"    

#Export Results to csv 
write.csv2(as.data.frame(res_I), file = "DESeq2_res_I.csv")
write.csv2(as.data.frame(res_X), file = "DESeq2_res_X.csv")





#################################################################################
#FIND MEAN COUNTS FOR EACH CLONE (MAPPED TO XINB3) FOR SUMMARY TABLE
#################################################################################
#need to subtract 0.5 to remove pseudocount that was added for plotting


#mean counts for Xinb3
xT_2392_x <- xT_2392[xT_2392$clone=="Xinb3",]
x2392 <- mean(xT_2392_x$count) - 0.5

xT_2384_x <- xT_2384[xT_2384$clone=="Xinb3",]
x2384 <- mean(xT_2384_x$count) - 0.5

xT_2320_x <- xT_2320[xT_2320$clone=="Xinb3",]
x2320 <- mean(xT_2320_x$count) - 0.5

xT_238_x <- xT_238[xT_238$clone=="Xinb3",]
x238 <- mean(xT_238_x$count) - 0.5

xT_2393_x <- xT_2393[xT_2393$clone=="Xinb3",]
x2393 <- mean(xT_2393_x$count) - 0.5

xT_2394_x <- xT_2394[xT_2394$clone=="Xinb3",]
x2394 <- mean(xT_2394_x$count) - 0.5

xT_2373_x <- xT_2373[xT_2373$clone=="Xinb3",]
x2373 <- mean(xT_2373_x$count) - 0.5

xT_239_x <- xT_239[xT_239$clone=="Xinb3",]
x239 <- mean(xT_239_x$count) - 0.5

xT_2323_x <- xT_2323[xT_2323$clone=="Xinb3",]
x2323 <- mean(xT_2323_x$count) - 0.5

xT_2310_x <- xT_2310[xT_2310$clone=="Xinb3",]
x2310 <- mean(xT_2310_x$count) - 0.5

xT_2398_x <- xT_2398[xT_2398$clone=="Xinb3",]
x2398 <- mean(xT_2398_x$count) - 0.5

xmeans <- c(x2392, x2384, x2320, x238, x2393, x2394, x2373, x239, x2323, x2310, x2398)



#mean counts for Iinb1
xT_2392_i <- xT_2392[xT_2392$clone=="Iinb1",]
i2392 <- mean(xT_2392_i$count) - 0.5

xT_2384_i <- xT_2384[xT_2384$clone=="Iinb1",]
i2384 <- mean(xT_2384_i$count) - 0.5

xT_2320_i <- xT_2320[xT_2320$clone=="Iinb1",]
i2320 <- mean(xT_2320_i$count) - 0.5

xT_238_i <- xT_238[xT_238$clone=="Iinb1",]
i238 <- mean(xT_238_i$count) - 0.5

xT_2393_i <- xT_2393[xT_2393$clone=="Iinb1",]
i2393 <- mean(xT_2393_i$count) - 0.5

xT_2394_i <- xT_2394[xT_2394$clone=="Iinb1",]
i2394 <- mean(xT_2394_i$count) - 0.5

xT_2373_i <- xT_2373[xT_2373$clone=="Iinb1",]
i2373 <- mean(xT_2373_i$count) - 0.5

xT_239_i <- xT_239[xT_239$clone=="Iinb1",]
i239 <- mean(xT_239_i$count) - 0.5

xT_2323_i <- xT_2323[xT_2323$clone=="Iinb1",]
i2323 <- mean(xT_2323_i$count) - 0.5

xT_2310_i <- xT_2310[xT_2310$clone=="Iinb1",]
i2310 <- mean(xT_2310_i$count) - 0.5

xT_2398_i <- xT_2398[xT_2398$clone=="Iinb1",]
i2398 <- mean(xT_2398_i$count) - 0.5

imeans <- c(i2392, i2384, i2320, i238, i2393, i2394, i2373, i239, i2323, i2310, i2398)
gene_number <- c("2392", "2384", "2320", "238", "2393", "2394", "2373", "239", "2323", "2310", "2398")

mean_counts <- data.frame(gene_number,xmeans, imeans)
