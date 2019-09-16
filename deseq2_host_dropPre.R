library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("Biobase") #for ExpressionSet
library("pheatmap") # for sample heatmap
library("VennDiagram") # for Venn diagram of #s of DEGs
library("tidyverse") # for wrangling, plotting
library("vegan") # for PCoA
library("ape") # for PCoA

setwd("~/Dropbox/FGB_Harvey_DESeq2/map2tran_JP/deseq2_host_dropPre/")

# Load data ---------------------------------------------------------------
countdata = read.table("../allcounts_Harvey_OfavB1_june2019.txt",header=TRUE,row.names = 1) 
head(countdata) 
length(countdata[,1])
# 51936 isogroups mapped

# make conditions table in excel
# write.csv(names(countdata),file="conds_prepost.csv")
conds <- read.csv("../conds_prepost.csv")
head(conds)

# change sample names
names(countdata) <- conds$sam
head(countdata)

# DROP FIRST TIME POINT PER REVIEWER SUGGESTION ------
countdata <- countdata %>% select(grep("recovery",names(countdata)))
# got rid of recovery?
length(grep("stress",names(countdata)))
length(grep("recovery",names(countdata)))
# yes

conds <- conds %>% filter(time=="recovery")
head(conds)

# split coral and symbiont ---------------
head(countdata)
tail(countdata)

sym.rows <- grep("sym", row.names(countdata))
coral.count <- countdata[-sym.rows,]
head(coral.count)
tail(coral.count)

sym.count <- countdata[sym.rows,]
head(sym.count)
tail(sym.count)

# rows add up to total? didn't miss anything?
nrow(coral.count) # 19937
nrow(sym.count) # 31999
nrow(coral.count)+nrow(sym.count) == nrow(countdata)
#yep

# CORAL ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts = colSums(coral.count)
min(totalCounts) # 1953
mean(totalCounts) # 2260258
max(totalCounts)  # 6961098

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = coral.count,
  colData = conds,
  design = ~ species + bank)

save(conds, coral.count, dds, file="ddsCoral_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(coral.count,1,mean)
table(means>3)
# FALSE  TRUE 
# 10275  9662 

means3 <- names(means[means>3])
head(means3)
length(means3)
#9662

coral.countFilt <- coral.count[row.names(coral.count) %in% means3,]
head(coral.countFilt)

totalCountsFilt <- colSums(coral.countFilt)
totalCountsFilt

min(totalCountsFilt) #1944
max(totalCountsFilt) #6933865
mean(totalCountsFilt) #2252948

# check sample order
table(names(coral.countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = coral.countFilt,
  colData = conds,
  design = ~ species + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("species", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out3 <- c("fav_recovery_west_21", "fav_recovery_west_30", "fav_recovery_east_35")
out2 <- c("frank_recovery_west_23")
out1 <- c("fav_recovery_east_17", "fav_recovery_east_18")

# remove out3 and out2
dim(coral.countFilt)
coral.count.out <- coral.countFilt %>% select(-one_of(c(out3, out2)))
dim(coral.count.out)

dim(conds)
head(conds)
conds.out <- conds %>% filter(!sam %in% c(out3,out2))
dim(conds.out)

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = coral.count.out,
  colData = conds.out,
  design = ~ species + bank)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)
#-- replacing outliers and refitting for 93 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 

#---Results
# log2 fold change (MLE): species fav vs frank 
resSpecies <- results(deds, independentFiltering = F, contrast=c("species","fav","frank"))
resSpecies

# log2 fold change (MLE): bank east vs west 
resBank <- results(deds, independentFiltering = F, contrast=c("bank", "east", "west"))
resBank

# Load and save
save(conds.out, coral.count.out, ddsFiltOut, deds,
     resSpecies, resBank, file="ddsCoralFiltOut.Rdata")
load("ddsCoralFiltOut.Rdata")

# how many genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resSpecies$padj < 0.05)
# FALSE  TRUE 
#  9408   183 

table(resBank$padj < 0.05)
# FALSE  TRUE 
#  9590     1

# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# FALSE  TRUE 
# 9043   619 
write.csv( as.data.frame(resSpecies), file="resSpecies.csv" ) 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE 
# 9287   375 
write.csv( as.data.frame(resBank), file="resBank.csv" ) 

# new vsd ------
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# differential expression as percentage of # of reads ----
isogroupsMapped <- 19937

# how many DEGs as percentage of transcriptome size?
round((table(resBank$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.01 %

round((table(resSpecies$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.92 %

# Save/Load Data ----------------------------------------------------------
save(conds.out, coral.count.out, ddsFiltOut, deds, 
     resSpecies, resBank, vsd, file="results.Rdata")
load("results.Rdata")

# annotated list of top DEGs -------
load("results.Rdata")

gg <- read.delim("~/Dropbox/genomes/ofav/orb_fav_iso2gene.tab", header=F)
head(gg)

annotsig <- merge(as.data.frame(resBank), gg, 
                  by.x=0, by.y=1) %>% select(Row.names, log2FoldChange, padj, V2)
head(annotsig)

# see top DE
annotsig %>% filter(padj<.05)
# Row.names         log2FoldChange    padj                                                                  V2
# isogroup11071      -2.301571 0.01634499 Thymosin beta-12 OS=Lateolabrax japonicus PE=1 SV=2 E(blastx)=2e-14