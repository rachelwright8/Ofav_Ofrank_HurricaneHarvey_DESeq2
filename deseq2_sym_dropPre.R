library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("Biobase") #for ExpressionSet
library("pheatmap") # for sample heatmap
library("VennDiagram") # for Venn diagram of #s of DEGs
library("tidyverse") # for wrangling, plotting
library("vegan") # for PCoA
library("ape") # for PCoA

setwd("~/Dropbox/FGB_Harvey_DESeq2/map2tran_JP/deseq2_sym_dropPre/")


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

# SYMBIONT ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts = colSums(sym.count)
min(totalCounts) # 906
mean(totalCounts) # 278759
max(totalCounts)  # 1177875

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = sym.count,
  colData = conds,
  design = ~ species + bank)

save(conds, sym.count, dds, file="ddsSym_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(sym.count,1,mean)
table(means>3)
# FALSE  TRUE 
# 18701 13298 

means3 <- names(means[means>3])
head(means3)
length(means3)
#13298

sym.countFilt <- sym.count[row.names(sym.count) %in% means3,]
head(sym.countFilt)

totalCountsFilt <- colSums(sym.countFilt)
totalCountsFilt

min(totalCountsFilt) #874
max(totalCountsFilt) #1129947
mean(totalCountsFilt) # 266070.5

# check sample order
table(names(sym.countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = sym.countFilt,
  colData = conds,
  design = ~ species + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("species", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out3 <- c("fav_recovery_west_21", "fav_recovery_east_35"	)
out2 <- c("frank_recovery_east_16", "frank_recovery_west_23", "fav_recovery_west_30")

# remove out3 and out2
dim(sym.countFilt)
sym.count.out <- sym.countFilt %>% select(-one_of(c(out3,out2)))
dim(sym.count.out)

dim(conds)
head(conds)
conds.out <- conds %>% filter(!sam %in% c(out3,out2))
dim(conds.out)

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = sym.count.out,
  colData = conds.out,
  design = ~ species + bank)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)
#-- replacing outliers and refitting for 71 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 

#---Results
# log2 fold change (MLE): species fav vs frank 
resSpecies <- results(deds, independentFiltering = F, contrast=c("species","fav","frank"))
resSpecies

# log2 fold change (MLE): bank east vs west 
resBank <- results(deds, independentFiltering = F, contrast=c("bank", "east", "west"))
resBank

# Load and save
save(conds.out, sym.count.out, ddsFiltOut, deds,
     resSpecies, resBank, file="ddsSymFiltOut.Rdata")
load("ddsSymFiltOut.Rdata")

# how many genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resSpecies$padj < 0.05)
# FALSE  TRUE 
# 13262     1  

table(resBank$padj < 0.05)
# FALSE  TRUE 
# 13250    13 

# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# FALSE  TRUE 
# 13015   283 
write.csv( as.data.frame(resSpecies), file="resSpecies.csv" ) 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE 
# 12890   408 
write.csv( as.data.frame(resBank), file="resBank.csv" ) 

# new vsd ------
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# differential expression as percentage of # of reads ----
isogroupsMapped <- 31999

# how many DEGs as percentage of transcriptome size?
round((table(resBank$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.04%

round((table(resSpecies$padj<0.05)[2]/isogroupsMapped*100),2)
#0.0 %

# Save/Load Data ----------------------------------------------------------
save(conds.out, sym.count.out, ddsFiltOut, deds, 
     resSpecies, resBank, vsd, file="results.Rdata")
load("results.Rdata")

#annotated list of top DEGs -------
load("results.Rdata")

gg <- read.delim("~/Dropbox/genomes/Bminutum_JParkinson/Bmin_iso2gene.tab", header=F)
head(gg)

annotsig <- merge(as.data.frame(resBank), gg, 
                  by.x=0, by.y=1) %>% select(Row.names, log2FoldChange, padj, V2)
head(annotsig)
# not annotated