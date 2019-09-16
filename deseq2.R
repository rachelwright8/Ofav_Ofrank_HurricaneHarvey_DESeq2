library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("Biobase") #for ExpressionSet
library("pheatmap") # for sample heatmap
library("VennDiagram") # for Venn diagram of #s of DEGs
library("tidyverse") # for wrangling, plotting
library("vegan") # for PCoA
library("ape") # for PCoA


setwd("~/Documents/davieslab/fgb_Harvey/map2tran_JP/deseq_host/")

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
mean(totalCounts) # 1,851,268
max(totalCounts)  # 6,961,098

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = coral.count,
  colData = conds,
  design = ~ species + time + bank)

save(conds, coral.count, dds, file="ddsCoral_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(coral.count,1,mean)
table(means>3)
# FALSE  TRUE 
# 10503  9434  

means3 <- names(means[means>3])
head(means3)
length(means3)
#9434

coral.countFilt <- coral.count[row.names(coral.count) %in% means3,]
head(coral.countFilt)

totalCountsFilt <- colSums(coral.countFilt)
totalCountsFilt

min(totalCountsFilt) #1944
max(totalCountsFilt) #6,929,957
mean(totalCountsFilt) #1,843,790

# check sample order
table(names(coral.countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = coral.countFilt,
  colData = conds,
  design = ~ species + time + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("species", "time", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out3 <- c("fav_recovery_west_21", "frank_recovery_west_23", "fav_recovery_west_30", 
          "fav_recovery_east_35", "frank_stress_east_39", "fav_stress_east_54")
# out2 <- c()
out1 <- c("fav_recovery_west_15", "fav_recovery_east_17", "fav_recovery_east_18", "fav_recovery_west_31",
          "frank_stress_east_38")

# remove out3 and out2
dim(coral.countFilt)
coral.count.out <- coral.countFilt %>% select(-one_of(c(out3)))
dim(coral.count.out)

dim(conds)
head(conds)
conds.out <- conds %>% filter(!sam %in% c(out3))
dim(conds.out)

# Reconstruct data object (filtered and outliers removed) ------------------------------

ddsFiltOut <- DESeqDataSetFromMatrix(
  countData = coral.count.out,
  colData = conds.out,
  design = ~ species + time + bank)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)
#-- replacing outliers and refitting for 98 genes
#-- DESeq argument 'minReplicatesForReplace' = 7 

#---Results
# log2 fold change (MLE): species fav vs frank 
resSpecies <- results(deds, independentFiltering = F, contrast=c("species","fav","frank"))
resSpecies

# log2 fold change (MLE): bank east vs west 
resBank <- results(deds, independentFiltering = F, contrast=c("bank", "east", "west"))
resBank

# log2 fold change (MLE): time recovery vs stress 
resHarvey <- results(deds, independentFiltering = F, contrast=c("time", "recovery", "stress"))
resHarvey

# Load and save
save(conds.out, coral.count.out, ddsFiltOut, deds,
     resSpecies, resBank, resHarvey, file="ddsCoralFiltOut.Rdata")
load("ddsCoralFiltOut.Rdata")

# how many genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resSpecies$padj < 0.05)
# FALSE  TRUE 
# 8604   769 

table(resBank$padj < 0.05)
# FALSE  TRUE 
# 9358    15

table(resHarvey$padj < 0.05)
# FALSE  TRUE 
# 9108   265 

# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# FALSE  TRUE 
# 9163   271
write.csv( as.data.frame(resSpecies), file="resSpecies.csv" ) 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE 
# 9108   326
write.csv( as.data.frame(resBank), file="resBank.csv" ) 


table(abs(resHarvey$log2FoldChange)>1.5)
# FALSE  TRUE 
# 9072   362  
write.csv( as.data.frame(resHarvey), file="resHarvey.csv" ) 

# new vsd ------
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# differential expression as percentage of # of reads ----

isogroupsMapped <- 19937

# how many DEGs as percentage of transcriptome size?
round((table(resBank$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.08 %
round((table(resHarvey$padj<0.05)[2]/isogroupsMapped*100),2)
# 1.33 %
round((table(resSpecies$padj<0.05)[2]/isogroupsMapped*100),2)
# 3.86 %

# Save/Load Data ----------------------------------------------------------
save(conds.out, coral.count.out, ddsFiltOut, deds, 
     resSpecies, resBank, resHarvey, vsd, file="results.Rdata")
load("results.Rdata")

# Explore with plots ------------------------------------------------------

# Load annotations
gg <- read.delim("~/Documents/genomes/ofav/orb_fav_iso2gene.tab", header=F)
head(gg)

#Sample distance heatmap
pheatmap(cor(assay(vsd)),border_color=NA, main="SampleHeatmap")

# Diagnostics -------------------------------------------------------------

#Dispersions plot
plotDispEsts(deds, main="Dispersion Plot Response")

#MA plot
plotMA(resSpecies, ylim = c(-1, 1), main="MA Plot Species") 
plotMA(resBank, ylim = c(-1, 1), main="MA Plot Bank") 
plotMA(resHarvey, ylim = c(-1, 1), main="MA Plot Pre Vs Post") 

# Write results for making heatmaps ---------------------------------------

###--------------Get pvals
head(resSpecies)
valsSpecies <- cbind(resSpecies$pvalue, resSpecies$padj)
head(valsSpecies)
colnames(valsSpecies)=c("pval.s", "padj.s")
length(valsSpecies[,1])
table(complete.cases(valsSpecies))

head(resBank)
valsBank <- cbind(resBank$pvalue, resBank$padj)
head(valsBank)
colnames(valsBank)=c("pval.b", "padj.b")
length(valsBank[,1])
table(complete.cases(valsBank))

head(resHarvey)
valsHarvey <- cbind(resHarvey$pvalue, resHarvey$padj)
head(valsHarvey)
colnames(valsHarvey)=c("pval.h", "padj.h")
length(valsHarvey[,1])
table(complete.cases(valsHarvey))

#Make rlogdata and pvals table
vsdpvals <- cbind(assay(vsd),valsSpecies, valsBank, valsHarvey)
head(vsdpvals)
dim(vsdpvals)
# 9434   47
table(complete.cases(vsdpvals))
  
write.csv(vsdpvals, "11june19_harvey_bm3_wald_RLDandPVALS.csv", quote=F)

# Venn diagram ---------
vsdpvals <- read.csv("11june19_harvey_bm3_wald_RLDandPVALS.csv")
row.names(vsdpvals) <- vsdpvals$X
vsdpvals$X <- NULL
vsdpvals <- as.data.frame(vsdpvals)
head(vsdpvals)

species <- row.names(vsdpvals[vsdpvals$padj.s<0.05 & !is.na(vsdpvals$padj.s),])
bank <- row.names(vsdpvals[vsdpvals$padj.b<0.05 & !is.na(vsdpvals$padj.b),])
harvey <- row.names(vsdpvals[vsdpvals$padj.h<0.05 & !is.na(vsdpvals$padj.h),])

candidates <- list("Species"=species,"Bank"=bank,"Harvey"=harvey)

prettyvenn <- venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral2", "forestgreen", "grey"),
  alpha = 0.5,
  label.col = c("black", "black", "black", "black", "black", "black", "black"),
  cex = 2.5,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("black", "black", "black"),
  cat.cex = 2.5,
  cat.fontfamily = "sans",
  cat.dist = c(0.08, 0.08, 0.08),
  cat.pos = 1
);
grid.draw(prettyvenn)

# what are the genes?
harvey_unique <- harvey[!harvey %in% species]
harvey_unique
harvey_unique_vld <- vsdpvals[harvey_unique,]
row.names(harvey_unique_vld) <- sub("c", "isogroup", row.names(harvey_unique_vld))

head(gg)

harvey_unique_annot <- gg[gg$V1 %in% row.names(harvey_unique_vld),]
# 125 annotated genes DE pre/post Harvey
harvey_unique_annot

# Make heatmap of unique genes
colnames(harvey_unique_vld)
exp <- harvey_unique_vld[c(1:41)]

gnames=c();expg=c()
for(i in row.names(exp)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp[i,])
  } 
}
row.names(expg)=gnames
expl=expg
# Calculate row means and subtract
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them
# Make colors
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1)(100)

# plot
svg(filename = "fig_HarveySpecificHeatmap.svg", width = 15, height = 12)
pheatmap(as.matrix(explc),color = heat.colors,cex=0.9,
         border_color=NA,clustering_distance_rows="correlation", cluster_cols=F)
dev.off()

# Write results for GO/KOG analysis -------------------------------------------

# by -log p-value
logs <- data.frame(cbind("gene"=row.names(resSpecies),"logP"=round(-log(resSpecies$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resSpecies$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
#  4081 5353 
logs$logP <- logs$logP*sign
logs$gene <- gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_Species_logP.csv",sep=",")

logs <- data.frame(cbind("gene"=row.names(resBank),"logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 4000 5434
logs$logP <- logs$logP*sign
logs$gene <- gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_Bank_logP.csv",sep=",")

logs <- data.frame(cbind("gene"=row.names(resHarvey),"logP"=round(-log(resHarvey$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resHarvey$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3689 5745
logs$logP <- logs$logP*sign
logs$gene <- gsub("c", "isogroup", logs$gene)
write.table(logs,quote=F,row.names=F,file="GO_Harvey_logP.csv",sep=",")

# heatmaps ------------
head(vsdpvals)
names(vsdpvals)
exp <- vsdpvals[c(1:41)]
head(exp)

#     Make p-value cut-offs
sig.s = row.names(vsdpvals[vsdpvals$padj.s<1e-3 & !is.na(vsdpvals$padj.s),])
length(sig.s)
# 40

sig.b = row.names(vsdpvals[vsdpvals$padj.b<0.05 & !is.na(vsdpvals$padj.b),])
length(sig.b)
# 15

sig.harvey = row.names(vsdpvals[vsdpvals$padj.h<0.01 & !is.na(vsdpvals$padj.h),])
length(sig.harvey)
# 59

# sig expression
exp.s <- exp[row.names(exp) %in% sig.s,]
nrow(exp.s)

exp.b <- exp[row.names(exp) %in% sig.b,]
nrow(exp.b)

exp.h <- exp[row.names(exp) %in% sig.harvey,]
nrow(exp.h)

# species heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.s)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.s[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them
explc_sort <- explc[c(frank,faveolata)]
names(explc_sort)

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.2)(100)

# cluster plot
svg(filename = "fig_heatmap_species_p-03.svg", width = 13, height = 5)
pheatmap(as.matrix(explc_sort),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=F)
dev.off()

# bank heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.b)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.b[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them
explc_sort <- explc[c(east,west)]
names(explc_sort)

# make colors 
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.2)(100)

# cluster plot
svg(filename = "fig_heatmap_bank_p.05.svg", width = 13, height = 2)
pheatmap(as.matrix(explc_sort),color = heat.colors,cex=0.9,border_color=NA,clustering_distance_rows="correlation", cluster_cols=F)
dev.off()

# pre/post Harvey heatmap ---------
# naming the rows by gene names
gnames=c();expg=c()
for(i in row.names(exp.h)){
  s=subset(gg,V1==i)
  if (length(s[,1])>0){
    gnames=append(gnames,paste(s$V2[1],i,sep="."))
    expg=rbind(expg,exp.h[i,])
  } 
}
row.names(expg)=gnames
expl=expg
means=apply(expl,1,mean) # means of rows
explc=expl-means # subtracting them
names(explc)

# sort it for plotting
explc_sort <- explc %>% t() %>% as.data.frame() %>% rownames_to_column("sam") %>%
  mutate(time = ifelse(grepl("recovery", sam)==T, "recovery", "stress"),
         species = ifelse(grepl("fav", sam)==T, "fav", "frank"),
         bank = ifelse(grepl("east", sam)==T, "east", "west")) %>%
  arrange(time, species, bank) %>%
  column_to_rownames("sam") %>%
  select(-time, -species, -bank) %>%
  t()
head(explc_sort)

# make colors 
heat.colors = colorRampPalette(rev(
  c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.2)(100)

# cluster plot
pdf(file = "fig_heatmap_harvey_p0.01.pdf", width = 15, height = 4)
pheatmap(explc_sort,color = heat.colors,cex=0.9,
         border_color=NA,
         clustering_distance_rows="correlation", 
         cluster_cols=F)
dev.off()

# PCoA ---------
load("results.Rdata")

# variance stabilized expression data
exp <- data.frame(assay(vsd))
head(exp)

# condition data
table(conds.out$sam == names(exp))
head(conds.out)

# compute dissimilarity indices
dd.veg <- vegdist(t(exp), "manhattan")
div.dd.veg <- dd.veg/1000
head(div.dd.veg)

# perform PERMANOVA  
set.seed(1)
adonisRes <- adonis(t(exp)~species+bank+time,data=conds.out,method="manhattan")
adonisRes
#           Df  SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# species    1   85715834 85715834  2.2590 0.05269  0.008 **
# bank       1   60671072 60671072  1.5990 0.03730  0.063 . 
# time       1   76352197 76352197  2.0123 0.04694  0.026 * 
# Residuals 37 1403907435 37943444         0.86307          
# Total     40 1626646538                  1.00000

# compute principal coordinate decomposition
dd.pcoa <- pcoa(div.dd.veg)
head(dd.pcoa)
scores <- dd.pcoa$vectors

# plotting PCoA----
margin <- 1.5

# play around with these numbers to see different axes
xaxis <- 1
yaxis <- 2

# PCoA for mid by site type
pdf(file = "fig_pca_host.pdf", width = 7, height = 6)
pca_plot <- plot(scores[,xaxis], scores[,yaxis],type="n", 
            main = "Host Gene Expression",
     xlim=c(min(scores[,xaxis])-margin,max(scores[,xaxis])+margin),
     ylim=c(min(scores[,2])-margin,max(scores[,2])+margin),
     mgp=c(2.3,1,0),
     xlab=paste("PCo", xaxis," (", 
                round(dd.pcoa$values$Relative_eig[xaxis]*100,1),"%)",sep=""),
     ylab=paste("PCo", yaxis," (", 
                round(dd.pcoa$values$Relative_eig[yaxis]*100,1),"%)",sep=""))
ordispider(scores,conds.out$species,label=F)
points(scores[conds.out$time=="recovery" & conds.out$species=="fav",xaxis],
       scores[conds.out$time=="recovery" & conds.out$species=="fav",yaxis],
       col="black", pch=15, cex=1.5) +
points(scores[conds.out$time=="recovery" & conds.out$species=="frank",xaxis],
         scores[conds.out$time=="recovery" & conds.out$species=="frank",yaxis],
         col="black", pch=16, cex=1.5) +
points(scores[conds.out$time=="stress" & conds.out$species=="fav",xaxis],
       scores[conds.out$time=="stress" & conds.out$species=="fav",yaxis],
       col="grey", pch=15, cex=1.5) +
points(scores[conds.out$time=="stress" & conds.out$species=="frank",xaxis],
         scores[conds.out$time=="stress" & conds.out$species=="frank",yaxis],
         col="grey", pch=16, cex=1.5)

# legend of sites 
legend("topright", 
       c("Sub-lethal stress (September)", "Recovery (October)"),
       pch=c(15,15), 
       col=c("grey","black"), cex=1, bty = "n")
legend("bottomright", 
       c(expression(italic("O. faveolata")),expression(italic("O. franksi"))),
       pch=c(15,16), 
       col=c("black","black"), cex=1, bty = "n")

#insert p value 
legend("topleft",
       paste("Host Species p = ",adonisRes$aov.tab$`Pr(>F)`[1], sep=" "), 
       cex=1, bty='n')  

legend("bottomleft", 
       paste("Collection Time p = ",adonisRes$aov.tab$`Pr(>F)`[3], sep=" "), 
       cex=1, bty='n')  

dev.off()

# annotated list of top DEGs -------
load("results.Rdata")

gg <- read.delim("~/Documents/genomes/ofav/orb_fav_iso2gene.tab", header=F)
head(gg)

annotsig <- merge(as.data.frame(resHarvey), gg, 
                  by.x=0, by.y=1) %>% select(Row.names, log2FoldChange, padj, V2)
head(annotsig)

# see top DE
annotsig %>% filter(padj<1e-5)


# look for your favorite genes
annotsig[grep("decay", annotsig$V2),]

annotsig[grep("ATF5", annotsig$V2),]
annotsig[grep("24084", annotsig$Row.names),]
