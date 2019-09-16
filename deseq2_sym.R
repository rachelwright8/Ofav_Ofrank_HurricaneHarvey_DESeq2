library("DESeq2") # for differntial gene expression analysis
library("arrayQualityMetrics") # to call outliers
library("Biobase") #for ExpressionSet
library("pheatmap") # for sample heatmap
library("VennDiagram") # for Venn diagram of #s of DEGs
library("tidyverse") # for wrangling, plotting
library("vegan") # for PCoA
library("ape") # for PCoA


setwd("~/Documents/davieslab/fgb_Harvey/map2tran_JP/deseq2_sym/")

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

# SYMBIONT ANALYSIS -------------------------------------------
#----------Total counts?
totalCounts = colSums(sym.count)
min(totalCounts) # 906
mean(totalCounts) # 216,259.2
max(totalCounts)  # 1,177,875

# Construct data object ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = sym.count,
  colData = conds,
  design = ~ species + time + bank)

save(conds, sym.count, dds, file="ddsSym_Start.Rdata")

# Set base mean minimum ---------------------------------------------------
means <- apply(sym.count,1,mean)
table(means>3)
# FALSE  TRUE 
# 20262 11737 

means3 <- names(means[means>3])
head(means3)
length(means3)
#11737

sym.countFilt <- sym.count[row.names(sym.count) %in% means3,]
head(sym.countFilt)

totalCountsFilt <- colSums(sym.countFilt)
totalCountsFilt

min(totalCountsFilt) #859
max(totalCountsFilt) #1,103,812
mean(totalCountsFilt) #202,569.3

# check sample order
table(names(sym.countFilt) == as.vector(conds$sam))

# Reconstruct data object (filtered) ------------------------------

ddsFilt <- DESeqDataSetFromMatrix(
  countData = sym.countFilt,
  colData = conds,
  design = ~ species + time + bank)

# Call outliers -----------------------------------------------------------
vsd <- varianceStabilizingTransformation(ddsFilt, blind=TRUE)
e <- ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
arrayQualityMetrics(e, intgroup=c("species", "time", "bank"), force=T)

# outliers that fail 3, 2, or 1 tests
out3 <- c("fav_recovery_west_21","fav_recovery_east_35")

out2 <- c("frank_recovery_west_23", "fav_recovery_west_30", "fav_stress_east_49", 
          "fav_stress_east_51", "fav_stress_east_54")

out1 <- c("frank_recovery_east_16")

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
  design = ~ species + time + bank)

# DESeq -------------------------------------------------------------------

#-------------DESeq pipeline in one step: makes large DESeqDataSet
deds <- DESeq(ddsFiltOut)
#-- replacing outliers and refitting for 20 genes
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
save(conds.out, sym.count.out, ddsFiltOut, deds,
     resSpecies, resBank, resHarvey, file="ddsSymFiltOut.Rdata")
load("ddsSymFiltOut.Rdata")

# how many genes pass multiplicity-corrected 0.05 FDR cutoff?
table(resSpecies$padj < 0.05)
# FALSE  TRUE 
# 11726     2  

table(resBank$padj < 0.05)
# FALSE  TRUE 
# 11707    21 

table(resHarvey$padj < 0.05)
# FALSE  TRUE 
# 10257  1471 

# Extract log-fold changes --------
table(abs(resSpecies$log2FoldChange)>1.5)
# FALSE  TRUE 
# 11673    64 
write.csv( as.data.frame(resSpecies), file="resSpecies.csv" ) 

table(abs(resBank$log2FoldChange)>1.5)
# FALSE  TRUE 
# 11430   307
write.csv( as.data.frame(resBank), file="resBank.csv" ) 


table(abs(resHarvey$log2FoldChange)>1.5)
# FALSE  TRUE 
# 10759   978  
write.csv( as.data.frame(resHarvey), file="resHarvey.csv" ) 

# new vsd ------
vsd <- varianceStabilizingTransformation(ddsFiltOut, blind=TRUE)

# differential expression as percentage of # of reads ----

isogroupsMapped <- 31999

# how many DEGs as percentage of transcriptome size?
round((table(resBank$padj<0.05)[2]/isogroupsMapped*100),2)
# 0.07 %
round((table(resHarvey$padj<0.05)[2]/isogroupsMapped*100),2)
# 4.6 %
round((table(resSpecies$padj<0.05)[2]/isogroupsMapped*100),2)
#0.01 %

# Save/Load Data ----------------------------------------------------------
save(conds.out, sym.count.out, ddsFiltOut, deds, 
     resSpecies, resBank, resHarvey, vsd, file="results.Rdata")
load("results.Rdata")

# Explore with plots ------------------------------------------------------

# Load annotations
gg0 <- read.delim("~/Documents/genomes/Bminutum_JParkinson/Bmin_iso2gene.tab", header=F)
head(gg0)

# get rid of rows that only have "-" (no annotation) by selecting rows that have any alphanumeric text
gg <- gg0 %>% filter(grepl('[A-Za-z0-9]',gg0$V2))
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
# 11737    46
table(complete.cases(vsdpvals))
  
write.csv(vsdpvals, "13june19_harvey_bm3_wald_RLDandPVALS.csv", quote=F)

# Venn diagram ---------
vsdpvals <- read.csv("13june19_harvey_bm3_wald_RLDandPVALS.csv")
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

# Write results for GO/KOG analysis -------------------------------------------

# by -log p-value
head(resSpecies)
row.names(resSpecies) <- gsub("sym", "", row.names(resSpecies))
logs <- data.frame(cbind("gene"=row.names(resSpecies),"logP"=round(-log(resSpecies$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resSpecies$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
#  5700 6037 
logs$logP <- logs$logP*sign
head(logs)
write.table(logs,quote=F,row.names=F,file="GO_Species_logP.csv",sep=",")

row.names(resBank) <- gsub("sym", "", row.names(resBank))
logs <- data.frame(cbind("gene"=row.names(resBank),"logP"=round(-log(resBank$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resBank$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 3966 7771
logs$logP <- logs$logP*sign
head(logs)
write.table(logs,quote=F,row.names=F,file="GO_Bank_logP.csv",sep=",")

row.names(resHarvey) <- gsub("sym", "", row.names(resHarvey))
logs <- data.frame(cbind("gene"=row.names(resHarvey),"logP"=round(-log(resHarvey$pvalue+1e-10,10),1)))
logs$logP <- as.numeric(as.character(logs$logP))
sign <- rep(1,nrow(logs))
sign[resHarvey$log2FoldChange<0]=-1  ##change to correct model
table(sign)
# -1    1 
# 2980 8757
logs$logP <- logs$logP*sign
head(logs)
write.table(logs,quote=F,row.names=F,file="GO_Harvey_logP.csv",sep=",")

# heatmaps ------------
head(vsdpvals)
names(vsdpvals)
exp <- vsdpvals[c(1:40)]
head(exp)

# categories for sorting
names(exp)
frank <- c(grep("frank",names(exp)))
faveolata <- c(grep("fav",names(exp)))
west <-  c(grep("west",names(exp)))
east <-  c(grep("east",names(exp)))
stress <-  c(grep("stress",names(exp)))
recovery <-  c(grep("recovery",names(exp)))

#     Make p-value cut-offs
sig.s = row.names(vsdpvals[vsdpvals$padj.s<0.05 & !is.na(vsdpvals$padj.s),])
length(sig.s)
# 2

sig.b = row.names(vsdpvals[vsdpvals$padj.b<0.05 & !is.na(vsdpvals$padj.b),])
length(sig.b)
# 21

sig.harvey = row.names(vsdpvals[vsdpvals$padj.h<0.0001 & !is.na(vsdpvals$padj.h),])
length(sig.harvey)
# 155
# 48

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

# only one gene is annotated:
# sp|Q8RXV3;Probable protein phosphatase 2C 59.comp29259

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
heat.colors = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=1.4)(100)

# cluster plot
pheatmap(as.matrix(explc_sort),color = heat.colors,cex=0.9,
         border_color=NA,clustering_distance_rows="correlation", cluster_cols=F)

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
  c("chocolate1","#FEE090","grey10", "cyan3","cyan")),bias=0.5)(100)

# cluster plot
pdf(file = "fig_heatmap_harvey_p0.0001.pdf", width = 10, height = 3.5)
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

#           Df  SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
# species    1   57214750  57214750 0.99076 0.02385  0.366   
# bank       1   81533746  81533746 1.41189 0.03399  0.125   
# time       1  181249751 181249751 3.13863 0.07555  0.007 **
# Residuals 36 2078931968  57748110         0.86661          
# Total     39 2398930214                   1.00000 

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
pdf(file = "fig_pca_sym.pdf", width = 7, height = 6)
plot(scores[,xaxis], scores[,yaxis],type="n", 
     main = expression(paste(italic("B. minutum "), "Gene Expression")),
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

gg <- read.delim("~/Documents/genomes/Bminutum_JParkinson/Bmin_iso2gene.tab", header=F)
head(gg)

row.names(resHarvey) <- gsub("sym", "", row.names(resHarvey))

annotsig <- merge(as.data.frame(resHarvey), gg, 
                  by.x=0, by.y=1) %>% select(Row.names, log2FoldChange, padj, V2)
head(annotsig)
annotsig[(annotsig$padj<0.01),]

# search by gene name
annotsig[grep("nonsense", annotsig$V2),]

# search by isogroup number
annotsig[grep("24084", annotsig$Row.names),]
