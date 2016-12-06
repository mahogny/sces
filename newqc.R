######################################################################
##### this is for miseq 13, containing both a crispr and ss2 plate
######################################################################


install.packages(c("gplots","scatterplot3d","org.Mm.eg.db","GO.db","plotrix"))
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")


library("pheatmap")
library("DESeq2")
library(DESeq)
library(plyr)
library(biomaRt)
library(gplots)
library(scatterplot3d)
#library(org.Mm.eg.db)
#library(GO.db)
library(plotrix)
library(Matrix)

######################################################################
### Load data ########################################################
######################################################################

listnongenes <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique","grna-numread")

# read normal gene counts
readcounttable <- function(){
  dat <- read.csv("counts/rawm13.csv",stringsAsFactors=FALSE)
  dim(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  dat <- dat[-grep("grna",rownames(dat)),]
  rownames(dat)[which(rownames(dat)=="CAS9")] <- "cas9"
  
  # read grna counts
  dat2 <- read.csv("counts/rawgrnam13.csv",stringsAsFactors=FALSE)
  dim(dat2)
  rownames(dat2) <- dat2[,1]
  dat2 <- dat2[,-1]
  rownames(dat2) <- sprintf("grna-%s",rownames(dat2))  #rename all to avoid clashes... not sure if this is a good idea
  colnames(dat2) <- colnames(dat)
  #would be better if this was done on the cluster side and written into one single file
  
  dat <- rbind(dat,dat2)
  dat
}
dat <- readcounttable()

sum(dat[,150])

### Construct plate design
getplatedesign <- function(){
  thecols <- rep(1:12,8)
  thecols <- c(thecols,thecols)
  therows <- as.double(sapply(1:8,function(i) rep(i,12)))
  therows <- c(therows,therows)
  theplate <- as.double(sapply(1:2,function(i) rep(i,96)))
  iswt  <-thecols<12 & therows %in% c(1,2)
  level1<-therows %in% c(3,4)
  level2<-therows %in% c(5,6)
  level3<-therows %in% c(7,8)
  levelnum <- iswt*1 + level1*2 + level2*3 + level3*4
  isgood <- levelnum>0
  cellcondition <- data.frame(
    col = thecols,
    row = therows,
    plate = theplate,
    iscr = (1:192)<=96,
    isss = (1:192)>96,
    iswt=iswt, level1=level1, level2=level2, level3=level3,
    levelnum=levelnum,
    isgood=isgood
  )
}
cellcondition <- getplatedesign()

######################################################################
### Read ensembl #####################################################
######################################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]

ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart=ensembl),stringsAsFactors=FALSE)
togenesym <- function(ensid){
  ensconvert$mgi_symbol[which(ensconvert$ensembl_gene_id %in% ensid)]
}

######################################################################
### QC for cells #####################################################
######################################################################

## mitochondrial content
count_cell <- dat
gene_count <- colSums(count_cell[grep("ENSMUSG",rownames(count_cell)),])
nogene_count <- count_cell[which(rownames(count_cell)=="no_feature"),]
nomapped_count <- count_cell[which(rownames(count_cell)=="not_aligned"),]
exon_prop <- rbind(gene_count, nogene_count, nomapped_count)
exon_prop <- t(t(exon_prop) / colSums(exon_prop))
mt_counts <- colSums(count_cell[mt_genes$ensembl_gene_id,])
mt_prop <- mt_counts / gene_count  

#pdf("qc_exonprop_mtprop_cell_2i_2.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(exon_prop[1,], mt_prop, pch=20, xlab="Proportion of reads mapped to exons",
     ylab="Proportion of reads mapped to mitochondrial genes", xlim=c(0,1), ylim=c(0,1))
abline(h=0.1, v=0.5, col="red", lty=2, lwd=2)
#dev.off()

#Plot result per plate
plot(mt_prop)

plot(colSums(dat>0))
plot(as.double(gene_count))
head(nomapped_count)

#estimate reads from the hiseq. on the 3 plate hiseq, 33% more reads
mean(gene_count[cellcondition$iscr])*1.28*2*300/20/2  #28% improvement. only got 50% of the reads. 150M vs 20M. twice num plates
mean(gene_count[cellcondition$isss])*2*150/20
#number of reads per cell
plot(as.double(colSums(dat[-which(rownames(dat3)=="grna-numread"),])))

#sum(as.double(colSums(dat3))[1:96])
#sum(as.double(colSums(dat3))[97:192])


plot(as.double(dat3[which(rownames(dat3)=="grna-numread"),])/
  as.double(colSums(dat3))) 



### TODO detected number of genes?

### Compute normalized counts
mynormalize <- function(dat){
  tokeep <- apply(dat,2,sum)>1e5  #don't try to normalize empty cells
  dat666 <- dat[-which(rownames(dat) %in% listnongenes),tokeep]
  sf <- estimateSizeFactorsForMatrix(dat666)
  sfnew <- rep(1,ncol(dat))
  sfnew[(1:ncol(dat))[tokeep]]<-sf
  for(i in 1:ncol(dat)){
    dat[,i] <- dat[,i]/sfnew[i]
  }
  dat
}
ncount_all <- mynormalize(dat)

### Comparing average gene expression among conditions

#ncount_all <- dat[-which(rownames(dat3)=="grna-numread"),]

gene_mean_cr <- rowMeans(ncount_all[,cellcondition$iscr]) 
gene_mean_ss <- rowMeans(ncount_all[,cellcondition$isss])
gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)
gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
#head(gene_mean_cr)

gene_mean_condition <- cbind(gene_mean_cr, gene_mean_ss)
colnames(gene_mean_condition) <- c("cr", "ss")
#pdf("comparison_gene_mean_conditions.pdf")
my_line <- function(x,y){
  points(x,y,pch=20, cex=.2, col="darkgray")
  abline(0, 1, col="black", lty=3, lwd=3)
}
pairs(gene_mean_condition, log="xy", panel=my_line, lower.panel=NULL)
#dev.off()

### Different ways of coloring the cells
colorbylevel <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$level1]<-"green"
  thecol[cellcondition$level2]<-"blue"
  thecol[cellcondition$level3]<-"red"
  thecol  
}
colorbyplate <- function(){
  thecol <- rep("black",nrow(cellcondition))
  #thecol[cellcondition$level3]<-"red"
  thecol[cellcondition$isss]<-"red"
  thecol
}

### PCA
length(which(gene_mean_ss>=0.00001))
#cor_ncount <- cor(ncount_all[gene_mean_ss>=0.00001,], method="spearman")
takecells <- cellcondition$isss #& !cellcondition$iswt 
cor_ncount <- cor(ncount_all[gene_mean_ss>=0.00001 ,takecells], method="spearman")  #olas way
#cor_ncount <- t(ncount_all[gene_mean_ss>=0.00001 ,takecells])  #ola does not scale
pca <- prcomp(cor_ncount, scale=FALSE)
#pdf("pca2d_cor_50.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(pca$x[,1], pca$x[,2], pch=20, col=colorbylevel()[takecells], xlab="PC1", ylab="PC2")
#dev.off()




######################################################################
### DEseq model testing ##############################################
######################################################################

# Only use normal genes, raw counts
datfordeseq  <- dat[grep("ENSMUSG",rownames(dat)),]

#compute average cas9 and ptma levels for each condition
calcrnacondition <- function(){
  forgenes <- c(sprintf("grna-ptma%s",1:5),"cas9")
  avglevel <- matrix(ncol = 4, nrow = length(forgenes))
  rownames(avglevel) <- c(sprintf("meanptma%s",1:5),"meancas9")
  for(i in 1:4){
    for(j in 1:length(forgenes)){
      avglevel[j,i] <- mean(as.double(ncount_all[which(rownames(ncount_all)==forgenes[j]),cellcondition$iscr & cellcondition$levelnum==i]))
      cellcondition$iscr 
    }
  }
  #and also average total ptma level
  avglevelptma <- colSums(avglevel[rownames(avglevel) %in% sprintf("meanptma%s",1:5),])
  
  #for each cell, pull out avg level
  avgforcell <- matrix(ncol=nrow(avglevel)+1, nrow=nrow(cellcondition))
  for(i in 1:ncol(avglevel)){
    for(j in 1:nrow(cellcondition)){
      if(cellcondition$levelnum[j]==i){
        avgforcell[j,] <- c(avglevel[,i], avglevelptma[cellcondition$levelnum[j]])
      }
    }
  }
  colnames(avgforcell) <- c(rownames(avglevel),"meanptma")
  
  #extract grna levels as conditions. use the normalized levels!
  cellconditiongrna <- t(ncount_all[
    c(which(rownames(ncount_all) %in% sprintf("grna-ptma%s",1:5)),
      which(rownames(ncount_all) %in% c("cas9"))),])
  colnames(cellconditiongrna) <- c(sprintf("ptma%s",1:5),"cas9")
  cellconditiongrna <- cbind(cellcondition,cellconditiongrna, avgforcell) 
  #head(cellconditiongrna)
  
  as.data.frame(cellconditiongrna,stringsAsFactors=FALSE)
}
cellconditiongrna <- calcrnacondition()

##### TODO why is the cas9 level lower on the highest level????


#######################################################################



############ Simplest model, just lowest infected cells vs highest infected cells

keepcells <- cellcondition$isss & cellcondition$levelnum%in%c(2,4)
cellcond <- data.frame(levelnum=factor(cellconditiongrna[keepcells,colnames(cellconditiongrna) %in% c("levelnum")]))
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells],
                              colData = cellcond,
                              design = ~ levelnum) #and plate? es cell type?
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds$condition <- relevel(dds$condition, ref="1")   #set the reference level to be all 0?
dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$padj),]
resOrdered
#resOrdered[1:5,]
togenesym(rownames(resOrdered)[1:5])
#which(resOrdered$padj>0.9997)
#for first miseq, ss plate only
#1 vs 4: "Ptma"   "Malat1" "Rab25"  "Nr5a2"  "Gnai3" 
#1 vs 3: "Narf"  "Scmh1" "Cdc45" "Gnai3" "Klf6" 
#2 vs 4: "Narf"    "Slc34a2" "Cdc45"   "H19"     "Gnai3"  


############ Model using average levels for each infection group

keepcells <- cellcondition$isss & cellcondition$levelnum%in%c(1,4)
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells],
                              colData = cellconditiongrna[keepcells,],
                              design = ~ meancas9 + meanptma) #and plate? es cell type?
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds$condition <- relevel(dds$condition, ref="1")   #set the reference level to be all 0?
dds <- DESeq(dds)
res <- results(dds)




select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



############ Bigger model, trying to use gRNA levels as contrasts

#### todo: do not use the dat with grna levels
### todo: 
keepcells <- cellcondition$isss
cellcond <- cellconditiongrna[keepcells,colnames(cellconditiongrna) %in% c("cas9",sprintf("ptma%s",1:5))]
rankMatrix(cellcond)
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells],
                              colData = cellcond,
                              design = ~ cas9 + ptma1 + ptma2 + ptma3 + ptma4 + ptma5) #and plate? es cell type?
#what about average level as a proxy to deal with noisyness?

dds <- dds[ rowSums(counts(dds)) > 1, ]  #suggested prefilter

dds$condition <- relevel(dds$condition, ref="untreated")   #set the reference level to be all 0?

dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$padj),]



select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)





######################################################################
### Comparison with olas previous data ###############################
######################################################################


