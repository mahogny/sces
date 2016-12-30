######################################################################
##### this is for the final 3 hiseq
######################################################################

#also need system libraries: libssl. libcurl?
installdeps <- function(){
  install.packages(c("gplots","scatterplot3d","plotrix","sqldf","Matrix"))
  install.packages("tsne")
  install.packages("Rtsne")
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("biomaRt","DESeq","DEseq2"))
  biocLite(c("monocle"))
  #biocLite("topGO")
  biocLite("limma")
}

library(monocle)
library(Rtsne)
library(tsne)
library(pheatmap)
library(DESeq2)
library(DESeq) #really need both?
library(plyr)
library(biomaRt)
library(gplots)
library(scatterplot3d)
library(plotrix)
library(Matrix)
library(sqldf)
#library(topGO)
library(limma)

######################################################################
### Load data ########################################################
######################################################################

listnongenes <- c(
  "no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique","grna-numread",
  "__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")

# read normal gene counts
readcounttablejhuma <- function(){
  dat <- read.csv("counts/henriksson_jhuma20161222.csv",stringsAsFactors=FALSE)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  rownames(dat)[which(rownames(dat)=="CAS9")] <- "cas9"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-1")] <- "grna-ptma1"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-2")] <- "grna-ptma2"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-3")] <- "grna-ptma3"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-4")] <- "grna-ptma4"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-5")] <- "grna-ptma5"
  dat
}
dat <- readcounttablejhuma()


### Construct plate design
getplatedesignjhuma <- function(){
  gene <- read.csv("jhumashrna.txt",stringsAsFactors = FALSE)[,1]
  gene <- rep(gene,3)
  mouse <- c(rep(1,13), rep(2,13), rep(3,13))
  isgood <- rep(TRUE,39)
  cellcondition <- data.frame(
    gene = gene,
    mouse = mouse,
    isgood = isgood
  )
  rownames(cellcondition) <- colnames(dat)
  cellcondition
}
cellcondition <- getplatedesignjhuma()




######################################################################
### Read ensembl #####################################################
######################################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]

ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart=ensembl),stringsAsFactors=FALSE)
togenesym <- function(ensid){
  out <- c()#terrible speed...
  for(n in ensid){
    ids <- which(ensconvert$ensembl_gene_id==n)
    if(length(ids)==0)
      more <- n
    else
      more <- ensconvert$mgi_symbol[ids]
    out <- c(out, more)    
  }
  #ensconvert$mgi_symbol[which(ensconvert$ensembl_gene_id %in% ensid)]
  out
}

toensid <- function(geneid){
  out <- c()#terrible speed...
  for(n in geneid){
    out <- c(out, ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol==n)])    
  }
  out
}


ensidXbp1 <- ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol=="Xbp1")] ##???
ensidErn1 <- ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol=="Ern1")] ##???

#hmmmm. TODO: isoform analysis!

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

par(cex.axis=1.5, cex.lab=1.5)
plot(exon_prop[1,], mt_prop, pch=20, xlab="Proportion of reads mapped to exons",
     ylab="Proportion of reads mapped to mitochondrial genes", xlim=c(0,1), ylim=c(0,1))
abline(h=0.1, v=0.5, col="red", lty=2, lwd=2)

#remove cells with too few reads (compare to mean for a plate instead?)
plot(gene_count)
plot(mt_prop)
cellcondition$isgood[gene_count<1e6] <- FALSE


#Plot result per plate
plot(mt_prop)


######################################################################
### DE analysis #####################################################
######################################################################


### Compute normalized counts
mynormalize <- function(dat){
  tokeep <- colSums(dat)>1e2 #& cellcondition$plate %in% c(1) #don't try to normalize empty cells
  dat666 <- dat[-which(rownames(dat) %in% listnongenes),tokeep]
  sf <- estimateSizeFactorsForMatrix(dat666)
  sfnew <- rep(1,ncol(dat))
  sfnew[(1:ncol(dat))[tokeep]]<-sf
  for(i in 1:ncol(dat)){
    dat[,i] <- dat[,i]/sfnew[i]
  }
  #print(sf)
  #print(sfnew)
  dat
}
ncount_all <- mynormalize(dat)
ncount_all <- ncount_all[-which(rownames(ncount_all) %in% listnongenes),]
log_ncount_all <- log(1+ncount_all)


mybatchremove <- function(dat){
  dat <- log_ncount_all
  tokeep <- cellcondition$isgood
  sf <- removeBatchEffect(dat[,tokeep], cellcondition$mouse[tokeep])
  sfnew <- dat
  sfnew[,(1:ncol(dat))[tokeep]]<-sf
  sfnew
}
log_newncount <- mybatchremove(log_ncount_all)
#removeBatchEffect(ncount_all, cellcondition$mouse)




#make plates
makeplatejhuma <- function(inp){
  inp <- as.double(inp)
  rbind(inp[1:13], inp[14:26], inp[27:39])
}
barplot(makeplatejhuma(ncount_all[ensidXbp1,])/mean(as.double(ncount_all[ensidXbp1,])),beside=TRUE)

  ncount_all[ensidXbp1,cellcondition$gene=="Xbp1" & cellcondition$mouse==1]
  ncount_all[ensidXbp1,cellcondition$gene=="ctrl" & cellcondition$mouse==1]
  

  ncount_all[ensidXbp1,cellcondition$gene=="Xbp1" & cellcondition$mouse==3]
  ncount_all[ensidXbp1,cellcondition$gene=="ctrl" & cellcondition$mouse==3]

mid<-3
ncount_all[ensidErn1,cellcondition$gene=="Ern1" & cellcondition$mouse==mid]
ncount_all[ensidErn1,cellcondition$gene=="ctrl" & cellcondition$mouse==mid]


    
colorbyjhumagene <- function(){
  thecol <- rep("gray",nrow(cellcondition))
  thecol[cellcondition$gene=="Xbp1"]<-"green"
  thecol[cellcondition$gene=="Ern1"]<-"blue"
  thecol[cellcondition$gene=="Irf4"]<-"black"
  thecol[cellcondition$gene=="ctrl"]<-"red"
  thecol  
}


colorbymouse <- function(){
  thecol <- rep("gray",nrow(cellcondition))
  thecol[cellcondition$mouse==1]<-"green"
  thecol[cellcondition$mouse==2]<-"blue"
  thecol[cellcondition$mouse==3]<-"red"
  thecol  
}


### PCA
#dopca <- function(){
  the_gene_mean <- rowMeans(ncount_all[,cellcondition$isgood])
  the_gene_mean <- the_gene_mean/sum(the_gene_mean)
  
  length(which(the_gene_mean>=0.00001))
  takecells <- cellcondition$isgood & cellcondition$mouse==3   & cellcondition$gene %in% c("ctrl","Xbp1")
  takecells <- cellcondition$isgood  & cellcondition$gene %in% c("ctrl","Xbp1")
  takecells <- cellcondition$isgood  & cellcondition$gene %in% c("ctrl","Ern1")
 
  # cor_ncount <- cor(log_ncount_all[the_gene_mean>=0.00001 ,takecells], method="spearman")  #better?
  cor_ncount <- cor(log_newncount[the_gene_mean>=0.00001 ,takecells], method="spearman")  #better?
  
  
  pca <- prcomp(cor_ncount, scale=FALSE)
  #pdf("pca2d_cor_50.pdf")
  par(cex.axis=1.5, cex.lab=1.5)
  #col <- colorbyjhumagene()
  col <- colorbymouse()
  plot(pca$x[,1], pca$x[,2], pch=20, col=col[takecells])#, xlab="PC1", ylab="PC2")
  
  #PC1 is mainly mouse. so is PC2
  
  #dev.off()
# }
# 
# dopca()

  

takecells <- cellcondition$isgood #& cellcondition$mouse==1   # %in% c(1) #cellcondition$isss #& !cellcondition$iswt 
geneids <- rownames(dat)[the_gene_mean>=0.00001]
runrtsne(geneids, takecells,colorbyjhumagene(),4)
runrtsne(geneids, takecells,colorbymouse(),4)



###################### deseq #########################

datfordeseq  <- dat[grep("ENSMUSG",rownames(dat)),]


keepcells <- cellcondition$isgood & cellcondition$mouse==2 & cellcondition$gene %in% c("ctrl","Ern1")
keepcells <- cellcondition$isgood & cellcondition$gene %in% c("ctrl","Ern1")
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells],
                              colData = data.frame(gene=factor(cellcondition$gene[keepcells])),
                              design = ~ gene)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
resOrdered[1:10,]
togenesym(rownames(resOrdered)[1:10])





