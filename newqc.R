##### this is for miseq 13, containing both a crispr and ss2 plate


install.packages(c("gplots","scatterplot3d","org.Mm.eg.db","GO.db","plotrix"))

source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(DESeq)
library(plyr)
library(biomaRt)
library(gplots)
library(scatterplot3d)
library(org.Mm.eg.db)
library(GO.db)
library(plotrix)

dat <- read.csv("counts/rawm13.csv",stringsAsFactors=FALSE)
dim(dat)
rownames(dat) <- dat[,1]
dat <- dat[,-1]
#dat <- dat[-which(rownames(dat) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")),]
dat <- dat[-grep("grna",rownames(dat)),]
#nongenes <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique",
#              rownames(dat)[grep("grna",rownames(dat))],rownames(dat)[grep("ERCC",rownames(dat))])


dat2 <- read.csv("counts/rawgrnam13.csv",stringsAsFactors=FALSE)
dim(dat2)
rownames(dat2) <- dat2[,1]
dat2 <- dat2[,-1]
#dat2 <- dat2[rownames(dat2) %in% c("ptma1","ptma2","ptma3","ptma4","ptma5"),]
#rownames(dat2)
rownames(dat2) <- sprintf("grna-%s",rownames(dat2))
colnames(dat2) <- colnames(dat)

dat3 <- rbind(dat,dat2)



thecols <- rep(1:12,8)
thecols <- c(thecols,thecols)
therows <- as.double(sapply(1:8,function(i) rep(i,12)))
therows <- c(therows,therows)
theplate <- as.double(sapply(1:2,function(i) rep(i,96)))



### Read ensembl #########################################################
# ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# ensembl_genes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_id',  'chromosome_name', 'gene_biotype'), mart=ensembl)
# rownames(ensembl_genes) <- ensembl_genes[,1]
# mt_genes <- ensembl_genes[ensembl_genes[,3]=="MT",]
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]
#which(ensembl_genes$gene_biotype=="Mt_rRNA")
#mt_genes

### QC for cells #####################################################

## mitochondrial content
count_cell <- dat
#dat <- dat[-which(rownames(dat) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")),]
gene_count <- colSums(count_cell[-which(rownames(count_cell) %in% c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")),])
#or pull ENS & MT?
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
#empty_debris = c(65, 71, 78, 82, 83, 84, 88, 89, 96)
#points(exon_prop[1,empty_debris],  mt_prop[empty_debris], pch=20, col="red")
abline(h=0.1, v=0.5, col="red", lty=2, lwd=2)
#dev.off()

#over the plate. higher mt content in ss2 plates. wtf?
plot(mt_prop)

plot(as.double(nomapped_count))
#<- count_cell[which(rownames(count_cell)=="not_aligned"),]

plot(colSums(dat>2))
nogene_count <- count_cell[which(rownames(count_cell)=="no_feature"),]
nomapped_count <- count_cell[which(rownames(count_cell)=="not_aligned"),]
plot(as.double(gene_count))
head(nomapped_count)

plot(as.double(colSums(dat))) ## all kind of genes. much fewer reads 
plot(as.double(colSums(dat2[-which(rownames(dat2)=="grna-numread"),]))) ##only grna. seems to make sense
plot(as.double(colSums(dat3[-which(rownames(dat3)=="grna-numread"),]))) ##combined 1 & 2 - a bit more for regular ss2

sum(as.double(colSums(dat3))[1:96])
sum(as.double(colSums(dat3))[97:192])


plot(as.double(dat3[which(rownames(dat3)=="grna-numread"),])/
  as.double(colSums(dat3))) 



### TODO detected number of genes?




### Comparing average gene expression among conditions

gene_mean_cr <- rowMeans(ncount_all[,theplate==1])
gene_mean_ss <- rowMeans(ncount_all[,theplate==2])

gene_mean_condition <- cbind(gene_mean_cr, gene_mean_ss)
colnames(gene_mean_condition) <- c("cr", "ss")
#pdf("comparison_gene_mean_conditions.pdf")
my_line <- function(x,y){
  points(x,y,pch=20, cex=.2, col="darkgray")
  abline(0, 1, col="black", lty=3, lwd=3)
}
pairs(gene_mean_condition, log="xy", panel=my_line, lower.panel=NULL)
#dev.off()


### Figure 1, PCA
cor_ncount <- cor(ncount_all[gene_mean>=50,], method="spearman")
pca <- prcomp(cor_ncount, scale=FALSE)

#pdf("pca2d_cor_50.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(pca$x[,1], pca$x[,2], pch=20, col="black", xlab="PC1", ylab="PC2")
points(pca$x[names(pca$x[,1])=="2i2",1], pca$x[names(pca$x[,2])=="2i2",2], pch=20, col="blue4")
points(pca$x[names(pca$x[,1])=="serum1",1], pca$x[names(pca$x[,2])=="serum1",2], pch=20, col="red")
points(pca$x[names(pca$x[,1])=="serum2",1], pca$x[names(pca$x[,2])=="serum2",2], pch=20, col="red4")
points(pca$x[names(pca$x[,1])=="a2i2",1], pca$x[names(pca$x[,2])=="a2i2",2], pch=20, col="darkgoldenrod2")
legend(0.9, 1, c("2i", "serum1"),
       pch=c(20, 20, 20, 20),
       col=c("blue4", "red", "red4", "darkgoldenrod2"), bty="n", cex=1.5)
#dev.off()
