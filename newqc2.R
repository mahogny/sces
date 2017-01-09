######################################################################
##### this is for the final 3 hiseq
######################################################################

#also need system libraries: libssl. libcurl?
installdeps <- function(){
  install.packages(c("gplots","scatterplot3d","plotrix","sqldf","Matrix","pi0"))
  install.packages("tsne")
  install.packages("Rtsne")
  install.packages("rPython")  #python-dev
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("biomaRt","DESeq","DEseq2"))
  biocLite(c("monocle"))
  biocLite(c("genefilter"))
  #biocLite("topGO")
  biocLite("limma")
  biocLite("scater")
  biocLite("scran")
  #mintlinux package: the mesa stuff. libx. freetype
}

library(RColorBrewer)
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
library(genefilter)
#library(topGO)
library(limma)
library(stringr)
library(reshape2)
library("scater")
library("scran")


######################################################################
### Load data ########################################################
######################################################################

listnongenes <- c(
  "no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique","grna-numread",
  "__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")

# read normal gene counts
readcounttable20161209one <- function(f){
  dat <- read.csv(f,stringsAsFactors=FALSE)
  #dim(dat)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  #  dat <- dat[-grep("grna",rownames(dat)),]
  rownames(dat)[which(rownames(dat)=="CAS9")] <- "cas9"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-1")] <- "grna-ptma1"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-2")] <- "grna-ptma2"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-3")] <- "grna-ptma3"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-4")] <- "grna-ptma4"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-5")] <- "grna-ptma5"

  dat
}
readcounttable20161209 <- function(){
  d1 <- readcounttable20161209one("counts/raw20161209_1.csv")
  d2 <- readcounttable20161209one("counts/raw20161209_2.csv")
  d3 <- readcounttable20161209one("counts/raw20161209_3.csv")
  dat <- cbind(d1,d2,d3)
  colnames(dat) <- sprintf("%s",1:ncol(dat))
  dat
}
dat <- readcounttable20161209()

### Regularized log
log1 <- function(x) as.double(log(1+x))


### Compute normalized counts
mynormalize <- function(dat){
  dat2 <- dat
  tokeep <- colSums(dat)>1e2 #& cellcondition$plate %in% c(1) #don't try to normalize empty cells
  #dat666 <- dat[-which(rownames(dat) %in% listnongenes),tokeep]  #for deseq
  dat666 <- dat[,tokeep] #for featurecount
  sf <- estimateSizeFactorsForMatrix(dat666)
  sfnew <- rep(1,ncol(dat))
  sfnew[(1:ncol(dat))[tokeep]]<-sf
  for(i in 1:ncol(dat)){
    dat2[,i] <- dat2[,i]/sfnew[i]
  }
  #print(sf)
  #print(sfnew)
  dat2
}
ncount_all <- mynormalize(dat)
#ncount_all <- ncount_all[-which(rownames(ncount_all) %in% listnongenes),] #for htseq output. not used anymore
log_ncount_all <- log1(ncount_all)


### Construct plate design
getplatedesign20161209 <- function(){
  #TODO get rid of one plate on hiseq3 completely! and one is empty as well but can ignore this
  
  #hiseq1 & hiseq2: swapped B,C,D. 
  
  numplates <- 4+4+4
  
  thecols <- rep(1:12,8)
  thecols <- rep(thecols,numplates)
  therows <- as.double(sapply(1:8,function(i) rep(i,12)))
  therows <- rep(therows,numplates)
  theplate <- as.double(sapply(1:numplates,function(i) rep(i,96)))
  level1<-thecols<12 & therows %in% c(1,2) & theplate<=10
  level2<-therows %in% c(3,4) & theplate<=10
  level3<-therows %in% c(5,6) & theplate<=10
  level4<-therows %in% c(7,8) & theplate<=10
  levelnum <- level1*1 + level2*2 + level3*3 + level4*4
  isgood <- levelnum>0
  isgood[theplate==1 & thecols>=7] <- FALSE #these are empty
  iscr <- theplate %in% c(4, 8, 9,10)
  isss <- theplate %in% c(1,2,3,  5,6,7)
  
  ### Batch groups for limma
  batchgroup <- rep(0,ncol(dat))
  batchgroup[theplate==1] <- 1
  batchgroup[theplate==2 & thecols<=6] <- 2
  batchgroup[theplate==2 & thecols>=6] <- 3
  batchgroup[theplate==3 & thecols<=6] <- 4
  batchgroup[theplate==3 & thecols>=6] <- 5
  batchgroup[theplate==4 & thecols<=6] <- 20
  batchgroup[theplate==4 & thecols>=6] <- 21
  batchgroup[theplate==5 & thecols<=6] <- 2
  batchgroup[theplate==5 & thecols>=6] <- 3
  batchgroup[theplate==6 & thecols<=6] <- 4
  batchgroup[theplate==6 & thecols>=6] <- 5
  batchgroup[theplate==7 & thecols<=6] <- 6
  batchgroup[theplate==7 & thecols>=6] <- 7
  batchgroup[theplate==8 & thecols<=6] <- 22
  batchgroup[theplate==8 & thecols>=6] <- 23
  batchgroup[theplate==9] <- 24
  batchgroup[theplate==10] <- 25

  fromseq <- c(rep(1,384),rep(2,384),rep(3,384)) 
  

    
  isfeeder <- rep(FALSE,ncol(dat))
  isfeeder[theplate %in% c(1,2,3,4,5,6,7,8,9,10) & thecols<=6 & isgood] <- TRUE
  isfeeder[theplate %in% c(1) & isgood] <- TRUE
  
  islif <- rep(FALSE,ncol(dat))
  islif[theplate %in% c(1,2,3,4,5,6,7,8,9,10) & thecols>=7 & isgood] <- TRUE
  islif[theplate %in% c(1) & isgood] <- TRUE
  
  cellcondition <- data.frame(
    col = thecols,
    row = therows,
    plate = theplate,
    iscr = iscr,
    isss = isss,
    level1=level1, level2=level2, level3=level3, level4=level4,
    levelnum=levelnum,
    isgood=isgood,
    batchgroup=batchgroup,
    fromseq=fromseq,
    isfeeder=isfeeder,
    islif=islif
  )
  rownames(cellcondition) <- colnames(dat)
  cellcondition
}
cellcondition <- getplatedesign20161209()

######################################################################
### Coloring methods #################################################
######################################################################

colorbygene <- function(ensid=ensidPtma,keep=rep(TRUE,ncol(dat))){
  thecol <- rep("red",nrow(cellcondition))
  v <- log_ncount_all[ensid,]
  rgb(0, v/max(v[keep]), 0, maxColorValue=1)
}
colorbylevel <- function(){
  thecol <- rep("red",nrow(cellcondition))
  thecol[cellcondition$level1]<-"#000000"
  thecol[cellcondition$level2]<-"#005500"
  thecol[cellcondition$level3]<-"#00AA00"
  thecol[cellcondition$level4]<-"#00FF00"
  thecol  
}
colorbychem <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$iscr]<-"red"
  thecol[cellcondition$isss]<-"blue"
  thecol  
}
colorbyplate <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$plate==1]<-"red"
  thecol[cellcondition$plate==2]<-"green"
  thecol[cellcondition$plate==3]<-"blue"
  #thecol[cellcondition$isss]<-"red"
  thecol
}
colorbyseq <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$fromseq==1]<-"green"
  thecol[cellcondition$fromseq==2]<-"blue"
  thecol[cellcondition$fromseq==3]<-"red"
  thecol  
}
colorbymedia <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$isfeeder==1]<-"green"
  thecol[cellcondition$islif==2]<-"blue"
  thecol  
}

######################################################################
### Read ensembl #####################################################
######################################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]

#Convert to gene symbols, or retain ID if no symbol
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

## Convert to ensembl IDs
toensid <- function(geneid){
  out <- c()#terrible speed...
  for(n in geneid){
    out <- c(out, ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol==n)])    
  }
  out
}

ensidPtma <- toensid("Ptma")#ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol=="Ptma")]


######################################################################
### Read olas data ###################################################
######################################################################
#ola figure 7d. genes repressed by PTMA
listolaptmagenes7d <- c("Vcan","Dst","Peg10","Parp12","Egr1","Nqo1","Sfrp1","Clptm1l","Ptp4a3","Psme3","Tcf3","Ung","Tomm70a","Eif2b5","Ywhah","Fgf4")
listolaptmagenes7d <- data.frame(
  genesym=listolaptmagenes7d, 
  ensid=rep("",length(listolaptmagenes7d)),stringsAsFactors=FALSE)
for(i in 1:nrow(listolaptmagenes7d)){
  listolaptmagenes7d$ensid[i] <- ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol==listolaptmagenes7d$genesym[i])]
}
colnames(listolaptmagenes7d) <- c("genesym","ensid")


#ola pluri-pot genes, taken from literature
listolapluripotency <- read.csv("ola_pluripotency_genes.csv",stringsAsFactors = FALSE)
colnames(listolapluripotency) <- c("ensid","genesym")


#also taken from literature I think
listolacellcycle <- read.csv("ola_cell_cycle_fig_genenames.csv",sep="\t",stringsAsFactors = FALSE)
colnames(listolacellcycle) <- c("ensid","genesym")
listola2c <- read.csv("ola_2C_genes_subset.csv",sep="\t",stringsAsFactors = FALSE)[,c(1,2)]
colnames(listola2c) <- c("ensid","genesym")


######################################################################
### QC for cells #####################################################
######################################################################

## mitochondrial content
count_cell <- dat
gene_count <- colSums(count_cell[grep("ENSMUSG",rownames(count_cell)),])
nogene_count <- count_cell[which(rownames(count_cell)=="no_feature"),]
nomapped_count <- count_cell[which(rownames(count_cell)=="not_aligned"),]
grna_count <- colSums(count_cell[grep("grna",rownames(count_cell)),])
exon_prop <- rbind(gene_count, nogene_count, nomapped_count)
exon_prop <- t(t(exon_prop) / colSums(exon_prop))
mt_counts <- colSums(count_cell[mt_genes$ensembl_gene_id,])
mt_prop <- mt_counts / gene_count  



#Plot result per plate
plot(grna_count)
plot(mt_prop)


#pdf("qc_exonprop_mtprop_cell_2i_2.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(exon_prop[1,], mt_prop, pch=20, xlab="Proportion of reads mapped to exons",
     ylab="Proportion of reads mapped to mitochondrial genes", xlim=c(0,1), ylim=c(0,1))
abline(h=0.1, v=0.5, col="red", lty=2, lwd=2)
#dev.off()



#remove cells with too few reads (compare to mean for a plate instead?)
png("plots/readcounts.png",width=800)
par(mfrow=c(1,2))
plot(gene_count[cellcondition$isss],main="SS read count")
plot(gene_count[cellcondition$iscr],main="DogSeq read count")
dev.off()

# plot(gene_count[cellcondition$plate %in% c(4,8,9,10)],ylim=c(0,200000))
# plot(gene_count[cellcondition$plate %in% c(4+4+4)],ylim=c(0,200000))
# length(which(gene_count[cellcondition$plate %in% c(4,8,9,10)]>100e3))
# plot(gene_count[cellcondition$plate %in% c(1,2,3,4)])
# plot(gene_count[cellcondition$plate %in% c(5,6,7,8)])
# plot(gene_count[cellcondition$plate %in% c(6)])
# mean(gene_count[cellcondition$isss])
# mean(gene_count[cellcondition$plate %in% c(4)])
# mean(gene_count[cellcondition$plate %in% c(8)])
#cellcondition$isss[which(order(gene_count)<96+10)] <- FALSE

#plot(gene_count[cellcondition$plate %in% c(4,8,9,10,11)],ylim=c(0,1e6))

length(which(gene_count[cellcondition$plate %in% c(4,8,9,10,11)]>20e3))

#remove bad cells
#cellcondition$isgood[cellcondition$isss & gene_count<1e5] <- FALSE

cellcondition$isgood[gene_count<100e3 & cellcondition$isss] <- FALSE #minimum exonic reads for regular chemistry


#mean(gene_count[cellcondition$isss])
#cellcondition$isss[which(order(gene_count)<96+10)] <- FALSE





### TODO detected number of genes?



### Comparing average gene expression among conditions


plotgenemean <- function(){
  gene_mean_cr <- rowMeans(ncount_all[,cellcondition$iscr]) 
  gene_mean_ss <- rowMeans(ncount_all[,cellcondition$isss])
  gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)
  gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
  
  gene_mean_condition <- cbind(gene_mean_cr, gene_mean_ss)
  colnames(gene_mean_condition) <- c("cr", "ss")
  #pdf("comparison_gene_mean_conditions.pdf")
  my_line <- function(x,y){
    points(x,y,pch=20, cex=.2, col="darkgray")
    abline(0, 1, col="black", lty=3, lwd=3)
  }
  pairs(gene_mean_condition, log="xy", panel=my_line, lower.panel=NULL)
  #dev.off()
}


### PCA for all expressed genes
dopca <- function(){
  gene_mean_cr <- rowMeans(ncount_all[,cellcondition$iscr]) 
  gene_mean_ss <- rowMeans(ncount_all[,cellcondition$isss])
  gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)
  gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
  
  length(which(gene_mean_ss>=0.00001))
  #cor_ncount <- cor(ncount_all[gene_mean_ss>=0.00001,], method="spearman")
#  takecells <- cellcondition$plate %in% c(1) & cellcondition$isss & cellcondition$isgood #& !cellcondition$iswt 
#  takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood #& !cellcondition$iswt 
  takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood &cellcondition$isfeeder  #& !cellcondition$iswt 
  takecells <- cellcondition$plate %in% c(1) & cellcondition$isss & cellcondition$isgood &cellcondition$isfeeder  #& !cellcondition$iswt 
  
  takecells[weird ] <- FALSE
  cor_ncount <- cor(ncount_all[gene_mean_ss>=0.00001 ,takecells], method="spearman")  #olas way
  #cor_ncount <- t(ncount_all[gene_mean_ss>=0.00001 ,takecells])  #ola does not scale
  pca <- prcomp(cor_ncount, scale=FALSE)
  #pdf("pca2d_cor_50.pdf")
  par(cex.axis=1.5, cex.lab=1.5)
  plot(pca$x[,3], pca$x[,2], pch=20, col=colorbylevel()[takecells], xlab="PC", ylab="PC")
  plot(pca$x[,1], pca$x[,2], pch=20, col=colorbyplate()[takecells], xlab="PC", ylab="PC")
#  weird <- (1:ncol(dat))[takecells][pca$x[,1]>2]
  #dev.off()
}



######################################################################
### LIMMA batch effect removal and QC ################################
######################################################################


#isgood <- cellcondition$isgood



### PCA
dopca2 <- function(){
  
  the_gene_mean <- rowMeans(ncount_all[,cellcondition$isgood & cellcondition$isss])
  the_gene_mean <- the_gene_mean/sum(the_gene_mean)
  
  length(which(the_gene_mean>=0.00001))
  #cor_ncount <- cor(ncount_all[gene_mean_ss>=0.00001,], method="spearman")
  takecells <- cellcondition$isss & cellcondition$isgood  # %in% c(1) #cellcondition$isss #& !cellcondition$iswt 
#  cor_ncount <- cor(ncount_all[gene_mean_ss>=0.000001 ,takecells], method="spearman")  #olas way
  cor_ncount <- cor(log_ncount_all[the_gene_mean>=0.00001 ,takecells], method="spearman")  #better?
  #cor_ncount <- t(ncount_all[gene_mean_ss>=0.00001 ,takecells])  #ola does not scale
  pca <- prcomp(cor_ncount, scale=FALSE)
  #pdf("pca2d_cor_50.pdf")
  par(cex.axis=1.5, cex.lab=1.5)
  plot(pca$x[,1], pca$x[,2], pch=20, cex=.4,col=colorbyseq()[takecells], xlab="PC1", ylab="PC2")
  plot(pca$x[,1], pca$x[,2], pch=20, cex=.4,col=colorbymedia()[takecells], xlab="PC1", ylab="PC2")
  plot(pca$x[,1], pca$x[,2], pch=20, cex=.8,col=colorbylevel()[takecells], xlab="PC1", ylab="PC2")
  #dev.off()
  
  clust <- kmeans(pca$x,centers=2) #or k-means on original data?
  ccol <- rep("blue",length(which(takecells)))
  ccol[clust$cluster==2] <- "red"
  plot(pca$x[,1], pca$x[,2], pch=20, cex=.8,col=ccol, xlab="PC1", ylab="PC2")
  
  g1<-log_ncount_all[ ,takecells][,clust$cluster==1]
  g2<-log_ncount_all[ ,takecells][,clust$cluster==2]

  df <- data.frame(p.value=rep(666,nrow(log_ncount_all)), x=rep(1,nrow(log_ncount_all)), y=rep(1,nrow(log_ncount_all)))
  for(i in (1:nrow(log_ncount_all))[the_gene_mean>=0.001]){
    t <- t.test(g1[i,],g2[i,])
    df[i,1] <- t$p.value
    df[i,2] <- t$estimate[1]
    df[i,3] <- t$estimate[2]
  }
#  df<-data.frame(pval=outt,fold=outf)  
  rownames(df)<-rownames(log_ncount_all)
  df <- df[order(df$p.value,decreasing = FALSE),]
  df[1:50,]
  data.frame(n=togenesym(rownames(df)[1:50]),p=df$pval[1:50],fold=df$fold[1:50])
}

######################################################################
### LIMMA model testing ##############################################
######################################################################
#need to define groups to normalize together. can make something called a subplate. and also define feeder vs nonfeeder
#limma: on log level, normalize . removebatch effect


#test genes different between media conditions

batchvector <- rep(0,ncol(dat))
batchvector[!cellcondition$isss] <- 1
#batch_log_ncount <- removeBatchEffect(log_ncount_all, batchvector)

# 
# TODO: order new TSO!
# 
#   QC according to vale: number of raeds mapped
# xi: write paper as much as possible
# TODO align reads again




######################################################################
### Comparison with olas previous data ###############################
######################################################################



### Average gene levels for cells instead of bulk
# genemeanperlevel <- cbind(
#   rowMeans(ncount_all[,cellcondition$level1]),
#   rowMeans(ncount_all[,cellcondition$level2]),
#   rowMeans(ncount_all[,cellcondition$level3]),
#   rowMeans(ncount_all[,cellcondition$level4]))
genemeanperlevel <- cbind(
  rowMeans(log_ncount_all[,cellcondition$level1]),
  rowMeans(log_ncount_all[,cellcondition$level2]),
  rowMeans(log_ncount_all[,cellcondition$level3]),
  rowMeans(log_ncount_all[,cellcondition$level4]))

#### Redoing fig7d. What is down-regulated by PTMA?
andid <- c(ensidPtma)#,toensid(c("Narf","Cdc45")))
b <- t(genemeanperlevel[
  rownames(genemeanperlevel) %in% c(listolaptmagenes7d$ensid,andid),
  c(1,2,3,4)])
#use order as in olas paper. ola sorts by fold change
b<-b[,c(listolaptmagenes7d$ensid, andid)] 
for(i in 1:ncol(b))
  b[,i] <- b[,i]/max(b[,i]) #normalize to highest level
colnames(b) <- togenesym(colnames(b)) #order preserving? not any more!
png("plots/ss2 olaprevtargets ptma.png",width = 1000,pointsize = 16)
par(cex.lab=0.2)
par(mfrow=c(1,1))
barplot(b, main="Genes per PTMA facs level",
        col=c("#000000","#880000","#AA0000","#FF0000"), cex.names=0.5,beside=TRUE) 
dev.off()

#TODO so correlation not great by FACS level. but better to reproduce plot by Ptma level then?
#TODO by ptma level. then by cherry level. just split into 4 bins. also only use SS?

#can also sort each gene exp by Ptma level and make a dot-plot?
dat2 <- ncount_all[,cellcondition$isgood & cellcondition$isss]
dat2 <- dat2[,order(dat2[ensidPtma,])]
plot(as.double(log1(dat2[listolaptmagenes7d$ensid[15],])))
#2,10,15 very bimodal  ... this is actually much clearer this way than below? would have to fit a sigmoidal.
#then two pieces of information; change in level and bimodality
plot(
  as.double(log1(dat2[ensidPtma,])),
  as.double(log1(dat2[listolaptmagenes7d$ensid[15],]))
  )
plot(
  as.double((dat2[ensidPtma,])),
  as.double(log1(dat2[listolaptmagenes7d$ensid[15],]))
)




#ptma goes down as it should
#Ptp4a3 looks nice! and Dst and Peg10 and Parp12*. but peg10 & dst the opposite way before
#Parp12 seems affected by NfkB or interferon in general. could be affeced by the transfection itself


### Pull out genes correlated with Ptma - by correlation
nametogenesym <- function(x){
  names(x) <- togenesym(names(x))
  x
}
#length(which(gene_mean_ss>=0.0001))
takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood  & cellcondition$levelnum>1
#takecells[weird ] <- FALSE
cor_genes <- cor(t(ncount_all[gene_mean_ss>=0.00001 | rownames(ncount_all)==ensidPtma ,takecells]), method="spearman")  #olas way
cor_ptma <- cor_genes[rownames(cor_genes)==ensidPtma]
names(cor_ptma) <- colnames(cor_genes)
nametogenesym(sort(cor_ptma)[1:30])
nametogenesym(sort(cor_ptma,decreasing = TRUE)[1:10])
plot(sort(cor_ptma,decreasing = TRUE))

##for ss 0.00001
# Pglyrp3    Gm10719       Lmo3    mcherry    Gm11168    Gm14762    Gm10801       Ccr1      Mss51    Gm26870 
# -0.4693054 -0.4404796 -0.4266743 -0.3945823 -0.3938720 -0.3822966 -0.3460647 -0.3308073 -0.3271836 -0.3176716 
# Ptma    Gm9800    Gm4617       Ncl       Erh     Hmgb2       Ran      Pcna     Psma5      Ppa1 
# 1.0000000 0.9623600 0.9419038 0.6859848 0.6818732 0.6742220 0.6652018 0.6650973 0.6635739 0.6607904 


######################################################################
### Plot gene vs gene ################################################
######################################################################

#Check some individual correlations
plotcor <- function(x,y){
  plot(x,y,cex.lab=1)
  cor(x,y,method="spearman")
}

keepcells <- cellcondition$isss & gene_count>300e3 & cellcondition$levelnum>0
#keepcells <- cellcondition$isss & cellcondition$levelnum>0
#keepcells <- cellcondition$iscr & gene_count>1e3 & cellcondition$levelnum>0

keepcells <- cellcondition$isss & gene_count>300e3 & cellcondition$levelnum>0

#keap1 supposed to interact with ptma. there could be a trend. keap1+ptma = degrade Nrf2? meah?
plot(log1(ncount_all[toensid("Nfe2l2"),keepcells]), #literature suggests interaction nfe2l2. https://en.wikipedia.org/wiki/NFE2L2
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount_all[toensid("Nfe2l1"),keepcells]), 
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)


plot(log1(ncount_all[toensid("Myc"),keepcells]),  #ptma in all tissues with myc. but clearly also without myc. 
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)

#note: if Gm4617 correlates so well with ptma then it can be used as a proxy for the cr2 data too

png("plots/ss2 ptma vs Gm9800.png",width = 800, height=600)
plot(log1(ncount_all[toensid("Gm9800"),keepcells]),  #super correlated
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 ptma vs Gm4617.png",width = 800, height=600)
plot(log1(ncount_all[toensid("Gm4617"),keepcells]),  #super correlated. predicted homolog of ptma http://www.ihop-net.org/UniPub/iHOP/gs/123843.html
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()
plot(log1(ncount_all[toensid("Ncl"),keepcells]), #https://en.wikipedia.org/wiki/Nucleolin  ribosomes
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount_all[toensid("Erh"),keepcells]), #suggested to be in cell cycle. and other things
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
png("plots/ss2 ptma vs Pglyrp3.png",width = 800, height=600)
plot(log1(ncount_all[toensid("Pglyrp3"),keepcells]),  #innate immune system. weak?
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()
plot(log1(ncount_all[toensid("Lmo3"),keepcells]),  #LIM domain something. weak?
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount_all[toensid("Prdx1"),keepcells]),  #https://en.wikipedia.org/wiki/Peroxiredoxin_1  strong. after release, causes cytokines to be released
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)

plot(log1(colSums(ncount_all[grep("grna-ptma",rownames(ncount_all)),keepcells])),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1,xlim=c(0,1.5))
plot(log1(ncount_all["grna-ptma2",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)



plot(log1(ncount_all[toensid("Npm1"),keepcells]),  #
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount_all[toensid("Hsp90ab1"),keepcells]),  #came out of deseq strongly
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)

#TODO get all the articles on ptma

#TODO deseq, high vs low ptma. and +- mcherry! level>1. and +- cas9

#ptma https://en.wikipedia.org/wiki/Thymosin_%CE%B11  - could it actually connect to pglyrp3

plot(log1(ncount_all[toensid("Egr1"),keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)


plot(log1(ncount_all[toensid("Ptp4a3"),keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)


######################################################################
### Plot gene vs gene, mcherry & cas9 ################################
######################################################################

keepcells <- cellcondition$isss & gene_count>300e3 & cellcondition$levelnum>0

# png("plots/ss2 ptma vs bfp.png",width = 800, height=600)
# plot(log1(ncount_all["bfp",keepcells]),
#      log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
# dev.off()
png("plots/ss2 ptma vs mcherry.png",width = 800, height=600)
plot(log1(ncount_all["mcherry",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 ptma vs cas9.png",width = 800, height=600)
plot(log1(ncount_all["cas9",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()

keepcells2 <- cellcondition$isss & gene_count>300e3 & cellcondition$levelnum>2
plot(log1(ncount_all["cas9",keepcells2]),
     log1(ncount_all[ensidPtma,keepcells2]),cex.lab=1)

#TODO a linear model with cas9 & mcherry vs ptma


# png("plots/ss2 bfp vs cas9.png",width = 800, height=600)
# plot(log1(ncount_all["bfp",keepcells]),
#      log1(ncount_all["cas9",keepcells]),cex.lab=1)
# dev.off()
png("plots/ss2 mcherry vs cas9.png",width = 800, height=600)
plot(log1(ncount_all["mcherry",keepcells]),
     log1(ncount_all["cas9",keepcells]),cex.lab=1) #good correlation
dev.off()

plot(log1(ncount_all["bfp",]))
plot(log1(ncount_all["mcherry",]))

plot(log1(ncount_all["cas9",keepcells]))
#plot(log1(ncount_all["bfp",keepcells]))
plot(log1(ncount_all["mcherry",keepcells]))

#check correlation bfp, cas9 and read count!
png("plots/ss2 cas9 vs readcount.png",width = 800, height=600)
plot(log1(ncount_all["cas9",keepcells]),
     log1(gene_count[keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 bfp vs readcount.png",width = 800, height=600)
plot(log1(ncount_all["bfp",keepcells]),
     log1(gene_count[keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 mcherry vs readcount.png",width = 800, height=600)
plot(log1(ncount_all["mcherry",keepcells]),
     log1(gene_count[keepcells]),cex.lab=1)
dev.off()





#####################################

keepcells <- cellcondition$iscr & gene_count>10e3 & cellcondition$levelnum>0
length(which(keepcells))


# png("plots/1kcr ptma vs bfp.png",width = 800, height=600)
# plot(log1(ncount_all["bfp",keepcells]),
#      log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
# dev.off()
png("plots/1kcr ptma vs mcherry.png",width = 800, height=600)
plot(log1(ncount_all["mcherry",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()
png("plots/1kcr ptma vs cas9.png",width = 800, height=600)
plot(log1(ncount_all["cas9",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
dev.off()

###############

keepcells <- cellcondition$iscr & gene_count>20e3 & cellcondition$levelnum>0
#keepcells <- cellcondition$iscr & gene_count>200e3 & cellcondition$levelnum>0

plot(log1(colSums(ncount_all[grep("grna-ptma",rownames(ncount_all)),keepcells])),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
plot(log1(colSums(ncount_all[grep("grna-ptma",rownames(ncount_all)),keepcells])) + log1(ncount_all["mcherry",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount_all["grna-ptma3",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount_all["mcherry",keepcells]),
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
##gRNA level is as good a predictor as mcherry for ptma KO. cells > 200k reads.

plot(log1(ncount_all[toensid("Gm9800"),keepcells]),  #super correlated
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount_all[toensid("Gm4617"),keepcells]),  #super correlated. predicted homolog of ptma http://www.ihop-net.org/UniPub/iHOP/gs/123843.html
     log1(ncount_all[ensidPtma,keepcells]),cex.lab=1)

log_metaptma <- 
  log1(ncount_all[toensid("Gm9800"),keepcells])*2 + 
  log1(ncount_all[toensid("Gm4617"),keepcells])*2 +
  log1(ncount_all[ensidPtma,keepcells])


#grna and mcherry as good predictors for PTMA level. but only when detecting the grna
plotcor(log1(colSums(ncount_all[grep("grna-ptma",rownames(ncount_all)),keepcells])), 
     log_metaptma)
plotcor(log1(ncount_all["mcherry",keepcells]),
     log_metaptma)


###############################################################################################

pheatmapbylevel <- function(geneids, keepcells=cellcondition$levelnum>0){
  #keepcells=cellcondition$levelnum>0
  #geneids <- listolaptmagenes7d$ensid
    
  lcounts <- log_ncount_all
  selectgenes <- which(rownames(lcounts) %in% geneids)
  
  lcounts <- lcounts[selectgenes,keepcells]
  rownames(lcounts) <- togenesym(rownames(lcounts))
  annotation_col <- data.frame(levelnum=factor(sprintf("level%s",cellcondition$levelnum[keepcells])))  
  #annotation_col
  annotation_colors <- list(levelnum=c(
    level1="#000000",level2="#005500",level3="#00AA00",level4="#00FF00"))
  #level1="#00FF00",level2="#000000",level3="#0000FF",level4="#AA0000"))
  pheatmap(
    #assay(dds)[select,],
    lcounts, 
    annotation_col=annotation_col,#gp = gpar(fill = "green")
    annotation_colors=annotation_colors,
    cluster_rows=TRUE, show_rownames=TRUE,
    cluster_cols=TRUE, show_colnames=TRUE)
  #todo: why is it not using the levels??
  #library(grid) #might be needed?
  #colnames(log2.norm.counts)
}

#  alsogenes <- c("Ptma")#sprintf("grna-ptma%s",1:5) #"cas9"
#                         c(ensidPtma,listolaptmagenes7d$ensid,alsogenes))
#"cas9"

# listannoyingcells <- c(126,109,165,144,134,120,131)
# cellcondition[listannoyingcells,]
# cellcondition$isss[listannoyingcells]<-FALSE

keepcells <- cellcondition$isss & gene_count>500e3 & cellcondition$levelnum>0
keepcells <- cellcondition$iscr & gene_count>10e3 & cellcondition$levelnum>0

pheatmapbylevel(c(ensidPtma,listolaptmagenes7d$ensid),keepcells)
pheatmapbylevel(c(ensidPtma,listolapluripotency$ensid),keepcells)

pheatmapbylevel(c(ensidPtma,listolacellcycle$ensid),keepcells)
pheatmapbylevel(c(listolacellcycle$ensid),keepcells)
pheatmapbylevel(c(ensidPtma,listola2c$ensid),keepcells)
#some really striking cells here for cell cycle!
#listolapluripotency$ensid

#which(order(gene_count)<96+40)


### todo: a k-means function as well. but can also sum up the gene exp for a list of them here
# isss: remove the 2 upper wells. same for iscr. actually combine with isgood?
# then plot. 9 ss cells show up low in every metric

#what is a good way of keeping track of cell definitions?

pheatmapcor <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount_all))){
  
  lcounts <- log_ncount_all
  selectgenes <- which(rownames(lcounts) %in% geneids)
  lcounts <- lcounts[selectgenes,keepcells]
  rownames(lcounts) <- togenesym(rownames(lcounts))
  
  thec <- cor(t(lcounts),method="spearman")
  
  #remove any N/A
  for(i in 1:ncol(thec)){
    ti <- is.na(thec[,i])
    if(length(ti)>0)
      thec[ti,i]  <-0
  }
  
  pheatmap(
    #assay(dds)[select,],
    thec, 
    cluster_rows=TRUE, show_rownames=TRUE,
    cluster_cols=TRUE, show_colnames=TRUE)
}
pheatmapcor(c(ensidPtma,listolapluripotency$ensid),cellcondition$isss)
pheatmapcor(c(ensidPtma,listolaptmagenes7d$ensid),cellcondition$isss)

#Ptma is definitely not involved in cell cycle
pheatmapcor(c(ensidPtma,listolacellcycle$ensid),cellcondition$isss)



######################################################################
### t-sne ###############################
######################################################################


#faraz:  6832
#cycles: 1ug: 


runrtsne <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount_all)),perplexity=30,max_iter=1000){
  #  geneids <- c(ensidPtma,listolacellcycle$ensid)
  #  keepcells <- cellcondition$isss & cellcondition$levelnum %in% c(1,4)
  lcounts <- log_ncount_all
  selectgenes <- which(rownames(lcounts) %in% geneids)
  lcounts <- lcounts[selectgenes,keepcells]
  #rownames(lcounts) <- togenesym(rownames(lcounts))
  
  d <- stats::dist(t(as.matrix(lcounts)))
  set.seed(0) # tsne has some stochastic steps (gradient descent) so need to set random 
  rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=perplexity, verbose = FALSE, max_iter=max_iter)
  
  #function to color by gene exp? ptma? grna? read count?
  
  # https://hms-dbmi.github.io/scw/analysis-of-heterogeneity-and-subpopulations.html
  rtsne_out$keepcells<-keepcells
  rtsne_out 
}
plottsne<-function(rtsne_out,col=colorbylevel(),title=""){
  plot(rtsne_out$Y, col=col[rtsne_out$keepcells], pch=16, main='',xlab=title,ylab="")
#  plot(rtsne_out$Y, col=col[keepcells], pch=16, main='')
}

#runrtsne(c(ensidPtma,listolacellcycle$ensid),cellcondition$isss)


#lots of blue in one group. all PTMA downreg in one group. controls 2c!
png("plots/tsne ss 2c.png")
x<-runrtsne(c(listola2c$ensid),cellcondition$isss, perplexity=30)  #3 clusters with perp30
#TODO: pull out cluster of cell for proper analysis? maybe no need. but can use this to define PTMA cutoff!
#ALSO: what is up with the stretch of cells going to the big cluster outside?
par(mfrow=c(1,2))
plottsne(x,col=colorbygene(),"Ptma")
plottsne(x,col=colorbygene("mcherry"),"Mcherry")
title("tSNE 2c",outer=TRUE)
dev.off()
#x<-runrtsne(c(listola2c$ensid),cellcondition$isgood, col=colorbychem())   #clearly separates by chemistry!

x<-runrtsne(c(listola2c$ensid),cellcondition$isss, col=colorbygene("mcherry"),perplexity = 30)  

#lots of blue in one group. all PTMA downreg in one group. controls pluripotency!
png("plots/tsne pluripot.png")
par(mfrow=c(1,2))
x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood & cellcondition$isss)  
plottsne(x,col=colorbygene(),"Ptma SS")
#plottsne(x,col=colorbygene("mcherry"),"mcherry SS")
x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood & cellcondition$iscr)
plottsne(x,col=colorbygene(),"Ptma DogSeq")
title("tSNE pluripotency",outer=TRUE)
dev.off()

png("plots/tsne all pluripot.png") #TODO - separates by chemistry
x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood)  
plottsne(x,col=colorbychem(),"Chemistry DogSeq")
#x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood, col=colorbychem())  
dev.off()


#lots of blue in one group
#PTMA level is spread over cell cycle. good! so not fluctuating then
png("plots/tsne ss cellcyc.png")
x<-runrtsne(c(listolacellcycle$ensid),cellcondition$isss)  
dev.off()

oddcells <- ((1:192)[cellcondition$isss])[which(x$Y[,1]>15)]
tokeep <- cellcondition$isss
tokeep[oddcells]<-FALSE
cellcondition[oddcells,] #mainly level2
runrtsne(c(listolacellcycle$ensid),tokeep)  #random
runrtsne(c(listola2c$ensid),tokeep)   #random
runrtsne(c(listola2c$ensid),tokeep)  #random
runrtsne(c(ensidPtma,listola2c$ensid),tokeep)  #structure!

#infection markers. what about interferon gamma?

colorbylevel


# runtsne <- function(){
#   # initialize counter to 0
#   x <- 0
#   epc <- function(x) {
#     x <<- x + 1
#     filename <- paste("d:\\plot", x, "jpg", sep=".")
#     cat("> Plotting TSNE to ", filename, " ")
#     
#     # plot to d:\\plot.x.jpg file of 2400x1800 dimension
#     jpeg(filename, width=2400, height=1800)
#     
#     plot(x, t='n', main="T-SNE")
#     text(x, labels=rownames(mydata))
#     dev.off()
#   }
#   
#   # run tsne (maximum iterations:500, callback every 100 epochs, target dimension k=5)
#   tsne_data <- tsne(mydata, k=5, epoch_callback=epc, max_iter=500, epoch=100)
# }



# monocle? can we find a progression?



################################################
##########################################################



#Show count "inp" for plate i
makeplate <- function(inp, platei){
  plateA <- matrix(0,ncol=12,nrow=8)
  for(i in 1:12){
    for(j in 1:8){
      cellcondition$col
      thei <- which(cellcondition$col==i & cellcondition$row==j & cellcondition$plate==platei)
      plateA[j,i] <- inp[thei]
    }
  }
  plateA
}
#plateA <- makeplate(colSums(dat[grep("grna-",rownames(dat)),]), 1)
#plateB <- makeplate(colSums(dat[grep("grna-ptma",rownames(dat)),]), 8)

        
#color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

#pheatmap(plateA,cluster_rows = FALSE, cluster_cols = FALSE, main = "DogSeq: gRNA count")

##############################
#TODO should normalize!!!! that is why there is a plate differenc. depth only.
#############

#for plate 1,2,3  5,6,7
plateA <- makeplate(colSums(dat[grep("grna-",rownames(dat)),]), 7)
breaksList = seq(0, 10, by = 1)
pheatmap(plateA,cluster_rows = FALSE, cluster_cols = FALSE, main = "ss2: gRNA count",
         colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(length(breaksList)),
         breaks = breaksList)
#plate 1. not so much. mainly left
#plate 2, left half lots. right very random
#plate 3. crazy!
#p5: fair, even
#p6: a lot on the right side. fits with p2! so it is from cleaning in the first round?
#p7: crazy. how does this make sense vs p3??

#Suitable for plate 4,8,9,10
platei <- 10
plateA <- makeplate(colSums(dat[grep("grna-",rownames(dat)),]), platei)
breaksList = seq(0, 200, by = 1)
pheatmap(plateA,cluster_rows = FALSE, cluster_cols = FALSE, main = sprintf("DogSeq: gRNA count, pool %s", platei),
         colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(length(breaksList)),
         breaks = breaksList)


#plate 4 looks ok! 
#plate 8 truncated to the left - but this means feeder or lif has it all, but not the other. why??
#cr: plate 9 looks great. plate 10 also great. one well has a ridiculous grna count(!).
#t-cell: limited number of grnas on plate.

platei <- 10
plateA <- makeplate(colSums(dat[grep("grna-ptma3",rownames(dat)),]), platei)
breaksList = seq(0, 40, by = 1)
pheatmap(plateA,cluster_rows = FALSE, cluster_cols = FALSE, main = sprintf("DogSeq: gRNA count, pool %s", platei),
         colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(length(breaksList)),
         breaks = breaksList)



################################# fluorescent proteins vs grna vs cas9 ####################################

#we have enough data to show grna and mcherry correlation
png("plots/cr20k mcherry vs total grna.png",w=800)
keepcells <- cellcondition$iscr & gene_count>20e3 & cellcondition$levelnum>0
plot(log1(ncount_all["mcherry",keepcells]), #good correlation but up to 10-100 fold difference
     log1(colSums(ncount_all[grep("grna-",rownames(ncount_all)),keepcells])),cex.lab=1)
dev.off()
#cas9 vs bfp is crap. this is a T2A construct. how can correlation be so bad?? likewise in their droplet data?
png("plots/cr20k bfp vs cas9.png",w=800)
keepcells <- cellcondition$iscr & gene_count>20e3 & cellcondition$levelnum>0
plot(log1(ncount_all["bfp",keepcells]), 
     log1(ncount_all["cas9",keepcells]),cex.lab=1)
dev.off()



#sort out the grna counts

#corr: what if we see 0-points as truncated points on the log scale?

### cr-data: what about trying to impute all genes which are zero?


######################################################################
### Differential expression (t-test) #################################
######################################################################

#t-test instead of deseq - valentin says no point doing deseq if having enough samples.
#this version of the t-test is faster because it doesn't use apply or other crap
myfastt <- function(groupA, groupB, keepgene = gene_mean_ss>=0.00001){
  fac <- rep(0,ncol(dat))
  fac[groupA] <- 1
  fac[groupB] <- 2
  tt <- log_ncount_all[keepgene,fac>0]
  mean1<- rowMeans(log_ncount_all[keepgene,fac==1])
  mean2<- rowMeans(log_ncount_all[keepgene,fac==2])
  out <- rowttests(as.matrix(tt), factor(fac[fac>0]))
  #out[order(out$p.value),]
  out <- cbind(rownames(out),mean1,mean2,out)
  colnames(out)[1] <- "ensembl_gene_id"
  colnames(out)[2] <- "mean1"
  colnames(out)[3] <- "mean2"
  sqldf("select * from out natural join ensconvert")
}

## differences between levels of grna (from facs)
orderttest <- function(out) out[order(out$p.value),]
t2vs4 <- myfastt(
  cellcondition$levelnum==2 & cellcondition$isgood & cellcondition$isss,
  cellcondition$levelnum==4 & cellcondition$isgood & cellcondition$isss)
t12vs34 <- myfastt(
  cellcondition$levelnum %in% c(1,2) & cellcondition$isgood & cellcondition$isss,
  cellcondition$levelnum %in% c(3,4) & cellcondition$isgood & cellcondition$isss)
orderttest(t2vs4)[1:20,]
orderttest(t12vs34)[1:20,]
intersect(orderttest(t2vs4)$mgi_symbol[1:30],  orderttest(t12vs34)$mgi_symbol[1:30])
#Eif4a2"  "Gm2830"  "Hacd2"   "Pmm1"    "Cebpzos" "Usp48"   "Gm15772"


## differences between chemistries
#mainly Gm genes different. And Rps
tchem <- myfastt(
  cellcondition$isgood & cellcondition$isss,
  cellcondition$isgood & cellcondition$iscr)
orderttest(tchem)[1:40,] 

### differences between growth conditions
tfeeder <- myfastt(
  cellcondition$isgood & cellcondition$isss & cellcondition$isfeeder,
  cellcondition$isgood & cellcondition$isss & cellcondition$islif)
orderttest(tfeeder)[1:40,] #
tfeeder2 <- myfastt(
  cellcondition$isgood & cellcondition$isss & cellcondition$isfeeder & cellcondition$levelnum %in% c(1,2),
  cellcondition$isgood & cellcondition$isss & cellcondition$islif & cellcondition$levelnum %in% c(1,2))
orderttest(tfeeder2)[1:40,] #




#what is the cut-off high and low Ptma? not obvious. can be at 1. or up to 2-3
hist(as.double(log_ncount_all[ensidPtma,]))

tlowhiPtma <- myfastt(
  log_ncount_all[ensidPtma,]<1  & cellcondition$isgood & cellcondition$isss,
  log_ncount_all[ensidPtma,]>=1 & cellcondition$isgood & cellcondition$isss)
orderttest(tlowhiPtma)[1:100,]


# ensembl_gene_id     mean1     mean2 statistic        dm      p.value mgi_symbol
# 7341 ENSMUSG00000023944 2.5861598 6.6578101 -22.02148 -4.071650 3.105824e-74   Hsp90ab1
# 256  ENSMUSG00000026238 0.2629102 4.1185336 -21.88734 -3.855623 1.327005e-73       Ptma
# 3632 ENSMUSG00000045128 2.2817501 5.5259196 -21.23189 -3.244170 1.603380e-70     Rpl18a
# 6569 ENSMUSG00000041841 1.1087969 3.7512533 -20.68481 -2.642456 5.980984e-68      Rpl37
# 4851 ENSMUSG00000034994 1.8546485 5.4083090 -19.91415 -3.553661 2.480590e-64       Eef2
# 255  ENSMUSG00000026234 1.8919222 5.3921190 -19.84674 -3.500197 5.135641e-64        Ncl
# 3028 ENSMUSG00000040952 2.1431352 5.1483719 -19.82328 -3.005237 6.615661e-64      Rps19
# 6688 ENSMUSG00000060036 1.4971144 4.6617431 -19.71159 -3.164629 2.207475e-63       Rpl3
# 4274 ENSMUSG00000025794 1.4819584 4.8052868 -19.69393 -3.323328 2.670854e-63      Rpl14
# 4109 ENSMUSG00000032294 2.1000371 5.9623721 -19.61676 -3.862335 6.138572e-63        Pkm
# 1321 ENSMUSG00000042244 2.3702680 0.3697308  19.59053  2.000537 8.144661e-63    Pglyrp3
# 6322 ENSMUSG00000025290 1.6900735 4.6523716 -19.58331 -2.962298 8.804820e-63      Rps24
# 5013 ENSMUSG00000057113 2.4047223 6.4469910 -19.53145 -4.042269 1.539902e-62       Npm1
# 2917 ENSMUSG00000030744 2.1451703 5.9539280 -19.48289 -3.808758 2.598734e-62       Rps3
# 7284 ENSMUSG00000008668 2.1690472 5.4090944 -19.42783 -3.240047 4.703013e-62      Rps18
# 7615 ENSMUSG00000024608 2.2678517 5.5787002 -19.41751 -3.310848 5.256138e-62      Rps14
# 1601 ENSMUSG00000028234 1.8633316 5.0450588 -19.31440 -3.181727 1.595495e-61      Rps20
# 3145 ENSMUSG00000003429 1.9525598 5.5563599 -19.09454 -3.603800 1.699048e-60      Rps11
# 7391 ENSMUSG00000057863 1.8298497 4.2493439 -19.01773 -2.419494 3.879634e-60      Rpl36
# 4271 ENSMUSG00000032518 2.6889902 6.1388280 -19.00779 -3.449838 4.316899e-60       Rpsa
# 2329 ENSMUSG00000072692 0.9238852 3.1301720 -19.00449 -2.206287 4.472890e-60    Rpl37rt
# 6600 ENSMUSG00000022283 1.5815919 4.7789735 -18.97371 -3.197382 6.225860e-60     Pabpc1
# 2713 ENSMUSG00000051695 1.2890594 4.0657307 -18.90767 -2.776671 1.265492e-59      Pcbp1
# 6337 ENSMUSG00000042354 1.0497911 4.0081902 -18.90123 -2.958399 1.356102e-59       Gnl3
# 144  ENSMUSG00000043716 1.5939318 5.3080389 -18.89842 -3.714107 1.397647e-59       Rpl7
# 2767 ENSMUSG00000057841 1.9547244 5.0888552 -18.89797 -3.134131 1.404459e-59      Rpl32
# 5213 ENSMUSG00000020372 2.3690931 5.8762672 -18.84604 -3.507174 2.452851e-59      Rack1
# 3058 ENSMUSG00000037563 0.9403649 3.3338979 -18.83467 -2.393533 2.771112e-59      Rps16
# 2650 ENSMUSG00000004980 1.5555389 4.7831787 -18.78346 -3.227640 4.801286e-59  Hnrnpa2b1
# 3315 ENSMUSG00000030695 2.2996978 6.3066492 -18.72457 -4.006951 9.030655e-59      Aldoa
# 1573 ENSMUSG00000059291 0.9910967 3.5094000 -18.64363 -2.518303 2.150865e-58      Rpl11
# 2396 ENSMUSG00000060419 1.0015533 3.2920716 -18.63495 -2.290518 2.360531e-58  Rps16-ps2
# 7242 ENSMUSG00000052146 1.3371594 4.0269179 -18.57931 -2.689759 4.284874e-58      Rps10
# 4508 ENSMUSG00000031320 1.6914519 5.2782513 -18.56692 -3.586799 4.893360e-58      Rps4x
# 3386 ENSMUSG00000025508 1.6547138 4.5107056 -18.56029 -2.855992 5.253527e-58      Rplp2
# 282  ENSMUSG00000081051 1.7186405 4.1275934 -18.52898 -2.408953 7.346186e-58    Gm15427
# 4960 ENSMUSG00000025362 1.7817472 4.5892126 -18.52746 -2.807465 7.466708e-58      Rps26
# 5492 ENSMUSG00000035530 1.2912423 4.3691622 -18.51184 -3.077920 8.826678e-58       Eif1
# 2357 ENSMUSG00000029614 1.4113413 4.7215349 -18.47814 -3.310194 1.266168e-57       Rpl6
# 1978 ENSMUSG00000028936 1.3270877 4.1371955 -18.45162 -2.810108 1.681829e-57      Rpl22
# 4827 ENSMUSG00000063457 1.9453491 4.7851291 -18.45081 -2.839780 1.696386e-57      Rps15
# 6594 ENSMUSG00000058600 1.3248627 4.2973426 -18.38619 -2.972480 3.386746e-57      Rpl30
# 2395 ENSMUSG00000029430 1.3810850 4.8474840 -18.38242 -3.466399 3.526342e-57        Ran
# 2484 ENSMUSG00000041453 1.7213421 4.5002102 -18.36356 -2.778868 4.314448e-57      Rpl21
# 4163 ENSMUSG00000037742 3.2139986 6.9286246 -18.31175 -3.714626 7.507810e-57     Eef1a1
# 1790 ENSMUSG00000028639 1.5677855 4.7459402 -18.27596 -3.178155 1.100682e-56       Ybx1
# 7026 ENSMUSG00000050299 1.1497729 3.3829488 -18.27147 -2.233176 1.154761e-56     Gm9843
# 3201 ENSMUSG00000061787 1.3696211 4.3724243 -18.25511 -3.002803 1.375267e-56      Rps17
# 4946 ENSMUSG00000025393 1.9794447 5.9055648 -18.23956 -3.926120 1.623883e-56      Atp5b
# 5895 ENSMUSG00000021270 2.5568564 6.1388565 -18.22793 -3.582000 1.838767e-56   Hsp90aa1
# 685  ENSMUSG00000062647 1.1586922 3.9550405 -18.20274 -2.796348 2.406430e-56      Rpl7a
# 7920 ENSMUSG00000024991 1.3879326 4.4087501 -18.19886 -3.020818 2.508122e-56      Eif3a
# 7033 ENSMUSG00000025613 1.3155670 4.9804116 -18.18419 -3.664845 2.933488e-56       Cct8
# 2670 ENSMUSG00000036371 1.7225652 4.8669042 -18.16396 -3.144339 3.640862e-56     Serbp1
# 3305 ENSMUSG00000032637 0.8840982 3.5294010 -18.15999 -2.645303 3.798376e-56     Atxn2l
# 5304 ENSMUSG00000060938 1.9073951 4.7983060 -18.15123 -2.890911 4.170630e-56      Rpl26
# 4944 ENSMUSG00000061315 1.1357732 4.1313042 -18.13357 -2.995531 5.036145e-56       Naca
# 4969 ENSMUSG00000020850 1.2251455 4.1335185 -18.12017 -2.908373 5.810113e-56      Prpf8
# 1964 ENSMUSG00000063524 2.3212197 6.0260431 -18.08519 -3.704823 8.438452e-56       Eno1
# 6663 ENSMUSG00000003970 2.5050358 5.7739992 -18.06424 -3.268963 1.055192e-55       Rpl8
# 2230 ENSMUSG00000047215 1.1443933 3.9434694 -18.03392 -2.799076 1.458035e-55       Rpl9
# 3160 ENSMUSG00000059070 1.9458407 5.0379259 -17.97248 -3.092085 2.806729e-55      Rpl18
# 5326 ENSMUSG00000078812 1.7443050 5.5418949 -17.94848 -3.797590 3.624710e-55      Eif5a
# 6384 ENSMUSG00000060373 1.1521287 4.6213351 -17.88681 -3.469206 6.990820e-55     Hnrnpc
# 4404 ENSMUSG00000084349 0.9729341 3.1247165 -17.88575 -2.151782 7.069901e-55   Rpl3-ps1
# 360  ENSMUSG00000040225 1.2753474 3.9699697 -17.87361 -2.694622 8.046027e-55     Prrc2c
# 2877 ENSMUSG00000036427 1.6383829 5.0783443 -17.87330 -3.439961 8.072541e-55       Gpi1
# 4120 ENSMUSG00000032399 2.9081844 6.2399268 -17.87213 -3.331742 8.173624e-55       Rpl4
# 5573 ENSMUSG00000016559 1.4068812 4.9210001 -17.79511 -3.514119 1.855058e-54      H3f3b
# 3983 ENSMUSG00000070319 0.9336509 3.7240742 -17.76253 -2.790423 2.623427e-54      Eif3g
# 356  ENSMUSG00000053332 1.3776553 4.2623159 -17.71000 -2.884661 4.585288e-54       Gas5
# 4139 ENSMUSG00000036781 1.2487346 4.1529012 -17.70973 -2.904167 4.598502e-54     Rps27l
# 3306 ENSMUSG00000030738 1.0608494 3.8502059 -17.69191 -2.789357 5.557344e-54      Eif3c
# 3273 ENSMUSG00000008683 1.4632612 4.3439967 -17.68462 -2.880735 6.004919e-54     Rps15a
# 3814 ENSMUSG00000000740 1.3296014 3.6861977 -17.68091 -2.356596 6.246704e-54      Rpl13
# 1264 ENSMUSG00000028081 1.9440640 5.2142934 -17.65704 -3.270229 8.049290e-54     Rps3a1
# 4958 ENSMUSG00000093674 2.6881315 5.3028108 -17.62838 -2.614679 1.091326e-53      Rpl41
# 2986 ENSMUSG00000012848 2.6252577 6.0462056 -17.58135 -3.420948 1.797886e-53       Rps5
# 4046 ENSMUSG00000009927 1.4100974 4.3563693 -17.56629 -2.946272 2.109573e-53      Rps25
# 7803 ENSMUSG00000071644 1.3873571 4.5243644 -17.55871 -3.137007 2.286292e-53      Eef1g
# 7645 ENSMUSG00000025428 1.6252902 5.1110340 -17.55539 -3.485744 2.368237e-53     Atp5a1
# 989  ENSMUSG00000027620 1.3122748 3.7966667 -17.55094 -2.484392 2.482596e-53      Rbm39
# 422  ENSMUSG00000026615 1.0347017 3.9994604 -17.52981 -2.964759 3.106467e-53       Eprs
# 4544 ENSMUSG00000080775 0.8085287 2.8093760 -17.49383 -2.000847 4.549332e-53     Gm6368
# 7276 ENSMUSG00000041881 1.1511759 4.1473322 -17.49031 -2.996156 4.722526e-53     Ndufa7
# 4226 ENSMUSG00000063856 1.3011855 4.3736443 -17.47479 -3.072459 5.567342e-53       Gpx1
# 1444 ENSMUSG00000047676 1.1021510 3.0060369 -17.44915 -1.903886 7.305302e-53  Rpsa-ps10
# 4890 ENSMUSG00000061904 1.4471015 4.7906957 -17.44334 -3.343594 7.769637e-53    Slc25a3
# 2795 ENSMUSG00000023456 1.5338994 4.9407072 -17.40207 -3.406808 1.203017e-52       Tpi1
# 5335 ENSMUSG00000018286 1.1450258 4.3445239 -17.39564 -3.199498 1.287801e-52      Psmb6
# 7566 ENSMUSG00000024359 1.2312163 4.4367710 -17.39210 -3.205555 1.336918e-52      Hspa9
# 4881 ENSMUSG00000020048 1.2491386 4.2138872 -17.38669 -2.964749 1.415790e-52    Hsp90b1
# 5294 ENSMUSG00000019505 1.7336193 5.2411844 -17.38402 -3.507565 1.456324e-52        Ubb
# 6591 ENSMUSG00000022234 1.3324907 4.7543684 -17.35662 -3.421878 1.946559e-52       Cct5
# 3410 ENSMUSG00000031490 1.3550136 4.6389157 -17.35641 -3.283902 1.950795e-52   Eif4ebp1
# 1513 ENSMUSG00000057280 1.5794983 0.1705183  17.34951  1.408980 2.098561e-52       Musk
# 7163 ENSMUSG00000014769 0.9894916 4.0205717 -17.33006 -3.031080 2.578282e-52      Psmb1
# 6974 ENSMUSG00000060636 1.2071299 3.6623542 -17.32144 -2.455224 2.824638e-52     Rpl35a
# 4034 ENSMUSG00000015656 1.7831191 5.3549311 -17.32016 -3.571812 2.862954e-52      Hspa8
# 2403 ENSMUSG00000070493 1.1091121 4.3558308 -17.31871 -3.246719 2.907320e-52     Chchd2

tlowhiPtmaNowt <- myfastt(
  log_ncount_all[ensidPtma,]<1  & cellcondition$isgood & cellcondition$isss & cellcondition$levelnum>1,
  log_ncount_all[ensidPtma,]>=1 & cellcondition$isgood & cellcondition$isss & cellcondition$levelnum>1)
orderttest(tlowhiPtmaNowt)[1:100,]


#TODO GO analysis?

write.table(tlowhiPtma$ensembl_gene_id[1:200], file="gotest_lowhiptma.csv")

#eif4a2?

#mean(out$mean1)

#TODO also use brennecke!


######################################################################
### brennecke, gozde is doing this ###############################################
######################################################################



######################################################################
### gRNA counts ###############################################
######################################################################

#TODO check again, which grna molecules do we have? can we look at the promoter in more detail?

### what are the relative fractions of different gRNAs? Here for smartseq
#grnatotc <- rowMeans(ncount_all[grep("grna-ptma",rownames(ncount_all)),cellcondition$iss & cellcondition$isgood])
#grnabg <- mean(rowMeans(ncount_all[setdiff(grep("grna-",rownames(ncount_all)),grep("grna-ptma",rownames(ncount_all))),]))
plotgrnapie <- function(ncount_all){
  grnatotc <- rowMeans(ncount_all[grep("grna-ptma",rownames(ncount_all)),])
  grnatotc <- grnatotc/sum(grnatotc)
#  png("plots/grna-ptmafraction SS2.png",w=800)
  #barplot(grnatotc, main="PTMA gRNA fractions", cex.names=0.5,beside=TRUE) 
  pie(grnatotc, labels = paste(paste(sprintf("grna-%s",1:5),round(grnatotc*100)),"%",sep=""),cex=1)
 # dev.off()
  grnatotc
}
png("plots/grna-ptmafraction pie.png",w=800)
par(mfrow=c(1,2),oma=c(0,0,2,0))
plotgrnapie(ncount_all[,cellcondition$isss & cellcondition$isgood])
plotgrnapie(ncount_all[,cellcondition$iscr & cellcondition$isgood])
title("Overall gRNA presence ss2 vs DogSeq",outer=TRUE)
dev.off()




# i<-2
# par(mfrow=c(1,2))
# plot(sort(as.double(log_ncount_all[grep(sprintf("grna-ptma%s",i),rownames(ncount_all)),cellcondition$isss & cellcondition$isgood])))
# plot(sort(as.double(log_ncount_all[grep(sprintf("grna-ptma%s",i),rownames(ncount_all)),cellcondition$iscr & cellcondition$isgood])))

#### Compare the dynamic range of gRNA detection
plotsortgrna <- function(log_ncount_all){
  out <- sort(as.double(log_ncount_all[grep(sprintf("grna-ptma%s",i),rownames(ncount_all)),]))
  for(i in 2:5){
    out <- cbind(out,sort(as.double(log_ncount_all[grep(sprintf("grna-ptma%s",i),rownames(ncount_all)),])))
  }  
  cols=c("black","red","green","blue","gray")
  plot(out[,1],type="l",ylim=c(0,max(out)),xlab="cell",ylab="log relative gRNA counts")
  for(i in 2:5){
    lines(out[,i],col=cols[i])
  }  
}
png("plots/grna dynamic range.png",w=800)
par(mfrow=c(1,2),oma=c(0,0,2,0))
plotsortgrna(log_ncount_all[,cellcondition$isss & cellcondition$isgood])
plotsortgrna(log_ncount_all[,cellcondition$iscr & cellcondition$isgood])
title("Dynamic range ss2 vs DogSeq",outer=TRUE)
dev.off()

#Result: Because dogseq seem quite consistent without any bumps, that's the one to go for. there is definitely a bias in gRNA presence

######################################################################
### Upstream screening ###############################################
######################################################################

takecells_ko <- log_ncount_all[ensidPtma,]<1  & cellcondition$isgood & cellcondition$isss
takecells_wt <- log_ncount_all[ensidPtma,]>=3 & cellcondition$isgood & cellcondition$isss
takecells_all <- log_ncount_all[ensidPtma,]>1 & cellcondition$isgood & cellcondition$isss

# takecells_wt <- cellcondition$levelnum %in% c(3,4) & cellcondition$isss & cellcondition$isgood 
# takecells_ko <- cellcondition$levelnum %in% c(1,2) & cellcondition$isss & cellcondition$isgood 

# takecells_wt <- cellcondition$levelnum %in% c(3,4) & cellcondition$isss & cellcondition$isgood 
# takecells_ko <- cellcondition$levelnum %in% c(1,2) & cellcondition$isss & cellcondition$isgood 
# takecells_all <- cellcondition$levelnum %in% c(1,2) & cellcondition$isss & cellcondition$isgood 
cor_wt <- cor(t(ncount_all[gene_mean_ss>=0.001 ,takecells_wt]), method="spearman")  
cor_ko <- cor(t(ncount_all[gene_mean_ss>=0.001 ,takecells_ko]), method="spearman") 
cor_all <- cor(t(ncount_all[gene_mean_ss>=0.001 ,takecells_all]), method="spearman") 
cor_wt <- melt(cor_wt)
cor_ko <- melt(cor_ko)
cor_all <- melt(cor_all)

cor_diff <- cor_wt
cor_diff[,3] <- cor_wt[,3] - cor_ko[,3]
cor_diff <- cor_diff[!(cor_diff$Var1 %in% c("cas9","mcherry")),]
cor_diff <- cor_diff[!(cor_diff$Var2 %in% c("cas9","mcherry")),]
cor_diff <- cor_diff[cor_diff$Var1!=cor_diff$Var2,]
cor_diff <- sqldf("select * from cor_diff natural join (select ensembl_gene_id as Var1, mgi_symbol as id1 from ensconvert) natural join (select ensembl_gene_id as Var2, mgi_symbol as id2 from ensconvert)")

cor_diff <- cor_diff[order(cor_diff[,3]),]
cor_diff[1:15,]
tail(cor_diff,n=15)
#TODO sqldf in genesym!
plot(cor_diff$value)

# 124  ENSMUSG00000092341 ENSMUSG00000026234 -0.8647368  Malat1     Ncl
# 3784 ENSMUSG00000026234 ENSMUSG00000092341 -0.8647368     Ncl  Malat1
# 992  ENSMUSG00000092341 ENSMUSG00000063229 -0.7985910  Malat1    Ldha
# 3798 ENSMUSG00000063229 ENSMUSG00000092341 -0.7985910    Ldha  Malat1
# 2657 ENSMUSG00000076258 ENSMUSG00000034994 -0.7533713 Gm23935    Eef2
# 3267 ENSMUSG00000034994 ENSMUSG00000076258 -0.7533713    Eef2 Gm23935
# 2781 ENSMUSG00000076258 ENSMUSG00000093674 -0.7513449 Gm23935   Rpl41
# 3269 ENSMUSG00000093674 ENSMUSG00000076258 -0.7513449   Rpl41 Gm23935
# 1285 ENSMUSG00000093674 ENSMUSG00000076281 -0.7512870   Rpl41 Gm24270
# 2749 ENSMUSG00000076281 ENSMUSG00000093674 -0.7512870 Gm24270   Rpl41
# 558  ENSMUSG00000092341 ENSMUSG00000030057 -0.7379396  Malat1    Cnbp
# 3791 ENSMUSG00000030057 ENSMUSG00000092341 -0.7379396    Cnbp  Malat1
# 765  ENSMUSG00000076281 ENSMUSG00000003429 -0.7346778 Gm24270   Rps11
# 1253 ENSMUSG00000003429 ENSMUSG00000076281 -0.7346778   Rps11 Gm24270
# 1240 ENSMUSG00000092341 ENSMUSG00000015656 -0.7163867  Malat1   Hspa8

#malat1 is a lncrna regulating splicing, located in nuclear specles http://onlinelibrary.wiley.com/doi/10.15252/embj.201592655/full  

downstreamPtma <- tlowhiPtma$ensembl_gene_id[which(tlowhiPtma$p.value<1e-50)]
length(downstreamPtma)
col <- rep("black",nrow(cor_diff))
col[cor_diff$Var1 %in% downstreamPtma | cor_diff$Var2 %in% downstreamPtma] <- "red"
col[cor_diff$Var1 %in% listola2c$ensid | cor_diff$Var2 %in% listola2c$ensid] <- "red"
col[cor_diff$Var1 %in% listolapluripotency$ensid | cor_diff$Var2 %in% listolapluripotency$ensid] <- "red"
plot(cor_diff$value, col=col,cex=0.5)

#things in high correlation with Ptma 
#cor_wt



######################################################################
### Genes changing in correlation depending on Ptma level ############
######################################################################


plotcor2 <- function(gene1,gene2,n1="g1",n2="g2"){
  takecells_ko <- log_ncount_all[ensidPtma,]<2  & cellcondition$isgood & cellcondition$isss
  takecells_wt <- log_ncount_all[ensidPtma,]>=3 & cellcondition$isgood & cellcondition$isss
  par(mfrow=c(1,2))
  plot(log1(ncount_all[gene1,takecells_wt]),
       log1(ncount_all[gene2,takecells_wt]),cex.lab=1,xlab=n1,ylab=n2)
  plot(log1(ncount_all[gene1,takecells_ko]),
       log1(ncount_all[gene2,takecells_ko]),cex.lab=1,xlab=n1,ylab=n2)
  print(length(which()))
  #cor(x,y,method="spearman")
}
#png("plots/test.png",width=800)
plotcor2(toensid("Hspa8"),toensid("Mir6236"))#cool. now broken??
#dev.off()
plotcor2(toensid("Malat1"),toensid("Ncl"))#cool
plotcor2(toensid("Malat1"),toensid("Ldha"))#cool
plotcor2(toensid("Mir6236"),toensid("Rps11"))#cool

for(i in 1:20){
  if(as.character(cor_diff$id1[i]) < as.character(cor_diff$id2[i])){
    png(sprintf("plots/cor/ptmaLow ptmaHigh %s %s %s.png",i,cor_diff$id1[i],cor_diff$id2[i]),width=800)   
    plotcor2(as.character(cor_diff$Var1[i]),as.character(cor_diff$Var2[i]),cor_diff$id1[i],cor_diff$id2[i])
    dev.off()
  }
}

for(i in (-20:0)+nrow(cor_diff)){
  if(as.character(cor_diff$id1[i]) < as.character(cor_diff$id2[i])){
    png(sprintf("plots/cor/ptmaLow ptmaHigh %s %s %s.png",i,cor_diff$id1[i],cor_diff$id2[i]),width=800)   
    plotcor2(as.character(cor_diff$Var1[i]),as.character(cor_diff$Var2[i]),cor_diff$id1[i],cor_diff$id2[i])
    dev.off()
  }
}

#Gm23935 is a known miRNA. correlation change vs 

######################################################################
### Analysis of cell stage and speed #################################
######################################################################


mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignmentsLow  <- cyclone(as.matrix(ncount_all[,ncount_all[ensidPtma,]<2  & cellcondition$isgood & cellcondition$isss]),mm.pairs)
assignmentsHigh <- cyclone(as.matrix(ncount_all[,ncount_all[ensidPtma,]>3  & cellcondition$isgood & cellcondition$isss]),mm.pairs)

plotphase <- function(assignments){
  col <- character(length(assignments$phases))
  col[assignments$phases=="G1"] <- "red"
  col[assignments$phases=="G2M"] <- "blue"
  col[assignments$phases=="S"] <- "darkgreen"
  #plot(assignments$score$G1, assignments$score$G2M, col=col, pch=16)
  cnt <- c(
    length(which(assignments$phases=="G1")),
    length(which(assignments$phases=="G2M")),
    length(which(assignments$phases=="S"))
  )
  pie(cnt,labels=c("G1","G2M","S"))
  
  #could sort and make a heatmap too?
}

#more cells in S-phase with PTMA high
png("plots/cycstage ptma.png",width=800)
par(mfrow=c(1,2))
plotphase(assignmentsLow)
plotphase(assignmentsHigh)
dev.off()

#get some measure of speed
#higher if ptma high. so ptma high means less cycling
mean(apply(assignmentsLow$normalized.scores,1,max)) 
mean(apply(assignmentsHigh$normalized.scores,1,max))   #higher confidence here
#TODO continous measure

#S score vs PTMA level?
keepcells <- cellcondition$isgood & cellcondition$isss
assignmentsAll <- cyclone(as.matrix(ncount_all[,keepcells]),mm.pairs)
png("plots/phase ptma vs S.png")
plot(assignmentsAll$normalized.scores$S,
     log1(ncount_all[ensidPtma,keepcells]))
dev.off()
png("plots/phase ptma vs G2M.png")
plot(assignmentsAll$normalized.scores$G2M,
     log1(ncount_all[ensidPtma,keepcells]))
dev.off()
png("plots/phase ptma vs G1.png")
plot(assignmentsAll$normalized.scores$G1,
     log1(ncount_all[ensidPtma,keepcells]))
dev.off()

#S and G2M scores lower in general for low PTMA. suggests faster cycling when PTMA low. again
#If S phase longer, does that mean high PTMA keep cells in synthesis? all ribosomal genes higher with high ptma.
#would this not suggest a shorter synth phase? unless the cell is struggling to up synthesis so proteins are not higher.
#or they are high but nothing happens

#Clearly, faster cycling with high Ptma!!!
png("plots/phase avscore sorted by PTMA.png", width=800)
#plot(apply(assignmentsAll$scores,1,mean)[order(ncount_all[ensidPtma,keepcells])])
plot(as.double(log_ncount_all[ensidPtma,keepcells]),
     apply(assignmentsAll$scores,1,mean),xlab="Ptma",ylab="Average cell cycling score")
dev.off()

#Would be nice to sort by Ptma level, then show the 3 classifications
#G1 vs G2M, color coded by Ptma. Black circles mainly in the mid=S
png("plots/phase G1 vs G2M colored by PTMA.png", width=800)
plot(assignmentsAll$normalized.scores$G1,
     assignmentsAll$normalized.scores$G2M, col=colorbygene()[keepcells])
dev.off()




######################################################################
### DEseq model testing NOT USED ##############################################
######################################################################

# Only use normal genes, raw counts. no grna or cas9 levels
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




keepcells <- cellcondition$isss & cellcondition$levelnum%in%c(1,4)
cellcond <- data.frame(levelnum=factor(cellconditiongrna[keepcells,colnames(cellconditiongrna) %in% c("levelnum")]))
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells], colData = cellcond,
                              design = ~ levelnum) #and plate? es cell type?
dds <- dds[ rowSums(counts(dds)) > 1, ]
#dds$condition <- relevel(dds$condition, ref="1")   #set the reference level to be all 0?
dds <- DESeq(dds)
res <- results(dds)

resOrdered <- res[order(res$padj),]
resOrdered
togenesym(rownames(resOrdered)[1:5])
#which(resOrdered$padj>0.9997)
#for first miseq, ss plate only
#1 vs 4: "Ptma"   "Malat1" "Rab25"  "Nr5a2"  "Gnai3" 
#1 vs 3: "Narf"  "Scmh1" "Cdc45" "Gnai3" "Klf6" 
#2 vs 4: "Narf"    "Slc34a2" "Cdc45"   "H19"     "Gnai3"  
### why is fgf4 not coming out????
#Narf,Cdc45??



############ Model using average levels for each infection group

keepcells <- cellcondition$isss & cellcondition$levelnum%in%c(1,4)
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells],
                              colData = cellconditiongrna[keepcells,],
                              design = ~ meancas9 + meanptma) #and plate? es cell type?
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)




select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)



############ Bigger model, trying to use gRNA levels as contrasts

keepcells <- cellcondition$isss
cellcond <- cellconditiongrna[keepcells,colnames(cellconditiongrna) %in% c("cas9",sprintf("ptma%s",1:5))]
rankMatrix(cellcond)
dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells],colData = cellcond,
                              design = ~ cas9 + ptma1 + ptma2 + ptma3 + ptma4 + ptma5) #and plate? es cell type?
dds <- dds[ rowSums(counts(dds)) > 1, ]  #suggested prefilter
#dds$condition <- relevel(dds$condition, ref="untreated")   #set the reference level to be all 0?
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]



######################################################################
### Monocle NOT USED ###############################
######################################################################


geneids <- listolapluripotency$ensid
keepcells <- cellcondition$isss & cellcondition$levelnum %in% c(1,4)
lcounts <- log_ncount_all
selectgenes <- which(rownames(lcounts) %in% geneids)
lcounts <- lcounts[selectgenes,keepcells]

pd <- new("AnnotatedDataFrame", data = cellcondition[keepcells,])
#fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)
#todo match row names?

HSMM <- newCellDataSet(as.matrix(lcounts),
                       phenoData = pd,
                       expressionFamily=negbinomial())
#filter
HSMM <- detectGenes(HSMM, min_expr = 0.1)
print(head(fData(HSMM)))



# diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
#                                       fullModelFormulaStr="~Media")
# ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

disp_table <- dispersionTable(HSMM)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.5 &
                           dispersion_empirical >= 2 * dispersion_fit)$gene_id
HSMM_myo <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components=2)
HSMM <- orderCells(HSMM, reverse=FALSE)

plot_cell_trajectory(HSMM, 1, 2, color_by="levelnum")


######################################################################
### clustering with CDR not used ###############################
######################################################################


#install.packages("devtools")
devtools::install_github("VCCRI/CIDR")
library(cidr)
sData <- scDataConstructor(tags)
#> 
#> cidr> sData <- determineDropoutCandidates(sData)
#> 
#> cidr> sData <- wThreshold(sData)
sData <- scCluster(sData)
plot(sData@PC[,c(1,2)], col=cols,pch=sData@clusters, main="CIDR", xlab="PC1", ylab="PC2")

#https://github.com/VCCRI/CIDR


#https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0805-z    

#https://github.com/epierson9/ZIFA

