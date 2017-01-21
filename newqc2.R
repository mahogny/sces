######################################################################
##### this is for the final 3 hiseq
######################################################################

#also need system libraries: libssl. libcurl?
installdeps <- function(){
  install.packages(c("gplots","scatterplot3d","plotrix","sqldf","Matrix","pi0"))
  install.packages("tsne")
  install.packages("gProfileR")
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
  biocLite("org.Mm.eg.db")
  #mintlinux package: the mesa stuff. libx. freetype. python-pip python-dev
  #run: pip install scipy
}

library(rgl)
library(gProfileR)
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
library(topGO)
library(org.Mm.eg.db)
library(limma)
library(stringr)
library(reshape2)
library(scater)
library(scran)
library(scLVM)
library(zoo)

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
  dat666 <- dat[grep("ENSMUSG",rownames(dat)),tokeep] #for featurecount
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
ncount <- mynormalize(dat)
ncount <- as.matrix(ncount)
#ncount <- ncount[-which(rownames(ncount) %in% listnongenes),] #for htseq output. not used anymore
log_ncount <- log(1+ncount)


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
  #isth <- theplate %in% c(12)

  
  ##th KO
  isgood[theplate==12]<-TRUE
  genelay <- as.matrix(read.csv("layouttcelltest.csv",stringsAsFactors=FALSE))
  genelay <- cbind(genelay,genelay,genelay,rep("f",8),rep("f",8),rep("f",8),rep("f",8)) #not really true
  mouserep <- c(1,1,1,2,2,2,3,3,3,0,0,0)  #TODO was changed
#  for(i in 1:3)
#    mouserep <- rbind(mouserep,mouserep)
  ko <- rep("",length(thecols))
  mouse <- rep(0,length(thecols))
  for(i in 1:length(thecols))  {
    if(theplate[i]==12){
      ko[i] <- genelay[therows[i],thecols[i]]
      mouse[i] <- mouserep[thecols[i]]
    }
  }
## TODO end change
  
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
    islif=islif,
    ko=ko,
    mouse=mouse,
    stringsAsFactors=FALSE
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
  v <- log_ncount[ensid,]
  rgb(0, v/max(v[keep]), 0, maxColorValue=1)
}
# colorbyptma <- function(keep=rep(TRUE,ncol(dat))){
#   thecol <- rep("red",nrow(cellcondition))
#   v <- log_ncount[ensidPtma,]
#   rgb(0, v/max(v[keep]), 0, maxColorValue=1)
#   col <- rep("#000000",length(v))
#   col[v>median(as.double(v))] <- "#00FF00"
#   col
# }
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

togenesymnames<-function(x){
  names(x)<-togenesym(names(x))
  x
}



######################################################################
### GO functions #####################################################
######################################################################

### Simplified use of topgo
stopgo <- function(all_data,DE_data,ID=c("genename","symbol","EnsemblID")){
  #all_data <- factor(all_data)#rownames(dat)
  #DE_data <-listolacellcycle$genesym
  #all_data <- rownames(ncount)
  #all_data <- unique(ensconvert$mgi_symbol)
  
  relevant.genes <- rep(1,length(all_data))
  relevant.genes[all_data %in% DE_data] <- 0
  names(relevant.genes) <- all_data
  relevant.genes <- as.factor(relevant.genes)
  
  GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes,nodeSize=5,annot = annFUN.org, 
                   mapping="org.Mm.eg.db", ID=ID, geneSel = function(p) p < 0.01)
  # apply Fisher's exact test with elimination mode:
  results <- runTest(GOdata.BP, algorithm = 'elim', statistic = 'fisher')
  GenTable(GOdata.BP, results, topNodes = 40)
}
## Even simpler use of topgo. background is all genes by default
# changed to gprofiler, much faster than topgo
stopgosym <- function(genelist,bg=unique(ensconvert$mgi_symbol)){
  x<-gprofiler(
    genelist,organism = "mmusculus", ordered_query = FALSE,
    correction_method = "fdr",
    hier_filtering = "moderate")
  x<-x[order(x$p.value),]
  x<-x[x$p.value < 1e-4,]
  x[,c("p.value","term.name")]
  
  #stopgo(bg,genelist,"symbol")
}
#mytopgo(rownames(ncount),orderttest(tlowhiPtma)[1:100,7],ID="genename")

  

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


listkedarcellcycle <- 
  rbind(cbind(read.table("kedar_g1.txt",sep="\t",stringsAsFactors = FALSE)[,1],"G1"),
        cbind(read.table("kedar_g2m.txt",sep="\t",stringsAsFactors = FALSE)[,1],"G2M"),
        cbind(read.table("kedar_s.txt",sep="\t",stringsAsFactors = FALSE)[,1],"S"))
listkedarcellcycle <- as.data.frame(listkedarcellcycle, stringsAsFactors = FALSE)
colnames(listkedarcellcycle) <- c("ensid","cellstage")
#listkedarcellcycle <- sqldf("select * from listkedarcellcycle natural join ensconvert")
listkedarcellcycle <- merge(listkedarcellcycle,ensconvert,by.x="ensid",by.y="ensembl_gene_id")
colnames(listkedarcellcycle)[3] <- "genesym"

######################################################################
### New QC for cells #################################################
######################################################################

#Numbef of grna reads
grna_count <- colSums(dat[grep("grna",rownames(dat)),])
#Number of exonic reads
gene_count <- colSums(dat[grep("ENSMUSG",rownames(dat)),])
#Number of detected genes
detected_genes <- colSums(dat>0)
#Proportion of mitochondrial reads
mt_counts <- colSums(dat[mt_genes$ensembl_gene_id,])
mt_prop <- mt_counts / gene_count  

findbadcells <- function(forcells,mt_cutoff=0.1){
  # Using scater to assign outliers
  libsize.drop <- isOutlier(gene_count[forcells], nmads=2, type="lower",log=TRUE)
  feature.drop <- isOutlier(detected_genes[forcells], nmads=2, type="lower", log=TRUE)
  mito.drop <- mt_prop[forcells]>mt_cutoff
  
  todrop <- libsize.drop | feature.drop | mito.drop
  print(dim(todrop))
  totc<-length(todrop)
  print(as.matrix(list(
    total=totc,
    kept=totc-length(which(todrop)),
    drop.by.size=length(which(libsize.drop)),
    drop.by.features=length(which(feature.drop)),
    drop.by.mito=length(which(mito.drop))
  )))
  
  todroptot <- rep(FALSE, ncol(dat))
  todroptot[(1:ncol(dat))[forcells][todrop]] <- TRUE
  
  par(mfrow=c(1,3))
  hist(log10(gene_count[forcells]), xlab="Mapped reads (log-scale)", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(gene_count[forcells][libsize.drop])
  abline(v=log10(v), col="red", lty=2, lwd=1.5) 
  hist(detected_genes[forcells], xlab="Number of detected genes", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(detected_genes[forcells][feature.drop])
  abline(v=v, col="red", lty=2, lwd=1.5) 
  hist(mt_prop[forcells], xlab="mtDNA %", main="",breaks=50, col="grey80", ylab="Number of cells")
  abline(v=mt_cutoff, col="red", lty=2, lwd=1.5)
  #dev.off()
  
  todroptot
}

png("plots/QC_ss.png",width=800)
x<-findbadcells(cellcondition$isss)
#as.double(gene_count[x]) #read counts of dropped genes
cellcondition$isgood[x]<-FALSE
dev.off()
png("plots/QC_dogseq.png",width=800)
x<-findbadcells(cellcondition$iscr) 
cellcondition$isgood[x]<-FALSE
dev.off()
png("plots/QC_th.png",width=800)
x<-findbadcells(cellcondition$plate==12)
cellcondition$isgood[x]<-FALSE
dev.off()
#plot(gene_count[cellcondition$plate==12]) #easily 2k reads. 30k cutoff would be ok


##better after QC!
#gene_mean_cr <- rowMeans(ncount[,cellcondition$iscr]) 
gene_mean_ss <- rowMeans(ncount[,cellcondition$isss])
gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
#gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)

### Comparing average gene expression among conditions
plotgenemean <- function(){
  gene_mean_cr <- rowMeans(ncount[,cellcondition$iscr]) 
  gene_mean_ss <- rowMeans(ncount[,cellcondition$isss])
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
  gene_mean_cr <- rowMeans(ncount[,cellcondition$iscr]) 
  gene_mean_ss <- rowMeans(ncount[,cellcondition$isss])
  gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)
  gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
  
  length(which(gene_mean_ss>=0.00001))
  #cor_ncount <- cor(ncount[gene_mean_ss>=0.00001,], method="spearman")
#  takecells <- cellcondition$plate %in% c(1) & cellcondition$isss & cellcondition$isgood #& !cellcondition$iswt 
#  takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood #& !cellcondition$iswt 
  takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood &cellcondition$isfeeder  #& !cellcondition$iswt 
  takecells <- cellcondition$plate %in% c(1) & cellcondition$isss & cellcondition$isgood &cellcondition$isfeeder  #& !cellcondition$iswt 
  
  takecells[weird ] <- FALSE
  cor_ncount <- cor(ncount[gene_mean_ss>=0.00001 ,takecells], method="spearman")  #olas way
  #cor_ncount <- t(ncount[gene_mean_ss>=0.00001 ,takecells])  #ola does not scale
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
  
  the_gene_mean <- rowMeans(ncount[,cellcondition$isgood & cellcondition$isss])
  the_gene_mean <- the_gene_mean/sum(the_gene_mean)
  
  length(which(the_gene_mean>=0.00001))
  #cor_ncount <- cor(ncount[gene_mean_ss>=0.00001,], method="spearman")
  takecells <- cellcondition$isss & cellcondition$isgood  # %in% c(1) #cellcondition$isss #& !cellcondition$iswt 
#  cor_ncount <- cor(ncount[gene_mean_ss>=0.000001 ,takecells], method="spearman")  #olas way
  cor_ncount <- cor(log_ncount[the_gene_mean>=0.00001 ,takecells], method="spearman")  #better?
  #cor_ncount <- t(ncount[gene_mean_ss>=0.00001 ,takecells])  #ola does not scale
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
  
  g1<-log_ncount[ ,takecells][,clust$cluster==1]
  g2<-log_ncount[ ,takecells][,clust$cluster==2]

  df <- data.frame(p.value=rep(666,nrow(log_ncount)), x=rep(1,nrow(log_ncount)), y=rep(1,nrow(log_ncount)))
  for(i in (1:nrow(log_ncount))[the_gene_mean>=0.001]){
    t <- t.test(g1[i,],g2[i,])
    df[i,1] <- t$p.value
    df[i,2] <- t$estimate[1]
    df[i,3] <- t$estimate[2]
  }
#  df<-data.frame(pval=outt,fold=outf)  
  rownames(df)<-rownames(log_ncount)
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
#batch_log_ncount <- removeBatchEffect(log_ncount, batchvector)

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
#   rowMeans(ncount[,cellcondition$level1]),
#   rowMeans(ncount[,cellcondition$level2]),
#   rowMeans(ncount[,cellcondition$level3]),
#   rowMeans(ncount[,cellcondition$level4]))
genemeanperlevel <- cbind(
  rowMeans(log_ncount[,cellcondition$level1]),
  rowMeans(log_ncount[,cellcondition$level2]),
  rowMeans(log_ncount[,cellcondition$level3]),
  rowMeans(log_ncount[,cellcondition$level4]))

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
png("plots/ss2 olaprevtargets ptma byfacs.png",width = 1000,pointsize = 16)
par(cex.lab=0.2)
par(mfrow=c(1,1))
barplot(b, main="Genes per PTMA facs level",
        col=c("#000000","#880000","#AA0000","#FF0000"), cex.names=0.5,beside=TRUE) 
dev.off()

################################# other way, by PTMA level

genemeanperlevel <- cbind(
  rowMeans(log_ncount[,cellcondition$isss & cellcondition$isgood & log_ncount[ensidPtma,]>3]),
  rowMeans(log_ncount[,cellcondition$isss & cellcondition$isgood & log_ncount[ensidPtma,]<3]))

#### Redoing fig7d. What is down-regulated by PTMA?
andid <- c(ensidPtma)#,toensid(c("Narf","Cdc45")))
b <- t(genemeanperlevel[
  rownames(genemeanperlevel) %in% c(listolaptmagenes7d$ensid,andid),
  c(1,2)])
#use order as in olas paper. ola sorts by fold change
b<-b[,c(listolaptmagenes7d$ensid, andid)] 
# for(i in 1:ncol(b))
#   b[,i] <- b[,i]/max(b[,i]) #normalize to highest level
colnames(b) <- togenesym(colnames(b)) #order preserving? not any more!
png("plots/ss2 olaprevtargets ptma byptma.png",width = 1000,pointsize = 16)
par(cex.lab=0.2)
par(mfrow=c(1,1))
barplot(b, main="Genes per PTMA facs level",
        col=c("#000000","#880000"), cex.names=0.5,beside=TRUE) 
dev.off()


######################


#TODO so correlation not great by FACS level. but better to reproduce plot by Ptma level then?
#TODO by ptma level. then by cherry level. just split into 4 bins. also only use SS?

#can also sort each gene exp by Ptma level and make a dot-plot?
dat2 <- ncount[,cellcondition$isgood & cellcondition$isss]
dat2 <- dat2[,order(dat2[ensidPtma,])]
plot(as.double(log1(dat2[listolaptmagenes7d$ensid[15],])))
#2,10,15 very bimodal  ... this is actually much clearer this way than below? would have to fit a sigmoidal.
#then two pieces of information; change in level and bimodality
plot(
  as.double(log1(dat2[ensidPtma,])),
  as.double(log1(dat2[listolaptmagenes7d$ensid[9],]))
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
##length(which(gene_mean_ss>=0.0001))
#takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood  & cellcondition$levelnum>1
takecells <- cellcondition$isss & cellcondition$isgood & log_ncount["cas9",]>3
#takecells[weird ] <- FALSE
cor_genes <- cor(t(ncount[gene_mean_ss>=0.00001 | rownames(ncount)==ensidPtma ,takecells]), method="spearman")  #olas way
cor_ptma <- cor_genes[rownames(cor_genes)==ensidPtma]
names(cor_ptma) <- colnames(cor_genes)
nametogenesym(sort(cor_ptma)[1:30])
nametogenesym(sort(cor_ptma,decreasing = TRUE)[1:10])
plot(sort(cor_ptma,decreasing = TRUE))




#update from full list:
# mcherry       Egr1      Ypel3    Aldh6a1        Fos      Ddit4       Sdc4       Phc3       Aire        Pcx 
# -0.5296585 -0.3515679 -0.3378253 -0.3210335 -0.3128279 -0.2884587 -0.2806149 -0.2745264 -0.2742435 -0.2649969 
# Tmem43    Pglyrp3    Nectin2    Slc4a11      Thap3       Pomk      Lpin1     Tm7sf3       Chd9     Fam35a 
# -0.2648565 -0.2614314 -0.2612546 -0.2587797 -0.2572245 -0.2563735 -0.2555303 -0.2513488 -0.2507193 -0.2504484 
# Pnrc1      Chac1     Cdkn1a     Pmaip1   Tmem184b       Hexb     Agpat4      Psrc1       Ctc1       Vmac 
# -0.2501441 -0.2494819 -0.2492637 -0.2469392 -0.2463481 -0.2459591 -0.2432356 -0.2410512 -0.2395969 -0.2395790 
# > nametogenesym(sort(cor_ptma,decreasing = TRUE)[1:10])
# Ptma    Gm9800    Gm4617       Erh    Ranbp1 Hspd1-ps3       Ncl     Hspd1   Gm21399     Hmgb2 
# 1.0000000 0.9501573 0.9229697 0.5796107 0.5685711 0.5679256 0.5585583 0.5567771 0.5448028 0.5417181 

##for ss 0.00001
# Pglyrp3    Gm10719       Lmo3    mcherry    Gm11168    Gm14762    Gm10801       Ccr1      Mss51    Gm26870 
# -0.4693054 -0.4404796 -0.4266743 -0.3945823 -0.3938720 -0.3822966 -0.3460647 -0.3308073 -0.3271836 -0.3176716 
# Ptma    Gm9800    Gm4617       Ncl       Erh     Hmgb2       Ran      Pcna     Psma5      Ppa1 
# 1.0000000 0.9623600 0.9419038 0.6859848 0.6818732 0.6742220 0.6652018 0.6650973 0.6635739 0.6607904 

#Gm9800    Gm4617  share a ridiculous part of the sequence, precisely



#cell cycle genes, subset
cor_cc <- cor_ptma
cor_cc <- cor_cc[names(cor_cc) %in% listkedarcellcycle$ensid]
nametogenesym(sort(cor_cc)[1:20])
nametogenesym(sort(cor_cc,decreasing = TRUE)[1:20])
plot(sort(cor_cc,decreasing = TRUE))
# Fos      Cdkn1a       Psrc1        Klf6        Plk2     Ccdc107        Xbp1        G2e3       Lrif1 
# -0.31282790 -0.24926366 -0.24105124 -0.23442641 -0.20045378 -0.18345354 -0.16242458 -0.14762064 -0.14209218 
# Plin3       Acsl3       Birc2        Mdc1     Hbb-bh1         Lpp       Zmym3     B4galt1        Ubr7 
# -0.13272100 -0.12543963 -0.11797901 -0.11356420 -0.11168566 -0.10243685 -0.09592027 -0.09224304 -0.08145435 
# Pkmyt1        Wwc1 
# -0.07723948 -0.06860648 
# Hmgb2       Pcna   Pcna-ps2       Lyar       Cbx3     Dynll1    Gm12355     Ube2d3      Hspa8       Nasp 
# 0.5417181  0.5368060  0.5218322  0.5028741  0.4889887  0.4843993  0.4818836  0.4776347  0.4689171  0.4628692 
# Slbp       Sod2       Sfpq     Srsf10       Rrm2      Nup37 Tubb4b-ps1     Tuba4a      Slmo2     Tubb4b 
# 0.4619750  0.4103535  0.4084943  0.4049469  0.4034397  0.3970951  0.3951333  0.3939280  0.3886303  0.3880144 
head(cor_cc)


######################################################################
### Plot gene vs gene ################################################
######################################################################

#Check some individual correlations
plotcor <- function(x,y,xlab="",ylab=""){
  plot(x,y,cex.lab=1,xlab=xlab, ylab=ylab)
  cor(x,y,method="spearman")
}

pgoneptma <- function(gene){
  gid <- toensid(gene)[1]
  plot(log1(ncount[gid,keepcells]),
       log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab=gene,ylab="Ptma")
}


keepcells <- cellcondition$isss & cellcondition$isgood
#>0
#keepcells <- cellcondition$isss & cellcondition$levelnum>0
#keepcells <- cellcondition$iscr & gene_count>1e3 & cellcondition$levelnum>0

#keepcells <- cellcondition$isss & gene_count>300e3 & cellcondition$levelnum>0

#2c markers accornding to ola
plot(log1(ncount[toensid("Fbxo15"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
#Fbxo15 slightly more correlation
#Rbm25  slightly more correlation
#Dst    little correlation
#

#ENSMUSG00000027712  Annexin V, apoptosis marker. no correlation at all
plot(log1(ncount[toensid("Anxa5"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
#https://www.ncbi.nlm.nih.gov/pubmed/24032683 ptma blocks caspase 3 activity. KO causes cell prolif defect
#caspase 3, not much happening
plot(log1(ncount[toensid("Casp3"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
#Keap1, nope. nor pgam5. possibly a positive correlation with bptf (http://www.sciencedirect.com/science/article/pii/S0891584911000700)
#https://en.wikipedia.org/wiki/BPTF
plot(log1(ncount[toensid("Bptf"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount[toensid("Dst"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)


plot(log1(ncount["ENSG00000141510",keepcells]), #this should be p53. why not here??
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount[toensid("Trp53"),keepcells]), #Trp53 positive correlation
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount[toensid("Greb1l"),keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

png("plots/ss2 some genes.png",width = 800)
par(mfrow=c(2,3))
plot(log1(ncount[toensid("Ypel3"),keepcells]), #known to induce cell arrest if inducted
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="Ypel3 (cell arrest)",ylab="Ptma")
plot(log1(ncount[toensid("Erh"),keepcells]), #suggested to be in cell cycle. and other things
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="Erh (most positive cor)",ylab="Ptma") #most positively correlated
plot(log1(ncount[toensid("Cdkn1a"),keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="Cdkn1a",ylab="Ptma")  #great correlation. negative
plot(log1(ncount[toensid("Hmgb2"),keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="Hmgb2",ylab="Ptma") #good correlation, positive
plot(log1(ncount[toensid("Pcna"),keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="Pcna",ylab="Ptma") #good correlation, positive
pgoneptma("Aire")
#plot(log1(ncount[toensid("Klf6"),keepcells]),
#     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="Klf6",ylab="Ptma")
dev.off()

##CC list 1
pgoneptma("Mapk13")
pgoneptma("Etv4")
pgoneptma("Ap3d1")
pgoneptma("Tsen54")
pgoneptma("Cdca7l")
pgoneptma("Casp2")
pgoneptma("Rad54b")
pgoneptma("Bora")

##CC het lis 2
pgoneptma("Ssr4")
pgoneptma("Adgre5")
pgoneptma("Ctcf")
pgoneptma("Ubxn11")
pgoneptma("Dnajb9")
pgoneptma("Hmgb3")
pgoneptma("Anp32e")
pgoneptma("Pttg1")

png("plots/ss2 cc by cor1.png",width = 800)
par(mfrow=c(2,3))
pgoneptma("Psrc1") #
pgoneptma("Ccdc107")#
pgoneptma("Plin3")#
pgoneptma("Acsl3")#
pgoneptma("Birc2")#
dev.off()
# Fos      Cdkn1a       Psrc1        Klf6        Plk2     Ccdc107        Xbp1        G2e3       Lrif1 
# Plin3       Acsl3       Birc2        Mdc1     Hbb-bh1         Lpp       Zmym3     B4galt1        Ubr7 
# Pkmyt1        Wwc1 

#pgoneptma("Casp2")

png("plots/ss2 cc by cor2.png",width = 800)
par(mfrow=c(2,3))
pgoneptma("Cbx3")#
pgoneptma("Gm12355")#
pgoneptma("Tubb4b-ps1")#
pgoneptma("Tuba4a")#
pgoneptma("Tubb4b")#
dev.off()

png("plots/ss2 cc by cor3.png",width = 800)
par(mfrow=c(2,3))
pgoneptma("Sdc4") ####
pgoneptma("Tmem43")##
pgoneptma("Nectin2") ##
pgoneptma("Slc4a11") #
pgoneptma("Pomk") ##
pgoneptma("Pmaip1") ##
dev.off()


### a copy
# mcherry       Egr1      Ypel3    Aldh6a1        Fos      Ddit4       Sdc4       Phc3       Aire        Pcx 
# -0.5296585 -0.3515679 -0.3378253 -0.3210335 -0.3128279 -0.2884587 -0.2806149 -0.2745264 -0.2742435 -0.2649969 
# Tmem43    Pglyrp3    Nectin2    Slc4a11      Thap3       Pomk      Lpin1     Tm7sf3       Chd9     Fam35a 
# -0.2648565 -0.2614314 -0.2612546 -0.2587797 -0.2572245 -0.2563735 -0.2555303 -0.2513488 -0.2507193 -0.2504484 
# Pnrc1      Chac1     Cdkn1a     Pmaip1   Tmem184b       Hexb     Agpat4      Psrc1       Ctc1       Vmac 


# Hmgb2       Pcna   Pcna-ps2       Lyar       Cbx3     Dynll1    Gm12355     Ube2d3      Hspa8       Nasp 
# Slbp       Sod2       Sfpq     Srsf10       Rrm2      Nup37 Tubb4b-ps1     Tuba4a      Slmo2     Tubb4b 


#positive cell cycle regulators
# Hmgb2       Pcna   Pcna-ps2       Lyar       Cbx3     Dynll1    Gm12355     Ube2d3      Hspa8       Nasp 

#negative cell cycle regulators. definitely Cdkn1a!
#Fos      Cdkn1a       Psrc1        Klf6 
#but seems better to do a t-test for this

#mayb weak positive correlation. Rela Eef2
#No correlation: Nfkb1. Casp3. Irf3
#Not expressed: Tlr4. Ifng. Ifnb1. 
plot(log1(ncount[toensid("Bax"),keepcells]), #maybe a positive correlation bax/ptma
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount[toensid("Bcl2"),keepcells]), #maybe a weak correlation to Bcl2. Bcl2 and Bax bind together
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
#p53 should upregulate Bax. Bax normally in cytosol but moving in for apoptosis; so expression might not matter

plot(log1(ncount[toensid("Trim25"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
#There seems to be a mild negative correlation to Trim25
#Correlate with innate immune response genes? 

#Map of toll like receptor pathway
#https://upload.wikimedia.org/wikipedia/commons/thumb/6/6b/Toll-like_receptor_pathways.svg/500px-Toll-like_receptor_pathways.svg.png
#Check if any Tlr expressed
rowMeans(log_ncount[ensconvert$ensembl_gene_id[grep("Tlr",ensconvert$mgi_symbol)],cellcondition$isgood & cellcondition$isss])
ensconvert$mgi_symbol[grep("Tlr",ensconvert$mgi_symbol)]

rowMeans(log_ncount[toensid("ENSMUSG00000063049"),cellcondition$isgood & cellcondition$isss])
rowMeans(log_ncount[toensid(c("Bax","Bad")),cellcondition$isgood & cellcondition$isss])
"Parp"
toensid(c("Ing2"))


#how is Ptma affecting mitochondria content? comparing to Ptma/total amount of non-mt reads
png("ptma vs mitochondria content ss.png")
keepcellsmt <- cellcondition$isss# & gene_count>800e3
plot(mt_prop[keepcellsmt], 
     log1(ncount[ensidPtma,keepcellsmt] / (1-mt_prop[keepcellsmt])),cex.lab=1,xlab="mtprop",ylab="Ptma/(1-mtprop)")
dev.off()

# plot(mt_prop[keepcellsmt], 
#      log1(ncount[toensid("Casp3"),keepcellsmt] / (1-mt_prop[keepcellsmt])),cex.lab=1)



#keap1 supposed to interact with ptma. there could be a trend. keap1+ptma = degrade Nrf2? meah?
plot(log1(ncount[toensid("Nfe2l2"),keepcells]), #literature suggests interaction nfe2l2. https://en.wikipedia.org/wiki/NFE2L2
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount[toensid("Nfe2l1"),keepcells]), 
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)


plot(log1(ncount[toensid("Ccnd1"),keepcells]),   
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

#note: if Gm4617 correlates so well with ptma then it can be used as a proxy for the cr2 data too

png("plots/ss2 ptma vs Gm9800.png",width = 800, height=600)
plot(log1(ncount[toensid("Gm9800"),keepcells]),  #super correlated
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 ptma vs Gm4617.png",width = 800, height=600)
plot(log1(ncount[toensid("Gm4617"),keepcells]),  #super correlated. predicted homolog of ptma http://www.ihop-net.org/UniPub/iHOP/gs/123843.html
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
dev.off()
plot(log1(ncount[toensid("Ncl"),keepcells]), #https://en.wikipedia.org/wiki/Nucleolin  ribosomes
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
png("plots/ss2 ptma vs Pglyrp3.png",width = 800, height=600)
plot(log1(ncount[toensid("Pglyrp3"),keepcells]),  #innate immune system. weak?
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
dev.off()
plot(log1(ncount[toensid("Lmo3"),keepcells]),  #LIM domain something. weak?
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount[toensid("Prdx1"),keepcells]),  #https://en.wikipedia.org/wiki/Peroxiredoxin_1  strong. after release, causes cytokines to be released
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

plot(log1(colSums(ncount[grep("grna-ptma",rownames(ncount)),keepcells])),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlim=c(0,1.5))
plot(log1(ncount["grna-ptma2",keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)



plot(log1(ncount[toensid("Npm1"),keepcells]),  #
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

plot(log1(ncount[toensid("Hsp90ab1"),keepcells]),  #came out of deseq strongly
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

#TODO get all the articles on ptma

#TODO deseq, high vs low ptma. and +- mcherry! level>1. and +- cas9

#ptma https://en.wikipedia.org/wiki/Thymosin_%CE%B11  - could it actually connect to pglyrp3

plot(log1(ncount[toensid("Egr1"),keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)


plot(log1(ncount[toensid("Ptp4a3"),keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)


#if I have which vary, I can use this for testing more genes
symCdk <- ensconvert$mgi_symbol[grep("Cdk",ensconvert$mgi_symbol)]
for(s in symCdk){
  png(sprintf("plots/gg/ptma vs %s.png",s))
  pgoneptma(s)
  dev.off()
}
#Pos: Cdk1 2(can be neg!) Cdk2ap1 Cdk4 Cdk7 20 Cdkn2aipnl   
#no idea: 2, 5 , 5rap3 9 12
#neg Cdk10 13 17 Cdkl13 Cdkn1a! n1b? n1c n2aip Cdkn2d Cdkn3
#independent Cdk11b 
##10 binds to Ets2 and modulates. 

#https://www.ncbi.nlm.nih.gov/pubmed/22547058
#cyclin K with cdk12,13, is important for self-renewing mES

symCyc <- ensconvert$mgi_symbol[grep("Ccn",ensconvert$mgi_symbol)]
for(s in symCyc){
  png(sprintf("plots/gg/ptma vs %s.png",s))
  pgoneptma(s)
  dev.off()
}

#pos: Ccnd1 Ccne1 Ccnk~  #D with 4/6. E with 2. both G1 and G1/S
#neg: Ccng1 Ccni  #G poorly studied??
#independent: Ccna2(mild pos) Ccnb1 Ccndb2 Ccnd3 Ccnd3-ps Ccne2 Ccni Ccnl1 Ccnl2 CCnt1 Ccnt2 Ccny

#
pgoneptmasave <- function(s){
  png(sprintf("plots/gg/ptma vs %s.png",s))
  pgoneptma(s)
  dev.off()
}

pgoneptma("Myc")
pgoneptma("p21") #=Cdkn1a = Ras!!

png(sprintf("plots/gg/ptma vs Raf1.png",s))
pgoneptma("Raf1")  #negative, cRaf. actually, could be positive
dev.off()

png(sprintf("plots/gg/ptma vs Rb1.png",s))
pgoneptma("Rb1")  #independent
dev.off()

Cdkn2aipnl

#myc. ras raf erk Rb 

#pgoneptmasave("Ink4a") #=Cdkn4a
pgoneptmasave("Arf")
pgoneptmasave("Anxa1")
pgoneptmasave("Cyr61")
pgoneptmasave("Edc3") #negative correlation
pgoneptmasave("Usp47")
pgoneptmasave("Glp2r") #no expression??


######################################################################
### Plot gene vs gene, mcherry & cas9 ################################
######################################################################

keepcells <- cellcondition$isss & cellcondition$isgood
keepcellsCas9high <- cellcondition$isss & cellcondition$isgood & log_ncount["cas9",]>3
keepcellsCas9low <- cellcondition$isss & cellcondition$isgood & log_ncount["cas9",]<3

# png("plots/ss2 ptma vs bfp.png",width = 800, height=600)
# plot(log1(ncount["bfp",keepcells]),
#      log1(ncount[ensidPtma,keepcells]),cex.lab=1)
# dev.off()
png("plots/ss2 ptma vs mcherry.png",width = 800, height=600)
plot(log1(ncount["mcherry",keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1,xlab="mcherry",ylab="Ptma")
dev.off()
png("plots/ss2 ptma vs mcherry vs cas9.png",width = 800, height=600)
par(mfrow=c(1,2))
plot(log1(ncount["mcherry",keepcellsCas9low]),
     log1(ncount[ensidPtma,keepcellsCas9low]),cex.lab=1,xlab="mcherry, cas9 low",ylab="Ptma",ylim=c(0,6),xlim=c(0,8))
plot(log1(ncount["mcherry",keepcellsCas9high]),
     log1(ncount[ensidPtma,keepcellsCas9high]),cex.lab=1,xlab="mcherry, cas9 high",ylab="Ptma",ylim=c(0,6),xlim=c(0,8))
dev.off()
png("plots/ss2 ptma vs cas9.png",width = 800, height=600)
plot(log1(ncount["cas9",keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1, xlab="cas9", ylab="Ptma")
dev.off()

keepcells2 <- cellcondition$isss & gene_count>300e3 & cellcondition$levelnum>2
plot(log1(ncount["cas9",keepcells2]),
     log1(ncount[ensidPtma,keepcells2]),cex.lab=1)

#TODO a linear model with cas9 & mcherry vs ptma


keepcells <- cellcondition$isss & cellcondition$isgood


# png("plots/ss2 bfp vs cas9.png",width = 800, height=600)
# plot(log1(ncount["bfp",keepcells]),
#      log1(ncount["cas9",keepcells]),cex.lab=1)
# dev.off()
png("plots/ss2 mcherry vs cas9.png",width = 800, height=600)
plot(log1(ncount["mcherry",keepcells]),
     log1(ncount["cas9",keepcells]),cex.lab=1) #good correlation
dev.off()

plot(log1(ncount["bfp",]))
plot(log1(ncount["mcherry",]))

plot(log1(ncount["cas9",keepcells]))
#plot(log1(ncount["bfp",keepcells]))
plot(log1(ncount["mcherry",keepcells]))

#check correlation bfp, cas9 and read count!
png("plots/ss2 cas9 vs readcount.png",width = 800, height=600)
plot(log1(ncount["cas9",keepcells]),
     log1(gene_count[keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 bfp vs readcount.png",width = 800, height=600)
plot(log1(ncount["bfp",keepcells]),
     log1(gene_count[keepcells]),cex.lab=1)
dev.off()
png("plots/ss2 mcherry vs readcount.png",width = 800, height=600)
plot(log1(ncount["mcherry",keepcells]),
     log1(gene_count[keepcells]),cex.lab=1)
dev.off()





#####################################

keepcells <- cellcondition$iscr & cellcondition$isgood
#keepcells <- cellcondition$iscr & gene_count>10e3 & cellcondition$levelnum>0
length(which(keepcells))


# png("plots/1kcr ptma vs bfp.png",width = 800, height=600)
# plot(log1(ncount["bfp",keepcells]),
#      log1(ncount[ensidPtma,keepcells]),cex.lab=1)
# dev.off()
png("plots/1kcr ptma vs mcherry.png",width = 800, height=600)
plot(log1(ncount["mcherry",keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
dev.off()
png("plots/1kcr ptma vs cas9.png",width = 800, height=600)
plot(log1(ncount["cas9",keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
dev.off()

###############

keepcells <- cellcondition$iscr & cellcondition$isgood & log1(ncount[ensidPtma,])<6.5
#keepcells <- cellcondition$iscr & gene_count>30e3 & cellcondition$isgood & log1(ncount[ensidPtma,])<6.5
#keepcells <- cellcondition$iscr & gene_count>20e3 & cellcondition$levelnum>0
#keepcells <- cellcondition$iscr & gene_count>200e3 & cellcondition$levelnum>0

png("plots/pairwise cr ptma vs grna.png",width = 800, height=600)
plotcor(log1(colSums(ncount[grep("grna-ptma",rownames(ncount)),keepcells])),
        log1(ncount[ensidPtma,keepcells]),xlab="log sum grna",ylab="ptma")
dev.off()
png("plots/pairwise cr ptma vs grna2.png",width = 800, height=600)
plotcor(colSums(log(1+ncount[grep("grna-ptma",rownames(ncount)),keepcells])),
        log1(ncount[ensidPtma,keepcells]),xlab="sum log grna",ylab="ptma")
dev.off()

png("plots/pairwise cr ptma vs mcherry.png",width = 800, height=600)
plotcor(log1(ncount["mcherry",keepcells]),
        log1(ncount[ensidPtma,keepcells]),xlab="mcherry",ylab="ptma")
dev.off()
#NOTE: exp(2)=7x more mcherry reads. so this means it might be easier to sequence and hence better at predicting

##gRNA level is as good a predictor as mcherry for ptmaKO. cells > 200k reads.

plot(log1(colSums(ncount[grep("grna-ptma",rownames(ncount)),keepcells])) + log1(ncount["mcherry",keepcells]),
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)

##All gRNAs are working .... no. wrong sign!
for(i in 1:5){
  png(sprintf("plots/pairwise cr ptma vs grna%s.png",i),width = 800, height=600)
  plot(log1(ncount[sprintf("grna-ptma%s",i),keepcells]),
       log1(ncount[ensidPtma,keepcells]),cex.lab=1)
  dev.off()
}




#TODO which grna is best predicting?

plot(log1(ncount[toensid("Gm9800"),keepcells]),  #super correlated
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)
plot(log1(ncount[toensid("Gm4617"),keepcells]),  #super correlated. predicted homolog of ptma http://www.ihop-net.org/UniPub/iHOP/gs/123843.html
     log1(ncount[ensidPtma,keepcells]),cex.lab=1)


#grna and mcherry as good predictors for PTMA level. but only when detecting the grna
plotcor(log1(colSums(ncount[grep("grna-ptma",rownames(ncount)),keepcells])), 
     log_metaptma)
plotcor(log1(ncount["mcherry",keepcells]),
     log_metaptma)



#an attempt at a linear model.
#but this is not quite the same as checking efficiency. does not take presence into account? or?
#keepcells <- cellcondition$iscr & gene_count>40e3 & cellcondition$isgood & log1(ncount[ensidPtma,])<6.5
keepcells <- cellcondition$isss & cellcondition$isgood 
#keepcells <- cellcondition$iscr & cellcondition$isgood  #now works but not as strongly as ss
lmdata <- t(ncount[c(ensidPtma, sprintf("grna-ptma%s",1:5)),keepcells])
lmdata <- log(1+lmdata)
colnames(lmdata) <- c("ptma",sprintf("grna%s",1:5))
mm <- lm(ptma ~ grna1 + grna2 + grna3 + grna4 + grna5 + 1, as.data.frame(lmdata))
mm
#annoyingly, it predicts better for the regular chemistry. likely because of the depth
#could we up ptma somehow?
#should divide by avg. amount of that grna


# For SS2:
# (Intercept)        grna1        grna2        grna3        grna4        grna5  
# 4.21398      1.78168     -0.07129     -1.46118     -1.40669     -0.71584  

sumgrna <- rowMeans(ncount[grep("grna-ptma",rownames(ncount)),cellcondition$isgood])
mm$coefficients[-1]/sumgrna
# grna1        grna2        grna3        grna4        grna5 
# 2.874554281 -0.008420633 -0.079980112 -0.021562446 -0.141865644 
# for ss2: normalized suggests grna5 strongest by factor 2x-7x. non-normalized, 3 and 4 strongest (2x of 5)
#grna1 has so little that we cannot assess strength
#note that model here is a bit messy, mixing log and non-log

######################################################################
### Most heterogenous genes ##########################################
######################################################################

## heterogenous, all cells (SS)
keepAll <- cellcondition$isgood & cellcondition$isss 
tnAll <- fitTechnicalNoise(
  ncount[,keepAll],
  fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 
hetAll  <- getVariableGenes(fit_type = 'log',ncount[,keepAll], tnAll$fit, plot=FALSE)
length(which(hetAll))

## heterogenious, Ptma + and -
keepLow  <- cellcondition$isgood & cellcondition$isss & log_ncount[ensidPtma,]<2
keepHigh <- cellcondition$isgood & cellcondition$isss & log_ncount[ensidPtma,]>4  
# keepLow  <- cellcondition$isgood & cellcondition$isss & log_ncount["mcherry",]<2
# keepHigh <- cellcondition$isgood & cellcondition$isss & log_ncount["mcherry",]>3
tnPtmaLow <- fitTechnicalNoise(
  ncount[,keepLow],
  fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 
tnPtmaHigh <- fitTechnicalNoise(
  ncount[,keepHigh],
  fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 

#Customized, only makes sense for log
getVariableGenes2 <- function(nCountsEndo, fit, method = "fit", threshold = 0.1, fit_type = NULL, plot = T) {
  LCountsEndo <- log10(nCountsEndo + 1)
  LmeansEndo <- rowMeans(LCountsEndo)
  Lcv2Endo = rowVars(LCountsEndo)/LmeansEndo^2
  score = fit$opts$offset * coefficients(fit)["a"] * 10^(-coefficients(fit)["k"] * LmeansEndo)
  is_het = score < Lcv2Endo & LmeansEndo > fit$opts$minmean
  ret = score
  ret[!is_het] <- 10000
  ret
}


hetPtmaLow  <- getVariableGenes2(fit_type = 'log',threshold=0.2, ncount[,keepLow],  tnPtmaLow$fit, plot=FALSE)
hetPtmaHigh <- getVariableGenes2(fit_type = 'log',threshold=0.2, ncount[,keepHigh], tnPtmaHigh$fit, plot=FALSE)

length(which(hetPtmaLow<1000))
length(which(hetPtmaHigh<1000))
diffHetHigh <- hetPtmaHigh[hetPtmaHigh<1000 & hetPtmaLow>1000]
diffHetLow <- hetPtmaLow[hetPtmaHigh>1000 & hetPtmaLow<1000]
length(diffHet)
#diffHetLog <- log(hetPtmaHigh) - log(hetPtmaLow)  #

hetdf <- function(diffHet, n=100){
  diffHet <- diffHet[order(diffHet)]
  diffHet <- diffHet[diffHet>=0]
  data.frame(sym=togenesym(names(diffHet)[1:n]),score=diffHet[1:n],me=gene_mean_ss[names(diffHet)[1:n]],stringsAsFactors=FALSE)
}

hetdf(hetPtmaLow)
hetdf(hetPtmaHigh)
hetdf(diffHet)
hetdf(diffHetLow)

#hm. does mitochondrial cutoff make sense?

length(which(gene_mean_ss>0.00001))

#Check enrichment
hetAlliffHet<-stopgosym(hetdf(diffHet)$sym) #crap
go_hetPtmaHigh<-stopgosym(hetdf(hetPtmaHigh)$sym)

go_diffHetLow<-stopgosym(hetdf(diffHetLow)$sym)
diffHetLow

listmosthet <- names(which(hetAll))
x<-stopgosym(listmosthet)
head(x)
x[1:30,]
#hetdf(diffHetLog) #method does not work well


length(listkedarcellcycle$ensid)
hetCChigh <- listkedarcellcycle$ensid[which(listkedarcellcycle$ensid %in% names(hetPtmaHigh[hetPtmaHigh<1000]))]
hetCClow <- listkedarcellcycle$ensid[which(listkedarcellcycle$ensid %in% names(hetPtmaLow[hetPtmaLow<1000]))]
togenesym(setdiff(hetCChigh, hetCClow))
togenesym(setdiff(hetCClow, hetCChigh))
length(togenesym(intersect(hetCClow, hetCChigh)))

length(setdiff(hetCChigh, hetCClow))
length(setdiff(hetCClow, hetCChigh))
hetCCall <- listkedarcellcycle$ensid[which(listkedarcellcycle$ensid %in% names(hetAll[hetAll<1000]))]
length(hetCCall)

hetA

hetPtmaHigh[hetPtmaHigh<1000]

######################################################################
### Plot heatmaps ####################################################
######################################################################

#Plot a heatmap, show levels on top
pheatmapbylevel <- function(geneids, keepcells=cellcondition$levelnum>0,normalize=FALSE){
  #keepcells=cellcondition$levelnum>0
  #geneids <- listolaptmagenes7d$ensid
    
  lcounts <- log_ncount
  selectgenes <- which(rownames(lcounts) %in% geneids)
  
  lcounts <- lcounts[selectgenes,keepcells]
  if(normalize){
    for(i in 1:nrow(lcounts)){
      s <- max(lcounts[i,])
      if(s==0) s<-1
      lcounts[i,] <- lcounts[i,]/s
    }
  }

  rownames(lcounts) <- togenesym(rownames(lcounts))
  annotation_col <- data.frame(levelnum=factor(sprintf("level%s",cellcondition$levelnum[keepcells])))  
  #annotation_col
  annotation_colors <- list(levelnum=c(
    level1="#000000",level2="#005500",level3="#00AA00",level4="#00FF00"))
  pheatmap(
    lcounts, 
    annotation_col=annotation_col,#gp = gpar(fill = "green")
    annotation_colors=annotation_colors,
    cluster_rows=TRUE, show_rownames=TRUE,
    cluster_cols=TRUE, show_colnames=FALSE)
}


pheatmapbygene <- function(geneids, keepcells=cellcondition$levelnum>0,normalize=FALSE,byid=ensidPtma){
  # keepcells=cellcondition$iss & cellcondition$isgood
  # geneids <- listkedarcellcycle$ensid
  # normalize <- FALSE
  
  lcounts <- log_ncount
  selectgenes <- which(rownames(lcounts) %in% geneids)
  
  lcounts <- lcounts[selectgenes,keepcells]
  if(normalize){
    for(i in 1:nrow(lcounts)){
      s <- max(lcounts[i,])
      if(s==0) s<-1
      lcounts[i,] <- lcounts[i,]/s
    }
  }
  
  cols <- greenred(100)
  rownames(lcounts) <- togenesym(rownames(lcounts))
  v<-log_ncount[byid,keepcells]
  scolor <- cols[ceiling(100*v/max(v))]  #rgb(0, v/max(v), 0, maxColorValue=1)

  heatmap.2(lcounts,trace="none",col=cols, ColSideColors = scolor, density.info="none", labCol = FALSE)
}


#  alsogenes <- c("Ptma")#sprintf("grna-ptma%s",1:5) #"cas9"
#                         c(ensidPtma,listolaptmagenes7d$ensid,alsogenes))
#"cas9"

# listannoyingcells <- c(126,109,165,144,134,120,131)
# cellcondition[listannoyingcells,]
# cellcondition$isss[listannoyingcells]<-FALSE

keepcells <- cellcondition$isss & cellcondition$isgood #cellcondition$levelnum>0
#keepcells <- cellcondition$iscr & cellcondition$isgood & cellcondition$levelnum>0

pheatmapbylevel(c(ensidPtma,listolaptmagenes7d$ensid),keepcells)
png("plots/heatmap pluripot ss.png",width=800)
pheatmapbygene(c(listolapluripotency$ensid),keepcells)
dev.off()

png("plots/heatmap 2c ss.png",width=800)
pheatmapbygene(c(listola2c$ensid),keepcells)
dev.off()

#pheatmapbylevel(c(names(which(hetAll))),keepcells)

png("plots/heatmap cellcycle ola ss.png",width=800)
#pheatmapbylevel(c(ensidPtma,listolacellcycle$ensid),keepcells)
pheatmapbygene(c(listolacellcycle$ensid),keepcells,normalize = TRUE)
dev.off()

png("plots/heatmap cellcycle kedar ss.png",width=800)
pheatmapbygene(c(listkedarcellcycle$ensid),keepcells,normalize = TRUE)
# pheatmapbygene(c(listkedarcellcycle$ensid[listkedarcellcycle$cellstage=="S"]),keepcells,normalize = TRUE)
# pheatmapbygene(c(listkedarcellcycle$ensid[listkedarcellcycle$cellstage=="G2M"]),keepcells,normalize = TRUE)
# pheatmapbygene(c(listkedarcellcycle$ensid[listkedarcellcycle$cellstage=="G1"]),keepcells,normalize = TRUE)
dev.off()

pheatmapbygene(c(listola2c$ensid),keepcells,normalize = TRUE)

#Ptma low -> most cell cycle markers low. G0?

pheatmapbylevel(c(listolacellcycle$ensid),keepcells)

pheatmapbylevel(c(ensidPtma,listola2c$ensid),keepcells,normalize = TRUE)
pheatmapbylevel(c(listola2c$ensid),keepcells)

#some really striking cells here for cell cycle!
#listolapluripotency$ensid

#which(order(gene_count)<96+40)


### todo: a k-means function as well. but can also sum up the gene exp for a list of them here
# isss: remove the 2 upper wells. same for iscr. actually combine with isgood?
# then plot. 9 ss cells show up low in every metric

#what is a good way of keeping track of cell definitions?

pheatmapcor <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount))){
  
  lcounts <- log_ncount
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
### t-sne ############################################################
######################################################################


#faraz:  6832
#cycles: 1ug: 


runrtsne <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount)),perplexity=30,max_iter=1000){
  #  geneids <- c(ensidPtma,listolacellcycle$ensid)
  #  keepcells <- cellcondition$isss & cellcondition$levelnum %in% c(1,4)
  lcounts <- log_ncount
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
plottsne<-function(rtsne_out,col=colorbylevel(),title="",labels=NULL){
  plot(rtsne_out$Y, col=col[rtsne_out$keepcells], pch=16, main='',xlab=title,ylab="")
  if(!is.null(labels)){
    text(rtsne_out$Y, labels=labels)
  }
#  plot(rtsne_out$Y, col=col[keepcells], pch=16, main='')
}

#runrtsne(c(ensidPtma,listolacellcycle$ensid),cellcondition$isss)


#lots of blue in one group. all PTMA downreg in one group. controls 2c!
png("plots/tsne ss 2c.png")
x<-runrtsne(c(listola2c$ensid),cellcondition$isss & cellcondition$isgood, perplexity=30)  #3 clusters with perp30
#TODO: pull out cluster of cell for proper analysis? maybe no need. but can use this to define PTMA cutoff!
#ALSO: what is up with the stretch of cells going to the big cluster outside?
par(mfrow=c(1,2))
plottsne(x,col=colorbygene(),"Ptma")
plottsne(x,col=colorbygene("mcherry"),"Mcherry")
title("tSNE 2c",outer=TRUE)
dev.off()
#x<-runrtsne(c(listola2c$ensid),cellcondition$isgood, col=colorbychem())   #clearly separates by chemistry!

#x<-runrtsne(c(listola2c$ensid),cellcondition$isss, col=colorbygene("mcherry"),perplexity = 30)  

#lots of blue in one group. all PTMA downreg in one group. controls pluripotency!
x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood & cellcondition$isss)  
y<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood & cellcondition$iscr)
png("plots/tsne pluripot.png")
par(mfrow=c(1,2))
plottsne(x,col=colorbygene(),"Ptma SS")
#plottsne(x,col=colorbymedia(),"Ptma SS")
#plottsne(x,col=colorbygene("mcherry"),"mcherry SS")
plottsne(y,col=colorbygene(),"Ptma DogSeq")
title("tSNE pluripotency",outer=TRUE)
dev.off()

png("plots/tsne 2c media.png")
par(mfrow=c(1,1))
x<-runrtsne(c(listola2c$ensid),cellcondition$isgood & cellcondition$isss)  
plottsne(x,col=colorbymedia(),"Ptma SS 2c colored by media")
dev.off()


png("plots/tsne all pluripot.png") #TODO - separates by chemistry
x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood)  
plottsne(x,col=colorbychem(),"Chemistry (pluripotency genes)")
#x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood, col=colorbychem())  
dev.off()

png("plots/tsne ss pluripot growth cond.png") #TODO - separates by chemistry
x<-runrtsne(c(listolapluripotency$ensid),cellcondition$isgood & cellcondition$isss)  
plottsne(x,col=colorbymedia(),"Growth condition (pluripotency genes)")
dev.off()



#lots of blue in one group
#PTMA level is spread over cell cycle. good! so not fluctuating then
#actually I see some correlation here
x<-runrtsne(c(listolacellcycle$ensid),cellcondition$isgood & cellcondition$isss)  
png("plots/tsne ss cellcyc ola.png")
par(mfrow=c(1,2))
plottsne(x,col=colorbygene(),"Ptma")
plottsne(x,col=colorbygene("mcherry"),"")
title("tSNE Ptma SS (ola cell cycle)",outer=TRUE)
dev.off()

#cell cycling, kedar style
x<-runrtsne(c(listkedarcellcycle$ensid),cellcondition$isgood & cellcondition$isss)  
png("plots/tsne ss cellcyc kedar.png")
par(mfrow=c(1,3))
plottsne(x,col=colorbycc(),"Cell cycle")
plottsne(x,col=colorbygene(),"Ptma")
plottsne(x,col=colorbygene("mcherry"),"")
title("tSNE Ptma SS (kedar cell cycle)",outer=TRUE)
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
plot(log1(ncount["mcherry",keepcells]), #good correlation but up to 10-100 fold difference
     log1(colSums(ncount[grep("grna-",rownames(ncount)),keepcells])),cex.lab=1)
dev.off()
#cas9 vs bfp is crap. this is a T2A construct. how can correlation be so bad?? likewise in their droplet data?
png("plots/cr20k bfp vs cas9.png",w=800)
keepcells <- cellcondition$iscr & gene_count>20e3 & cellcondition$levelnum>0
plot(log1(ncount["bfp",keepcells]), 
     log1(ncount["cas9",keepcells]),cex.lab=1)
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
  tt <- log_ncount[keepgene,fac>0]
  mean1<- rowMeans(log_ncount[keepgene,fac==1])
  mean2<- rowMeans(log_ncount[keepgene,fac==2])
  out <- rowttests(as.matrix(tt), factor(fac[fac>0]))
  #out[order(out$p.value),]
  out <- cbind(rownames(out),mean1,mean2,out)
  colnames(out)[1] <- "ensembl_gene_id"
  colnames(out)[2] <- "mean1"
  colnames(out)[3] <- "mean2"
  sqldf("select * from out natural join ensconvert")
}
orderttest <- function(out, pcutoff=1) {
  out<-out[order(out$p.value),]
  out[out$p.value<pcutoff,]
}
## differences between levels of grna (from facs)
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




#what is the cut-off high and low Ptma? not obvious. can be at 1. or up to 2-3. likely 3-4
hist(as.double(log_ncount[ensidPtma,]))

#Compare ptma low/high
tlowhiPtma <- myfastt(
  log_ncount[ensidPtma,]<3  & cellcondition$isgood & cellcondition$isss,
  log_ncount[ensidPtma,]>4  & cellcondition$isgood & cellcondition$isss)  #was cut off 1, for all mito
orderttest(tlowhiPtma)[1:100,]



# tlowhiPtmaNowt <- myfastt(
#   log_ncount[ensidPtma,]<1  & cellcondition$isgood & cellcondition$isss & cellcondition$levelnum>1,
#   log_ncount[ensidPtma,]>=1 & cellcondition$isgood & cellcondition$isss & cellcondition$levelnum>1)
# orderttest(tlowhiPtmaNowt)[1:100,]


#TODO GO analysis?

write.table(tlowhiPtma$ensembl_gene_id[1:200], file="gotest_lowhiptma.csv")

#eif4a2?

#mean(out$mean1)

#x<-mytopgosym(listolacellcycle$genesym)

go_lowhiptma<-stopgosym(orderttest(tlowhiPtma)[1:100,7])
go_lowhiptma

orderttest(tlowhiPtma,1e-9)


######################################################################
### gRNA counts ######################################################
######################################################################

#TODO check again, which grna molecules do we have? can we look at the promoter in more detail?

### what are the relative fractions of different gRNAs? Here for smartseq
#grnatotc <- rowMeans(ncount[grep("grna-ptma",rownames(ncount)),cellcondition$iss & cellcondition$isgood])
#grnabg <- mean(rowMeans(ncount[setdiff(grep("grna-",rownames(ncount)),grep("grna-ptma",rownames(ncount))),]))
plotgrnapie <- function(ncount){
  grnatotc <- rowMeans(ncount[grep("grna-ptma",rownames(ncount)),])
  grnatotc <- grnatotc/sum(grnatotc)
#  png("plots/grna-ptmafraction SS2.png",w=800)
  #barplot(grnatotc, main="PTMA gRNA fractions", cex.names=0.5,beside=TRUE) 
  pie(grnatotc, labels = paste(paste(sprintf("grna-%s",1:5),round(grnatotc*100)),"%",sep=""),cex=1)
 # dev.off()
  grnatotc
}
plotgrnabar <- function(ncount){
  grnatotc <- rowMeans(ncount[grep("grna-ptma",rownames(ncount)),])
  grnatotc <- grnatotc/sum(grnatotc)
  names(grnatotc) <- paste(paste(sprintf("grna-%s",1:5),round(grnatotc*100)),"%",sep="")
  barplot(grnatotc, main="", cex.names=0.5,beside=TRUE) 
#  pie(grnatotc, labels = paste(paste(sprintf("grna-%s",1:5),round(grnatotc*100)),"%",sep=""),cex=1)
#  grnatotc
}
png("plots/grnaptmafraction pie.png",w=800)
par(mfrow=c(1,2),oma=c(0,0,2,0))
plotgrnapie(ncount[,cellcondition$isss & cellcondition$isgood])
plotgrnapie(ncount[,cellcondition$iscr & cellcondition$isgood & cellcondition$plate %in% c(9,10)])
#plotgrnapie(ncount[,cellcondition$iscr & cellcondition$isgood])
title("Overall gRNA presence ss2 vs DogSeq",outer=TRUE)
dev.off()

png("plots/grnaptmafraction barplot dogseq.png",w=800)
par(mfrow=c(1,1),oma=c(0,0,2,0))
plotgrnabar(ncount[,cellcondition$iscr & cellcondition$isgood & cellcondition$plate %in% c(9,10)])
dev.off()



# i<-2
# par(mfrow=c(1,2))
# plot(sort(as.double(log_ncount[grep(sprintf("grna-ptma%s",i),rownames(ncount)),cellcondition$isss & cellcondition$isgood])))
# plot(sort(as.double(log_ncount[grep(sprintf("grna-ptma%s",i),rownames(ncount)),cellcondition$iscr & cellcondition$isgood])))

#### Compare the dynamic range of gRNA detection
plotsortgrna <- function(log_ncount){
  out <- sort(as.double(log_ncount[grep(sprintf("grna-ptma%s",1),rownames(ncount)),]))
  for(i in 2:5){
    out <- cbind(out,sort(as.double(log_ncount[grep(sprintf("grna-ptma%s",i),rownames(ncount)),])))
  }  
  cols=c("black","red","green","blue","gray")
  plot(out[,1],type="l",ylim=c(0,max(out)),xlab="cell",ylab="log relative gRNA counts")
  for(i in 2:5){
    lines(out[,i],col=cols[i])
  }  
}
png("plots/grna dynamic range.png",w=800)
par(mfrow=c(1,2),oma=c(0,0,2,0))
plotsortgrna(log_ncount[,cellcondition$isss & cellcondition$isgood])
plotsortgrna(log_ncount[,cellcondition$iscr & cellcondition$isgood])
title("Dynamic range ss2 vs DogSeq",outer=TRUE)
dev.off()

#Result: Because dogseq seem quite consistent without any bumps, that's the one to go for. there is definitely a bias in gRNA presence

######################################################################
### Upstream screening ###############################################
######################################################################

#or constrain by mcherry! 
takecells_ko <- log_ncount[ensidPtma,]<3.5  & cellcondition$isgood & cellcondition$isss
takecells_wt <- log_ncount[ensidPtma,]>=3.5 & cellcondition$isgood & cellcondition$isss
takecells_all <- log_ncount[ensidPtma,]>1 & cellcondition$isgood & cellcondition$isss

# takecells_wt <- cellcondition$levelnum %in% c(3,4) & cellcondition$isss & cellcondition$isgood 
# takecells_ko <- cellcondition$levelnum %in% c(1,2) & cellcondition$isss & cellcondition$isgood 

# takecells_wt <- cellcondition$levelnum %in% c(3,4) & cellcondition$isss & cellcondition$isgood 
# takecells_ko <- cellcondition$levelnum %in% c(1,2) & cellcondition$isss & cellcondition$isgood 
# takecells_all <- cellcondition$levelnum %in% c(1,2) & cellcondition$isss & cellcondition$isgood 
cor_wt <- cor(t(ncount[gene_mean_ss>=0.0001 & rownames(ncount) %in% listmosthet,takecells_wt]), method="spearman")  
cor_ko <- cor(t(ncount[gene_mean_ss>=0.001 ,takecells_ko]), method="spearman") 
cor_all <- cor(t(ncount[gene_mean_ss>=0.001 ,takecells_all]), method="spearman") 
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
plot(cor_diff$value)
#### note. dangerous if we don't deal with homologs. should use sailfish or filter out duplicates in final count table
### but should not filter out the grnas comparing their counts! 

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

#TOO get rid of Gm9800 and Gm4617

downstreamPtma <- tlowhiPtma$ensembl_gene_id[which(tlowhiPtma$p.value<1e-10)]
downstreamPtma
length(downstreamPtma)

col <- rep("black",nrow(cor_diff))
col[cor_diff$Var1 %in% downstreamPtma | cor_diff$Var2 %in% downstreamPtma] <- "red"
col[cor_diff$Var1 %in% listola2c$ensid | cor_diff$Var2 %in% listola2c$ensid] <- "red"
col[cor_diff$Var1 %in% listolapluripotency$ensid | cor_diff$Var2 %in% listolapluripotency$ensid] <- "red"
plot(cor_diff$value, col=col,cex=0.5)

#instead of cor, highly varying genes in mcherry or ptma +? not quite the same but...

#things in high correlation with Ptma 
#cor_wt


#https://en.wikipedia.org/wiki/MALAT1-associated_small_cytoplasmic_RNA
#in nuclear speckle domains

######################################################################
### Genes changing in correlation depending on Ptma level ############
######################################################################


plotcor2 <- function(gene1,gene2,n1="g1",n2="g2"){
  takecells_ko <- log_ncount[ensidPtma,]<2  & cellcondition$isgood & cellcondition$isss
  takecells_wt <- log_ncount[ensidPtma,]>=3 & cellcondition$isgood & cellcondition$isss
  par(mfrow=c(1,2))
  plot(log1(ncount[gene1,takecells_wt]),
       log1(ncount[gene2,takecells_wt]),cex.lab=1,xlab=n1,ylab=n2)
  plot(log1(ncount[gene1,takecells_ko]),
       log1(ncount[gene2,takecells_ko]),cex.lab=1,xlab=n1,ylab=n2)
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
# assignmentsLow  <- cyclone(as.matrix(ncount[,ncount[ensidPtma,]<3.5  & cellcondition$isgood & cellcondition$isss]),mm.pairs)
# assignmentsHigh <- cyclone(as.matrix(ncount[,ncount[ensidPtma,]>3.5  & cellcondition$isgood & cellcondition$isss]),mm.pairs)
#PTMA 2 and 3 - original cut-offs
cycass <- cyclone(as.matrix(ncount),mm.pairs)

plotphase <- function(phases){
  cnt <- c(
    length(which(phases=="G1")),
    length(which(phases=="G2M")),
    length(which(phases=="S"))
  )
  pie(cnt,labels=c("G1","G2M","S"))
  
  #could sort and make a heatmap too?
  cnt
}

#more cells in S-phase with PTMA high
png("plots/cycstage ptma.png",width=800)
par(mfrow=c(1,2))
# plotphase(assignmentsLow$phases)  #difference huge at minimum to ptma 3.5
# plotphase(assignmentsHigh$phases)  
plotphase(cycass$phases[ncount[ensidPtma,]<3.5  & cellcondition$isgood & cellcondition$isss])
plotphase(cycass$phases[ncount[ensidPtma,]>3.5  & cellcondition$isgood & cellcondition$isss])
dev.off()


##################### In each group, Ptma level?
png("plots/ptma in different phases.png")
v1 <-ncount[ensidPtma,cycass$phases=="S" & cellcondition$isss & cellcondition$isgood]
v2 <- ncount[ensidPtma,cycass$phases=="G2M" & cellcondition$isss & cellcondition$isgood]
v3 <- ncount[ensidPtma,cycass$phases=="G1" & cellcondition$isss & cellcondition$isgood]
ox<-c(rep(1,length(v1)),rep(2,length(v2)),rep(3,length(v3)))
plot(ox+runif(length(ox))*0.4,log1(c(v1,v2,v3)),xlab="S                                 G2M                                G1",ylab="Ptma")
dev.off() 
mean(v1)
mean(v2)
mean(v3)
t.test(log1(v1),log1(v3))

#get some measure of speed
#higher if ptma high. so ptma high means less cycling
# mean(apply(assignmentsLow$normalized.scores,1,max)) 
# mean(apply(assignmentsHigh$normalized.scores,1,max))   #higher confidence here
#TODO continous measure

#S score vs PTMA level?
# png("plots/phase ptma vs S.png")
# plot(cycass$normalized.scores$S,
#      log1(ncount[ensidPtma,keepcells]))
# dev.off()
# png("plots/phase ptma vs G2M.png")
# plot(cycass$normalized.scores$G2M,
#      log1(ncount[ensidPtma,keepcells]))
# dev.off()
# png("plots/phase ptma vs G1.png")
# plot(cycass$normalized.scores$G1,
#      log1(ncount[ensidPtma,keepcells]))
# dev.off()

#S and G2M scores lower in general for low PTMA. suggests faster cycling when PTMA low. again
#If S phase longer, does that mean high PTMA keep cells in synthesis? all ribosomal genes higher with high ptma.
#would this not suggest a shorter synth phase? unless the cell is struggling to up synthesis so proteins are not higher.
#or they are high but nothing happens

keepcells <- cellcondition$isgood & cellcondition$isss

#Cell cycling does not seem to change anymore
png("plots/phase avscore sorted by PTMA.png", width=800)
par(mfrow=c(1,1))
avscore <- apply(cycass$scores[keepcells,],1,mean)[order(ncount[ensidPtma,keepcells])]
k<-60
avscoreroll <- rollmean(avscore,k)
plot(avscore,ylab="Average cell cycling score",pch=16,col="green")
lines((1:length(avscoreroll)+k/2),avscoreroll,col="blue")
dev.off()

#Would be nice to sort by Ptma level, then show the 3 classifications
#G1 vs G2M, color coded by Ptma. Black circles mainly in the mid=S
png("plots/phase G1 vs G2M colored by PTMA.png", width=800)
col<-rep("black",ncol(ncount))
col[ncount[ensidPtma,]>median(ncount[ensidPtma,keepcells])] <- "green"
plot(cycass$normalized.scores$G1[keepcells],
     cycass$normalized.scores$G2M[keepcells], col=col, pch=16,cex=.5,xlab="G1",ylab="G2M")
dev.off()

#Continous diagram of cell cycle stage vs ptma level
plotphaseimage <- function(k=4){
  o<-order(ncount[ensidPtma,keepcells])
  v <- cycass$normalized.scores[keepcells,][o,]
  #base <- cumsum(rep(0.2,))
  #plot(cumsum(v$G1),type="l",ylim=c(0,200))  
  #lines(cumsum(v$G2M),col="red")
  #lines(cumsum(v$S),col="green")
  #print(dim(v))
  v <- cbind(
    rollmean(v[,1],k),
    rollmean(v[,2],k),
    rollmean(v[,3],k))
  m <- matrix(3,nrow=nrow(v),ncol=100)
  for(i in 1:nrow(v)){
    cut1<-round(v[i,1]*100) #G1
    cut2<-round((v[i,1]+v[i,3])*100) #G2M
    m[i,1:cut1]<-1
    m[i,cut1:cut2]<-2
  }
  image(1:nrow(m),1:ncol(m),m,col=c("#000000","#FF0000","#0000FF"),xlab="Cell index (smoothened)",ylab="G1/G2M/S  Log(Ptma)")
}
png("plots/phase image vs Ptma.png", width=800)
plotphaseimage(8)
lval <- sort(log(1+ncount[ensidPtma,keepcells]))
lines(lval*100/max(lval),col="yellow")
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
      avglevel[j,i] <- mean(as.double(ncount[which(rownames(ncount)==forgenes[j]),cellcondition$iscr & cellcondition$levelnum==i]))
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
  cellconditiongrna <- t(ncount[
    c(which(rownames(ncount) %in% sprintf("grna-ptma%s",1:5)),
      which(rownames(ncount) %in% c("cas9"))),])
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
lcounts <- log_ncount
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




##############################################

#rename a few variables


#construct and initialize new scLVM object
Y <- log_ncount[,cellcondition$isgood & cellcondition$plate==1]
genes_het_bool <- as.vector(hetAll) 
geneID <- rownames(ncount) 
tech_noise <- as.vector(tnAll$techNoiseLog) #technical noise
sclvm <- new("scLVM")
sclvm <- init(sclvm,Y=Y,tech_noise = tech_noise)

ens_ids_cc <- getEnsembl('GO:0007049')

CellCycleARD = fitFactor(sclvm,geneSet = ens_ids_cc, k=20,use_ard = TRUE)
#Get cell-cycle factor
Kcc = CellCycle$K
Xcc = CellCycle$X
image(Kcc,xaxt = "n", yaxt = "n", xlab = 'cells', ylab = 'cells')
title('Similarity matrix based on cell cycle')

sclvm = varianceDecomposition(sclvm, K=Kcc, idx = idx_het)

# get variance components
results_var = getVarianceComponents(sclvm)
var_filtered = results_var$var[results_var$conv,] # filter out genes for which vd has not converged
head(var_filtered)

# get corrected expression levels
Ycorr = getCorrectedExpression(sclvm)
dim(Ycorr)

var_mean = apply(var_filtered,2,mean)
colors = c('Green','Blue','Gray')
pie(var_mean, , col = colors)


#################################################################



