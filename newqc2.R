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

######################################################################
### Load data ########################################################
######################################################################

listnongenes <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique","grna-numread")

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
  #rownames(dat)[grep("grna",rownames(dat))]
  
  #there is the -F- in grna-Ptma-F-4 ... and it is capital
  
  # # read grna counts
  # dat2 <- read.csv("counts/rawgrnam13.csv",stringsAsFactors=FALSE)
  # dim(dat2)
  # rownames(dat2) <- dat2[,1]
  # dat2 <- dat2[,-1]
  # rownames(dat2) <- sprintf("grna-%s",rownames(dat2))  #rename all to avoid clashes... not sure if this is a good idea
  # colnames(dat2) <- colnames(dat)
  # #would be better if this was done on the cluster side and written into one single file
  # 
  # dat <- rbind(dat,dat2)
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

sum(dat[,150])

### Construct plate design
getplatedesign20161209 <- function(){
  #TODO get rid of one plate on hiseq3 completely! and one is empty as well but can ignore this
  
  #hiseq1 & hiseq2: swapped B,C,D. 
  
  numplates <- 4+4+4

  fromhiseq <- c(rep(1,384),rep(2,384),rep(3,384)) #todo
    
  thecols <- rep(1:12,8)
  thecols <- rep(thecols,numplates)
  therows <- as.double(sapply(1:8,function(i) rep(i,12)))
  therows <- rep(therows,numplates)
  theplate <- as.double(sapply(1:numplates,function(i) rep(i,96)))
  level1<-thecols<12 & therows %in% c(1,2)
  level2<-therows %in% c(3,4)
  level3<-therows %in% c(5,6)
  level4<-therows %in% c(7,8)
  levelnum <- level1*1 + level2*2 + level3*3 + level4*4
  isgood <- levelnum>0
  iscr <- theplate %in% c(4,8)
  isss <- !iscr #not quite!
  cellcondition <- data.frame(
    col = thecols,
    row = therows,
    plate = theplate,
    iscr = iscr,
    isss = isss,
    level1=level1, level2=level2, level3=level3, level4=level4,
    levelnum=levelnum,
    isgood=isgood
  )
  rownames(cellcondition) <- colnames(dat)
  cellcondition
}
cellcondition <- getplatedesign20161209()


######################################################################
### Read olas data #####################################################
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


ensidPtma <- ensconvert$ensembl_gene_id[which(ensconvert$mgi_symbol=="Ptma")]

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

#pdf("qc_exonprop_mtprop_cell_2i_2.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(exon_prop[1,], mt_prop, pch=20, xlab="Proportion of reads mapped to exons",
     ylab="Proportion of reads mapped to mitochondrial genes", xlim=c(0,1), ylim=c(0,1))
abline(h=0.1, v=0.5, col="red", lty=2, lwd=2)
#dev.off()



#remove cells with too few reads (compare to mean for a plate instead?)
plot(gene_count)
plot(gene_count[cellcondition$plate %in% c(4,8,9,10)],ylim=c(0,200000))
length(which(gene_count[cellcondition$plate %in% c(4,8,9,10)]>100e3))
plot(gene_count[cellcondition$plate %in% c(1,2,3,4)])
plot(gene_count[cellcondition$plate %in% c(5,6,7,8)])
plot(gene_count[cellcondition$plate %in% c(6)])
mean(gene_count[cellcondition$isss])
mean(gene_count[cellcondition$plate %in% c(4)])
mean(gene_count[cellcondition$plate %in% c(8)])
cellcondition$isss[which(order(gene_count)<96+10)] <- FALSE

plot(grna_count)
plot(mt_prop)
#mean(gene_count[cellcondition$isss])
#cellcondition$isss[which(order(gene_count)<96+10)] <- FALSE


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
log_ncount_all <- log(1+ncount_all)


### Comparing average gene expression among conditions

#ncount_all <- dat[-which(rownames(dat3)=="grna-numread"),]
gene_mean_cr <- rowMeans(ncount_all[,cellcondition$iscr]) 
gene_mean_ss <- rowMeans(ncount_all[,cellcondition$isss])
gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)
gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
#head(gene_mean_cr)

plotgenemean <- function(){
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
dopca <- function(){
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
}




######################################################################
### DEseq model testing ##############################################
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
for(i in 1:ncol(b)) b[,i] <- b[,i]/max(b[,i]) #normalize to highest level
colnames(b) <- togenesym(colnames(b)) #order preserving? not any more!
png("plots/ssptmadown.png",width = 800,pointsize = 16)
par(cex.lab=0.2)
barplot(b, main="PTMA downreg: wt vs high",
        col=c("#000000","#880000","#AA0000","#FF0000"), cex.names=0.5,beside=TRUE) 
dev.off()

#ptma goes down as it should
#Ptp4a3 looks nice! and Dst and Peg10 and Parp12*. but peg10 & dst the opposite way before
#Parp12 seems affected by NfkB or interferon in general. could be affeced by the transfection itself


#Check some individual correlations
log1 <- function(x) as.double(log(1+x))
plot(log1(ncount_all[toensid("Ptp4a3"),]),
     log1(ncount_all[ensidPtma,]),cex.lab=1)




pheatmapbylevel <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount_all))){
  
  lcounts <- log_ncount_all
  selectgenes <- which(rownames(lcounts) %in% geneids)
  
  lcounts <- lcounts[selectgenes,keepcells]
  rownames(lcounts) <- togenesym(rownames(lcounts))
  annotation_col <- data.frame(levelnum=factor(sprintf("level%s",cellcondition$levelnum[keepcells])))  
  #annotation_col
  annotation_colors <- list(levelnum=c(
    level1="#00FF00",level2="#000000",level3="#0000FF",level4="#AA0000"))
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

listannoyingcells <- c(126,109,165,144,134,120,131)
cellcondition[listannoyingcells,]
cellcondition$isss[listannoyingcells]<-FALSE


pheatmapbylevel(c(ensidPtma,listolaptmagenes7d$ensid),cellcondition$isss)
pheatmapbylevel(c(ensidPtma,listolapluripotency$ensid),cellcondition$isss)
pheatmapbylevel(c(ensidPtma,listolacellcycle$ensid),cellcondition$isss)
pheatmapbylevel(c(listolacellcycle$ensid),cellcondition$isss)
pheatmapbylevel(c(ensidPtma,listola2c$ensid),cellcondition$isss)
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


runrtsne <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount_all))){
  #  geneids <- c(ensidPtma,listolacellcycle$ensid)
  #  keepcells <- cellcondition$isss & cellcondition$levelnum %in% c(1,4)
  lcounts <- log_ncount_all
  selectgenes <- which(rownames(lcounts) %in% geneids)
  lcounts <- lcounts[selectgenes,keepcells]
  #rownames(lcounts) <- togenesym(rownames(lcounts))
  
  d <- stats::dist(t(as.matrix(lcounts)))
  set.seed(0) # tsne has some stochastic steps (gradient descent) so need to set random 
  rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=10, verbose = FALSE)
  plot(rtsne_out$Y, col=colorbylevel()[keepcells], pch=16, main='tSNE')
  
  #function to color by gene exp? ptma? grna? read count?
  
  # https://hms-dbmi.github.io/scw/analysis-of-heterogeneity-and-subpopulations.html
  rtsne_out 
}
runrtsne(c(ensidPtma,listolacellcycle$ensid),cellcondition$isss)
runrtsne(c(listolacellcycle$ensid),cellcondition$isss)  #lots of blue in one group
runrtsne(c(listola2c$ensid),cellcondition$isss)  #lots of blue in one group
runrtsne(c(listolapluripotency$ensid),cellcondition$isss)  #lots of blue in one group
x<-runrtsne(c(listolacellcycle$ensid),cellcondition$isss)  #lots of blue in one group

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

######################################################################
### Monocle ###############################
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




################################################
##########################################################

#
#dat$S1_L001.count


#make plates
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
plateA <- makeplate(colSums(dat[grep("grna-",rownames(dat)),]), 4)
#plateB <- makeplate(colSums(dat[grep("grna-ptma",rownames(dat)),]), 8)

pheatmap(plateA,cluster_rows = FALSE, cluster_cols = FALSE, main = "DogSeq: gRNA count")
