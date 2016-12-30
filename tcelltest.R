#from outset2b

genelay <- as.matrix(read.csv("layouttcelltest.csv"))
genelay <- cbind(genelay,genelay,genelay,rep("f",8),rep("f",8),rep("f",8),rep("f",8)) #not really true
genelay
mouserep <- c(1,1,1,2,2,2,3,3,3,0,0)
for(i in 1:3)
  mouserep <- rbind(mouserep,mouserep)
genelay
as.array(genelay)
#[13]



tocelli <- function(row,col){
  out <- c()
  for(ro in row){
    for(co in col){
      out <- c(out, 11*96 +  (ro-1)*12 + (co-1) + 1)
    }
  }
  out
}

acelli <- matrix(nrow = 8, ncol=12)
#arow <- matrix(nrow = 8, ncol=12)
for(i in 1:12){
  for(j in 1:8){
    acelli[j,i] <- tocelli(j,i)
  }
}
cellgene <- rep("",12*96)
for(i in 1:96){
  cellgene[(4+4+3)*96+i-1] <- genelay[cellcondition$row[i], cellcondition$col[i]]
}
#cellcondition

# 365 AGAGAGAAGTCCCCGCGC    fits lag3 as it should!

x <- 365
x <- x - 3*96

i <- floor((x-1)/12)
j <- x - i*12
j # is the column
i+1 # is the row

i*12+j




#indexCheck <- indexErn

indexCtrl <- which(cellgene=="Thy1control")
#only 303 and 321 have any grnas detected
#indexCtrl <- c(303,321)
indexCtrl <- c(306,309,315,318)


# Only use normal genes, raw counts. no grna or cas9 levels
#datfordeseq  <- dat[grep("ENSMUSG",rownames(dat)),]

indexCheck <- indexIl4

sum(dat[,indexIl4[1]])
colSums(dat[,indexIl4])[1:10]
add <- cbind(rowSums(dat[,indexIl4]),rowSums(dat[,indexCtrl]))[1:10,]
cellcond <- data.frame(levelnum=factor(c(1,2)))
dds <- DESeqDataSetFromMatrix(countData = add, colData = cellcond,design = ~ levelnum) 
dds <- dds[ rowSums(counts(dds)) > 0.1, ]
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
data.frame(sym=togenesym(rownames(resOrdered)),padj=resOrdered$padj,stringsAsFactors=FALSE)




compare2 <- function(indexCheck, indexC=indexCtrl){
  keepcells <- 1:ncol(dat) %in% c(indexC,indexCheck)
  cellcond <- data.frame(levelnum=1:ncol(dat) %in% c(indexC))
  cellcond <- cellcond[keepcells,,drop=FALSE]
  dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells], colData = cellcond,
                                design = ~ levelnum) 
  dds <- dds[ rowSums(counts(dds)) > 0.1, ]
  #dds$condition <- relevel(dds$condition, ref="1")   #set the reference level to be all 0?
  dds <- DESeq(dds)
  res <- results(dds)
  res
}

showres <- function(res,N=min(50,nrow(res))){
  resOrdered <- res[order(res$padj),]
  #resOrdered
  data.frame(sym=togenesym(rownames(resOrdered[1:N,])),padj=resOrdered$padj[1:N],stringsAsFactors=FALSE)
}



indexXbp <- which(cellgene=="Xbp1")
indexErn <- which(cellgene=="Ern1")
indexIl2 <- which(cellgene=="Il2")
indexIl4 <- which(cellgene=="Il4")
indexEtv2 <- which(cellgene=="Etv2")

dat[grep(toensid("Il4"),rownames(dat)),120:384]


resXbp <- compare2(indexXbp)
resErn <- compare2(indexErn)
resIl4 <- compare2(indexIl4)
resIl2 <- compare2(indexIl2)
resEtv2 <- compare2(indexEtv2)
#resXbp <- res
#ddsXbp <- dds

# resOrdered <- resXbp[order(resXbp$padj),]
# sampleGOdata <- new("topGOdata",description = "Simple session", ontology = "BP",
#                     allGenes = rownames(dat), geneSel = rownames(resOrdered)[1:20],
#                       nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
# resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

res2 <- compare2(c(313,325),c(303,315))
showres(res2)

res <- resIl2

showres(resXbp,2)
cat(paste(showres(resErn,100)$sym),sep="\n")
intersect(showres(resXbp,20)$sym,showres(resErn,20)$sym)
#Cdc45 and Gnai common??
cat(ensconvert$mgi_symbol,file = "background.txt")

showres(resEtv2,20)

showres(resIl4,20)
showres(resIl2,20)

v <- showres(resIl4)
v$sym[1:20]
  
#cbind(resIl2,resIl4,resErn,resXbp)
#["Il4",]

#51 is xbp1!
#327 empty
#330 is xbp1

which(rownames(resIl4)==toensid("Il4"))

toensid("Il4")

#poor mans normalization
# ncount_all <- dat
# for(i in 1:ncol(ncount_all)){
#   ncount_all[,i] <- ncount_all[,i]/sum(ncount_all[,i])
# }
# #colMeans(ncount_all)

### PCA
#dopca <- function(){
  #cor_ncount <- cor(ncount_all[gene_mean_ss>=0.00001,], method="spearman")
  takecells <- cellcondition$plate==4+4+4 & cellgene!="f" & cellcondition$row<=5 & colMeans(dat)>1 & !((1:ncol(dat)) %in% c(1075,1084))    #cellgene!="Ern1" & (1:n)

  takecells <- cellcondition$plate==4+4+4 & cellgene!="f" #& cellcondition$row<=5 & colMeans(dat)>1 #& !((1:ncol(dat)) %in% c(1075,1084))    #cellgene!="Ern1" & (1:n)
  
  colorbymouse <- function(){
    thecol <- rep("black",nrow(cellcondition))
    thecol[cellcondition$col %in% c(1,2,3)]<-"green"
    thecol[cellcondition$col %in% c(4,5,6)]<-"blue"
    thecol[cellcondition$col %in% c(7,8,9)]<-"red"
    thecol  
  }
  colorbygene <- function(){
    thecol <- rep("black",nrow(cellcondition))
    thecol[cellgene=="Thy1control"]<-"green"
#    thecol[cellgene=="Xbp1"]<-"red"
#    thecol[cellgene=="F2rl1"]<-"red"
    thecol[cellgene=="Il2"]<-"red"
    thecol  
  }
  
  the_mean_ss <- rowMeans(dat[,takecells])
  takegenes <- the_mean_ss>=0.00001 && !(1:nrow(dat) %in% grep("grna",rownames(dat)))
    #colMeans(dat[,takecells])
  cor_ncount <- cor(ncount_all[takegenes ,takecells], method="spearman")  #olas way
  #cor_ncount <- cor(log(1+dat)[the_mean_ss>=0.001 ,takecells], method="spearman")  #olas way
  #cor_ncount <- t(ncount_all[gene_mean_ss>=0.00001 ,takecells])  #ola does not scale
  pca <- prcomp(cor_ncount, scale=FALSE)
  
  
#  pdf("tcellPCA.pdf")
  ax1 <- 3
  ax2 <- 2
  par(cex.axis=1.5, cex.lab=1.5)
  plot(pca$x[,ax1], pca$x[,ax2], pch=5, col="white",xlab="PC1", ylab="PC2")
  text(pca$x[,ax1], pca$x[,ax2], cellgene[takecells], cex=0.5, col = colorbygene()[takecells])
  text(pca$x[,ax1], pca$x[,ax2], (1:ncol(dat))[takecells], cex=1, col = colorbymouse()[takecells])
  
  
  par(cex.axis=1.5, cex.lab=1.5)
  plot(pca$x[,3], pca$x[,4], pch=5, col="white",xlab="PC3", ylab="PC4")
  text(pca$x[,3], pca$x[,4], cellgene[takecells], cex=0.5)
  text(pca$x[,3], pca$x[,4], (1:ncol(dat))[takecells], cex=1)
  #  dev.off()

  pca$rotation[,1]
  
#  takegenes <- the_mean_ss>=0.00001 && !(1:nrow(dat) %in% grep("grna",rownames(dat)))
  runrtsne(rownames(dat)[takegenes],takecells,colorbygene())
  