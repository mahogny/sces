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
      out <- c(out, 3*96 +  (ro-1)*12 + (co-1) + 1)
    }
  }
  out
}

acelli <- matrix(nrow = 8, ncol=12)
#arow <- matrix(nrow = 8, ncol=12)
for(i in 1:12){
  for(j in 1:8){
    acelli[j,i] <- tocelli(i,j)
  }
}
cellgene <- rep("",384)
for(i in (3*96):(4*96)){
  cellgene[i] <- genelay[cellcondition$row[i], cellcondition$col[i]]
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




compare2 <- function(indexCheck){
  keepcells <- 1:ncol(dat) %in% c(indexCtrl,indexCheck)
  cellcond <- data.frame(levelnum=1:ncol(dat) %in% c(indexCtrl))
  cellcond <- cellcond[keepcells,,drop=FALSE]
  dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells], colData = cellcond,
                                design = ~ levelnum) 
  dds <- dds[ rowSums(counts(dds)) > 0.1, ]
  #dds$condition <- relevel(dds$condition, ref="1")   #set the reference level to be all 0?
  dds <- DESeq(dds)
  res <- results(dds)
  res
}

showres <- function(res,N=min(50,nrow(resOrdered))){
  resOrdered <- res[order(res$padj),]
  #resOrdered
  data.frame(sym=togenesym(rownames(resOrdered[1:N])),padj=resOrdered$padj[1:N],stringsAsFactors=FALSE)
}


indexCtrl <- which(cellgene=="Thy1control")
#only 303 and 321 have any grnas detected
#indexCtrl <- c(303,321)
indexCtrl <- c(306,309,315,318)

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

res <- resIl2

showres(resXbp,2)
showres(resErn,20)
intersect(showres(resXbp,20)$sym,showres(resErn,20)$sym)
#Cdc45 and Gnai common??

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


