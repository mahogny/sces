colorbyko <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$ko=="Thy1control"]<-"green"
  #    thecol[cellgene=="Xbp1"]<-"red"
  #    thecol[cellgene=="F2rl1"]<-"red"
  thecol[cellcondition$ko=="Il2"]<-"red"
  thecol  
}

thcompare2 <- function(indexCheck, indexC=which(cellcondition$ko=="Thy1control" & cellcondition$isgood)){
  datfordeseq  <- dat[grep("ENSMUSG",rownames(dat)),]
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

thshowres <- function(res,N=min(50,nrow(res))){
  resOrdered <- res[order(res$padj),]
  #resOrdered
  data.frame(sym=togenesym(rownames(resOrdered[1:N,])),padj=resOrdered$padj[1:N],stringsAsFactors=FALSE)
}

#list of all KOed genes
unique(cellcondition$ko)




resXbp <- thcompare2(which(cellcondition$ko=="Xbp1" & cellcondition$isgood))
thshowres(resXbp)

resIl4 <- thcompare2(which(cellcondition$ko=="Il4" & cellcondition$isgood))
thshowres(resIl4)
# 1            Mt1 0.003299927
# 2        Tnfrsf9 0.035100147
# 3           Gzme 0.066055273
# 4          Il1r2 0.149377093
# 5           Glrx 0.182959546

resIl2 <- thcompare2(which(cellcondition$ko=="Il2" & cellcondition$isgood))
thshowres(resIl2)
# 1            Mt1 0.0001656131
# 2        Tnfrsf9 0.0004647096
# 3            Mt2 0.1625141649





resErn1 <- thcompare2(which(cellcondition$ko=="Ern1" & cellcondition$isgood))
thshowres(resErn1)
# 1        Gm10800 0.07426989
# 2          Rn7s1 0.07426989
# 3        Gm10801 0.23049494
# 4        Gm26870 0.23049494
# 5        Gm10717 0.23049494
# 6        Gm21738 0.23049494
# 7        Gm10720 0.23319601
# 8           Ifng 0.23319601
# 9        Gm10715 0.39286539
# 10       Fam107b 0.52934105
# 11       Gm10718 0.53291564


resIl13 <- thcompare2(which(cellcondition$ko=="Il13" & cellcondition$isgood))
thshowres(resIl13)
#crap


resEtv2 <- thcompare2(which(cellcondition$ko=="Etv2" & cellcondition$isgood))
thshowres(resEtv2)
#crap




resLag3 <- thcompare2(which(cellcondition$ko=="Lag3" & cellcondition$isgood))
thshowres(resLag3)
# 1   Gm10720 0.004177696
# 2   Gm17535 0.005059304
# 3   Gm10800 0.031526989
# 4    Slc2a3 0.031526989
# 5   Gm26870 0.031526989
# 6   Gm21738 0.031526989
# 7      Gzme 0.031526989
# 8    Lmbrd1 0.032973689
# 9   Gm11168 0.032973689
# 10  Gm10715 0.032973689
# 11  Gm10717 0.039047484
# 12  Gm10801 0.041399058
# 13  Gm10718 0.041399058
# 14     Ly6a 0.041399058
# 15   Dnajc8 0.056890698
# 16  Gm10719 0.062021002
# 17 Ebna1bp2 0.062794883
# 18     Cycs 0.062794883
# 19   Eif4a1 0.062794883
# 20 Mapkapk2 0.067296311
# 21  Gm10123 0.067296311
# 22  Gm10722 0.069528974
# 23     Tma7 0.092053055
# 24    H2-K1 0.092053055
# 25   Gm4149 0.094170099
# 26    Aldoa 0.099306955
# 27  Gm10721 0.099306955
# 28      Pkm 0.099306955
# 29     Eef2 0.099306955

resF2 <- thcompare2(which(cellcondition$ko=="F2rl1" & cellcondition$isgood))
thshowres(resF2)
#crap


#indexEtv2 <- which(cellgene=="Zc3h12a") #??

resC <- thcompare2(which(cellcondition$ko=="Cxcr7" & cellcondition$isgood))
thshowres(resC)
#crap


resZ <- thcompare2(which(cellcondition$ko=="Zc3h12a" & cellcondition$isgood))
thshowres(resZ)
#

#####################



thsne <- function(){
  takecells <- cellcondition$isgood & cellcondition$ko!="" &cellcondition$ko %in% c("Thy1control","Il4","Lag3")
  
  the_gene_mean <- rowMeans(ncount[,takecells])
  the_gene_mean <- the_gene_mean/sum(the_gene_mean)
  
  takegenes <- the_gene_mean>=0.00001 && !(1:nrow(dat) %in% grep("grna",rownames(dat)))
  x <- runrtsne(rownames(dat)[takegenes],takecells,perplexity = 10)
  plottsne(x,colorbyko())
}

