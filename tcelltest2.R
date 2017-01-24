#289 and 290 empty
which(cellcondition$ko=="Xbp1")-96*8+1

d <- function(...){
  v <- unlist(list(...))-96*3
  print(v %% 12)
  print((v %/% 12)+1)
  invisible(v+96*11)
}


#xbp1
d(334)             # row 4. col 10      [3]
d(339,342,345,346) # row 5. col 3 6 9 10

#thy1ctrl
d(303,310) #row 2, col 3 and 10 [3]
d(322)     #row 3, col 10
cellcondition$ko[d(303,310,322)]
#here is one set of thycontrol, inferred

# 303, 306, 309
# 315, 318, 321

#il4. only one type found
d(353, 356) # row 6, column 5 and 8 [2].   and row 6?
#infer...
# 350, 353, 356

#lag3
d(365)           # row 7  col 5      [2]
d(373,374,380)   # row 8  col 1 2 8  [2] ????
#inferring...
# 362,365,368   
# 
#B: 367, 370,373,374

#il13
d(296)         # row 1  col 8
d(302,305,308) # row 2  col 2 5 8
cellcondition$ko[d(302,305,308)]
#inferring...
#
#

#Il2
d(314,320)
d(329)

#Ern1
d(328,331)
cellcondition$ko[d(328,331)]
which(cellcondition$ko=="Ern1")


th_fastt <- function(groupA, groupB, keepgene = gene_mean_ss>=0.00001){
  fac <- rep(0,ncol(dat))
  fac[groupA] <- 1
  fac[groupB] <- 2
  tt <- th_ncount[keepgene,fac>0]
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
# orderttest <- function(out, pcutoff=1) {
#   out<-out[order(out$p.value),]
#   out[out$p.value<pcutoff,]
# }


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



#xbpA In 334,   #most likely then also 
#xbpB In 339, 342, 345 + 346  #346... extra test? yes!

#xbpA to be expected in: 327 330 333 + 334
339-12
339-12-12

#note!!! different order of col/rows for plate 12 vs the other plates!! ... or not

cellcondition$ko<-""
cellcondition$ko[96*8-1+c(327,330,333,   339,342,345)]<-"Xbp1"

cellcondition$ko[96*8-1+c(303, 306, 309,   315,318,321)]<-"Thy1control"

cellcondition$ko[96*8-1+c(303, 306, 309,   315,318,321)] #what is it??


#96*8-1+c(327,330,333,   339,342,345)
resXbp <- thcompare2(which(cellcondition$ko=="Xbp1" & cellcondition$isgood))
thshowres(resXbp)

resIl4 <- thcompare2(which(cellcondition$ko=="Il4" & cellcondition$isgood))
thshowres(resIl4)
# 1            Mt1 0.003299927
# 2        Tnfrsf9 0.035100147
# 3           Gzme 0.066055273
# 4          Il1r2 0.149377093
# 5           Glrx 0.182959546
summary(gene_count[which(cellcondition$ko=="Il4" & cellcondition$isgood)])

resIl2 <- thcompare2(which(cellcondition$ko=="Il2" & cellcondition$isgood))
thshowres(resIl2)
# 1            Mt1 0.0001656131
# 2        Tnfrsf9 0.0004647096
# 3            Mt2 0.1625141649
summary(gene_count[which(cellcondition$ko=="Il2" & cellcondition$isgood)])





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
summary(gene_count[which(cellcondition$ko=="Il13" & cellcondition$isgood)])
gene_count[which(cellcondition$ko=="Il13" & cellcondition$isgood)]


resEtv2 <- thcompare2(which(cellcondition$ko=="Etv2" & cellcondition$isgood))
thshowres(resEtv2)
#crap



# 365 AGAGAGAAGTCCCCGCGC    fits lag3 as it should!

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
gene_count[which(cellcondition$ko=="Lag3" & cellcondition$isgood)]

resF2 <- thcompare2(which(cellcondition$ko=="F2rl1" & cellcondition$isgood))
thshowres(resF2)
#crap
gene_count[which(cellcondition$ko=="F2rl1" & cellcondition$isgood)]


#indexEtv2 <- which(cellgene=="Zc3h12a") #??

resC <- thcompare2(which(cellcondition$ko=="Cxcr7" & cellcondition$isgood))
thshowres(resC)
#crap


resZ <- thcompare2(which(cellcondition$ko=="Zc3h12a" & cellcondition$isgood))
thshowres(resZ)
#

#####################



tcellgenes <- c(
  "Il17a","Runx3","Ets1","Hlx","Junb","Gfi1","Ikzf1","Ikzf3","Rora","Stat3","Ahr","Ep300","Hif1a",
  "Gata3","Il1a","Il1b","Il2","Il4","Il5","Il12a","Il12b","Il13","Il10",
  
  "Atf4","Atf6",
  
  "Ilf2", #=Nf45
  "Ilf3", #=Nf90
  "Il2ra","Il2rb","Il2rg", #these form a complex with Il2
  "Spi1",
  "Runx1",
  "Nfil3", #bind to IL13 promoter  
  "Bcl2l1","Bcl6","Sp1",
  "Ap1","Maf",  #these two sit close to the IL4 locus. Ap1 bindings sites are all over in xi screen
  "Fos","Jun", #wikipeda, cFos and cJun these generally interact with Ap1 for cell proliferation
  
  "Gli2",
  #  "Tcrb",
  "Cd3d","Cd3e","Cd3g","Cd4","Cd44","Sell", #Sell=Cd62l
  "Cd40",
  "Cyp11a1",
  "Tbx21",
  "Ly6c1","Ly6c2",  "Ly9", #michal mentioned Ly9
  "Ly6e", # from marioni
  "Xbp1",
  "Ern1",
  "Jak1",
  "Nfatc1", "Nfatc2", "Nfatc3", "Nfatc4", "Nfat5",
  "Rel",           "Rela",          "Relb","Lck",
  
  #  outneg$s9_stg[grep("Map",outneg$s9_stg)],
  "Ppp2r2d", #from another paper. TNF related
  "Ppp2r1a", #from regev... pops up with xbp1? maybe one of the il4
  "Zfp36", #from regev
  
  #  "Ern1", #xbp1
  
  "Npm1", "Rab7", "Fus","Notch1","Notch3",#ALL
  "Nfkb1",
  
  "Vgf", "Lag3","Dusp4", "Perp",#my own addition
  "Myd88","Tlr4","Tnfaip3","Irak4","Sharpin","Nedd8","Ube2f", "Gpatch8",    #DC paper   . Gpath8 and Myd88 interesting   
  "Rag1",
  "Bach2",
  "Foxp3",    #gata3 high in tregs
  "Irf4",
  "Batf",
  "Ddit3",
  "Ikbke", "Xcl1",#http://www.cell.com/cell-reports/fulltext/S2211-1247%2816%2930696-9?elsca1=etoc&amp;elsca2=email&amp;elsca3=2211-1247_20160712_16_2_&amp;elsca4=Cell%20Press|Webinar
  "Ifng",           "Ifngr1",         "Ifngr2",
  "Stat1","Stat2","Stat4","Stat5a","Stat5b","Stat6",
  "Scara3",
  "Trim25")


tcellgenes <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol %in% tcellgenes]




thsne <- function(){
  takecells <- cellcondition$isgood & cellcondition$ko!="" &cellcondition$ko %in% c("Thy1control","Il4","Lag3")
  
  the_gene_mean <- rowMeans(ncount[,takecells])
  the_gene_mean <- the_gene_mean/sum(the_gene_mean)
  
  takegenes <- the_gene_mean>=0.00001 && !(1:nrow(dat) %in% grep("grna",rownames(dat)))
  x <- runrtsne(rownames(dat)[takegenes],takecells,perplexity = 10)
  plottsne(x,colorbyko())
}


#??colorbrew

colorgenes <- unique(cellcondition$ko)#c("Thy1control","Xbp1","Il2","Il4","Il13","Ern1")
kogenepal <- function(){
  rainbow(length(colorgenes))
}
colorbyko <- function(){
  thecol <- rep("black",nrow(cellcondition))
  
  cols <- kogenepal()
  for(i in 1:length(colorgenes)){
    thecol[cellcondition$ko==colorgenes[i]]<-cols[i]
    
  }
  #    thecol[cellgene=="Xbp1"]<-"red"
  #    thecol[cellgene=="F2rl1"]<-"red"
#  thecol[cellcondition$ko=="Il2"]<-"red"
  thecol  
}

colorbymouse <- function(){
  thecol <- rep("black",nrow(cellcondition))
  thecol[cellcondition$col %in% c(1,2,3)]<-"green"
  thecol[cellcondition$col %in% c(4,5,6)]<-"blue"
  thecol[cellcondition$col %in% c(7,8,9)]<-"red"
  thecol  
}
#rainbow(6)

runrtsne2 <- function(geneids, keepcells=rep(TRUE,ncol(log_ncount)),perplexity=30,max_iter=1000,log_ncount=log_ncount,dims=2){
  #  geneids <- c(ensidPtma,listolacellcycle$ensid)
  #  keepcells <- cellcondition$isss & cellcondition$levelnum %in% c(1,4)
  lcounts <- log_ncount
  selectgenes <- which(rownames(lcounts) %in% geneids)
  lcounts <- lcounts[selectgenes,keepcells]
  #rownames(lcounts) <- togenesym(rownames(lcounts))
  print(dims)
  d <- stats::dist(t(as.matrix(lcounts)))
  set.seed(0) # tsne has some stochastic steps (gradient descent) so need to set random 
  rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=perplexity, verbose = FALSE, max_iter=max_iter,dims = dims)
  
  #function to color by gene exp? ptma? grna? read count?
  
  # https://hms-dbmi.github.io/scw/analysis-of-heterogeneity-and-subpopulations.html
  rtsne_out$keepcells<-keepcells
  rtsne_out 
}
plottsne2<-function(rtsne_out,col=colorbylevel(),title="",labels=NULL){
#  plot(rtsne_out$Y, col=col[rtsne_out$keepcells], pch=16, main='',xlab=title,ylab="")
  print(ncol(rtsne_out$Y))
  if(ncol(rtsne_out$Y)==2){
    plot(rtsne_out$Y, col=col[rtsne_out$keepcells], pch=16, main='',xlab=title,ylab="",cex=2)
    text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=labels,cex=0.5) #,col=col[rtsne_out$keepcells]
  } else {
    plot3d(rtsne_out$Y[,1], rtsne_out$Y[,2],rtsne_out$Y[,3],col=col[rtsne_out$keepcells], radius=0.3,type="s")
  }  
}

thsne2 <- function(perplexity=10,col=colorbyko(),dims=2){
  genes<-unique(cellcondition$ko)
  #  genes<-c("Thy1control","Il4","Lag3","Xbp1")
  takecells <- cellcondition$isgood & cellcondition$ko!="" & cellcondition$ko %in% genes & cellcondition$ko != "f"
  print(dims)
  # 
  # the_gene_mean <- rowMeans(ncount[,takecells])
  # the_gene_mean <- the_gene_mean/sum(the_gene_mean)
  # the_gene_mean>=0.00001 && (1:nrow(dat) %in% grep("grna",rownames(dat)))
  # takegenes <- the_gene_mean>=0.00001 && !(1:nrow(dat) %in% grep("grna",rownames(dat)))
  
  takegenes <- rownames(ncount) %in% tcellgenes
  
  x <- runrtsne2(rownames(dat)[takegenes],takecells,perplexity = perplexity, log_ncount=th_ncount,dims=dims)
  plottsne2(x,col,labels = cellcondition$ko[takecells])
}

## Batch Normalize
th_ncount <- log_ncount

th_gene_mean <- rowMeans(th_ncount[,cellcondition$isgood & cellcondition$ko!=""])
th_gene_var <- rowVar(th_ncount[,cellcondition$isgood & cellcondition$ko!=""])
for(i in 1:3){
#  ctrli <- cellcondition$mouse==i & cellcondition$isgood & cellcondition$ko=="Thy1control"
  ctrli <- cellcondition$mouse==i & cellcondition$isgood 
  which(ctrli)
  th_ctrl_mean <- rowMeans(th_ncount[,ctrli])
  th_corr <- th_ctrl_mean - th_gene_mean
  th_ncount[,cellcondition$mouse==i & cellcondition$isgood] <- th_ncount[,cellcondition$mouse==i & cellcondition$isgood] - th_corr
}
#th_gene_mean <- the_gene_mean/sum(the_gene_mean)
th_gene_mean[tcellgenes]

rowMeans(th_ncount[tcellgenes,cellcondition$mouse==1 & cellcondition$isgood & cellcondition$ko=="Thy1control"])
rowMeans(th_ncount[tcellgenes,cellcondition$mouse==2 & cellcondition$isgood & cellcondition$ko=="Thy1control"])




pdf("plots/tcellSNE.pdf")
par(mfrow=c(1,1))
thsne2(20)
dev.off()


pdf("plots/tcellSNE mouse.pdf")
par(mfrow=c(1,1))
thsne2(20,colorbymouse())
dev.off()


thsne2(20,dims=3)

th_readcount <- rowMeans(dat[tcellgenes,cellcondition$plate==12])
plot(sort(as.double(th_readcount)))

togenesymnames(sort(th_readcount))


mean(colSums(dat[,cellcondition$plate==12]))
4e6*0.75/40e3 #75x more reads

mean(colSums(dat[,cellcondition$plate==12 & cellcondition$col==1]))
mean(colSums(dat[,cellcondition$plate==12 & cellcondition$col==3]))
mean(colSums(dat[,cellcondition$plate==12 & cellcondition$col==4]))
mean(colSums(dat[,cellcondition$plate==12 & cellcondition$col==7]))
mean(colSums(dat[,cellcondition$plate==12 & cellcondition$col==8]))
mean(colSums(dat[,cellcondition$plate==12 & cellcondition$col==9]))


