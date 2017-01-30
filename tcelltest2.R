#289 and 290 empty
which(cellcondition$ko=="Xbp1")-96*8+1

d <- function(...){
  v <- unlist(list(...))-96*3
  print(v %% 12)
  print((v %/% 12)+1)
  invisible(v+96*11)
}


which(cellcondition$ko %in% c("Thy1control","Xbp1","Ern1") & cellcondition$isgood)-96*8

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


factor(cellcondition$mouse)

thcompare2 <- function(indexCheck, indexC=which(cellcondition$ko=="Thy1control" & cellcondition$isgood)){
  datfordeseq  <- dat[grep("ENSMUSG",rownames(dat)),]
  
  keepcells <- 1:ncol(dat) %in% c(indexC,indexCheck)
  cellcond <- data.frame(levelnum=1:ncol(dat) %in% c(indexC),mouse=factor(cellcondition$mouse) )
#  cellcond <- data.frame(levelnum=1:ncol(dat) %in% c(indexC))
  cellcond <- cellcond[keepcells,,drop=FALSE]
  dds <- DESeqDataSetFromMatrix(countData = datfordeseq[,keepcells], colData = cellcond,
                                design = ~ mouse + levelnum) 
  dds <- dds[ rowSums(counts(dds)) > 0.1, ]
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
#only crap


resIl4 <- thcompare2(which(cellcondition$ko=="Il4" & cellcondition$isgood))
thshowres(resIl4)
# 1   Tnfrsf9 2.780979e-05
# 2       Mt1 7.534897e-05
# 3     Il1r2 1.513342e-03
# 4    Ifitm1 2.989957e-03
# 5      Gzme 3.898597e-03
# 6      Gas5 9.266760e-03
# 7      Ecm1 9.266760e-03
# 8     Ccnd2 9.266760e-03
# 9      Glrx 9.266760e-03
# 10      Id2 1.469670e-02
# 11      Mt2 1.880963e-02
# 12    Bnip3 4.766777e-02
# 13     Gzmf 8.043582e-02
# 14    Psma1 8.855692e-02
# 15      Grp 9.638755e-02
# 16  Akr1c18 1.071125e-01
# 17 Tnfrsf18 1.296007e-01
# 18     Thy1 1.296007e-01
# 19     Sat1 1.296007e-01
# 20   Tmsb4x 1.296007e-01
# 21    Rpl17 1.296007e-01
# 22  Fam107b 1.467275e-01
# 23    Il2ra 1.467275e-01


resIl2 <- thcompare2(which(cellcondition$ko=="Il2" & cellcondition$isgood))
thshowres(resIl2)
# 1        Tnfrsf9 1.639168e-09 **
# 2            Mt1 7.221683e-07
# 3        Akr1c18 2.881618e-04
# 4           Ccl3 5.681071e-04 ***
# 5            Mt2 7.663344e-03
# 6         Eef1a1 3.464967e-02
# 7         S100a6 5.384586e-02
# 8            Pth 6.356927e-02
# 9          Il2ra 7.718310e-02 ***
# 10         Rpl14 8.426667e-02
# 11          Sub1 8.426667e-02
# 12         Bnip3 1.185127e-01
# 13          Ccl4 1.185127e-01 ***
# 14          Spp1 1.520285e-01 ***
# 15       Fam107b 1.740263e-01
# 16      Serpinf1 1.740263e-01
# 17        Ifitm1 1.787332e-01
# 18          Cd3d 1.787332e-01
# 19         Ccnd2 1.874422e-01
# 20         Psme2 1.874422e-01
# 21        Tmsb4x 1.955266e-01
# 22        Glipr2 2.395757e-01
# 23        mt-Nd1 2.825438e-01
# 24       Gm13394 2.957604e-01
# 25       Tnfrsf4 2.957604e-01
# 26       Slc25a5 2.957604e-01
# 27          Il13 2.957604e-01 ****
# 28        Tmsb10 3.081991e-01
# 29          Calr 3.198068e-01
# 30           Dbi 3.281118e-01
# 31         Park7 3.281118e-01
# 32         Fxyd5 3.281118e-01
# 33         Serp1 3.319988e-01
# 34           Ak2 3.319988e-01
# 35       Sec61a1 3.319988e-01
# 36          Gzmb 3.610896e-01
# 37       Gm24270 3.658631e-01
# 38 1810037I17Rik 3.695285e-01
# 39        Ranbp1 3.695285e-01
# 40        Akr1a1 3.741557e-01
# 41          Il10 4.250913e-01 ***






resErn1 <- thcompare2(which(cellcondition$ko=="Ern1" & cellcondition$isgood))
thshowres(resErn1)
# 1    Tnfrsf9 0.003533888
# 2       Spp1 0.020893123
# 3    Akr1c18 0.020893123
# 4        Mt1 0.029896605
# 5       Ccl4 0.029896605
# 6       Calr 0.046149673
# 7     Glipr2 0.049004346
# 8     Ptp4a2 0.049518675
# 9        Pth 0.051586069
# 10      Sub1 0.059269604
# 11       Pkm 0.061632140
# 12      Gas5 0.062331036
# 13       Ak2 0.087282525
# 14   Tnfrsf4 0.109819708
# 15      Ifng 0.109819708
# 16   Gm10715 0.118020368
# 17     Serp1 0.145122599
# 18   S100a13 0.145122599
# 19 Serpinb6b 0.145122599
# 20    Cdkn1a 0.145122599
# 21       Dbi 0.146209423
# 22   Fam107b 0.146209423
# 23      Cst7 0.146209423
# 24   Gm11942 0.146209423
# 25     Bnip3 0.146209423
# 26       Mt2 0.146209423
# 27    Tmsb10 0.167297396
# 28     Atp5b 0.179501688
# 29   Gm10288 0.200421650
# 30    Rps27a 0.230835321
# 31     Il2ra 0.235025020
# 32   Arl6ip5 0.235025020
# 33   Gm10717 0.235025020
# 34  Mapkapk2 0.245513993
# 35     Cytip 0.245513993
# 36   Gm10801 0.245513993
# 37     Rpl22 0.245513993
# 38   Lilrb4a 0.245513993
# 39     Pdia6 0.245513993
# 40    Ctla2a 0.245513993
# 41      Gzmf 0.245513993
# 42    Gm8172 0.254452064
# 43 Rpl14-ps1 0.265161093
# 44     Hspe1 0.268700437
# 45      Dad1 0.268700437
# 46      Rac2 0.268700437
# 47    Cox6a1 0.287526738
# 48     Psma1 0.304830708
# 49       Id2 0.304830708
# 50   Gm26870 0.308671208



resIl13 <- thcompare2(which(cellcondition$ko=="Il13" & cellcondition$isgood))
thshowres(resIl13)
#crap
#summary(gene_count[which(cellcondition$ko=="Il13" & cellcondition$isgood)])
#gene_count[which(cellcondition$ko=="Il13" & cellcondition$isgood)]


resEtv2 <- thcompare2(which(cellcondition$ko=="Etv2" & cellcondition$isgood))
thshowres(resEtv2)
# 1    Tnfrsf4 0.0001727499
# 2    Akr1c18 0.0001727499
# 3      Ccnd2 0.0002996489
# 4       Gzma 0.0004633158
# 5    Gm13456 0.0045901485
# 6       Ccl3 0.0213945656
# 7    Tnfrsf9 0.0214297938
# 8     Ifitm1 0.0214297938
# 9       Prf1 0.0214297938
# 10    Tmsb10 0.0342955549
# 11       Pkm 0.0376420099
# 12      Il13 0.0464462778
# 13      Gzme 0.0464462778
# 14      Spp1 0.0577689011
# 15       Mt1 0.0577689011
# 16      Igkc 0.0634564120
# 17    Lgals1 0.0634564120
# 18     Ctla4 0.0812496389
# 19     Il2ra 0.0812496389
# 20   Gm10720 0.0812496389
# 21       Son 0.0812496389
# 22    Glipr2 0.0902940891
# 23   Gm12715 0.0902940891
# 24    Slc2a3 0.0902940891
# 25      Gzmd 0.0902940891
# 26       Vim 0.0961470547
# 27   Gm10801 0.0961470547
# 28    Serbp1 0.0961470547
# 29   Gm10715 0.0961470547
# 30      Ccl4 0.0961470547


# 365 AGAGAGAAGTCCCCGCGC    fits lag3 as it should!

resLag3 <- thcompare2(which(cellcondition$ko=="Lag3" & cellcondition$isgood))
thshowres(resLag3)
# 1        Gm10720 0.01793612
# 2         Lmbrd1 0.04052486
# 3           Xrn2 0.04052486
# 4        Gm10800 0.04052486
# 5         Slc2a3 0.04052486
# 6          Aldoa 0.04052486
# 7        Gm26870 0.04052486
# 8        Gm11168 0.04052486
# 9        Gm17535 0.04052486
# 10       Gm10715 0.04052486
# 11          Eef2 0.04052486
# 12          Mdp1 0.04052486
# 13          Ly6a 0.04052486


gene_count[which(cellcondition$ko=="Lag3" & cellcondition$isgood)]

resF2 <- thcompare2(which(cellcondition$ko=="F2rl1" & cellcondition$isgood))
thshowres(resF2)
# 1         Ctla2a 0.01621271
# 2           Ccl3 0.07058116
# 3      Hnrnpa2b1 0.09969749
# 4         Tmsb4x 0.09969749
# 5           Gzma 0.09969749
# 6          Ccnd2 0.12941277
# 7          Rps18 0.30676669
# 8          Rps19 0.60920445
# 9          Il2ra 0.92324472




#indexEtv2 <- which(cellgene=="Zc3h12a") #??

resC <- thcompare2(which(cellcondition$ko=="Cxcr7" & cellcondition$isgood))
thshowres(resC)
#crap




resZ <- thcompare2(which(cellcondition$ko=="Zc3h12a" & cellcondition$isgood))
thshowres(resZ)
#crap

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


