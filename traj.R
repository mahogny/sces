#trajectory analysis

#TODO: maybe only 

#G1 has lowest Ptma. So it goes to G0 from G1

#G1 -> S -> G2 -> M ->
#black -> blue -> red > 
#actually ptma is low between G1 and S

colorbycc <- function(){
  thecol <- rep("green",nrow(cellcondition))
  thecol[cycass$phases=="S"] <- "#0000FF"
  thecol[cycass$phases=="G1"] <- "#000000"
  thecol[cycass$phases=="G2M"] <- "#FF0000"
  thecol
}

### PCA for all expressed genes
dopca <- function(){
  # gene_mean_cr <- rowMeans(ncount[,cellcondition$iscr]) 
  # gene_mean_ss <- rowMeans(ncount[,cellcondition$isss])
  # gene_mean_cr <- gene_mean_cr/sum(gene_mean_cr)
  # gene_mean_ss <- gene_mean_ss/sum(gene_mean_ss)
  
  takegene <- rownames(dat) %in% listkedarcellcycle$ensid
  
  #length(which(gene_mean_ss>=0.00001))
  #cor_ncount <- cor(ncount[gene_mean_ss>=0.00001,], method="spearman")
  #  takecells <- cellcondition$plate %in% c(1) & cellcondition$isss & cellcondition$isgood #& !cellcondition$iswt 
  #  takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood #& !cellcondition$iswt 
  #takecells <- cellcondition$plate %in% c(1,2,3) & cellcondition$isss & cellcondition$isgood &cellcondition$isfeeder  #& !cellcondition$iswt 
  takecells <-  cellcondition$isss & cellcondition$isgood 
  
  
  #and mcherry
         
  #olas way
  cor_ncount <- cor(ncount[takegene ,takecells], method="spearman")  #olas way
  pca_cor <- prcomp(cor_ncount, scale=FALSE)
  pca <- pca_cor
  
  #normal way
  pca_log <- prcomp(t(log_ncount[takegene, takecells]), scale=FALSE)
  pca_nolog <- prcomp(t(ncount[takegene, takecells]), scale=FALSE)
  
  pca <- pca_nolog
  pca <- pca_log
#  pca <- cor_ncount <-   #normal way
  
  
  pdf("plots/traj1.pdf")
  par(cex.axis=1.5, cex.lab=1.5)
  col<-colorbygene()
  col<-colorbycc()
  col[ncount[ensidPtma,]<6 & ncount["mcherry",]>2]<-"green"
  summary(factor(col[takecells]))
  #length(pca$x[,1])
  
  plot3d(pca$x[,1], pca$x[,2], pca$x[,3], col=col[takecells], type="s", radius=0.4, size=10)
  plot(pca$x[,1], pca$x[,2], pch=20, col=col[takecells], xlab="PC", ylab="PC")
  dev.off()

}


#can use WT samples to define the normal trajectory. for each point then, find closest mapping (ref sample), subtract, what is left?
### can find closest 3 points and take average of these! which genes are most variable in the residual?
## gaussian kernel over distance? plot all distances

#easily up to 5 points

#todo tune
pointwt <- c(ncount[ensidPtma,]>6 & ncount["mcherry",]<5 & cellcondition$isss & cellcondition$isss)
pointarrest <- c(ncount[ensidPtma,]<3 & ncount["mcherry",]>3 & cellcondition$isss & cellcondition$isss)
sum(pointwt)
sum(pointarrest)

distgenes <- rownames(dat) %in% listkedarcellcycle$ensid
distcount <- log_ncount[distgenes,pointwt]
diff_log <- matrix(0,nrow=nrow(log_ncount),ncol=sum(pointarrest))
rownames(diff_log)<-rownames(log_ncount)
rpi <- 1
for(pi in which(pointarrest)){
  #careful not to mix log and non-log here
  distances <- apply(distcount,2,function(x) sum((x-log_ncount[distgenes,pi])^2) )
  #distcount[1:5,] - t(log_ncount[distgenes[1:5],pi])
  print(rpi)
  
  #nd <- exp(-distances*distances/mean(distances*distances/30)) #100 for 5 points 
  #nd <- nd/sum(nd)
  #round(tail(sort(nd,decreasing = FALSE),n=10),digits = 2)
  #localav <- rowSums(log_ncount[1:4,]*nd)
  localav <- rowMeans(log_ncount[,order(distances)[1:5]])
  
  #localav could be improved by eliminating 0s
  
  diff_log[,rpi] <- log_ncount[,pi]-localav
  rpi <- rpi+1
  
  #mean(distances)
  
  # pdf("plots/test.pdf")
  # hist(nd)
  # dev.off()
}

#can also sort these outliers by Ptma level!


#of these genes, which are % most off 
pointwt_mean <- rowMeans(log_ncount[,pointwt])
diff_mean <- rowMeans(abs(diff_log))
sum(pointwt_mean>1)
difflist <- sort((diff_mean/pointwt_mean)[pointwt_mean>0.001]) #0.001 for grna. 1 for cas9. 0.001 anxa1 high, apoptosis!
tail(difflist,n=40)
togenesymnames((tail(difflist,n=40)))

difflist[toensid("Gprasp2")] #12% difference

# Cryz   Tmem134  Slc25a20    Ctdsp2   Gm15800      Asb8    Armcx1     Decr2    Shkbp1     N4bp3     Usp47     Cyr61      Edc3 
# 1.493463  1.507454  1.519667  1.522504  1.564482  1.618232  1.639756  1.655886  1.661381  1.749266  1.790783  1.871028  2.028713 

#screen for genes involved in murine hema stem cell repopulation
#eg.  Arhgef5, Armcx1, Cadps2, Crispld1, Emcn, Foxa3, Fstl1, Glis2, Gprasp2, Gpr56, Myct1, Nbea, P2ry14, Smarca2, Sox4, Stat4, and Zfp251
#https://www.ncbi.nlm.nih.gov/pubmed/26880577

#one GO gruop for apoptosis
#http://www.informatics.jax.org/go/term/GO:1900117

#not well studied. 63% diff, >1
#KO of this gene causes improved repopulation of stem cells    same for Gprasp2
#https://www.ncbi.nlm.nih.gov/pubmed/?term=Armcx1

#getting nothing out
stopgosym(togenesym(names(tail(difflist,n=100))))

#Anxa1 0.001, definitely cell death
#Glp2r also apoptosis? https://en.wikipedia.org/wiki/Glucagon-like_peptide_2_receptor
#Edc3 1 mrna capping removal
#1 https://en.wikipedia.org/wiki/CYR61   many processes, relevant!

#in particular, use cell cycle genes for the proximity calculation? or all but mt?

#PCA on cell cycle genes. maybe only WT, then try to map in the other cells?



#WT: high ptma and mcherry


#DEseq on smaller quantiles

#markers for apoptosis. GO!

