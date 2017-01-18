
seqs <- c(
  "CTATAAGTATCCCTTGGAGAACCACcttggcgccgcgtgagtcccccacGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttgcaatagcgccgggactagggGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttgctgcgctcagccaatagcgcGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttgttcggaatcgagccaatgagGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttggcgcagcgcgcgccaagccgGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA")


d <- read.table("/ztuff/ebi/sces.git/counts/alignPtma1.stat")
d <- cbind(d,read.table("/ztuff/ebi/sces.git/counts/alignPtma2.stat"))
d <- cbind(d,read.table("/ztuff/ebi/sces.git/counts/alignPtma3.stat"))
d <- cbind(d,read.table("/ztuff/ebi/sces.git/counts/alignPtma4.stat"))
d <- cbind(d,read.table("/ztuff/ebi/sces.git/counts/alignPtma5.stat"))
atcg <- d[1:4,]
startat <- d[-(1:4),]


par(mfrow=c(1,ncol(d)))
for(i in 1:ncol(d))
  pie(atcg[,i],labels=c("A","T","C","G"),main=sprintf("Ptma-%s first base",i))


plotpile2 <- function(n){
  plot(26:48,startat[26:48,n],type="l",main=sprintf("Ptma-%s startpos",n),xlab="",ylab="")  
  for(i in 29:48){
    text(labels = str_sub(seqs[n],i,i),x=i,y=startat[i,n])
  }
}

png("plots/grna startpos.png",width=800)
#par(mfrow=c(1,ncol(d)))
par(mfrow=c(2,3))
for(i in 1:ncol(d))
  plotpile2(i)
dev.off()

#does not seem to align text right. opposite way!
#2 should peak on CAAT... no. on the T. but it is a G, modified. 

#could make a logo. with a matrix based on the starting position