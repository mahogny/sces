#this approach does not seem to work well. mainly precision seem crap and there is no way to look at what happen to the sequence

bed <- read.table("counts/grna.bed",stringsAsFactors=FALSE,sep="\t")

#max(bed[,3])
#min(bed[,2])

bedcount <- matrix(0,5,120)
r<-1
dat2 <- bed[bed[,1] %in% sprintf("grna-Ptma-F-%s",1:5),]
for(i in 1:nrow(dat2))
  bedcount[r,dat2[i,2]:dat2[i,3]] <- bedcount[r,dat2[i,2]:dat2[i,3]] + 1
bgbed <- bedcount[1,]

bedcount <- matrix(0,5,120)
for(r in 1:5){
  dat2 <- bed[bed[,1]==sprintf("grna-Ptma-F-%s",r),]
  for(i in 1:nrow(dat2)){
    bedcount[r,dat2[i,2]] <- bedcount[r,dat2[i,2]] + 1
    #    bedcount[r,dat2[i,2]:dat2[i,3]] <- bedcount[r,dat2[i,2]:dat2[i,3]] + 1
  }
}
for(i in 1:5){
  bedcount[i,] <- bedcount[i,]/max(bedcount[i,])
}
# bgbed <- bgbed/sum(bgbed)
#lines(t(bedcount[1:5,]),type="l")  


plot(bedcount[5,],type="l")  
plot(bedcount[2,30:50],type="l")  
plot(bgbed[26:50],type="l")  
plot(bgbed,type="l")  

seqs <- c(
  "CTATAAGTATCCCTTGGAGAACCACcttggcgccgcgtgagtcccccacGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttgcaatagcgccgggactagggGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttgctgcgctcagccaatagcgcGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttgttcggaatcgagccaatgagGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",
  "CTATAAGTATCCCTTGGAGAACCACcttggcgcagcgcgcgccaagccgGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA")

plotpile <- function(n){
  #base<-26
  plot(26:48,bedcount[n,26:48],type="l")  
  for(i in 29:48){
    text(labels = str_sub(seqs[n],i+1,i+1),x=i,y=bedcount[n,i])
  }
}
plotpile(2)

#pileups are of poor precision. should do a multiple alignment. and then from there pull out starting position

#note for #2, a G is forced onto it in the beginning
#grep -h --color GGGGTTTTAGAGCTAGAAATAGCAA S5*clip2*  | sort

#grep --color GAGGTTTTAGAGCTAGAAATAGCAA S5*clip2*

#should subtract the other grnas as background!

sqldf("select V1,count(*) from bed group by V1")
apply

#are many short mapping to 4, but not to 1?

#26 to 50?
#CTATAAGTATCCCTTGGAGAACCAC%sGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA",dat$seq[i])  

#>grna-Ptma-F-1
#CTATAAGTATCCCTTGGAGAACCACcttggcgccgcgtgagtcccccacGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGA


#might be best to remap but with only -Ptma genes, to handle incomplete ones better

# > grnatotc
# grna-ptma1 grna-ptma2 grna-ptma3 grna-ptma4 grna-ptma5 
# 0.01031705 0.09118875 0.21959086 0.61636088 0.06254246 