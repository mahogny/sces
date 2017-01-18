######################################################################
### QC for cells #####################################################
######################################################################

## mitochondrial content
count_cell <- dat
gene_count <- colSums(count_cell[grep("ENSMUSG",rownames(count_cell)),])
nogene_count <- count_cell[which(rownames(count_cell)=="no_feature"),]
nomapped_count <- count_cell[which(rownames(count_cell)=="not_aligned"),]
exon_prop <- rbind(gene_count, nogene_count, nomapped_count)
exon_prop <- t(t(exon_prop) / colSums(exon_prop))
mt_counts <- colSums(count_cell[mt_genes$ensembl_gene_id,])
mt_prop <- mt_counts / gene_count  



#Plot result per plate
plot(grna_count)
plot(mt_prop)


#pdf("qc_exonprop_mtprop_cell_2i_2.pdf")
par(cex.axis=1.5, cex.lab=1.5)
plot(exon_prop[1,], mt_prop, pch=20, xlab="Proportion of reads mapped to exons",
     ylab="Proportion of reads mapped to mitochondrial genes", xlim=c(0,1), ylim=c(0,1))
abline(h=0.1, v=0.5, col="red", lty=2, lwd=2)
#dev.off()



#remove cells with too few reads (compare to mean for a plate instead?)
png("plots/readcounts.png",width=800)
par(mfrow=c(1,2))
plot(gene_count[cellcondition$isss],main="SS read count")
plot(gene_count[cellcondition$iscr],main="DogSeq read count")
dev.off()

# plot(gene_count[cellcondition$plate %in% c(4,8,9,10)],ylim=c(0,200000))
# plot(gene_count[cellcondition$plate %in% c(4+4+4)],ylim=c(0,200000))
# length(which(gene_count[cellcondition$plate %in% c(4,8,9,10)]>100e3))
# plot(gene_count[cellcondition$plate %in% c(1,2,3,4)])
# plot(gene_count[cellcondition$plate %in% c(5,6,7,8)])
# plot(gene_count[cellcondition$plate %in% c(6)])
# mean(gene_count[cellcondition$isss])
# mean(gene_count[cellcondition$plate %in% c(4)])
# mean(gene_count[cellcondition$plate %in% c(8)])
#cellcondition$isss[which(order(gene_count)<96+10)] <- FALSE

#plot(gene_count[cellcondition$plate %in% c(4,8,9,10,11)],ylim=c(0,1e6))

length(which(gene_count[cellcondition$plate %in% c(4,8,9,10,11)]>20e3))

#remove bad cells
#cellcondition$isgood[cellcondition$isss & gene_count<1e5] <- FALSE

cellcondition$isgood[gene_count<100e3 & cellcondition$isss] <- FALSE #minimum exonic reads for regular chemistry


#mean(gene_count[cellcondition$isss])
#cellcondition$isss[which(order(gene_count)<96+10)] <- FALSE





### TODO detected number of genes?



