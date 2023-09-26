library(reshape2)

### Load files ##
lfsr=read.table(file='LFSR.tsv',header = T,row.names = 1)
effsize=read.table(file='effect_size.tsv',header = T,row.names = 1)
effsize_se=read.table(file='effect_size_se.tsv',header = T,row.names = 1)

### Calculate CI(95%) of effect size
ciu=effsize+1.96*effsize_se
cil=effsize-1.96*effsize_se

### Set effect size overlapping 0 to NA.
effsize[ciu>0 & cil<0] <- NA

### Set effect size with LFSR > 0.05 to NA.
effsize[lfsr>0.05] <- NA

### Set negative LFSR values to 0.
lfsr[lfsr<0] <- 0 

effsize=cbind(rownames(effsize),effsize)
mat=na.omit(cbind(melt(effsize),melt(effsize_se)[,2],melt(lfsr)[,2]))
colnames(mat)=c('gene','tissue','effsize','effsize_se','lfsr')
mat=mat[with(mat, order(tissue, lfsr)),]

write.table(file='signif.sbgenes.txt',mat,sep='\t',col.names=T,row.names = F,quote=F)
