#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
fname <- args[1]
width <- as.numeric(args[2])

library(wccsom)
dat <- read.table(fname,sep='\t')
l = length(dat)
out <- array(rep(0,l*l), dim=c(l,l))
for (i in 1:l) {
    if ( i %% 100 == 1 ) print(i)
    for (j in 1:l) {
        out[i,j] <- wccsom::wcc(dat[,i],dat[,j],width)
    }
}
fout = paste(fname,'.wcc_width',width,sep='')
write.table(out, fout, quote=FALSE, sep = '\t', row.names=FALSE, col.names=FALSE)
