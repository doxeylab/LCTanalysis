#!/usr/bin/Rscript

tb <- read.delim('effectors.tb', sep = "\t")

pdf('effectors_auto.pdf')

plot(tb$Cluster_Identity*100, tb$Num_Seqs, main='Number of effector sequence clusters by toxin family', ylab='Number of clusters', xlab='Percent identity', ylim=c(5,800), xlim=c(100, 50), pch=c(21,22,23,24,25)[factor(tb$Type)], cex=1.75, bg=c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f')[factor(tb$Type)])
legend('topright', pch=c(21,22,23,24,25), col=c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f'), legend=unique(tb$Type))

lines(tb[c(1:11), 2]*100, tb[c(1:11), 3])
lines(tb[c(12:22), 2]*100, tb[c(12:22), 3])
lines(tb[c(23:33), 2]*100, tb[c(23:33), 3])
lines(tb[c(34:44), 2]*100, tb[c(34:44), 3])
lines(tb[c(45:55), 2]*100, tb[c(45:55), 3])

dev.off()

