## Figure S5 - Calculate pphub stats

consensus.prom.pphub <- read.table('consensus_promoters_clusters_lps141_853.bed')

lps141.promoters <- gsub(' ','',apply(consensus.prom.pphub,1,function(x){paste(x[1],x[2],x[3], sep='_')}))
consensus.prom.pphub <- cbind(lps141.promoters,consensus.prom.pphub)

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/DESeq_results.RData')
lps141.baseline.gene.expression <- res2vs0$baseMean
names(lps141.baseline.gene.expression) <- rownames(res2vs0)

a <- pphub.metrics.fun(consensus.prom.pphub)                                                                                                                                                             


###################################
## Distances and gene number histograms

distances <- c(); gene.nums <- c()

for (cluster in unique(consensus.prom.pphub$V4)){
  print(cluster)
  cluster.matrix <- consensus.prom.pphub[consensus.prom.pphub$V4==cluster,]
  distances <- c(distances,max(cluster.matrix$V3) - min(cluster.matrix$V2))
  gene.nums <- c(gene.nums,length(unique(cluster.matrix$V9)))
}

pdf(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/paperfigures/cluster_descrip_histograms.pdf', width = 8, height = 9)
par(mfrow=c(2,1))
hist(distances,breaks=50, xlab = 'Cluster sizes', ylab= 'Frequency', main='Cluster sizes')
hist(gene.nums,breaks=50, xlab = 'Number of genes', ylab= 'Frequency', main='Number of genes per cluster')
dev.off()


## Genes in hubs expression comparison

hub.genes <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/geneLists/genesinhubs.txt')

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/DESeq_results.RData')



expr.in.hubs <- res2vs0[rownames(res2vs0)%in%hub.genes$V1,'baseMean']
expr.not.in.hubs <- res2vs0[!rownames(res2vs0)%in%hub.genes$V1,'baseMean']

pdf(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/paperfigures/genesinhubs_expressionBoxplot.pdf',
    width = 5,height = 6)
boxplot(log2(expr.not.in.hubs),log2(expr.in.hubs), frame.plot=F,
        names=c('Genes not in hubs','Genes in hubs'),ylab='Log2 Expression')
dev.off()
t.test(log2(expr.not.in.hubs),log2(expr.in.hubs))

