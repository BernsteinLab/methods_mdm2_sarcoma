## Expression of looped genes vs non-looped genes

hub.genes <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/geneLists/genesinhubs.txt')
looped.genes <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/geneLists/looped_genes_LPS141_EP_PP_merged.txt')

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/DESeq_results.RData')



boxplot(log2(res2vs0[rownames(res2vs0)%in%looped.genes$V1,'baseMean']), log2(res2vs0[!rownames(res2vs0)%in%looped.genes$V1,'baseMean']))
t.test(log2(res2vs0[rownames(res2vs0)%in%looped.genes$V1,'baseMean']), log2(res2vs0[!rownames(res2vs0)%in%looped.genes$V1,'baseMean']))
