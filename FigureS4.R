###############################################################################
##
## Use the gene sets for Jun/Runx, p53 and MDM2 and compare with the 
## gene expression at the drug treatment data
##
###############################################################################

library(scales)

## Load gene sets
runx.targets <- read.table('runx_geneset.txt')$V1
cjun.targets <- read.table('cjun_geneset.txt')$V1
mdm2.targets <- read.table('mdm2_geneset.txt')$V1
all.mdm2.targets <- read.table('all_mdm2_geneset.txt')$V1
p53.targets <- read.table('p53_geneset.txt')$V1

runx.cjun.common <- intersect(cjun.targets,runx.targets)


## Load lps141 and lps853 drug treatment

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/DESeq_results.RData')
load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/timecourse_count_matrix_design.RData')
res2h.up <- rownames(res2vs0)[which(res2vs0$padj<0.05&res2vs0$log2FoldChange>0)]
res2h.down <- rownames(res2vs0)[which(res2vs0$padj<0.05&res2vs0$log2FoldChange<0)]

res4h.up <- setdiff(rownames(res4vs0)[which(res4vs0$padj<0.05&res4vs0$log2FoldChange>0)],res2h.up)
res4h.down <- setdiff(rownames(res4vs0)[which(res4vs0$padj<0.05&res4vs0$log2FoldChange<0)],res2h.down)

res6h.up <- setdiff(setdiff(rownames(res6vs0)[which(res6vs0$padj<0.05&res6vs0$log2FoldChange>0)],res2h.up), res4h.up)
res6h.down <- setdiff(setdiff(rownames(res6vs0)[which(res6vs0$padj<0.05&res6vs0$log2FoldChange<0)],res2h.down), res4h.down)

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/lps853_diff_genes.RData')
load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/lps853_DESeq_results.RData')

lps853.res4hvs0.up <- setdiff(lps853.res4hvs0.up,lps853.res2hvs0.up)
lps853.res4hvs0.down <- setdiff(lps853.res4hvs0.down,lps853.res2hvs0.down)

lps853.res6hvs0.up <- setdiff(setdiff(lps853.res6hvs0.up,lps853.res2hvs0.up),lps853.res4hvs0.up)
lps853.res6hvs0.down <- setdiff(setdiff(lps853.res6hvs0.down,lps853.res2hvs0.down),lps853.res4hvs0.down)


## Load knockdown data
load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/knowckdown_cohort2/shRNA_1_4_DEGsvsNT.RData')

shRNA4.up <- rownames(results(dds, name="Group_shRNA4_vs_NT"))[which(results(dds, name="Group_shRNA4_vs_NT")$padj<0.05 & results(dds, name="Group_shRNA4_vs_NT")$log2FoldChange>0)]
shRNA4.down <- rownames(results(dds, name="Group_shRNA4_vs_NT"))[which(results(dds, name="Group_shRNA4_vs_NT")$padj<0.05 & results(dds, name="Group_shRNA4_vs_NT")$log2FoldChange<0)]

shRNA1.up <- rownames(results(dds, name="Group_shRNA1_vs_NT"))[which(results(dds, name="Group_shRNA1_vs_NT")$padj<0.05 & results(dds, name="Group_shRNA1_vs_NT")$log2FoldChange>0)]
shRNA1.down <- rownames(results(dds, name="Group_shRNA1_vs_NT"))[which(results(dds, name="Group_shRNA1_vs_NT")$padj<0.05 & results(dds, name="Group_shRNA1_vs_NT")$log2FoldChange<0)]


################################################################################
## Data analysis

################################################################################
## Heatmaps

library(ComplexHeatmap)
library(circlize)

LPS141.FC.matrix <- cbind(res2vs0$log2FoldChange, res4vs0$log2FoldChange, res6vs0$log2FoldChange)
rownames(LPS141.FC.matrix) <- rownames(res2vs0)

pdf(file = 'hdm201_upGenes_1ps141.pdf', width = 4, height = 5)
Heatmap(LPS141.FC.matrix[unique(c(res2h.up,res4h.up,res6h.up)),], 
        row_order = c(res2h.up,res4h.up,res6h.up), cluster_columns = F,
        col = colorRamp2(c(-2, 0, 2), c("blue4", "white", "red4")), row_labels = rep('',length(c(res2h.up,res4h.up,res6h.up))))
dev.off()

high.genes <- c('GMNN','H1-1','H4C1','H4C2','H1-2','H1-4','H2BC5','H4C5','H1-3','H2BC9','H2BC10','H4C8',
                'H4C9','ZNF184','H2BC13','H2AC14','H2AC16','H3C12','H2AC17')
pdf(file = 'hdm201_downGenes_1ps141.pdf', width = 4, height = 5)
Heatmap(LPS141.FC.matrix[unique(c(res2h.down,res4h.down,res6h.down)),], 
        row_order = c(res2h.down,res4h.down,high.genes,setdiff(res6h.down,high.genes)), cluster_columns = F,
        col = colorRamp2(c(-2, 0, 2), c("blue4", "white", "red4")), 
        row_labels = rep('',length(unique(c(res2h.down,res4h.down,high.genes,setdiff(res6h.down,high.genes))))))
dev.off()


LPS853.FC.matrix <- cbind(lps853.res2vs0$log2FoldChange, lps853.res4vs0$log2FoldChange, lps853.res6vs0$log2FoldChange)
rownames(LPS853.FC.matrix) <- rownames(lps853.res2vs0)

pdf(file = 'hdm201_upGenes_1ps853.pdf', width = 4, height = 5)
Heatmap(LPS853.FC.matrix[unique(c(lps853.res2hvs0.up,lps853.res4hvs0.up,lps853.res6hvs0.up)),], 
        row_order = unique(c(lps853.res2hvs0.up,lps853.res4hvs0.up,lps853.res6hvs0.up)), cluster_columns = F, 
        col = colorRamp2(c(-2, 0, 2), c("blue4", "white", "red4")), row_labels = rep('',length(unique(c(lps853.res2hvs0.up,lps853.res4hvs0.up,lps853.res6hvs0.up)))))
dev.off()

pdf(file = 'hdm201_downGenes_1ps853.pdf', width = 4, height = 5)
Heatmap(LPS853.FC.matrix[unique(c(lps853.res2hvs0.down,lps853.res4hvs0.down,lps853.res6hvs0.down)),], 
        row_order = unique(c(lps853.res2hvs0.down,lps853.res4hvs0.down,lps853.res6hvs0.down)), cluster_columns = F, 
        col = colorRamp2(c(-2, 0, 2), c("blue4", "white", "red4")), row_labels = rep('',length(unique(c(lps853.res2hvs0.down,lps853.res4hvs0.down,lps853.res6hvs0.down)))))
dev.off()
################################################################################
## Drug Treatment

sum(p53.targets%in%res2h.up); sum(p53.targets%in%res4h.up); sum(p53.targets%in%res6h.up)

pdf(file = 'hdm201_p53_enrichment_1ps141.pdf', width = 5, height = 8)
barplot(c(length(p53.targets)/14928,sum(p53.targets%in%res2h.up)/length(res2h.up), 
                sum(p53.targets%in%res4h.up)/length(res4h.up),
                sum(p53.targets%in%res6h.up)/length(res6h.up),
        sum(p53.targets%in%res2h.down)/length(res2h.down), 
        sum(p53.targets%in%res4h.down)/length(res4h.down),
        sum(p53.targets%in%res6h.down)/length(res6h.down)),
        col=c('#009052'), main='LPS141 P53 enrichment', ylim = c(0,0.5))
dev.off()

pdf(file = 'hdm201_mdm2_enrichment_1ps141.pdf', width = 5, height = 8)
barplot(c(length(mdm2.targets)/14928,sum(mdm2.targets%in%res2h.up)/length(res2h.up), 
          sum(mdm2.targets%in%res4h.up)/length(res4h.up),
          sum(mdm2.targets%in%res6h.up)/length(res6h.up),
          sum(mdm2.targets%in%res2h.down)/length(res2h.down), 
          sum(mdm2.targets%in%res4h.down)/length(res4h.down),
          sum(mdm2.targets%in%res6h.down)/length(res6h.down)),
        col=c('#F15A29'), main='LPS141 MDM2 enrichment', ylim = c(0,0.5))
dev.off()

pdf(file = 'hdm201_jun_runx_enrichment_1ps141.pdf', width = 5, height = 8)
barplot(c(length(runx.cjun.common)/14928,sum(runx.cjun.common%in%res2h.up)/length(res2h.up), 
          sum(runx.cjun.common%in%res4h.up)/length(res4h.up),
          sum(runx.cjun.common%in%res6h.up)/length(res6h.up),
          sum(runx.cjun.common%in%res2h.down)/length(res2h.down), 
          sum(runx.cjun.common%in%res4h.down)/length(res4h.down),
          sum(runx.cjun.common%in%res6h.down)/length(res6h.down)),
        col=c('#00AEEF'), main='LPS141 Jun/RUNX enrichment', ylim = c(0,0.5))
dev.off()


## Enrichments

#P53
fisher.test(cbind(c(sum(p53.targets%in%res2h.up),length(res2h.up)),c(length(p53.targets),14928)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%res4h.up),length(res4h.up)),c(length(p53.targets),14928)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%res6h.up),length(res6h.up)),c(length(p53.targets),14928)))$p.value

fisher.test(cbind(c(sum(p53.targets%in%res2h.down),length(res2h.down)),c(length(p53.targets),14928)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%res4h.down),length(res4h.down)),c(length(p53.targets),14928)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%res6h.down),length(res6h.down)),c(length(p53.targets),14928)))$p.value

prop.test(cbind(c(sum(p53.targets%in%res2h.up),length(res2h.up)),c(length(p53.targets),14928)))$p.value
prop.test(cbind(c(sum(p53.targets%in%res4h.up),length(res4h.up)),c(length(p53.targets),14928)))$p.value
prop.test(cbind(c(sum(p53.targets%in%res6h.up),length(res6h.up)),c(length(p53.targets),14928)))$p.value

prop.test(cbind(c(sum(p53.targets%in%res2h.down),length(res2h.down)),c(length(p53.targets),14928)))$p.value
prop.test(cbind(c(sum(p53.targets%in%res4h.down),length(res4h.down)),c(length(p53.targets),14928)))$p.value
prop.test(cbind(c(sum(p53.targets%in%res6h.down),length(res6h.down)),c(length(p53.targets),14928)))$p.value

#MDM2
fisher.test(cbind(c(sum(mdm2.targets%in%res2h.up),length(res2h.up)),c(length(mdm2.targets),14928)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%res4h.up),length(res4h.up)),c(length(mdm2.targets),14928)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%res6h.up),length(res6h.up)),c(length(mdm2.targets),14928)))$p.value

fisher.test(cbind(c(sum(mdm2.targets%in%res2h.down),length(res2h.down)),c(length(mdm2.targets),14928)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%res4h.down),length(res4h.down)),c(length(mdm2.targets),14928)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%res6h.down),length(res6h.down)),c(length(mdm2.targets),14928)))$p.value

prop.test(cbind(c(sum(mdm2.targets%in%res2h.up),length(res2h.up)),c(length(mdm2.targets),14928)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%res4h.up),length(res4h.up)),c(length(mdm2.targets),14928)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%res6h.up),length(res6h.up)),c(length(mdm2.targets),14928)))$p.value

prop.test(cbind(c(sum(mdm2.targets%in%res2h.down),length(res2h.down)),c(length(mdm2.targets),14928)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%res4h.down),length(res4h.down)),c(length(mdm2.targets),14928)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%res6h.down),length(res6h.down)),c(length(mdm2.targets),14928)))$p.value

#Jun/RUNX
fisher.test(cbind(c(sum(runx.cjun.common%in%res2h.up),length(res2h.up)),c(length(runx.cjun.common),14928)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%res4h.up),length(res4h.up)),c(length(runx.cjun.common),14928)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%res6h.up),length(res6h.up)),c(length(runx.cjun.common),14928)))$p.value

fisher.test(cbind(c(sum(runx.cjun.common%in%res2h.down),length(res2h.down)),c(length(runx.cjun.common),14928)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%res4h.down),length(res4h.down)),c(length(runx.cjun.common),14928)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%res6h.down),length(res6h.down)),c(length(runx.cjun.common),14928)))$p.value

prop.test(cbind(c(sum(runx.cjun.common%in%res2h.up),length(res2h.up)),c(length(runx.cjun.common),14928)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%res4h.up),length(res4h.up)),c(length(runx.cjun.common),14928)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%res6h.up),length(res6h.up)),c(length(runx.cjun.common),14928)))$p.value

prop.test(cbind(c(sum(runx.cjun.common%in%res2h.down),length(res2h.down)),c(length(runx.cjun.common),14928)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%res4h.down),length(res4h.down)),c(length(runx.cjun.common),14928)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%res6h.down),length(res6h.down)),c(length(runx.cjun.common),14928)))$p.value


'''
LPS141
P53: 
2h up 6.903125e-27; 4h up 3.367244e-13; 6h up 0.1501788
2h down 0.0001178147; 4h down 0.7015377; 6h down 0.2149655

MDM2: 
2h up 0.09024407; 4h up 0.0002260856; 6h up 2.591865e-06
2h down 0.000158838; 4h down 0.002117753; 6h down 0.001116289
 
Jun/Runx: 
2h up 0.9002664; 4h up 2.817619e-13; 6h up 9.510888e-05
2h down 0.06001509; 4h down 3.412751e-07; 6h down 0.003656574
'''


'''
LPS141 - Adj pVal 
P53: 
2h up 1.242562e-25; 4h up 5.387591e-12; 6h up 6.007151e-01
2h down 1.413776e-03; 4h down 1.000000e+00; 6h down 6.448964e-01

MDM2: 
2h up 4.512203e-01; 4h up 2.260856e-03; 6h up 3.628611e-05
2h down 1.747218e-03; 4h down 1.694203e-0; 6h down 1.004661e-02
 
Jun/Runx: 
2h up 1.000000e+00; 4h up 4.789952e-12; 6h up 1.236415e-03
2h down 3.600905e-01; 4h down 5.119127e-06; 6h down 2.559602e-02
'''



##########################################################################
### Add nutlin data to the figure

##/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/nutlin_rna_deseq

sum(rownames(lps141.DE.up.genes)%in%mdm2.genes$V1)/length(rownames(lps141.DE.up.genes))
## DEGs 76; 6 MDM2-bound // 1214   14928 // 0.07894737
sum(rownames(lps141.DE.up.genes)%in%p53.genes$V1)/length(rownames(lps141.DE.up.genes))
## DEGs 76; 19 P53-bound // 478    14928 // 0.25

sum(rownames(lps141.DE.down.genes)%in%mdm2.genes$V1)/length(rownames(lps141.DE.down.genes))
## DEGs 265; 54 MDM2-bound // 1214   14928 // 0.2037736
sum(rownames(lps141.DE.down.genes)%in%p53.genes$V1)/length(rownames(lps141.DE.down.genes))
## DEGs 265; 7 P53-bound // 478    14928 // 0.02641509


pdf(file = 'hdm201_p53_enrichment_1ps141_withNutlin.pdf', width = 6.5, height = 8)
barplot(c(length(p53.targets)/14928,sum(p53.targets%in%res2h.up)/length(res2h.up), 
          sum(p53.targets%in%res4h.up)/length(res4h.up),
          sum(p53.targets%in%res6h.up)/length(res6h.up),
          sum(p53.targets%in%res2h.down)/length(res2h.down), 
          sum(p53.targets%in%res4h.down)/length(res4h.down),
          sum(p53.targets%in%res6h.down)/length(res6h.down),
          19/76, 7/265),
        col=c('#009052'), main='LPS141 P53 enrichment', ylim = c(0,0.5))
dev.off()

prop.test(rbind(c(19,76),c(478, 14112)))$p.value
prop.test(rbind(c(7,265),c(478, 14112)))$p.value

pdf(file = 'hdm201_mdm2_enrichment_1ps141_withNutlin.pdf', width = 6.5, height = 8)
barplot(c(length(mdm2.targets)/14928,sum(mdm2.targets%in%res2h.up)/length(res2h.up), 
          sum(mdm2.targets%in%res4h.up)/length(res4h.up),
          sum(mdm2.targets%in%res6h.up)/length(res6h.up),
          sum(mdm2.targets%in%res2h.down)/length(res2h.down), 
          sum(mdm2.targets%in%res4h.down)/length(res4h.down),
          sum(mdm2.targets%in%res6h.down)/length(res6h.down),
          6/76, 54/265),
        col=c('#F15A29'), main='LPS141 MDM2 enrichment', ylim = c(0,0.5))
dev.off()

prop.test(rbind(c(6,76),c(1214, 14112)))$p.value
prop.test(rbind(c(54,265),c(1214, 14112)))$p.value

pdf(file = 'hdm201_jun_runx_enrichment_1ps141_withNutlin.pdf', width = 6.5, height = 8)
barplot(c(length(runx.cjun.common)/14928,sum(runx.cjun.common%in%res2h.up)/length(res2h.up), 
          sum(runx.cjun.common%in%res4h.up)/length(res4h.up),
          sum(runx.cjun.common%in%res6h.up)/length(res6h.up),
          sum(runx.cjun.common%in%res2h.down)/length(res2h.down), 
          sum(runx.cjun.common%in%res4h.down)/length(res4h.down),
          sum(runx.cjun.common%in%res6h.down)/length(res6h.down),
          19/76, 36/265),
        col=c('#00AEEF'), main='LPS141 Jun/RUNX enrichment', ylim = c(0,0.5))
dev.off()
prop.test(rbind(c(19,76),c(2051, 14112)))$p.value
prop.test(rbind(c(36,265),c(2051, 14112)))$p.value

prop.test(c(19,2051),c(76, 14112))$p.value
prop.test(c(36,2051),c(265,14112))$p.value


fisher.test(rbind(c(19,76),c(2051, 14112)))$p.value
fisher.test(rbind(c(36,265),c(2051, 14112)))$p.value
