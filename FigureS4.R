library(scales)

## Load gene sets
runx.targets <- read.table('runx_geneset.txt')$V1
cjun.targets <- read.table('cjun_geneset.txt')$V1
mdm2.targets <- read.table('mdm2_geneset.txt')$V1
all.mdm2.targets <- read.table('all_mdm2_geneset.txt')$V1
p53.targets <- read.table('p53_geneset.txt')$V1

runx.cjun.common <- intersect(cjun.targets,runx.targets)


load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/lps853_diff_genes.RData')
load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/lps853_DESeq_results.RData')

lps853.res4hvs0.up <- setdiff(lps853.res4hvs0.up,lps853.res2hvs0.up)
lps853.res4hvs0.down <- setdiff(lps853.res4hvs0.down,lps853.res2hvs0.down)

lps853.res6hvs0.up <- setdiff(setdiff(lps853.res6hvs0.up,lps853.res2hvs0.up),lps853.res4hvs0.up)
lps853.res6hvs0.down <- setdiff(setdiff(lps853.res6hvs0.down,lps853.res2hvs0.down),lps853.res4hvs0.down)



pdf(file = 'hdm201_p53_enrichment_1ps853.pdf', width = 5, height = 8)
barplot(c(length(p53.targets)/15853,sum(p53.targets%in%lps853.res2hvs0.up)/length(lps853.res2hvs0.up), 
          sum(p53.targets%in%lps853.res4hvs0.up)/length(lps853.res4hvs0.up),
          sum(p53.targets%in%lps853.res6hvs0.up)/length(lps853.res6hvs0.up),
          sum(p53.targets%in%lps853.res2hvs0.down)/length(lps853.res2hvs0.down), 
          sum(p53.targets%in%lps853.res4hvs0.down)/length(lps853.res4hvs0.down),
          sum(p53.targets%in%lps853.res6hvs0.down)/length(lps853.res6hvs0.down)),
        col=c('#009052'), main='LPS853 P53 enrichment', ylim = c(0,0.6))
dev.off()

pdf(file = 'hdm201_mdm2_enrichment_1ps853.pdf', width = 5, height = 8)
barplot(c(length(mdm2.targets)/15853,sum(mdm2.targets%in%lps853.res2hvs0.up)/length(lps853.res2hvs0.up), 
          sum(mdm2.targets%in%lps853.res4hvs0.up)/length(lps853.res4hvs0.up),
          sum(mdm2.targets%in%lps853.res6hvs0.up)/length(lps853.res6hvs0.up),
          sum(mdm2.targets%in%lps853.res2hvs0.down)/length(lps853.res2hvs0.down), 
          sum(mdm2.targets%in%lps853.res4hvs0.down)/length(lps853.res4hvs0.down),
          sum(mdm2.targets%in%lps853.res6hvs0.down)/length(lps853.res6hvs0.down)),
        col=c('#F15A29'), main='LPS853 MDM2 enrichment', ylim = c(0,0.6))
dev.off()

pdf(file = 'hdm201_jun_runx_enrichment_1ps853.pdf', width = 5, height = 8)
barplot(c(length(runx.cjun.common)/15853,sum(runx.cjun.common%in%lps853.res2hvs0.up)/length(lps853.res2hvs0.up), 
          sum(runx.cjun.common%in%lps853.res4hvs0.up)/length(lps853.res4hvs0.up),
          sum(runx.cjun.common%in%lps853.res6hvs0.up)/length(lps853.res6hvs0.up),
          sum(runx.cjun.common%in%lps853.res2hvs0.down)/length(lps853.res2hvs0.down), 
          sum(runx.cjun.common%in%lps853.res4hvs0.down)/length(lps853.res4hvs0.down),
          sum(runx.cjun.common%in%lps853.res6hvs0.down)/length(lps853.res6hvs0.down)),
        col=c('#00AEEF'), main='LPS853 Jun/RUNX enrichment', ylim = c(0,0.6))
dev.off()


## Enrichments

#P53
fisher.test(cbind(c(sum(p53.targets%in%lps853.res2hvs0.up),length(lps853.res2hvs0.up)),c(length(p53.targets),15853)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%lps853.res4hvs0.up),length(lps853.res4hvs0.up)),c(length(p53.targets),15853)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%lps853.res6hvs0.up),length(lps853.res6hvs0.up)),c(length(p53.targets),15853)))$p.value

fisher.test(cbind(c(sum(p53.targets%in%lps853.res2hvs0.down),length(lps853.res2hvs0.down)),c(length(p53.targets),15853)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%lps853.res4hvs0.down),length(lps853.res4hvs0.down)),c(length(p53.targets),15853)))$p.value
fisher.test(cbind(c(sum(p53.targets%in%lps853.res6hvs0.down),length(lps853.res6hvs0.down)),c(length(p53.targets),15853)))$p.value

#MDM2
fisher.test(cbind(c(sum(mdm2.targets%in%lps853.res2hvs0.up),length(lps853.res2hvs0.up)),c(length(mdm2.targets),15853)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%lps853.res4hvs0.up),length(lps853.res4hvs0.up)),c(length(mdm2.targets),15853)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%lps853.res6hvs0.up),length(lps853.res6hvs0.up)),c(length(mdm2.targets),15853)))$p.value

fisher.test(cbind(c(sum(mdm2.targets%in%lps853.res2hvs0.down),length(lps853.res2hvs0.down)),c(length(mdm2.targets),15853)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%lps853.res4hvs0.down),length(lps853.res4hvs0.down)),c(length(mdm2.targets),15853)))$p.value
fisher.test(cbind(c(sum(mdm2.targets%in%lps853.res6hvs0.down),length(lps853.res6hvs0.down)),c(length(mdm2.targets),15853)))$p.value

#Jun/RUNX
fisher.test(cbind(c(sum(runx.cjun.common%in%lps853.res2hvs0.up),length(lps853.res2hvs0.up)),c(length(runx.cjun.common),15853)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%lps853.res4hvs0.up),length(lps853.res4hvs0.up)),c(length(runx.cjun.common),15853)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%lps853.res6hvs0.up),length(lps853.res6hvs0.up)),c(length(runx.cjun.common),15853)))$p.value

fisher.test(cbind(c(sum(runx.cjun.common%in%lps853.res2hvs0.down),length(lps853.res2hvs0.down)),c(length(runx.cjun.common),14928)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%lps853.res4hvs0.down),length(lps853.res4hvs0.down)),c(length(runx.cjun.common),14928)))$p.value
fisher.test(cbind(c(sum(runx.cjun.common%in%lps853.res6hvs0.down),length(lps853.res6hvs0.down)),c(length(runx.cjun.common),14928)))$p.value

'''
LPS853
P53: 
2h up 1.129179e-20; 4h up 2.771314e-15; 6h up 1.92008e-07
2h down 0.4152864; 4h down 0.0606269; 6h down 0.0628328

MDM2: 
2h up 1; 4h up 0.02755491; 6h up  0.06411499
2h down 1; 4h down 7.088267e-05; 6h down 0.1601807

Jun/Runx: 
2h up 0.05948488; 4h up 1; 6h up 0.04495789
2h down 0.2923259; 4h down 9.373959e-16; 6h down 0.0143339
'''

'''
LPS853 - adj Pvals
P53: 
2h up 2.032523e-19; 4h up 4.434102e-14; 6h up 2.880120e-06
2h down 1.000000e+00; 4h down 5.948488e-01; 6h down 5.948488e-01

MDM2: 
2h up 1.000000e+00; 4h up 3.306589e-01; 6h up  5.948488e-01
2h down 1.000000e+00; 4h down 9.923574e-04; 6h down 9.610840e-01

Jun/Runx: 
2h up 5.948488e-01; 4h up 1.000000e+00; 6h up 4.945367e-01
2h down 1.000000e+00; 4h down 1.593573e-14; 6h down 1.863407e-01
'''

## Proptests

#P53
prop.test(cbind(c(sum(p53.targets%in%lps853.res2hvs0.up),length(lps853.res2hvs0.up)),c(length(p53.targets),15853)))$p.value
prop.test(cbind(c(sum(p53.targets%in%lps853.res4hvs0.up),length(lps853.res4hvs0.up)),c(length(p53.targets),15853)))$p.value
prop.test(cbind(c(sum(p53.targets%in%lps853.res6hvs0.up),length(lps853.res6hvs0.up)),c(length(p53.targets),15853)))$p.value

prop.test(cbind(c(sum(p53.targets%in%lps853.res2hvs0.down),length(lps853.res2hvs0.down)),c(length(p53.targets),15853)))$p.value
prop.test(cbind(c(sum(p53.targets%in%lps853.res4hvs0.down),length(lps853.res4hvs0.down)),c(length(p53.targets),15853)))$p.value
prop.test(cbind(c(sum(p53.targets%in%lps853.res6hvs0.down),length(lps853.res6hvs0.down)),c(length(p53.targets),15853)))$p.value

#MDM2
prop.test(cbind(c(sum(mdm2.targets%in%lps853.res2hvs0.up),length(lps853.res2hvs0.up)),c(length(mdm2.targets),15853)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%lps853.res4hvs0.up),length(lps853.res4hvs0.up)),c(length(mdm2.targets),15853)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%lps853.res6hvs0.up),length(lps853.res6hvs0.up)),c(length(mdm2.targets),15853)))$p.value

prop.test(cbind(c(sum(mdm2.targets%in%lps853.res2hvs0.down),length(lps853.res2hvs0.down)),c(length(mdm2.targets),15853)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%lps853.res4hvs0.down),length(lps853.res4hvs0.down)),c(length(mdm2.targets),15853)))$p.value
prop.test(cbind(c(sum(mdm2.targets%in%lps853.res6hvs0.down),length(lps853.res6hvs0.down)),c(length(mdm2.targets),15853)))$p.value

#Jun/RUNX
prop.test(cbind(c(sum(runx.cjun.common%in%lps853.res2hvs0.up),length(lps853.res2hvs0.up)),c(length(runx.cjun.common),15853)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%lps853.res4hvs0.up),length(lps853.res4hvs0.up)),c(length(runx.cjun.common),15853)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%lps853.res6hvs0.up),length(lps853.res6hvs0.up)),c(length(runx.cjun.common),15853)))$p.value

prop.test(cbind(c(sum(runx.cjun.common%in%lps853.res2hvs0.down),length(lps853.res2hvs0.down)),c(length(runx.cjun.common),15853)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%lps853.res4hvs0.down),length(lps853.res4hvs0.down)),c(length(runx.cjun.common),15853)))$p.value
prop.test(cbind(c(sum(runx.cjun.common%in%lps853.res6hvs0.down),length(lps853.res6hvs0.down)),c(length(runx.cjun.common),15853)))$p.value

