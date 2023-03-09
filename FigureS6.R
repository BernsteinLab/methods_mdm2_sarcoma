

## Total mdm2 peaks per cell line

pdf('barplots/totalMDM2binding_1.pdf', width = 2.4, height = 4)
barplot(rbind(142,86), beside=T,ylim=c(0,1450), names=c('hct116','u2os'), las=2)
dev.off()

pdf('barplots/totalMDM2binding_2.pdf', width = 2.4, height = 4)
barplot(rbind(678,1414), beside=T,ylim=c(0,1450), names=c('lps853','lps141'), las=2)
dev.off()


##################################################################
## Generate barplots on MDM2 Promoter-Enhancer categories, 
## and overlaps with hubs, p53 and YY1


lps141.mdm2.total <- 1414
lps141.mdm2.promoters <- 711
lps141.mdm2.enhancers <- 606


lps853.mdm2.total <- 678
lps853.mdm2.promoters <- 371
lps853.mdm2.enhancers <- 271


hct116.mdm2.total <- 142
hct116.mdm2.promoters <- 49
hct116.mdm2.enhancers <- 67
  
  
u2os.mdm2.total <- 86
u2os.mdm2.promoters <- 32
u2os.mdm2.enhancers <- 46




bar.mat <- cbind(c(lps141.mdm2.promoters/lps141.mdm2.total,lps141.mdm2.enhancers/lps141.mdm2.total),
      c(lps853.mdm2.promoters/lps853.mdm2.total,lps853.mdm2.enhancers/lps853.mdm2.total),
      c(hct116.mdm2.promoters/hct116.mdm2.total,hct116.mdm2.enhancers/hct116.mdm2.total),
      c(u2os.mdm2.promoters/u2os.mdm2.total,u2os.mdm2.enhancers/u2os.mdm2.total))

colnames(bar.mat) <- c('LPS141','LPS853','HCT116','U2OS')
rownames(bar.mat) <- c('Promoters','Enhancers')

pdf(file = 'promoter_enhancer_mdm2_barplots/stacked_barplot_enhProm.pdf', height = 5, width = 5)
barplot(bar.mat,  names=c('LPS141','LPS853','HCT116','U2OS'))
dev.off()

pdf(file = 'promoter_enhancer_mdm2_barplots/besides_barplot_enhProm.pdf', height = 5, width = 5)
barplot(bar.mat, beside=T,  names=c('LPS141','LPS853','HCT116','U2OS'))
dev.off()
