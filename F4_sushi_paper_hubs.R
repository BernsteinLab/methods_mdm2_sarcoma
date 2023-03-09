#######################################################
## Create Sushi plots


library('Sushi')

path='/Volumes/seq_epiprod02/'

gene.body <- rbind(c('chr6',27780000,27920000),c('chr2',27300000,27440000),c('chr14',24180000,24320000))
gene.body <- data.frame('chr'=c('chr6','chr2','chr14'), 'start'=c(27780000,27300000,24180000), 'end'=c(27920000,27440000,24320000))

gene.bed <- read.table(paste0(path,'sgaldon/liposarcoma/time-course/gene_clusters/all_genes.bed.txt'), sep = '\t', header = F)
gene.bed <- gene.bed[,c(3,5,6,13,12,4)]
colnames(gene.bed) <- c('chrom','start','stop','gene','score','strand')
gene.bed$score <- '.'
gene.bed <- gene.bed[-c(grep(pattern = 'alt', gene.bed$chrom),grep(pattern = 'fix', gene.bed$chrom),grep(pattern = 'random', gene.bed$chrom)),]

gene.file <- read.table(gzfile('~/Downloads/hg38.refGene.gtf.gz'), sep = '\t', header = F)
gene.file <- gene.file[grep('exon',gene.file$V3),]
gene.name <- unlist(lapply(strsplit(gene.file$V9,';'), function(x){gsub(" ","",gsub(" gene_name","",x[grep('gene_name',x)]))}))
gene.bed <- cbind(gene.file[,c(1,4,5)], gene.name, gene.file[,c(6,7)])
colnames(gene.bed) <- c('chrom','start','stop','gene','score','strand')
gene.bed <- gene.bed[-c(grep(pattern = 'alt', gene.bed$chrom),grep(pattern = 'fix', gene.bed$chrom),grep(pattern = 'random', gene.bed$chrom)),]


#Load bedgraphs
print('Loading MDM2 genome files.....')
mdm2.bedGraph <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/sushi/210223_LPS141_DMSO_MDM2.ucsc.bedGraph'), sep = '\t', header = F)
mdm2.2.bedGraph <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/sushi/LPS141_200104_DMSO_MDM2.ucsc.bedGraph'), sep = '\t', header = F)
print('Loading H3K27ac genome files.....')
k27.bedGraph <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/sushi/LPS141_200910_201204_DMSO_201208_K27ac.ucsc.bedGraph'), sep = '\t', header = F)
print('Loading Input genome files.....')
Input.bedGraph <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/sushi/LPS141_200910_201204_DMSO_201208.ucsc.bedGraph'), sep = '\t', header = F)
print('Loading YY1 genome files.....')
yy1.bedGraph <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/sushi/220317_LPS853_DMSO_cell_signaling_YY1_S23.ucsc.bedgraph'), sep = '\t', header = F)


mdm2.bed <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_mdm2.bed'), sep = '\t', header = F)

#Load loops
LPS141_1.EP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/LPS141_1_E_P_filter9_loops.bedpe'), sep = '\t', header = F)
LPS141_1.EP.bedpe <- cbind(LPS141_1.EP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(LPS141_1.EP.bedpe)[1]),LPS141_1.EP.bedpe[,c(8,7,7)],rep(1,dim(LPS141_1.EP.bedpe)[1]))
colnames(LPS141_1.EP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

LPS141_1.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/LPS141_1_P_P_filter9_loops.bedpe'), sep = '\t', header = F)
LPS141_1.PP.bedpe <- cbind(LPS141_1.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(LPS141_1.PP.bedpe)[1]),LPS141_1.PP.bedpe[,c(8,7,7)],rep(1,dim(LPS141_1.PP.bedpe)[1]))
colnames(LPS141_1.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

LPS853_1.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/LPS853_1_P_P_filter9_loops.bedpe'), sep = '\t', header = F)
LPS853_1.PP.bedpe <- cbind(LPS853_1.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(LPS853_1.PP.bedpe)[1]),LPS853_1.PP.bedpe[,c(8,7,7)],rep(1,dim(LPS853_1.PP.bedpe)[1]))
colnames(LPS853_1.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

DD10.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/loop_annotation/DD10/DD10_P_P_filter9_loops.bedpe'), sep = '\t', header = F)
DD10.PP.bedpe <- cbind(DD10.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(DD10.PP.bedpe)[1]),DD10.PP.bedpe[,c(8,7,7)],rep(1,dim(DD10.PP.bedpe)[1]))
colnames(DD10.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

DD20.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/loop_annotation/DD20/DD20_P_P_filter9_loops.bedpe'), sep = '\t', header = F)
DD20.PP.bedpe <- cbind(DD20.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(DD20.PP.bedpe)[1]),DD20.PP.bedpe[,c(8,7,7)],rep(1,dim(DD20.PP.bedpe)[1]))
colnames(DD20.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

DD31.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/chip_master_datasets/loop_annotation/DD31/DD31_P_P_filter3_loops.bedpe'), sep = '\t', header = F)
DD31.PP.bedpe <- cbind(DD31.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(DD31.PP.bedpe)[1]),DD31.PP.bedpe[,c(8,7,7)],rep(1,dim(DD31.PP.bedpe)[1]))
colnames(DD31.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831496.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831496_P_P_loops_filter5.bedpe'), sep = '\t', header = F)
trim_SRR5831496.PP.bedpe <- cbind(trim_SRR5831496.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831496.PP.bedpe)[1]),trim_SRR5831496.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831496.PP.bedpe)[1]))
colnames(trim_SRR5831496.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831497.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831497_P_P_loops_filter14.bedpe'), sep = '\t', header = F)
trim_SRR5831497.PP.bedpe <- cbind(trim_SRR5831497.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831497.PP.bedpe)[1]),trim_SRR5831497.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831497.PP.bedpe)[1]))
colnames(trim_SRR5831497.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831498.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831498_P_P_loops_filter14.bedpe'), sep = '\t', header = F)
trim_SRR5831498.PP.bedpe <- cbind(trim_SRR5831498.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831498.PP.bedpe)[1]),trim_SRR5831498.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831498.PP.bedpe)[1]))
colnames(trim_SRR5831498.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831499.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831499_P_P_loops_filter14.bedpe'), sep = '\t', header = F)
trim_SRR5831499.PP.bedpe <- cbind(trim_SRR5831499.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831499.PP.bedpe)[1]),trim_SRR5831499.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831499.PP.bedpe)[1]))
colnames(trim_SRR5831499.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831500.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831500_P_P_loops_filter19.bedpe'), sep = '\t', header = F)
trim_SRR5831500.PP.bedpe <- cbind(trim_SRR5831500.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831500.PP.bedpe)[1]),trim_SRR5831500.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831500.PP.bedpe)[1]))
colnames(trim_SRR5831500.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831501.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831501_P_P_loops_filter9.bedpe'), sep = '\t', header = F)
trim_SRR5831501.PP.bedpe <- cbind(trim_SRR5831501.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831501.PP.bedpe)[1]),trim_SRR5831501.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831501.PP.bedpe)[1]))
colnames(trim_SRR5831501.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831502.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831502_P_P_loops_filter9.bedpe'), sep = '\t', header = F)
trim_SRR5831502.PP.bedpe <- cbind(trim_SRR5831502.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831502.PP.bedpe)[1]),trim_SRR5831502.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831502.PP.bedpe)[1]))
colnames(trim_SRR5831502.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831503.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831503_P_P_loops_filter19.bedpe'), sep = '\t', header = F)
trim_SRR5831503.PP.bedpe <- cbind(trim_SRR5831503.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831503.PP.bedpe)[1]),trim_SRR5831503.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831503.PP.bedpe)[1]))
colnames(trim_SRR5831503.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

trim_SRR5831504.PP.bedpe <- read.csv(paste0(path,'sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/loop_annotation/trim_SRR5831504_P_P_loops_filter19.bedpe'), sep = '\t', header = F)
trim_SRR5831504.PP.bedpe <- cbind(trim_SRR5831504.PP.bedpe[,c(1,2,3,4,5,6)],rep(NA,dim(trim_SRR5831504.PP.bedpe)[1]),trim_SRR5831504.PP.bedpe[,c(8,7,7)],rep(1,dim(trim_SRR5831504.PP.bedpe)[1]))
colnames(trim_SRR5831504.PP.bedpe) <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name' ,'score', 'strand1','strand2', 'samplenumber')

pdf(file = paste0('/Volumes/seq_epiprod02/sgaldon/liposarcoma/hichip_datasets/terra_output/hichipper_datasets/mdm2_high_hubs.pdf'),width = 9, height = 12)
for (n in seq(1:dim(gene.body)[1])){
  par(mfrow=c(12,1),mar = c(0.4, 3.5, 0.1, 0.1))
  
  #Gene Track Plot
  pg = plotGenes(geneinfo = gene.bed[which(gene.bed$chrom==gene.body[n,1]),], chrom = gene.body[n,1], chromstart = gene.body[n,2] ,chromend = gene.body[n,3] ,
                 bheight=0.2,plotgenetype="arrow",bentline=FALSE,
                 labeloffset=.4,fontsize=0.9,arrowlength = 0.025,
                 labeltext=TRUE, labelat = 'middle')
  labelgenome(gene.body[n,1], gene.body[n,2],gene.body[n,3],n=3,scale="Mb")
  
  #Plot ChIP tracks and peaks
  ##MDM2
  plotBedgraph(signal = mdm2.bedGraph,gene.body[n,1],gene.body[n,2],gene.body[n,3], transparency=.50, color = '#F15A29')
  axis(side=2,las=2,tcl=.2)
  mtext("MDM2",side=2,line=2.2,cex=0.6,font=1)
  plotBedgraph(signal = mdm2.2.bedGraph,gene.body[n,1],gene.body[n,2],gene.body[n,3], transparency=.50, color = '#F15A29')
  axis(side=2,las=2,tcl=.2)
  mtext("MDM2",side=2,line=2.2,cex=0.6,font=1)
  
  tryCatch(
    {
      plotBed(mdm2.bed, gene.body[n,1],gene.body[n,2],gene.body[n,3], type = "region", color = '#F15A29', height = 0.4)
    },
    error=function(cond) {
      message(paste0('No MDM2 peaks for region ',n))
    }
  )
  
  ##YY1
  plotBedgraph(signal = yy1.bedGraph,gene.body[n,1],gene.body[n,2],gene.body[n,3], transparency=.50, color = '#FFC800')
  axis(side=2,las=2,tcl=.2)
  mtext("YY1",side=2,line=2.2,cex=0.6,font=1)
  
  #H3K27ac
  plotBedgraph(signal = k27.bedGraph,gene.body[n,1],gene.body[n,2],gene.body[n,3], transparency=.50, color = 'gray55')
  axis(side=2,las=2,tcl=.2)
  mtext("H3k27ac",side=2,line=2.2,cex=0.6,font=1)
  
  #Plot loops
  plotBedpe(bedpedata = LPS141_1.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = LPS141_1.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("LPS141",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = LPS853_1.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = LPS853_1.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("LPS853",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = DD10.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = DD10.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("DD10",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = DD20.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = DD20.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("DD20",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = DD31.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = DD31.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("DD31",side=2,line=2.2,cex=.75,font=1)

  # plotBedpe(bedpedata = trim_SRR5831496.PP.bedpe,chrom = gene.body[n,1],
  #           chromstart = gene.body[n,2],chromend = gene.body[n,3],
  #           heights = trim_SRR5831496.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  # labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
  #             chromend = gene.body[n,3],n=4,scale="Mb")
  # axis(side=2,las=2,tcl=.2)
  # mtext("Naive",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831497.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831497.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Naive",side=2,line=2.2,cex=.75,font=1)

  plotBedpe(bedpedata = trim_SRR5831498.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831498.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Naive",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831499.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831499.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Naive",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831500.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831500.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Naive",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831501.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831501.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Th17",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831502.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831502.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Th17",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831503.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831503.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Th17",side=2,line=2.2,cex=.75,font=1)
  
  plotBedpe(bedpedata = trim_SRR5831504.PP.bedpe,chrom = gene.body[n,1],
            chromstart = gene.body[n,2],chromend = gene.body[n,3],
            heights = trim_SRR5831504.PP.bedpe$score,plottype="loops", color = 'black', flip = TRUE)
  labelgenome(chrom = gene.body[n,1], chromstart = gene.body[n,2],
              chromend = gene.body[n,3],n=4,scale="Mb")
  axis(side=2,las=2,tcl=.2)
  mtext("Th17",side=2,line=2.2,cex=.75,font=1)
}

dev.off()




