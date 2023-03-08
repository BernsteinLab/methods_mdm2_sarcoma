#############################################################
## Use promoter and enhancer files to create a loop table
## with gene annotations

sample.name=commandArgs(TRUE)[1]


enhancers <- read.table(paste0(sample.name,'/enhancers_totalLoops.tmp'))
enhancers_paste <- apply(enhancers[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
enhancers_paste <-gsub(" ","",enhancers_paste)
promoters <- read.table(paste0(sample.name,'/promoters_loopAnnot_totalLoops.tmp'))
promoters_paste <- apply(promoters[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
promoters_paste <-gsub(" ","",promoters_paste)
total.loops <- read.table(paste0(sample.name,'/totalLoops.tmp'))
total.loops1_paste <- apply(total.loops[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
total.loops1_paste <-gsub(" ","",total.loops1_paste)
total.loops2_paste <- apply(total.loops[,c('V4','V5','V6')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
total.loops2_paste <-gsub(" ","",total.loops2_paste)

total.loops_paste <- cbind(total.loops1_paste,total.loops2_paste,total.loops$V8)

promoter.annotation <- read.csv(paste0(sample.name,'/promoter_annotation/EP_3_promoter_annotation.txt'), sep = '\t')
promoters_paste_annot <- apply(promoter.annotation[,c('Chr','Start','End')], 1, function(x){paste(x[1], as.numeric(x[2])-1, x[3], sep = '_')})
promoters_paste_annot <-gsub(" ","",promoters_paste_annot)
promoter.annotation <- cbind(promoters_paste_annot,promoter.annotation)

## MDM2-bound
if (file.exists(paste0(sample.name,'/promoter_totalLoops_mdm2Bound.tmp'))){
  promoter.mdm2 <- read.table(paste0(sample.name,'/promoter_totalLoops_mdm2Bound.tmp'))
  promoter.mdm2_paste <- apply(promoter.mdm2[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
  promoter.mdm2_paste <-gsub(" ","",promoter.mdm2_paste)
  
  enhancers.mdm2 <- read.table(paste0(sample.name,'/enhancers_totalLoops_mdm2Bound.tmp'))
  enhancers.mdm2_paste <- apply(enhancers.mdm2[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
  enhancers.mdm2_paste <-gsub(" ","",enhancers.mdm2_paste)
  
  #Promoter
  mdm2.counts.matrix <- read.csv('/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/quant_DD10_mdm2.tmp',sep='\t')
  mdm2.counts.matrix <- cbind(gsub(" ","",apply(mdm2.counts.matrix, 1, function(x){paste(x[1], x[2], x[3], sep = '_')})),mdm2.counts.matrix)
  anchor.regions <- gsub(" ","",apply(promoter.mdm2, 1, function(x){paste(x[1], x[2], x[3], sep = '_')}))
  mdm2.binding.regions <- cbind(gsub(" ","",apply(promoter.mdm2, 1, function(x){paste(x[4], x[5], x[6], sep = '_')})))
  
  count.1 <- c()
  for (region in mdm2.binding.regions){
    mdm2.region <- mdm2.counts.matrix[which(mdm2.counts.matrix[,1]==region),]
    count.1 <- c(count.1,mean(mdm2.region[,5]+1, na.rm = T)/mean(mdm2.region[,6]+1, na.rm = T))
  }
  
  mdm2.counts.lvl <- cbind(mdm2.binding.regions,count.1)
  rownames(mdm2.counts.lvl) <- anchor.regions
  
  #Enhancer
  anchor.regions <- gsub(" ","",apply(enhancers.mdm2, 1, function(x){paste(x[1], x[2], x[3], sep = '_')}))
  mdm2.binding.regions <- cbind(gsub(" ","",apply(enhancers.mdm2, 1, function(x){paste(x[4], x[5], x[6], sep = '_')})))
  
  count.1 <- c()
  for (region in mdm2.binding.regions){
    mdm2.region <- mdm2.counts.matrix[which(mdm2.counts.matrix[,1]==region),]
    count.1 <- c(count.1,mean(mdm2.region[,5]+1, na.rm = T)/mean(mdm2.region[,6]+1, na.rm = T))
  }
  
  mdm2.enhancers.lvl <- cbind(mdm2.binding.regions,count.1)
  rownames(mdm2.enhancers.lvl) <- anchor.regions
  
} else{
  promoter.mdm2_paste <- NULL
  enhancers.mdm2_paste <- NULL
  mdm2.counts.lvl <- NULL
  mdm2.enhancers.lvl <- NULL
}

##################################################
promoters.list <- c(); enhancers.list <- c(); genes.list <- c()
mdm2.promoters.list <- c(); mdm2.enhancers.list <- c()

mdm2.promoters.count1 <- c(); mdm2.enhancers.count1 <- c()

for (loop in c(1:dim(total.loops_paste)[1])){
  #print(total.loops_paste[loop,])
  if (length(promoters_paste[which(promoters_paste%in%total.loops_paste[loop,1])])>0){
    promoter.region <- promoters_paste[which(promoters_paste%in%total.loops_paste[loop,1])]
    promoters.list <- c(promoters.list,promoter.region)
    gene.name <- promoter.annotation[which(promoter.annotation$promoters_paste_annot%in%promoter.region),17]
    genes.list <- c(genes.list, as.character(gene.name))
    
    mdm2.promoters.list <- c(mdm2.promoters.list,promoter.region%in%promoter.mdm2_paste)
    if (promoter.region%in%promoter.mdm2_paste){
      mdm2.promoters.count1 <- c(mdm2.promoters.count1,max(as.numeric(mdm2.counts.lvl[rownames(mdm2.counts.lvl)%in%promoter.region,2],na.rm = T)))
    }else{
      mdm2.promoters.count1 <- c(mdm2.promoters.count1, NA)
    }
  }
  if (length(promoters_paste[which(promoters_paste%in%total.loops_paste[loop,2])])>0){
    promoter.region <- promoters_paste[which(promoters_paste%in%total.loops_paste[loop,2])]
    promoters.list <- c(promoters.list,promoter.region)
    gene.name <- promoter.annotation[which(promoter.annotation$promoters_paste_annot%in%promoter.region),17]
    genes.list <- c(genes.list, as.character(gene.name))
    
    mdm2.promoters.list <- c(mdm2.promoters.list,promoter.region%in%promoter.mdm2_paste)
    if (promoter.region%in%promoter.mdm2_paste){
      mdm2.promoters.count1 <- c(mdm2.promoters.count1,max(as.numeric(mdm2.counts.lvl[rownames(mdm2.counts.lvl)%in%promoter.region,2],na.rm = T)))
    }else{
      mdm2.promoters.count1 <- c(mdm2.promoters.count1, NA)
    }
  }
  if (length(enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,1])])>0){
    enhancer.region <- enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,1])]
    enhancers.list <- c(enhancers.list,enhancer.region)
    
    mdm2.enhancers.list <- c(mdm2.enhancers.list,enhancer.region%in%enhancers.mdm2_paste)
    if (enhancer.region%in%enhancers.mdm2_paste){
      mdm2.enhancers.count1 <- c(mdm2.enhancers.count1,max(as.numeric(mdm2.enhancers.lvl[rownames(mdm2.enhancers.lvl)%in%enhancer.region,2],na.rm = T)))
    }else{
      mdm2.enhancers.count1 <- c(mdm2.enhancers.count1, NA)
    }
  }
  if (length(enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,2])])>0){
    enhancer.region <- enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,2])]
    enhancers.list <- c(enhancers.list,enhancer.region)
    
    mdm2.enhancers.list <- c(mdm2.enhancers.list,enhancer.region%in%enhancers.mdm2_paste)
    if (enhancer.region%in%enhancers.mdm2_paste){
      mdm2.enhancers.count1 <- c(mdm2.enhancers.count1,max(as.numeric(mdm2.enhancers.lvl[rownames(mdm2.enhancers.lvl)%in%enhancer.region,2],na.rm = T)))
    }else{
      mdm2.enhancers.count1 <- c(mdm2.enhancers.count1, NA)
    }
  }
}


table_loop_info <- cbind(paste(total.loops_paste[,1],total.loops_paste[,2],sep='_'),total.loops_paste[,3],genes.list,
                         promoters.list,enhancers.list,mdm2.promoters.list,mdm2.enhancers.list,mdm2.promoters.count1,
                         mdm2.enhancers.count1)
colnames(table_loop_info) <- c('Loop','Score','Gene Name','Promoter','Enhancer','MDM2 Promoter','MDM2 Enhancer','210114_MDM2_Promoter',
                               '210114_MDM2_Enhancer')


write.table(table_loop_info,file = paste0(sample.name,'/EP_loops/EP_3_loop_annotation_file.tsv'), sep='\t',quote = F,row.names = F,col.names = T)



############################################################################
## P-P loops

promoters <- read.table(paste0(sample.name,'/PP_loops/PP_promoters_loopAnnot_totalLoops.tmp'))
promoters_paste <- apply(promoters[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
promoters_paste <-gsub(" ","",promoters_paste)

total.loops <- read.table(paste0(sample.name,'/PP_loops/PP_totalLoops.tmp'))
total.loops1_paste <- apply(total.loops[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
total.loops1_paste <-gsub(" ","",total.loops1_paste)
total.loops2_paste <- apply(total.loops[,c('V4','V5','V6')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
total.loops2_paste <-gsub(" ","",total.loops2_paste)

total.loops_paste <- cbind(total.loops1_paste,total.loops2_paste,total.loops$V8)

promoter.annotation <- read.csv(paste0(sample.name,'/PP_loops/promoter_annotation/PP_3_promoter_annotation.txt'), sep = '\t')
promoters_paste_annot <- apply(promoter.annotation[,c('Chr','Start','End')], 1, function(x){paste(x[1], as.numeric(x[2])-1, x[3], sep = '_')})
promoters_paste_annot <-gsub(" ","",promoters_paste_annot)
promoter.annotation <- cbind(promoters_paste_annot,promoter.annotation)


## MDM2-bound
promoter.mdm2 <- read.table(paste0(sample.name,'/PP_loops/PP_promoter_totalLoops_mdm2Bound.tmp'))
promoter.mdm2_paste <- apply(promoter.mdm2[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
promoter.mdm2_paste <-gsub(" ","",promoter.mdm2_paste)

#Promoter
mdm2.counts.matrix <- read.csv('/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/quant_DD10_mdm2.tmp',sep='\t')
mdm2.counts.matrix <- cbind(gsub(" ","",apply(mdm2.counts.matrix, 1, function(x){paste(x[1], x[2], x[3], sep = '_')})),mdm2.counts.matrix)
anchor.regions <- gsub(" ","",apply(promoter.mdm2, 1, function(x){paste(x[1], x[2], x[3], sep = '_')}))
mdm2.binding.regions <- cbind(gsub(" ","",apply(promoter.mdm2, 1, function(x){paste(x[4], x[5], x[6], sep = '_')})))

count.1 <- c()
for (region in mdm2.binding.regions){
  mdm2.region <- mdm2.counts.matrix[which(mdm2.counts.matrix[,1]==region),]
  count.1 <- c(count.1,mean(mdm2.region[,5]+1, na.rm = T)/mean(mdm2.region[,6]+1, na.rm = T))
}

mdm2.counts.lvl <- cbind(mdm2.binding.regions,count.1)
rownames(mdm2.counts.lvl) <- anchor.regions

##########################################

promoters.list1 <- c();promoters.list2 <- c(); genes.list1 <- c();genes.list2 <- c()
mdm2.promoters.list1 <- c();mdm2.promoters.list2 <- c()

mdm2.promoters1.count1 <- c(); mdm2.promoters2.count1 <- c()

for (loop in c(1:dim(total.loops_paste)[1])){
  if (length(promoters_paste[which(promoters_paste%in%total.loops_paste[loop,1])])>0){
    promoter.region <- promoters_paste[which(promoters_paste%in%total.loops_paste[loop,1])]
    promoters.list1 <- c(promoters.list1,promoter.region)
    gene.name <- promoter.annotation[which(promoter.annotation$promoters_paste_annot%in%promoter.region),17]
    genes.list1 <- c(genes.list1, as.character(gene.name))
    
    mdm2.promoters.list1 <- c(mdm2.promoters.list1,promoter.region%in%promoter.mdm2_paste)
    if (promoter.region%in%promoter.mdm2_paste){
      mdm2.promoters1.count1 <- c(mdm2.promoters1.count1,max(as.numeric(mdm2.counts.lvl[rownames(mdm2.counts.lvl)%in%promoter.region,2],na.rm = T)))
    }else{
      mdm2.promoters1.count1 <- c(mdm2.promoters1.count1, NA)
    }
  }
  if (length(promoters_paste[which(promoters_paste%in%total.loops_paste[loop,2])])>0){
    promoter.region <- promoters_paste[which(promoters_paste%in%total.loops_paste[loop,2])]
    promoters.list2 <- c(promoters.list2,promoter.region)
    gene.name <- promoter.annotation[which(promoter.annotation$promoters_paste_annot%in%promoter.region),17]
    genes.list2 <- c(genes.list2, as.character(gene.name))
    
    mdm2.promoters.list2 <- c(mdm2.promoters.list2,promoter.region%in%promoter.mdm2_paste)
    if (promoter.region%in%promoter.mdm2_paste){
      mdm2.promoters2.count1 <- c(mdm2.promoters2.count1,max(as.numeric(mdm2.counts.lvl[rownames(mdm2.counts.lvl)%in%promoter.region,2],na.rm = T)))
    }else{
      mdm2.promoters2.count1 <- c(mdm2.promoters2.count1, NA)
    }
  }
}



table_loop_info <- cbind(paste(total.loops_paste[,1],total.loops_paste[,2],sep='_'),total.loops_paste[,3],genes.list1, genes.list2,
                         promoters.list1,promoters.list2,mdm2.promoters.list1,mdm2.promoters.list2,mdm2.promoters1.count1,
                         mdm2.promoters2.count1)

colnames(table_loop_info) <- c('Loop','Score','Gene.Name.1', 'Gene.Name.2',
                               'Promoter1','Promoter2','MDM2.Promoter1','MDM2.Promoter2','210114_MDM2_Promoter1','210114_MDM2_Promoter2')


write.table(table_loop_info,file = paste0(sample.name,'/PP_loops/PP_3_loop_annotation_file.tsv'), sep='\t',quote = F,row.names = F,col.names = T)



############################################################################
## E-E loops


enhancers <- read.table(paste0(sample.name,'/EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp'))
enhancers_paste <- apply(enhancers[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
enhancers_paste <-gsub(" ","",enhancers_paste)

total.loops <- read.table(paste0(sample.name,'/EE_loops/EE_totalLoops.tmp'))
total.loops1_paste <- apply(total.loops[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
total.loops1_paste <-gsub(" ","",total.loops1_paste)
total.loops2_paste <- apply(total.loops[,c('V4','V5','V6')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
total.loops2_paste <-gsub(" ","",total.loops2_paste)

total.loops_paste <- cbind(total.loops1_paste,total.loops2_paste,total.loops$V8)

enhancers.annotation <- read.csv(paste0(sample.name,'/EE_loops/enhancers_annotation/EE_3_enhancers_annotation.txt'), sep = '\t')
enhancers_paste_annot <- apply(enhancers.annotation[,c('Chr','Start','End')], 1, function(x){paste(x[1], as.numeric(x[2])-1, x[3], sep = '_')})
enhancers_paste_annot <-gsub(" ","",enhancers_paste_annot)
enhancers.annotation <- cbind(enhancers_paste_annot,enhancers.annotation)


## MDM2-bound
enhancers.mdm2 <- read.table(paste0(sample.name,'/EE_loops/EE_enhancers_totalLoops_mdm2Bound.tmp'))
enhancers.mdm2_paste <- apply(enhancers.mdm2[,c('V1','V2','V3')], 1, function(x){paste(x[1], x[2], x[3], sep = '_')})
enhancers.mdm2_paste <-gsub(" ","",enhancers.mdm2_paste)

#Enhancer
anchor.regions <- gsub(" ","",apply(enhancers.mdm2, 1, function(x){paste(x[1], x[2], x[3], sep = '_')}))
mdm2.binding.regions <- cbind(gsub(" ","",apply(enhancers.mdm2, 1, function(x){paste(x[4], x[5], x[6], sep = '_')})))

count.1 <- c()
for (region in mdm2.binding.regions){
  mdm2.region <- mdm2.counts.matrix[which(mdm2.counts.matrix[,1]==region),]
  count.1 <- c(count.1,mean(mdm2.region[,5]+1, na.rm = T)/mean(mdm2.region[,6]+1, na.rm = T))
}

mdm2.enhancers.lvl <- cbind(mdm2.binding.regions,count.1)
rownames(mdm2.enhancers.lvl) <- anchor.regions


##########################################

enhancers.list1 <- c();enhancers.list2 <- c()
mdm2.enhancers.list1 <- c();mdm2.enhancers.list2 <- c()

mdm2.enhancers1.count1 <- c(); mdm2.enhancers2.count1 <- c()


for (loop in c(1:dim(total.loops_paste)[1])){
  if (length(enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,1])])>0){
    enhancer.region <- enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,1])]
    enhancers.list1 <- c(enhancers.list1,enhancer.region)
    
    mdm2.enhancers.list1 <- c(mdm2.enhancers.list1,enhancer.region%in%enhancers.mdm2_paste)
    if (enhancer.region%in%enhancers.mdm2_paste){
      mdm2.enhancers1.count1 <- c(mdm2.enhancers1.count1,max(as.numeric(mdm2.enhancers.lvl[rownames(mdm2.enhancers.lvl)%in%enhancer.region,2],na.rm = T)))
    }else{
      mdm2.enhancers1.count1 <- c(mdm2.enhancers1.count1, NA)
    }
  }
  if (length(enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,2])])>0){
    enhancer.region <- enhancers_paste[which(enhancers_paste%in%total.loops_paste[loop,2])]
    enhancers.list2 <- c(enhancers.list2,enhancer.region)
    
    mdm2.enhancers.list2 <- c(mdm2.enhancers.list2,enhancer.region%in%enhancers.mdm2_paste)
    if (enhancer.region%in%enhancers.mdm2_paste){
      mdm2.enhancers2.count1 <- c(mdm2.enhancers2.count1,max(as.numeric(mdm2.enhancers.lvl[rownames(mdm2.enhancers.lvl)%in%enhancer.region,2],na.rm = T)))
    }else{
      mdm2.enhancers2.count1 <- c(mdm2.enhancers2.count1, NA)
    }
  }
}




table_loop_info <- cbind(paste(total.loops_paste[,1],total.loops_paste[,2],sep='_'),total.loops_paste[,3],
                         enhancers.list1,enhancers.list2,mdm2.enhancers.list1,mdm2.enhancers.list2,mdm2.enhancers1.count1,
                         mdm2.enhancers2.count1)

colnames(table_loop_info) <- c('Loop','Score','enhancer1','enhancer2','MDM2.enhancer1','MDM2.enhancer2','210114_MDM2_Enhancer1',
                               '210114_MDM2_Enhancer2')


write.table(table_loop_info,file = paste0(sample.name,'/EE_loops/EE_3_loop_annotation_file.tsv'), sep='\t',quote = F,row.names = F,col.names = T)

stop("Finishing the script here")

#######################################################################################
## Cytoscape files
## P-P loops

setwd("/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/")

sample.name='DD10'
a <- read.table(file = paste0(sample.name,'/PP_loops/PP_3_loop_annotation_file.tsv'), sep='\t', header=T)

interactions.file = matrix(ncol=3)
colnames(interactions.file) <- c('loop_1','loop_2','score')
descriptive.file = matrix(ncol=4)
colnames(descriptive.file) <- c('Promoter','Gene','MDM2','MDM2.quant.1')

for (loop in c(1:dim(a)[1])){
  loop.split <- unlist(strsplit(a[loop,1], split = '_'))
  interactions.file <- rbind(interactions.file, c(paste(loop.split[1],loop.split[2],loop.split[3], sep = '_'), 
                                                  paste(loop.split[4],loop.split[5],loop.split[6], sep = '_'), a[loop,2]))
  if (!a[loop,'Promoter1']%in%descriptive.file[,'Promoter']){
    descriptive.file <- rbind(descriptive.file, c(a[loop,'Promoter1'], a[loop,'Gene.Name.1'], a[loop,'MDM2.Promoter1'],a[loop,'X210114_MDM2_Promoter1']))
  }
  if (!a[loop,'Promoter2']%in%descriptive.file[,'Promoter']){
    descriptive.file <- rbind(descriptive.file, c(a[loop,'Promoter2'], a[loop,'Gene.Name.2'], a[loop,'MDM2.Promoter2'],a[loop,'X210114_MDM2_Promoter2']))
  }
}


interactions.file <- interactions.file[2:dim(interactions.file)[1],]
descriptive.file <- descriptive.file[2:dim(descriptive.file)[1],]


write.table(interactions.file,file = paste0(sample.name,'/PP_loops/PP_3_loop_cytoscape.tsv'), quote = F, sep = '\t', row.names = F, col.names = T)
write.table(descriptive.file,file = paste0(sample.name,'/PP_loops/PP_3_promDescription_cytoscape.tsv'), quote = F, sep = '\t', row.names = F, col.names = T)

