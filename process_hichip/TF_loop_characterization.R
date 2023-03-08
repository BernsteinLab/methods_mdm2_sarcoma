###########################################################
##
## Use the tsv files from the hiChIP datasets to extract
## the regions associated to each of the TFs


lps141.1.EE <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/EE_loops/EE_3_loop_annotation_file.tsv', 
           header=TRUE, sep = '\t')

lps141.1.EP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/EP_loops/EP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

lps141.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')


lps141.p53.enhancers<- unique(c(lps141.1.EE[which(lps141.1.EE$p53.enhancer1),'enhancer1'], 
                                lps141.1.EE[which(lps141.1.EE$p53.enhancer2),'enhancer2'],
                                lps141.1.EP[which(lps141.1.EP$p53.Enhancer),'Enhancer']))

lps141.p53.promoters<- unique(c(lps141.1.PP[which(lps141.1.PP$p53.Promoter1),'promoter1'], 
                                lps141.1.PP[which(lps141.1.PP$p53.Promoter2),'promoter2'],
                                lps141.1.EP[which(lps141.1.EP$p53.Promoter),'Promoter']))

lps141.mdm2.enhancers<- unique(c(lps141.1.EE[which(lps141.1.EE$MDM2.enhancer1),'enhancer1'],
                                lps141.1.EE[which(lps141.1.EE$MDM2.enhancer2),'enhancer2'],
                                lps141.1.EP[which(lps141.1.EP$MDM2.Enhancer),'Enhancer']))

lps141.mdm2.promoters<- unique(c(lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1),'promoter1'], 
                                lps141.1.PP[which(lps141.1.PP$MDM2.Promoter2),'promoter2'],
                                lps141.1.EP[which(lps141.1.EP$MDM2.Promoter),'Promoter']))

lps141.cjun.enhancers<- unique(c(lps141.1.EE[which(lps141.1.EE$cJun.enhancer1),'enhancer1'],
                                 lps141.1.EE[which(lps141.1.EE$cJun.enhancer2),'enhancer2'],
                                 lps141.1.EP[which(lps141.1.EP$cJun.Enhancer),'Enhancer']))

lps141.cjun.promoters<- unique(c(lps141.1.PP[which(lps141.1.PP$cJun.Promoter1),'promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$cJun.Promoter2),'promoter2'],
                                 lps141.1.EP[which(lps141.1.EP$cJun.Promoter),'Promoter']))

lps141.runx.enhancers<- unique(c(lps141.1.EE[which(lps141.1.EE$runx.enhancer1),'enhancer1'],
                                 lps141.1.EE[which(lps141.1.EE$runx.enhancer2),'enhancer2'],
                                 lps141.1.EP[which(lps141.1.EP$RUNX.Enhancer),'Enhancer']))

lps141.runx.promoters<- unique(c(lps141.1.PP[which(lps141.1.PP$runx.Promoter1),'promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$runx.Promoter2),'promoter2'],
                                 lps141.1.EP[which(lps141.1.EP$RUNX.Promoter),'Promoter']))

dev.off()
VennDiagram::draw.pairwise.venn(length(lps141.p53.enhancers), length(lps141.mdm2.enhancers), 
                                length(intersect(lps141.p53.enhancers,lps141.mdm2.enhancers)))
dev.off()
VennDiagram::draw.pairwise.venn(length(lps141.p53.promoters), length(lps141.mdm2.promoters), 
                                length(intersect(lps141.p53.promoters,lps141.mdm2.promoters)))
dev.off()
VennDiagram::draw.pairwise.venn(length(lps141.cjun.promoters), length(lps141.mdm2.promoters), 
                                length(intersect(lps141.cjun.promoters,lps141.mdm2.promoters)))
dev.off()
VennDiagram::draw.pairwise.venn(length(lps141.cjun.enhancers), length(lps141.mdm2.enhancers), 
                                length(intersect(lps141.cjun.enhancers,lps141.mdm2.enhancers)))
dev.off()
VennDiagram::draw.pairwise.venn(length(lps141.runx.promoters), length(lps141.mdm2.promoters), 
                                length(intersect(lps141.runx.promoters,lps141.mdm2.promoters)))
dev.off()
VennDiagram::draw.pairwise.venn(length(lps141.runx.enhancers), length(lps141.mdm2.enhancers), 
                                length(intersect(lps141.runx.enhancers,lps141.mdm2.enhancers)))

##################################################################################################
### MDM2 at both promoter sides


promoter.both <- unique(c(lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==T),'Promoter1'],
         lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==T),'Promoter2']))


promoter.onlyOne <- unique(c(lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==F),'Promoter1'],
         lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==F&lps141.1.PP$MDM2.Promoter2==T),'Promoter2']))


promoter.toenh <- unique(c(lps141.1.EP[which(lps141.1.EP$MDM2.Promoter==T&lps141.1.EP$MDM2.Enhancer==F),'Promoter']))

mdm2.total.promoters <- unique(c(lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==T),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==T),'Promoter2'],
                               lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==F),'Promoter1'],
                               lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==F&lps141.1.PP$MDM2.Promoter2==T),'Promoter2'],
                               lps141.1.EP[which(lps141.1.EP$MDM2.Promoter==T&lps141.1.EP$MDM2.Enhancer==F),'Promoter']))
p53.total.promoters <- unique(c(lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==T),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==T),'Promoter2'],
                                 lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==F),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$p53.Promoter1==F&lps141.1.PP$p53.Promoter2==T),'Promoter2'],
                                 lps141.1.EP[which(lps141.1.EP$p53.Promoter==T&lps141.1.EP$p53.Enhancer==F),'Promoter']))
cjun.total.promoters <- unique(c(lps141.1.PP[which(lps141.1.PP$cJun.Promoter1==T&lps141.1.PP$cJun.Promoter2==T),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$cJun.Promoter1==T&lps141.1.PP$cJun.Promoter2==T),'Promoter2'],
                                 lps141.1.PP[which(lps141.1.PP$cJun.Promoter1==T&lps141.1.PP$cJun.Promoter2==F),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$cJun.Promoter1==F&lps141.1.PP$cJun.Promoter2==T),'Promoter2'],
                                 lps141.1.EP[which(lps141.1.EP$cJun.Promoter==T&lps141.1.EP$cJun.Enhancer==F),'Promoter']))
runx.total.promoters <- unique(c(lps141.1.PP[which(lps141.1.PP$runx.Promoter1==T&lps141.1.PP$runx.Promoter2==T),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$runx.Promoter1==T&lps141.1.PP$runx.Promoter2==T),'Promoter2'],
                                 lps141.1.PP[which(lps141.1.PP$runx.Promoter1==T&lps141.1.PP$runx.Promoter2==F),'Promoter1'],
                                 lps141.1.PP[which(lps141.1.PP$runx.Promoter1==F&lps141.1.PP$runx.Promoter2==T),'Promoter2'],
                                 lps141.1.EP[which(lps141.1.EP$runx.Promoter==T&lps141.1.EP$runx.Enhancer==F),'Promoter']))

setdiff(unique(c(promoter.onlyOne,promoter.toenh)),promoter.both)

barplot(c(length(promoter.both), length(setdiff(promoter.onlyOne,promoter.both)), length(setdiff(promoter.toenh,c(promoter.onlyOne,promoter.both)))), main = 'mdm2',
        col=c(alpha('red3', 0.2), alpha('yellow3', 0.2),alpha('purple2', 0.2)))


mdm2.PP.both <- list()
a <- lps141.1.PP[which(lps141.1.PP$MDM2.Promoter1==T&lps141.1.PP$MDM2.Promoter2==T),]
for (loop in c(1:dim(a)[1])){
        if (a[loop,'Promoter1']%in%names(mdm2.PP.both)){
                mdm2.PP.both[a[loop,'Promoter1']] <- unlist(mdm2.PP.both[a[loop,'Promoter1']])+1
        } else{
                mdm2.PP.both[a[loop,'Promoter1']] = 1
        }
        if (a[loop,'Promoter2']%in%names(mdm2.PP.both)){
                mdm2.PP.both[a[loop,'Promoter2']] <- unlist(mdm2.PP.both[a[loop,'Promoter2']])+1
        } else{
                mdm2.PP.both[a[loop,'Promoter2']] = 1
        }
}

p53.PP.both <- list()
a <- lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==T),]
for (loop in c(1:dim(a)[1])){
        if (a[loop,'Promoter1']%in%names(p53.PP.both)){
                p53.PP.both[a[loop,'Promoter1']] <- unlist(p53.PP.both[a[loop,'Promoter1']])+1
        } else{
                p53.PP.both[a[loop,'Promoter1']] = 1
        }
        if (a[loop,'Promoter2']%in%names(p53.PP.both)){
                p53.PP.both[a[loop,'Promoter2']] <- unlist(p53.PP.both[a[loop,'Promoter2']])+1
        } else{
                p53.PP.both[a[loop,'Promoter2']] = 1
        }
}

cjun.PP.both <- list()
a <- lps141.1.PP[which(lps141.1.PP$cJun.Promoter1==T&lps141.1.PP$cJun.Promoter2==T),]
for (loop in c(1:dim(a)[1])){
        if (a[loop,'Promoter1']%in%names(cjun.PP.both)){
                cjun.PP.both[a[loop,'Promoter1']] <- unlist(cjun.PP.both[a[loop,'Promoter1']])+1
        } else{
                cjun.PP.both[a[loop,'Promoter1']] = 1
        }
        if (a[loop,'Promoter2']%in%names(cjun.PP.both)){
                cjun.PP.both[a[loop,'Promoter2']] <- unlist(cjun.PP.both[a[loop,'Promoter2']])+1
        } else{
                cjun.PP.both[a[loop,'Promoter2']] = 1
        }
}

runx.PP.both <- list()
a <- lps141.1.PP[which(lps141.1.PP$runx.Promoter1==T&lps141.1.PP$runx.Promoter2==T),]
for (loop in c(1:dim(a)[1])){
        if (a[loop,'Promoter1']%in%names(runx.PP.both)){
                runx.PP.both[a[loop,'Promoter1']] <- unlist(runx.PP.both[a[loop,'Promoter1']])+1
        } else{
                runx.PP.both[a[loop,'Promoter1']] = 1
        }
        if (a[loop,'Promoter2']%in%names(runx.PP.both)){
                runx.PP.both[a[loop,'Promoter2']] <- unlist(runx.PP.both[a[loop,'Promoter2']])+1
        } else{
                runx.PP.both[a[loop,'Promoter2']] = 1
        }
}

promoter.both <- unique(c(lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==T),'Promoter1'],
                          lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==T),'Promoter2']))


promoter.onlyOne <- unique(c(lps141.1.PP[which(lps141.1.PP$p53.Promoter1==T&lps141.1.PP$p53.Promoter2==F),'Promoter1'],
                             lps141.1.PP[which(lps141.1.PP$p53.Promoter1==F&lps141.1.PP$p53.Promoter2==T),'Promoter2']))


promoter.toenh <- unique(c(lps141.1.EP[which(lps141.1.EP$p53.Promoter==T&lps141.1.EP$p53.Enhancer==F),'Promoter']))


setdiff(unique(c(promoter.onlyOne,promoter.toenh)),promoter.both)

barplot(c(length(promoter.both), length(setdiff(promoter.onlyOne,promoter.both)), length(setdiff(promoter.toenh,c(promoter.onlyOne,promoter.both)))), main = 'p53',
        col=c(alpha('red3', 0.2), alpha('yellow3', 0.2),alpha('purple2', 0.2)))

n=1 #Number of anchors it has to be connected to

barplot(c(sum(mdm2.PP.both>=n)/length(mdm2.total.promoters), sum(p53.PP.both>=n)/length(p53.total.promoters),
          sum(cjun.PP.both>=n)/length(cjun.total.promoters),sum(runx.PP.both>=n)/length(runx.total.promoters)),
        col=c('purple3'), names.arg = c('MDM2','p53','cJun','RUNX'), las=2)


#LPS141

barplot(c(3917/(3917+4133),450/(450+567),5938/(5938+17977),7396/(7396+15485)),col='blue3', names.arg = c('MDM2','p53','cJun','RUNX'), las=2, ylim=c(0,0.6))
abline(h = 0.3338672, lty=2, lwd=2)

barplot(c(sum(mdm2.PP.both>=n)/(3917+4133), sum(p53.PP.both>=n)/(450+567),
          sum(cjun.PP.both>=n)/(5938+17977),sum(runx.PP.both>=n)/(7396+15485)),
        col=c('purple3'), names.arg = c('MDM2','p53','cJun','RUNX'), las=2, ylim=c(0,0.35))


################################################################################
################################################################################

n=1

a <- LPS141.1.PP[which(LPS141.1.PP$MDM2.Promoter1==T&LPS141.1.PP$MDM2.Promoter2==T),]
mdm2.PP.both <- calculate.PP.frequency(a)


a <- LPS141.1.PP[which(LPS141.1.PP$p53.Promoter1==T&LPS141.1.PP$p53.Promoter2==T),]
p53.PP.both <- calculate.PP.frequency(a)


a <- LPS141.1.PP[which(LPS141.1.PP$cJun.Promoter1==T&LPS141.1.PP$cJun.Promoter2==T),]
cjun.PP.both <- calculate.PP.frequency(a)


a <- LPS141.1.PP[which(LPS141.1.PP$runx.Promoter1==T&LPS141.1.PP$runx.Promoter2==T),]
runx.PP.both <- calculate.PP.frequency(a)



n=1 #Number of P-P loops for a given promoter

pp.sums <- list(mdm2.number.promoters=c(),p53.number.promoters=c(),cjun.number.promoters=c(),runx.number.promoters=c())
for (z in c(1:1000)){
        
        LPS141.random <- randomize.TF.binding(lps141.1.PP,lps141.1.EP, 3917, 450, 5938,7396)
        print(z)
        a <- LPS141.random$PP.matrix[which(LPS141.random$PP.matrix$MDM2.Promoter1==T&LPS141.random$PP.matrix$MDM2.Promoter2==T),]
        mdm2.PP.random <- calculate.PP.frequency(a)
        pp.sums[['mdm2.number.promoters']] <- c(pp.sums[['mdm2.number.promoters']],sum(mdm2.PP.random>n))
        a <- LPS141.random$PP.matrix[which(LPS141.random$PP.matrix$p53.Promoter1==T&LPS141.random$PP.matrix$p53.Promoter2==T),]
        p53.PP.random <- calculate.PP.frequency(a)
        pp.sums[['p53.number.promoters']] <- c(pp.sums[['p53.number.promoters']],sum(p53.PP.random>n))
        a <- LPS141.random$PP.matrix[which(LPS141.random$PP.matrix$cJun.Promoter1==T&LPS141.random$PP.matrix$cJun.Promoter2==T),]
        cjun.PP.random <- calculate.PP.frequency(a)
        pp.sums[['cjun.number.promoters']] <- c(pp.sums[['cjun.number.promoters']],sum(cjun.PP.random>n))
        a <- LPS141.random$PP.matrix[which(LPS141.random$PP.matrix$runx.Promoter1==T&LPS141.random$PP.matrix$runx.Promoter2==T),]
        runx.PP.random <- calculate.PP.frequency(a)
        pp.sums[['runx.number.promoters']] <- c(pp.sums[['runx.number.promoters']],sum(runx.PP.random>n))
}

#load('lps141_pp_loops_sum_n1.RData')
n=1
plot(density(pp.sums$mdm2.number.promoters), main='MDM2 randomized PP loops (>1)', frame.plot=F)
abline(v = sum(mdm2.PP.both>n), lty=2)
sum(mdm2.PP.both>n)/mean(pp.sums$mdm2.number.promoters)

plot(density(pp.sums$p53.number.promoters), main='p53 randomized PP loops (>1)', frame.plot=F)
abline(v = sum(p53.PP.both>n), lty=2)
sum(p53.PP.both>n)/mean(pp.sums$p53.number.promoters)

plot(density(pp.sums$cjun.number.promoters), main='cJun randomized PP loops (>1)', frame.plot=F)
abline(v = sum(cjun.PP.both>n), lty=2)
sum(cjun.PP.both>n)/mean(pp.sums$cjun.number.promoters)

plot(density(pp.sums$runx.number.promoters), main='RUNX randomized PP loops (>1)', frame.plot=F)
abline(v = sum(runx.PP.both>n), lty=2)
sum(runx.PP.both>n)/mean(pp.sums$runx.number.promoters)

################################################################################
## Check overall degree of connectivity

promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_3_promDescription_cytoscape.tsv', sep='\t', header=T)
LPS141.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')
n=3
n=9
## Select a threshold to the number of PETs in the loop
LPS141.1.PP <- LPS141.1.PP[which(LPS141.1.PP$Score>n),]

write.table(LPS141.1.PP[,c(5,6,2)], file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_9_loops_cytoscape.tsv', quote = F, col.names = T, row.names = F, sep='\t')

promoter.description.file <- promoter.description.file[which(promoter.description.file$Promoter%in%c(LPS141.1.PP$Promoter1,LPS141.1.PP$Promoter2)),]

total.edges.per.promoter <- c()
for (promoter in promoter.description.file$Promoter){
        total.edges.per.promoter <- c(total.edges.per.promoter,length(which(LPS141.1.PP$Promoter1==promoter | LPS141.1.PP$Promoter2==promoter)))
}

description.with.degree <- cbind(promoter.description.file, total.edges.per.promoter)

boxplot(description.with.degree[which(description.with.degree$MDM2==T),'total.edges.per.promoter'],
        description.with.degree[which(description.with.degree$MDM2==F),'total.edges.per.promoter'])

################################################################################
## Define clusters

#LPS141.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_3_loop_annotation_file.tsv', 
#                          header=TRUE, sep = '\t')


promoter.list <- unique(c(LPS141.1.PP[,'Promoter1'],LPS141.1.PP[,'Promoter2']))


clusters.info <- calculate.clusters(LPS141.1.PP, promoter.list)
names(clusters.info) <- paste0('cluster_',c(1:length(clusters.info)))
#save(clusters.info, file='lps141_clusters_info.RData')

#load(file='lps141_clusters_info.RData')
#promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_3_promDescription_cytoscape.tsv', sep='\t', header=T)

cluster.order <- c()
for(promoter in promoter.description.file$Promoter){
        found=F
        for (cluster in names(clusters.info)){
                if (promoter%in%clusters.info[[cluster]]){
                        cluster.order <- c(cluster.order, cluster)
                        found=T
                }
        }
        if (found==F){
                print(promoter)
        }
}

promoter.description.file <- cbind(promoter.description.file,cluster.order)
promoter.degree.description <- cbind(description.with.degree,cluster.order)
#write.table(promoter.degree.description, file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', quote = F, col.names = T, row.names = F, sep='\t')
promoter.description.file <- read.table(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', sep = '\t', header=T)
promoter.degree.description <- read.table(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', sep = '\t', header=T)


total.mdm2 <- sum(promoter.description.file$MDM2); total.promoters <-length(unique(promoter.description.file$Promoter))
mdm2.per.cluster <- c(); elements.per.cluster <- c(); clusters.p.values <- c()
for (cluster in unique(promoter.description.file$cluster.order)){
        mdm2.binding <- sum(promoter.description.file[which(promoter.description.file$cluster.order==cluster),'MDM2'])
        elems.in.cluster <- length(promoter.description.file[which(promoter.description.file$cluster.order==cluster),'Promoter'])
        mdm2.per.cluster <- c(mdm2.per.cluster,mdm2.binding)
        elements.per.cluster <- c(elements.per.cluster,elems.in.cluster)
        clusters.p.values <- c(clusters.p.values,fisher.test(rbind(c(mdm2.binding,total.mdm2),c(elems.in.cluster-mdm2.binding,total.promoters-total.mdm2)), alternative='greater')$p.value)
}

names(mdm2.per.cluster) <- unique(promoter.description.file$cluster.order)
names(elements.per.cluster) <- unique(promoter.description.file$cluster.order)
names(clusters.p.values) <- unique(promoter.description.file$cluster.order)

mdm2.per.cluster.filt <- mdm2.per.cluster[which(elements.per.cluster>4)]
elements.per.cluster.filt <- elements.per.cluster[which(elements.per.cluster>4)]
clusters.p.values.filt <- clusters.p.values[which(elements.per.cluster>4)]

which(p.adjust(clusters.p.values.filt, method = 'BH')<0.05)
#write.table(unique(promoter.description.file[which(promoter.description.file$cluster.order%in%names(which(clusters.p.values.filt<0.05&mdm2.per.cluster.filt>4))&promoter.description.file$MDM2==T),'Gene']), file = 'aaaa.tsv',quote = F, row.names = F, col.names = F)

significant.clusters <- names(which(clusters.p.values<0.5&elements.per.cluster>4))
non.significant.clusters <- names(which(clusters.p.values>0.5&elements.per.cluster>4))
small.non.significant.clusters <- names(which(elements.per.cluster<4))

## "cluster_1"    "cluster_9"    "cluster_28"   "cluster_29"   "cluster_32"   "cluster_564"  "cluster_725"  "cluster_770"  "cluster_1128" "cluster_1288"
## "cluster_1309" "cluster_1383" "cluster_1394" "cluster_1511" "cluster_1517"

write.table(x = (promoter.description.file[which(promoter.description.file$cluster.order=='cluster_35'),'Gene']), file = 'aaaa.tsv',quote = F, row.names = F, col.names = F)

mdm2.bound.clusters <- unique(promoter.description.file[which(promoter.description.file$cluster.order%in%names(which(clusters.p.values.filt<0.05&mdm2.per.cluster.filt>4))&promoter.description.file$MDM2==T),'Gene'])

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/up_down_clusters_timecourse.RData')
mdm2.bound.clusters[which(mdm2.bound.clusters%in%names(down.genes))]

library(ComplexHeatmap)
library(circlize)

aaa <- cbind(promoter.description.file[which(promoter.description.file$cluster.order%in%names(which(clusters.p.values.filt<0.05&mdm2.per.cluster.filt>4))&promoter.description.file$MDM2==T),'Gene'],
      promoter.description.file[which(promoter.description.file$cluster.order%in%names(which(clusters.p.values.filt<0.05&mdm2.per.cluster.filt>4))&promoter.description.file$MDM2==T),'cluster.order'])
ht1 = ComplexHeatmap::Heatmap(FC.matrix[aaa[,1][aaa[,1]%in%rownames(FC.matrix)],],cluster_rows = F, cluster_columns = F)
ha_row = ComplexHeatmap::rowAnnotation(df = data.frame(aaa[,2][aaa[,1]%in%rownames(FC.matrix)]), 
                       col = list("cluster_35" =  "green", "cluster_233" = "orange","cluster_239" =  "red", "cluster_496" = "yellow",
                                  "cluster_1015" =  "gray", "cluster_1022" = "pink","cluster_1197" =  "blue"), width = unit(1, "cm"))

ht_list = ht1 + ha_row

draw(ht_list)

##################
## Add FC expression at Kd

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/knowckdown_cohort2/shRNA_1_4_DEGsvsNT.RData')
load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/DESeq_results.RData')

promoter.degree.description <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', sep = '\t', header=T)

expression.values.kd <- c();expression.values.tc <- c()
for (gene in promoter.degree.description$Gene){
  if (gene%in%rownames(results(dds, name='Group_shRNA4_vs_NT'))){
    gene.expr <- results(dds, name='Group_shRNA4_vs_NT')[gene,'log2FoldChange']
    expression.values.kd <- c(expression.values.kd,gene.expr)
  }else{
    expression.values.kd <- c(expression.values.kd, 0)
  }
  if (gene%in%rownames(res2vs0)){
    gene.expr <- res2vs0[gene,'log2FoldChange']
    expression.values.tc <- c(expression.values.tc,gene.expr)
  }else{
    expression.values.tc <- c(expression.values.tc, 0)
  }
}


promoter.degree.description.FCkd <- cbind(promoter.degree.description,expression.values.kd,expression.values.tc)

write.table(promoter.degree.description.FCkd,'/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_1/PP_loops/PP_9_promDescription_cyto_cluster_KDexpr.tsv', sep = '\t', col.names = T, row.names = F, quote=F)

##############################################################################
## Check expression at clusters

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/timecourse_count_matrix_design.RData')
load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/DESeq_results.RData')
a = promoter.description.file[which(promoter.description.file$cluster.order=='cluster_1288'&promoter.description.file$MDM2==T),'Gene']
b = promoter.description.file[which(promoter.description.file$cluster.order=='cluster_1288'&promoter.description.file$MDM2==F),'Gene']

boxplot(log2(tpm.mat[which(rownames(tpm.mat)%in%a),1]), log2(tpm.mat[which(rownames(tpm.mat)%in%b),1]))

FC.matrix <- cbind(res2vs0$log2FoldChange, res4vs0$log2FoldChange, res6vs0$log2FoldChange)
rownames(FC.matrix) <- rownames(res2vs0)

#1517
ComplexHeatmap::Heatmap(FC.matrix[c('ZNF250','ZNF16','ZNF251','ZNF517','ZNF7','COMMD5','ZNF34'),],)
ComplexHeatmap::Heatmap(FC.matrix[c('TSTA3','TIGD5','ZNF623','PYCR3','ZC3H3','MROH6','ZNF707','PUF60'),],)
ComplexHeatmap::Heatmap(FC.matrix[c('LNCOC1','LY6K','DENND3','DENND3','LY6E'),],)

#1288
ComplexHeatmap::Heatmap(FC.matrix[c('FAF2','ARL10','NOP16','CLTB','RNF44','GPRIN1','TSPAN17'),],)
ComplexHeatmap::Heatmap(FC.matrix[c('PDLIM7','B4GALT7','TMED9','DDX41','LMAN2','NSD1'),],)


## Compare with k27ac

k27ac <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_countPeaks_files.txt', sep = '\t', header = T) #24 and 25

cluster_x <- c('TLXNA','TMEM39B','KHDRBS1','PTP4A2','COL16A1','PEF1')

significant.clusters; non.significant.clusters

k27ac.exprs.gene <- calc.exprs.k27ac(promoter.description.file[which(promoter.description.file$cluster.order=='cluster_35'),'Gene'])
k27ac.exprs.gene <- calc.exprs.k27ac(promoter.description.file[which(promoter.description.file$cluster.order%in%significant.clusters),'Gene'])
non.sig.k27ac.exprs.gene <- calc.exprs.k27ac(promoter.description.file[which(promoter.description.file$cluster.order%in%non.significant.clusters),'Gene'])

plot(log2(k27ac.exprs.gene[[2]]),k27ac.exprs.gene[[1]], frame.plot=F, pch=19)
text(log2(k27ac.exprs.gene[[2]]),k27ac.exprs.gene[[1]], labels=k27ac.exprs.gene[[3]], pos=3)

k27ac$Gene.Name[which(k27ac$Gene.Name%in%cluster_x)]
tpm.mat[rownames(tpm.mat)%in%cluster_x,]



###############################################################################
## Test genes from significant clusters, genes from non-significant clusters and 
## Genes not belonging to clustrers at all

mdm2.significant.clusters <- promoter.description.file[which(promoter.description.file$cluster.order%in%significant.clusters&promoter.description.file$MDM2==T),'Gene']
mdm2Unbound.significant.clusters <- promoter.description.file[which(promoter.description.file$cluster.order%in%significant.clusters&promoter.description.file$MDM2==F),'Gene']
mdm2.nonSignificant.clusters <- promoter.description.file[which(promoter.description.file$cluster.order%in%non.significant.clusters&promoter.description.file$MDM2==T),'Gene']
mdm2Unbound.nonSignificant.clusters <- promoter.description.file[which(promoter.description.file$cluster.order%in%non.significant.clusters&promoter.description.file$MDM2==F),'Gene']
mdm2.small.nonSignificant.clusters <- promoter.description.file[which(promoter.description.file$cluster.order%in%small.non.significant.clusters&promoter.description.file$MDM2==T),'Gene']
mdm2Unbound.small.nonSignificant.clusters <- promoter.description.file[which(promoter.description.file$cluster.order%in%small.non.significant.clusters&promoter.description.file$MDM2==F),'Gene']

mdm2.significant.exprs <- calc.exprs.k27ac(mdm2.significant.clusters)
mdm2Unbound.significant.exprs <- calc.exprs.k27ac(mdm2Unbound.significant.clusters)
mdm2.nonSignificant.exprs <- calc.exprs.k27ac(mdm2.nonSignificant.clusters)
mdm2Unbound.nonSignificant.exprs <- calc.exprs.k27ac(mdm2Unbound.nonSignificant.clusters)
mdm2.small.nonSignificant.exprs <- calc.exprs.k27ac(mdm2.small.nonSignificant.clusters)
mdm2Unbound.small.nonSignificant.exprs <- calc.exprs.k27ac(mdm2Unbound.small.nonSignificant.clusters)

## Two categories
boxplot(c(mdm2.significant.exprs[[1]],mdm2.nonSignificant.exprs[[1]],mdm2Unbound.significant.exprs[[1]],mdm2Unbound.nonSignificant.exprs[[1]]), 
        c(mdm2.small.nonSignificant.exprs[[1]], mdm2Unbound.small.nonSignificant.exprs[[1]]), frame.plot=F, main='k27ac', names=c('Cluster','non-cluster'))
boxplot(c(log2(mdm2.significant.exprs[[2]]),log2(mdm2.nonSignificant.exprs[[2]]),log2(mdm2Unbound.significant.exprs[[2]]),log2(mdm2Unbound.nonSignificant.exprs[[2]])), 
        c(log2(mdm2.small.nonSignificant.exprs[[2]]), log2(mdm2Unbound.small.nonSignificant.exprs[[2]])), frame.plot=F, main='Gene expression', names=c('Cluster','non-cluster'))

#####
boxplot(c(mdm2.significant.exprs[[1]],mdm2.nonSignificant.exprs[[1]]),c(mdm2Unbound.significant.exprs[[1]],mdm2Unbound.nonSignificant.exprs[[1]]), 
        mdm2.small.nonSignificant.exprs[[1]], mdm2Unbound.small.nonSignificant.exprs[[1]], frame.plot=F, main='k27ac', names=c('MDM2','non-MDM2','MDM2','non-MDM2'),las=2)
boxplot(c(log2(mdm2.significant.exprs[[2]]),log2(mdm2.nonSignificant.exprs[[2]])),c(log2(mdm2Unbound.significant.exprs[[2]]),log2(mdm2Unbound.nonSignificant.exprs[[2]])), 
        log2(mdm2.small.nonSignificant.exprs[[2]]), log2(mdm2Unbound.small.nonSignificant.exprs[[2]]), frame.plot=F, main='Gene expression', names=c('MDM2','non-MDM2','MDM2','non-MDM2'),las=2)

#####
boxplot(mdm2.significant.exprs[[1]],mdm2Unbound.significant.exprs[[1]],mdm2.nonSignificant.exprs[[1]],mdm2Unbound.nonSignificant.exprs[[1]], 
        mdm2.small.nonSignificant.exprs[[1]], mdm2Unbound.small.nonSignificant.exprs[[1]], frame.plot=F, main='k27ac', names=c('MDM2','non-MDM2','MDM2','non-MDM2','MDM2','non-MDM2'),las=2)
boxplot(log2(mdm2.significant.exprs[[2]]),log2(mdm2Unbound.significant.exprs[[2]]),log2(mdm2.nonSignificant.exprs[[2]]),log2(mdm2Unbound.nonSignificant.exprs[[2]]), 
        log2(mdm2.small.nonSignificant.exprs[[2]]), log2(mdm2Unbound.small.nonSignificant.exprs[[2]]), frame.plot=F, main='Gene expression', names=c('MDM2','non-MDM2','MDM2','non-MDM2','MDM2','non-MDM2'),las=2)


###############################################################################
## Add expression in table

expression.values <- c()
for (gene in promoter.degree.description$Gene){
        if (gene%in%rownames(tpm.mat)){
                expression.values <- c(expression.values,mean(tpm.mat[gene,c(1,2)]))
        }else{
                expression.values <- c(expression.values,NA)
        }
}

promoter.degree.description.expression <- cbind(promoter.degree.description, expression.values)
which(!is.na(promoter.degree.description.expression$expression.values))

plot(log2(promoter.degree.description.expression$expression.values[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters)]),
     log2(promoter.degree.description.expression$total.edges.per.promoter[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters)]),
     pch=19, col=alpha('black',0.4),frame.plot=F, xlab = 'log2 Expression', ylab= ' Log2 Degree of connectivity')

points(log2(promoter.degree.description.expression$expression.values[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]),
       log2(promoter.degree.description.expression$total.edges.per.promoter[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]), 
       pch=19, col=alpha('orange3',1))

text(log2(promoter.degree.description.expression$expression.values[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]),
     log2(promoter.degree.description.expression$total.edges.per.promoter[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]),
     labels = promoter.degree.description.expression$cluster.order[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)])

#############################################################################################
#############################################################################################
#############################################################################################

lps853.1.EE <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/EE_loops/EE_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

lps853.1.EP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/EP_loops/EP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

lps853.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

lps853.p53.enhancers<- unique(c(lps853.1.EE[which(lps853.1.EE$p53.enhancer1==T),'enhancer1'], 
                                lps853.1.EE[which(lps853.1.EE$p53.enhancer2==T),'enhancer2'],
                                lps853.1.EP[which(lps853.1.EP$p53.Enhancer==T),'Enhancer']))

lps853.p53.promoters<- unique(c(lps853.1.PP[which(lps853.1.PP$p53.Promoter1),'promoter1'], 
                                lps853.1.PP[which(lps853.1.PP$p53.Promoter2),'promoter2'],
                                lps853.1.EP[which(lps853.1.EP$p53.Promoter),'Promoter']))

lps853.mdm2.enhancers<- unique(c(lps853.1.EE[which(lps853.1.EE$MDM2.enhancer1),'enhancer1'],
                                 lps853.1.EE[which(lps853.1.EE$MDM2.enhancer2),'enhancer2'],
                                 lps853.1.EP[which(lps853.1.EP$MDM2.Enhancer),'Enhancer']))

lps853.mdm2.promoters<- unique(c(lps853.1.PP[which(lps853.1.PP$MDM2.Promoter1),'promoter1'], 
                                 lps853.1.PP[which(lps853.1.PP$MDM2.Promoter2),'promoter2'],
                                 lps853.1.EP[which(lps853.1.EP$MDM2.Promoter),'Promoter']))


write.table(t(matrix(unlist(strsplit(lps853.p53.enhancers, split = '_')),nrow=3)),
            '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/anchors_bed/lps853_p53_enhancers.bed', sep = '\t',
            col.names = F, row.names = F, quote = F)

write.table(t(matrix(unlist(strsplit(lps853.p53.promoters, split = '_')),nrow=3)),
            '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/anchors_bed/lps853_p53_promoters.bed', sep = '\t',
            col.names = F, row.names = F, quote = F)

dev.off()

VennDiagram::draw.pairwise.venn(length(lps853.p53.enhancers), length(lps853.mdm2.enhancers), 
                                length(intersect(lps853.p53.enhancers,lps853.mdm2.enhancers)))
dev.off()
VennDiagram::draw.pairwise.venn(length(lps853.p53.promoters), length(lps853.mdm2.promoters), 
                                length(intersect(lps853.p53.promoters,lps853.mdm2.promoters)))



n=1

a <- LPS853.1.PP[which(LPS853.1.PP$MDM2.Promoter1==T&LPS853.1.PP$MDM2.Promoter2==T),]
mdm2.PP.both <- calculate.PP.frequency(a)


a <- LPS853.1.PP[which(LPS853.1.PP$p53.Promoter1==T&LPS853.1.PP$p53.Promoter2==T),]
p53.PP.both <- calculate.PP.frequency(a)


a <- LPS853.1.PP[which(LPS853.1.PP$cJun.Promoter1==T&LPS853.1.PP$cJun.Promoter2==T),]
cjun.PP.both <- calculate.PP.frequency(a)


a <- LPS853.1.PP[which(LPS853.1.PP$runx.Promoter1==T&LPS853.1.PP$runx.Promoter2==T),]
runx.PP.both <- calculate.PP.frequency(a)



n=3 #Number of P-P loops for a given promoter

pp.sums <- list(mdm2.number.promoters=c(),p53.number.promoters=c(),cjun.number.promoters=c(),runx.number.promoters=c())
for (z in c(1:1000)){
        
        LPS853.random <- randomize.TF.binding(LPS853.1.PP,LPS853.1.EP, 1996, 873, 6458, 8023)
        print(z)
        a <- LPS853.random$PP.matrix[which(LPS853.random$PP.matrix$MDM2.Promoter1==T&LPS853.random$PP.matrix$MDM2.Promoter2==T),]
        mdm2.PP.random <- calculate.PP.frequency(a)
        pp.sums[['mdm2.number.promoters']] <- c(pp.sums[['mdm2.number.promoters']],sum(mdm2.PP.random>n))
        a <- LPS853.random$PP.matrix[which(LPS853.random$PP.matrix$p53.Promoter1==T&LPS853.random$PP.matrix$p53.Promoter2==T),]
        p53.PP.random <- calculate.PP.frequency(a)
        pp.sums[['p53.number.promoters']] <- c(pp.sums[['p53.number.promoters']],sum(p53.PP.random>n))
        a <- LPS853.random$PP.matrix[which(LPS853.random$PP.matrix$cJun.Promoter1==T&LPS853.random$PP.matrix$cJun.Promoter2==T),]
        cjun.PP.random <- calculate.PP.frequency(a)
        pp.sums[['cjun.number.promoters']] <- c(pp.sums[['cjun.number.promoters']],sum(cjun.PP.random>n))
        a <- LPS853.random$PP.matrix[which(LPS853.random$PP.matrix$runx.Promoter1==T&LPS853.random$PP.matrix$runx.Promoter2==T),]
        runx.PP.random <- calculate.PP.frequency(a)
        pp.sums[['runx.number.promoters']] <- c(pp.sums[['runx.number.promoters']],sum(runx.PP.random>n))
}

#load('lps853_pp_loops_sum_n1.RData')
n=3
plot(density(pp.sums$mdm2.number.promoters), main='MDM2 randomized PP loops (>1)', frame.plot=F)
abline(v = sum(mdm2.PP.both>n), lty=2)
sum(mdm2.PP.both>n)/mean(pp.sums$mdm2.number.promoters)

plot(density(pp.sums$p53.number.promoters), main='p53 randomized PP loops (>1)', frame.plot=F)
abline(v = sum(p53.PP.both>n), lty=2)
sum(p53.PP.both>n)/mean(pp.sums$p53.number.promoters)

plot(density(pp.sums$cjun.number.promoters), main='cJun randomized PP loops (>1)', frame.plot=F)
abline(v = sum(cjun.PP.both>n), lty=2)
sum(cjun.PP.both>n)/mean(pp.sums$cjun.number.promoters)

plot(density(pp.sums$runx.number.promoters), main='RUNX randomized PP loops (>1)', frame.plot=F)
abline(v = sum(runx.PP.both>n), lty=2)
sum(runx.PP.both>n)/mean(pp.sums$runx.number.promoters)


#LPS853

barplot(rbind(c(2237/(2237+1012),913/(913+1514),10324/(10324+26469),12029/(12029+19511)),
              c(1012/(2237+1012),1514/(913+1514),26469/(10324+26469),19511/(12029+19511))), 
        col=c('blue3', alpha('blue3', 0.2)), names.arg = c('MDM2','p53','cJun','RUNX'), las=2)

barplot(c(1996/(1996+951),873/(873+1462),6458/(6458+18542),8023/(8023+14466)),col='blue3', names.arg = c('MDM2','p53','cJun','RUNX'), las=2, ylim=c(0,0.7))
abline(h = 0.3497655, lty=2, lwd=2)

n=3

barplot(c(sum(mdm2.PP.both>=n)/(1996+951), sum(p53.PP.both>=n)/(873+1462),
          sum(cjun.PP.both>=n)/(6458+18542),sum(runx.PP.both>=n)/(8023+14466)),
        col=c('purple3'), names.arg = c('MDM2','p53','cJun','RUNX'), las=2, ylim=c(0,0.5))


sum(mdm2.PP.both>=n)/sum(mdm2.PP.random>=n)
sum(p53.PP.both>=n)/sum(p53.PP.random>=n)
sum(cjun.PP.both>=n)/sum(cjun.PP.random>=n)
sum(runx.PP.both>=n)/sum(runx.PP.random>=n)


################################################################################
## Check overall degree of connectivity

promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_promDescription_cytoscape.tsv', sep='\t', header=T)
LPS853.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

total.edges.per.promoter <- c()
for (promoter in promoter.description.file$Promoter){
        total.edges.per.promoter <- c(total.edges.per.promoter,length(which(LPS853.1.PP$Promoter1==promoter | LPS853.1.PP$Promoter2==promoter)))
}

description.with.degree <- cbind(promoter.description.file, total.edges.per.promoter)

boxplot(log2(description.with.degree[which(description.with.degree$MDM2==T),'total.edges.per.promoter']),
        log2(description.with.degree[which(description.with.degree$MDM2==F),'total.edges.per.promoter']))

################################################################################
## Define clusters

LPS853.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

promoter.list <- unique(c(LPS853.1.PP[,'Promoter1'],LPS853.1.PP[,'Promoter2']))


clusters.info <- calculate.clusters(LPS853.1.PP, promoter.list)
names(clusters.info) <- paste0('cluster_',c(1:length(clusters.info)))
#save(clusters.info, file='clusters_info.RData')

load(file='clusters_info.RData')
promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_promDescription_cytoscape.tsv', sep='\t', header=T)

cluster.order <- c()
for(promoter in promoter.description.file$Promoter){
        found=F
        for (cluster in names(clusters.info)){
                if (promoter%in%clusters.info[[cluster]]){
                        cluster.order <- c(cluster.order, cluster)
                        found=T
                }
        }
        if (found==F){
                print(promoter)
        }
}

promoter.description.file <- cbind(promoter.description.file,cluster.order)
promoter.degree.description <- cbind(description.with.degree,cluster.order)

total.mdm2 <- sum(promoter.description.file$MDM2); total.promoters <-length(unique(promoter.description.file$Promoter))
mdm2.per.cluster <- c(); elements.per.cluster <- c(); clusters.p.values <- c()
for (cluster in unique(promoter.description.file$cluster.order)){
        mdm2.binding <- sum(promoter.description.file[which(promoter.description.file$cluster.order==cluster),'MDM2'])
        elems.in.cluster <- length(promoter.description.file[which(promoter.description.file$cluster.order==cluster),'Promoter'])
        mdm2.per.cluster <- c(mdm2.per.cluster,mdm2.binding)
        elements.per.cluster <- c(elements.per.cluster,elems.in.cluster)
        clusters.p.values <- c(clusters.p.values,fisher.test(rbind(c(mdm2.binding,total.mdm2),c(elems.in.cluster-mdm2.binding,total.promoters-total.mdm2)))$p.value)
}

names(mdm2.per.cluster) <- unique(promoter.description.file$cluster.order)
names(elements.per.cluster) <- unique(promoter.description.file$cluster.order)
names(clusters.p.values) <- unique(promoter.description.file$cluster.order)

mdm2.per.cluster.filt <- mdm2.per.cluster[which(elements.per.cluster>20)]
elements.per.cluster.filt <- elements.per.cluster[which(elements.per.cluster>20)]
clusters.p.values.filt <- clusters.p.values[which(elements.per.cluster>20)]

which(p.adjust(clusters.p.values.filt, method = 'BH')<0.05)
write.table(unique(promoter.description.file[which(promoter.description.file$cluster.order%in%names(which(clusters.p.values.filt<0.05&mdm2.per.cluster.filt>5))),'Gene']), file = 'aaaa.tsv',quote = F, row.names = F, col.names = F)

significant.clusters <- names(which(clusters.p.values.filt<0.05&mdm2.per.cluster.filt>5))
                                                             
#write.table(x = (promoter.description.file[which(promoter.description.file$cluster.order=='cluster_621'),'Gene']), file = 'aaaa.tsv',quote = F, row.names = F, col.names = F)

##############################################################################
## Check expression at clusters

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/matrices/timecourse_count_matrix_design.RData')

a = promoter.description.file[which(promoter.description.file$cluster.order=='cluster_539'&promoter.description.file$MDM2==T),'Gene']
b = promoter.description.file[which(promoter.description.file$cluster.order=='cluster_539'&promoter.description.file$MDM2==F),'Gene']

boxplot(log2(tpm.mat[which(rownames(tpm.mat)%in%a),1]), log2(tpm.mat[which(rownames(tpm.mat)%in%b),1]))


###############################################################################
## Add expression in table

expression.values <- c()
for (gene in promoter.degree.description$Gene){
        if (gene%in%rownames(tpm.mat)){
                expression.values <- c(expression.values,mean(tpm.mat[gene,c(1,2)]))
        }else{
                expression.values <- c(expression.values,NA)
        }
}

promoter.degree.description.expression <- cbind(promoter.degree.description, expression.values)
which(!is.na(promoter.degree.description.expression$expression.values))

plot(log2(promoter.degree.description.expression$expression.values[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters)]),
     log2(promoter.degree.description.expression$total.edges.per.promoter[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters)]),
     pch=19, col=alpha('black',0.4),frame.plot=F, xlab = 'log2 Expression', ylab= ' Log2 Degree of connectivity')

points(log2(promoter.degree.description.expression$expression.values[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]),
      log2(promoter.degree.description.expression$total.edges.per.promoter[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]), 
      pch=19, col=alpha('orange3',1))

text(log2(promoter.degree.description.expression$expression.values[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]),
     log2(promoter.degree.description.expression$total.edges.per.promoter[which(!is.na(promoter.degree.description.expression$expression.values)&promoter.degree.description.expression$cluster.order%in%significant.clusters&promoter.degree.description.expression$MDM2==T)]),
     labels = c())

###############################################################################
##
promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_2/PP_loops/PP_3_promDescription_cytoscape.tsv', sep='\t', header=T)
LPS853.2.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_2/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')
n=3
n=9
## Select a threshold to the number of PETs in the loop
LPS853.2.PP <- LPS853.2.PP[which(LPS853.2.PP$Score>n),]

write.table(LPS853.2.PP[,c(5,6,2)], file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_2/PP_loops/PP_9_loops_cytoscape.tsv', quote = F, col.names = T, row.names = F, sep='\t')

promoter.description.file <- promoter.description.file[which(promoter.description.file$Promoter%in%c(LPS853.2.PP$Promoter1,LPS853.2.PP$Promoter2)),]

write.table(promoter.description.file, file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_2/PP_loops/PP_9_promDescription_cytoscape.tsv', quote = F, col.names = T, row.names = F, sep='\t')

###############################################################################
## Create clusters

LPS853.1.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

LPS853.1.PP <- LPS853.1.PP[LPS853.1.PP$Score>9,]

promoter.list <- unique(c(LPS853.1.PP[,'Promoter1'],LPS853.1.PP[,'Promoter2']))


clusters.info <- calculate.clusters(LPS853.1.PP, promoter.list)
names(clusters.info) <- paste0('cluster_',c(1:length(clusters.info)))
#save(clusters.info, file='lps853_clusters_info.RData')

#load(file='lps853_clusters_info.RData')
#promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cytoscape.tsv', sep='\t', header=T)
#write.table(promoter.description.file[which(promoter.description.file$Promoter%in%promoter.list),], file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cytoscape.tsv', quote=F, row.names = F, sep = '\t')

cluster.order <- c()
for(promoter in promoter.description.file$Promoter){
  found=F
  for (cluster in names(clusters.info)){
    if (promoter%in%clusters.info[[cluster]]){
      cluster.order <- c(cluster.order, cluster)
      found=T
    }
  }
  if (found==F){
    print(promoter)
  }
}

load('/Volumes/seq_epiprod02/sgaldon/liposarcoma/rna_master_datasets/knowckdown_cohort2/shRNA_1_4_DEGsvsNT.RData')

results(dds, name='Group_shRNA4_vs_NT')

promoter.description.file <- cbind(promoter.description.file,cluster.order)
promoter.degree.description <- cbind(description.with.degree,cluster.order)
#write.table(promoter.description.file, file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', quote = F, col.names = T, row.names = F, sep='\t')
promoter.description.file <- read.table(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_promDescription_cyto_cluster.tsv', sep = '\t', header=T)
promoter.degree.description <- read.table(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_promDescription_cyto_cluster.tsv', sep = '\t', header=T)

#promoter.degree.description <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', sep = '\t', header=T)


###############################################################################
## Functions


randomize.TF.binding <- function(PP.loop.matrix=PP.loop.matrix, EP.loop.matrix=EP.loop.matrix,mdm2.number=mdm2.number,p53.number=p53.number,cjun.number=cjun.number,runx.number=runx.number){
        #LPS853: MDM2(1996), p53(873), cJun(6458), RUNX(8023)
        #LPS141: MDM2(3917), p53(450), cJun(5938), RUNX(7396)
        promoters.list <- unique(c(PP.loop.matrix$Promoter1,PP.loop.matrix$Promoter2,EP.loop.matrix$Promoter))
        
        mdm2.sample <- sample(promoters.list,mdm2.number); p53.sample <- sample(promoters.list,p53.number)
        cjun.sample <- sample(promoters.list,cjun.number); runx.sample <- sample(promoters.list,runx.number)
        
        new.EP.matrix <- EP.loop.matrix[,c(1:13)]; new.EP.matrix[,c(6:13)] <- F
        new.EP.matrix[which(new.EP.matrix[,'Promoter']%in%mdm2.sample),'MDM2.Promoter'] <- T
        new.EP.matrix[which(new.EP.matrix[,'Promoter']%in%p53.sample),'p53.Promoter'] <- T
        new.EP.matrix[which(new.EP.matrix[,'Promoter']%in%cjun.sample),'cJun.Promoter'] <- T
        new.EP.matrix[which(new.EP.matrix[,'Promoter']%in%runx.sample),'RUNX.Promoter'] <- T
        
        new.PP.matrix <- PP.loop.matrix[,c(1:14)]; new.PP.matrix[,c(7:14)] <- F
        new.PP.matrix[which(new.PP.matrix[,'Promoter1']%in%mdm2.sample),'MDM2.Promoter1'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter2']%in%mdm2.sample),'MDM2.Promoter2'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter1']%in%p53.sample),'p53.Promoter1'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter2']%in%p53.sample),'p53.Promoter2'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter1']%in%cjun.sample),'cJun.Promoter1'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter2']%in%cjun.sample),'cJun.Promoter2'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter1']%in%runx.sample),'runx.Promoter1'] <- T
        new.PP.matrix[which(new.PP.matrix[,'Promoter2']%in%runx.sample),'runx.Promoter2'] <- T
        
        return(list('PP.matrix'=new.PP.matrix,'EP.matrix'=new.EP.matrix))
}




calculate.PP.frequency <- function(PP.matrix){
        TF.PP.random <- list()
        for (loop in c(1:dim(PP.matrix)[1])){
                if (PP.matrix[loop,'Promoter1']%in%names(TF.PP.random)){
                        TF.PP.random[PP.matrix[loop,'Promoter1']] <- unlist(TF.PP.random[PP.matrix[loop,'Promoter1']])+1
                } else{
                        TF.PP.random[PP.matrix[loop,'Promoter1']] = 1
                }
                if (PP.matrix[loop,'Promoter2']%in%names(TF.PP.random)){
                        TF.PP.random[PP.matrix[loop,'Promoter2']] <- unlist(TF.PP.random[PP.matrix[loop,'Promoter2']])+1
                } else{
                        TF.PP.random[PP.matrix[loop,'Promoter2']] = 1
                }
        }
        return(TF.PP.random)
}


calculate.clusters <- function(PP.loop.matrix, promoters.list){
        cluster.list <- list()
        n = 0
        for (promoter in promoter.list){
                if(promoter%in%PP.loop.matrix$Promoter1){
                        connected.promoter <- PP.loop.matrix[which(PP.loop.matrix$Promoter1==promoter),'Promoter2']
                        if(promoter%in%unlist(cluster.list)){
                                for (elem in c(1:n)){
                                        if (promoter%in%cluster.list[[elem]]){
                                                cluster.list[[elem]] <- unique(c(cluster.list[[elem]],promoter, connected.promoter))}
                                }
                        } else if(any(connected.promoter%in%unlist(cluster.list))){
                                for (elem in c(1:n)){
                                        if (any(c(promoter,connected.promoter)%in%cluster.list[[elem]])){
                                                cluster.list[[elem]] <- unique(c(cluster.list[[elem]],promoter, connected.promoter))}
                                }}else{
                                        n = n+1
                                        cluster.list[[n]] <- c(promoter, connected.promoter)
                                }
                }else if(promoter%in%PP.loop.matrix$Promoter2){
                        connected.promoter <- PP.loop.matrix[which(PP.loop.matrix$Promoter2==promoter),'Promoter1']
                        if(promoter%in%unlist(cluster.list)){
                                for (elem in c(1:n)){
                                        if (promoter%in%cluster.list[[elem]]){
                                                cluster.list[[elem]] <- unique(c(cluster.list[[elem]],promoter, connected.promoter))}
                                }
                        } else if(any(connected.promoter%in%unlist(cluster.list))){
                                for (elem in c(1:n)){
                                        if (any(c(promoter,connected.promoter)%in%cluster.list[[elem]])){
                                                cluster.list[[elem]] <- unique(c(cluster.list[[elem]],promoter, connected.promoter))}
                                }}else{
                                        n = n+1
                                        cluster.list[[n]] <- c(promoter, connected.promoter)
                                }
                }
        }
        names(cluster.list) <- paste0('cluster_',c(1:length(cluster.list)))
        for (cluster1 in names(cluster.list)){
                for (cluster2 in names(cluster.list)){
                        if (!cluster1==cluster2){
                                if(any(cluster.list[[cluster1]]%in%cluster.list[[cluster2]])){
                                        cluster.list[[cluster2]] <- unique(c(cluster.list[[cluster1]],cluster.list[[cluster2]]))
                                        cluster.list <- cluster.list[names(cluster.list) != cluster1]
                                }
                        }
                }
        }
        return(cluster.list)
}


calc.exprs.k27ac <- function(gene.ids){
  gene.h3k27ac <- c(); gene.exprs <- c(); gene.name <- c()
  for (gene in gene.ids){
    if (gene%in%k27ac$Gene.Name){
      if (gene%in%rownames(tpm.mat)){
        if(length(grep('promoter',k27ac[k27ac$Gene.Name==gene,c(8)]))>0|any(abs(k27ac[k27ac$Gene.Name==gene,c(10)])<5000)){
          to.use <- which(apply((k27ac[k27ac$Gene.Name==gene,c(24,25)]+1)/(k27ac[k27ac$Gene.Name==gene,c(36,37)]+1), 1, mean)>2&(abs(k27ac[k27ac$Gene.Name==gene,c(10)])<5000))
          if(length(to.use)>0){
            levels <- mean(apply((k27ac[k27ac$Gene.Name==gene,c(24,25)]+1)/(k27ac[k27ac$Gene.Name==gene,c(36,37)]+1), 1, mean)[to.use])
          }else{levels = 0}
        }else{levels = 0}
        gene.h3k27ac <- c(gene.h3k27ac,levels)
        gene.exprs <- c(gene.exprs, mean(tpm.mat[gene,c(1,2)]))
        gene.name <- c(gene.name,gene)
      }
    }
  }
  return(list(gene.h3k27ac,gene.exprs,gene.name))
}




###############################################################################
## Tumors

DD10.PP <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/DD10/PP_loops/PP_3_loop_annotation_file.tsv', 
                          header=TRUE, sep = '\t')

DD10.PP <- DD10.PP[DD10.PP$Score>9,]

promoter.list <- unique(c(DD10.PP[,'Promoter1'],DD10.PP[,'Promoter2']))


clusters.info <- calculate.clusters(DD10.PP, promoter.list)
names(clusters.info) <- paste0('cluster_',c(1:length(clusters.info)))
#save(clusters.info, file='lps853_clusters_info.RData')

#load(file='lps853_clusters_info.RData')
#promoter.description.file <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cytoscape.tsv', sep='\t', header=T)
#write.table(promoter.description.file[which(promoter.description.file$Promoter%in%promoter.list),], file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cytoscape.tsv', quote=F, row.names = F, sep = '\t')

cluster.order <- c()
for(promoter in promoter.description.file$Promoter){
  found=F
  for (cluster in names(clusters.info)){
    if (promoter%in%clusters.info[[cluster]]){
      cluster.order <- c(cluster.order, cluster)
      found=T
    }
  }
  if (found==F){
    print(promoter)
  }
}

promoter.description.file <- cbind(promoter.description.file,cluster.order)
promoter.degree.description <- cbind(description.with.degree,cluster.order)
#write.table(promoter.description.file, file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_9_promDescription_cyto_cluster.tsv', quote = F, col.names = T, row.names = F, sep='\t')
promoter.description.file <- read.table(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_promDescription_cyto_cluster.tsv', sep = '\t', header=T)
promoter.degree.description <- read.table(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS853_1/PP_loops/PP_3_promDescription_cyto_cluster.tsv', sep = '\t', header=T)

