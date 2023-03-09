####################################
## Test expression estability in hubs


## Load LPS TMP data and promoter hubs info

hub.genes <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/geneLists/genesinhubs.txt')
mdm2.hub.genes <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/geneLists/mdm2_bound_genesinhubs.txt')

hub.clusters <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/pphub_metrics/consensus_promoters_clusters_lps141_853.bed')
hub.clusters <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/pphub_metrics/mdm2Enriched_consensus_promoters_clusters_lps141_853.bed')


x <- read.csv("~/Dropbox (Partners HealthCare)/Liposarcoma/FINAL_COHORT_SAMP_ID.csv")

a <- read.delim("~/Dropbox (Partners HealthCare)/Liposarcoma/Computational analyses/Chadi_analysis/RNA-Seq_All_Cohorts/Liposarcoma_Cohort1_TPM_Matrix.txt", sep = "\t", row.names = 1, header = TRUE, dec = '.',skip = 1)

b <- read.delim("~/Dropbox (Partners HealthCare)/Liposarcoma/Computational analyses/Chadi_analysis/RNA-Seq_All_Cohorts/Liposarcoma_Cohort2_TPM_Matrix.txt", sep = "\t", row.names = 1, header = TRUE, dec = '.',skip = 1)

DD.N.DEseq <- read.table('~/Dropbox (Partners HealthCare)/Liposarcoma/Computational analyses/Sam_analysis/RNAseq/primary tumors/210330_DEseq/DD_N_deseq.txt')

y <- cbind(a,b)[c(2:dim(a)[1]),]

colnames(y) <-c('DD10', 'DD15', 'DD1', 'DD24', 'DD27', 'DD2', 'DD3', 'DD_SS',  'N10',
                'N12',  'N14',  'N27',  'N6',  'WD10',  'WD12',  'WD14', 'WD1', 'WD21', 'WD4',
                'WD5', 'DD15_24', 'DD15_27', "DD20",    "DD24",    "DD30",    "DD31",    "DD35",
                "DD36",    "N15" ,    "N19" ,    "N22" ,    "N30" ,    "N31"  ,   "N33"  ,   "N34",
                "WD13"  ,"WD14",    "WD15",    "WD19"  ,  "WD1" ,    "WD22" ,   "WD23",
                "WD24"  ,  "WD28" ,   "WD30"  ,  "WD33"  ,  "WD34" ,   "WD36")
IDs <- x$Sample_ID

high_adiposity <- c('WD36','WD34','WD24','WD33','WD12','WD15','WD23','WD22','WD4','WD10','WD14')

low_adiposity <- c('WD19','WD13','WD28','WD30','WD5','DD30','DD20','DD10','DD2','DD1','DD15','DD35','DD36','DD31','DD_SS')

exp <- y[,colnames(y) %in% IDs]
row.names(exp) <- row.names(y)
exp <- exp[,order(names(exp))]


DD_N <- exp[,1:22]
DD_WD <- exp[,c(1:10,23:39)]
N_WD <- exp[,11:39]
N_DD_WD <- exp[,c(11:22,1:10,23:39)]

high.adip <- exp[,high_adiposity]
low.adip <- exp[,low_adiposity]

DD_N <- DD_N[!rowSums(DD_N)<5,]
DD_N <- DD_N[which(rowMeans(!is.na(DD_N)) > 0.5), ]

DD_WD <- DD_WD[!rowSums(DD_WD)<5,]
DD_WD <- DD_WD[which(rowMeans(!is.na(DD_WD)) > 0.5), ]

N_DD_WD <- N_DD_WD[!rowSums(N_DD_WD)<5,]
N_DD_WD <- N_DD_WD[which(rowMeans(!is.na(N_DD_WD)) > 0.5), ]

high.adip <- high.adip[!rowSums(high.adip)<5,]
high.adip <- high.adip[which(rowMeans(!is.na(high.adip)) > 0.5), ]

low.adip <- low.adip[!rowSums(low.adip)<5,]
low.adip <- low.adip[which(rowMeans(!is.na(low.adip)) > 0.5), ]

DEGenes <- rownames(DD.N.DEseq)[DD.N.DEseq$pvalue<0.05&(!is.na(DD.N.DEseq$pvalue))]
DEGenes<-DEGenes[DEGenes%in%rownames(N_DD_WD)]

subset <- DEGenes
subset <- rownames(N_DD_WD)

p.values <- c()

for (n in c(1:1000)){

  DD.mean.cor <- c(); WD.mean.cor <- c(); N.mean.cor <- c()
  DD.mean.cor.random <- c(); WD.mean.cor.random <- c(); N.mean.cor.random <- c()
  
  high.mean.cor <- c(); low.mean.cor <- c()
  
  for (cluster in unique(hub.clusters$V4)){
    subcluster <- hub.clusters[hub.clusters$V4==cluster,]
    subcluster.genes <- unique(subcluster$V9)
    if (sum(subcluster.genes%in%rownames(N_DD_WD))>2){
      DD.cormat <- cor(t(N_DD_WD[rownames(N_DD_WD)%in%subcluster.genes,c(13:22)]))
      DD.cormat.random <- cor(t(N_DD_WD[sample(subset,size = sum(subcluster.genes%in%rownames(N_DD_WD))),c(13:22)]))
      WD.cormat <- cor(t(N_DD_WD[rownames(N_DD_WD)%in%subcluster.genes,c(23:39)]))
      WD.cormat.random <- cor(t(N_DD_WD[sample(subset,size = sum(subcluster.genes%in%rownames(N_DD_WD))),c(23:39)]))
      N.cormat <- cor(t(N_DD_WD[rownames(N_DD_WD)%in%subcluster.genes,c(1:12)]))
      N.cormat.random <- cor(t(N_DD_WD[sample(subset,size = sum(subcluster.genes%in%rownames(N_DD_WD))),c(1:12)]))
      
      high.cormat <- cor(t(high.adip[rownames(high.adip)%in%subcluster.genes,]))
      low.cormat <- cor(t(low.adip[rownames(low.adip)%in%subcluster.genes,]))
      
      if (dim(DD.cormat)[1]>3){
        corlowtri <- DD.cormat[lower.tri(DD.cormat)]
        DD.mean.cor <- c(DD.mean.cor,mean(abs(corlowtri)))
        
        corlowtri <- DD.cormat.random[lower.tri(DD.cormat.random)]
        DD.mean.cor.random <- c(DD.mean.cor.random,mean(abs(corlowtri)))
      }
      
      if (dim(WD.cormat)[1]>3){
        corlowtri <- WD.cormat[lower.tri(WD.cormat)]
        WD.mean.cor <- c(WD.mean.cor,mean(abs(corlowtri)))
        
        corlowtri <- WD.cormat.random[lower.tri(WD.cormat.random)]
        WD.mean.cor.random <- c(WD.mean.cor.random,mean(abs(corlowtri)))
      }
      if (dim(N.cormat)[1]>3){
        corlowtri <- N.cormat[lower.tri(N.cormat)]
        N.mean.cor <- c(N.mean.cor,mean(abs(corlowtri)))
        
        corlowtri <- N.cormat.random[lower.tri(N.cormat.random)]
        N.mean.cor.random <- c(N.mean.cor.random,mean(abs(corlowtri)))
      }
      
      if (dim(high.cormat)[1]>3){
        corlowtri <- high.cormat[lower.tri(high.cormat)]
        high.mean.cor <- c(high.mean.cor,mean(abs(corlowtri)))
      }
      
      if (dim(low.cormat)[1]>3){
        corlowtri <- low.cormat[lower.tri(low.cormat)]
        low.mean.cor <- c(low.mean.cor,mean(abs(corlowtri)))
      }
    }
  }
  p.values <- c(p.values,t.test(DD.mean.cor-DD.mean.cor.random,N.mean.cor-N.mean.cor.random)$p.value)
}

pdf(file = 'pphub_correlation_dd_wd_n.pdf', width=3.7,height = 5)
boxplot(DD.mean.cor,WD.mean.cor,N.mean.cor, names=c('DD','WD','N'),main='PPhub genes',
        ylab='Correlation', frame.plot=F,ylim=c(0,1))
dev.off()
t.test(DD.mean.cor,N.mean.cor)

pdf(file = 'pphub_correlation_ddlps.pdf', width=3.7,height = 5)
boxplot(DD.mean.cor,DD.mean.cor.random, names=c('pphub genes','random genes'),las=2,main='PPhub corr in DDLPS',
        ylab='Correlation', frame.plot=F,cex.axis=0.7,ylim=c(0,1))
dev.off()

