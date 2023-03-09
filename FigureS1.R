## Plot the TCGA mdm2 expression distribution in ddlps vs normal

xena.study <- read.table('EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena', sep='\t', header=T)
mdm2.expr <- as.numeric(xena.study[10744,2:dim(xena.study)[2]])

sample.id <- colnames(xena.study)
sample.id <- sample.id[2:length(sample.id)]

names(mdm2.expr) <- sample.id

id.code <- unlist(lapply(sample.id, function(x){strsplit(x, split='.', fixed=T)[[1]][2]}))
sample.type <- as.numeric(unlist(lapply(sample.id, function(x){strsplit(x, split='.', fixed=T)[[1]][4]})))

tissue.source <- read.table('tcga_code_tables/tissueSourceSite.tsv', sep='\t', header=T)
code.ids <- c(); code.name <- c()
for (code in unique(id.code)){
  if (code%in%tissue.source$TSS.Code){
    code.name <- c(code.name, tissue.source[which(tissue.source$TSS.Code==code),"Study.Name"])
    code.ids <- c(code.ids, code)
  }
}
names(code.name) <- code.ids
sample.ids.code <- code.name[id.code]
#names(mdm2.expr) <- code.name[id.code] <- we dont want to change the names, they have the same order

sample.type.data <- read.table('tcga_code_tables/sampleType.tsv', sep='\t', header=T, row.names=1)
sample.type <- sample.type.data[sample.type,]
rownames(sample.type) <- names(mdm2.expr)


mdm2.cnv <- read.table('MDM2_cnv.tsv', header=T, sep='\t', row.names = 1)
mdm2.cnv.names <- gsub("-", ".", rownames(mdm2.cnv), fixed = TRUE)
mdm2.cnv.names <- substr(mdm2.cnv.names,1,nchar(mdm2.cnv.names)-1)

mdm2.cnv <- mdm2.cnv[!duplicated(mdm2.cnv.names),]

names(mdm2.cnv) <- mdm2.cnv.names[!duplicated(mdm2.cnv.names)]

ddlps.samples <- read.table('DDLPS_identifiers.txt')$V1
ddlps.samples <- gsub("-", ".", ddlps.samples, fixed = TRUE)[c(1:length(ddlps.samples)-1)]


##########################################################
## We have defined:
## mdm2 expression -> with TCGA names
## sample.type -> same order as mdm2 expression, only containing the sample type
## sample.ids.code -> same order as mdm2 expression, only containing the tumor type
## mdm2.cnv -> mdm2 amplification
## ddlps.samples -> samples corresponding to ddlps tumors

mdm2.cnv[names(mdm2.cnv)%in%rownames(sample.type)[sample.type$Short.Letter.Code=='NT']]

avg.mdm2.exprs <- c()
for (tumor in unique(names(mdm2.expr))){
  avg.mdm2.exprs <- c(avg.mdm2.exprs, mean(mdm2.expr[which(names(mdm2.expr)==tumor)], na.rm=T))
}

names(avg.mdm2.exprs) <- unique(names(mdm2.expr))

normal.d <- density(mdm2.expr[rownames(sample.type)[sample.type$Short.Letter.Code=='NT']])
ddlps.d <- density(mdm2.expr[ddlps.samples], na.rm=T)

pdf('ddlps_normal_mdm2_expression.pdf', height = 5, width = 5)
plot(normal.d, xlim=c(8,17), 
     frame.plot=F, main='Normal - DDLPS MDM2 expression', xlab='MDM2 expression', col=alpha('darkblue', 0.4))
polygon(normal.d, col=alpha("darkblue", 0.4), border="darkblue")
lines(ddlps.d, col=alpha('darkorange', 0.4))
polygon(ddlps.d, col=alpha("darkorange", 0.4), border="darkorange")

dev.off()


## Plot the adipocity score as calculated from GSEA dataset

library(GSVA)

adipogenesis <- read.table('/Volumes/seq_epiprod02/sgaldon/liposarcoma/geneset_signatures/HALLMARK_ADIPOGENESIS.v2022.1.Hs.grp', 
                       sep='\t', header = T)[,1]

GeneLists <- list(adipogenesis)

# Tumors
x <- read.table("../Cohort2/Liposarcoma_Cohort2_TPM_Matrix.txt", head=TRUE, sep="\t", row.names=1)
x.2 <- read.table("../Cohort1/Liposarcoma_Cohort1_TPM_Matrix.txt", head=TRUE, sep="\t", row.names=1)

data.order <- c('N34','N15','N22','N33','N31','N30','N14','N19','N10','N12',
                'WD36','WD34','WD24','WD33','WD12','WD15','WD23','WD22','WD4','WD10','WD14','WD19','WD13','WD28','WD30','WD5',
                'DD30','DD20','DD10','DD2','DD1','DD15','DD35','DD36','DD31','DD_SS')
data <- cbind(x, x.2)[,data.order]

M <- as.matrix(data)

G <- gsva(M, GeneLists)


pdf(file = '/Volumes/seq_epiprod02/sgaldon/liposarcoma/geneset_signatures/Figures/mesenchymal_adipocytic_signatures.pdf')

par(mfrow=c(2, 2), mar=c(5.1, 4.1, 4.1, 2.1))
plot(c(1:36),t(G[1,]), col = c(rep('darkred', 10), rep('darkblue', 16), rep('grey', 10)), pch=19, 
     frame.plot=F, main='ADIPOGENESIS', ylab = 'GSVA score', xlab='')
dev.off()
