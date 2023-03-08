
## Generate the bed files using the metagenes script for figure 3

###############################################################################
## P53 binding on P53 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/p53_merged_geneset_lps141_220505_p53_SInorm.gz', 
                     sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",1800)),dec=".") 

rownames <- myData[,4]

p53.matrix <- myData[,-c(1:6)]
rownames(p53.matrix) <- rownames


x <- 1:300

y <- apply(p53.matrix[,c(1:300)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(301:600)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(601:900)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(901:1200)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(1201:1500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(1501:1800)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_p53Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#42BB6E', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='P53 regions - P53 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#42BB6E",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#00741F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#00741F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_p53Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_p53Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_p53Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()




###############################################################################
## MDM2 binding on P53 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/p53_merged_geneset_lps141_LPS141_220505_MDM2_SInorm.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",1800)),dec=".") 

rownames <- myData[,4]

mdm2.matrix <- myData[,-c(1:6)]
rownames(mdm2.matrix) <- rownames


x <- 1:300

y <- apply(mdm2.matrix[-4,c(1:300)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(301:600)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(601:900)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(901:1200)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(1201:1500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(1501:1800)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_mdm2Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='palevioletred2', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='P53 regions - MDM2 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("palevioletred2",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#B1210F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#B1210F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_mdm2Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_mdm2Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_mdm2Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()


###############################################################################
## P53 binding on MDM2 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/mdm2_merged_geneset_lps141_220505_p53_SInorm.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",1800)),dec=".") 

rownames <- myData[,4]

p53.matrix <- myData[,-c(1:6)]
rownames(p53.matrix) <- rownames


x <- 1:300

y <- apply(p53.matrix[,c(1:300)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(301:600)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(601:900)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(901:1200)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(1201:1500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(1501:1800)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_p53Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#42BB6E', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='MDM2 regions - P53 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#42BB6E",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#00741F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#00741F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_p53Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_p53Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_p53Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2200), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()



###############################################################################
## MDM2 binding on MDM2 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/mdm2_merged_geneset_lps141_LPS141_220505_MDM2_SInorm.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",1800)),dec=".") 

rownames <- myData[,4]

mdm2.matrix <- myData[,-c(1:6)]
rownames(mdm2.matrix) <- rownames


x <- 1:300

y <- apply(mdm2.matrix[,c(1:300)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(301:600)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(601:900)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(901:1200)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(1201:1500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(1501:1800)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_mdm2Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='palevioletred2', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='P53 regions - MDM2 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("palevioletred2",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#B1210F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#B1210F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_mdm2Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_mdm2Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/MDM2regions_mdm2Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()




###############################################################################
## YY1 binding on P53 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/p53_merged_geneset_lps141_LPS853_220330_YY1_SInorm.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",900)),dec=".") 

rownames <- myData[,4]

yy1.matrix <- myData[,-c(1:6)]
rownames(yy1.matrix) <- rownames


x <- 1:300

y <- apply(yy1.matrix[-4,c(1:300)],2,mean)
NT = smooth.spline(x, y, spar=0.6)

y <- apply(yy1.matrix[-4,c(301:600)],2,mean)
T2h = smooth.spline(x, y, spar=0.6)

y <- apply(yy1.matrix[-4,c(601:900)],2,mean)
T24h = smooth.spline(x, y, spar=0.6)



library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_yy1Binding.pdf', width = 5, height = 7)
plot(predict(NT)$y, col='palevioletred2', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='P53 regions - YY1 binding')



lines(predict(T2h)$y, col='#F15A29', lwd=2, type='l')



lines(predict(T24h)$y, col='#B1210F', lwd=2, type='l')
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_yy1Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_yy1Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/P53regions_yy1Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(0,45), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()





################################################################################
## Generate the metagene plots using the tag SI norm without the RPGC norm. 
## correct by background levels



###############################################################################
## P53 binding on P53 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/metagene_normalization_iterations/p53_merged_geneset_lps141_220505_p53_SInorm_tagnoNorm_10kb.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",12600)),dec=".") 

rownames <- myData[,4]

p53.matrix <- myData[,-c(1:6)]
rownames(p53.matrix) <- rownames


x <- 1:2100

y <- apply(p53.matrix[,c(1:2100)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(2101:4200)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(4201:6300)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(6301:8400)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(8401:10500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(10501:12600)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_p53Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#42BB6E', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='P53 regions - P53 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#42BB6E",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#00741F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#00741F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_p53Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_p53Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_p53Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()



###############################################################################
## MDM2 binding on P53 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/metagene_normalization_iterations/p53_merged_geneset_lps141_LPS141_220505_MDM2_SInorm_tagnoNorm_10kb.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",12600)),dec=".") 

rownames <- myData[,4]

mdm2.matrix <- myData[,-c(1:6)]
rownames(mdm2.matrix) <- rownames


x <- 1:2100

y <- apply(mdm2.matrix[-4,c(1:2100)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(2101:4200)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(4201:6300)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(6301:8400)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(8401:10500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[-4,c(10501:12600)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_mdm2Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='palevioletred2', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='P53 regions - MDM2 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("palevioletred2",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#B1210F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#B1210F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_mdm2Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_mdm2Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/P53regions_mdm2Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()


###############################################################################
## P53 binding on MDM2 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/metagene_normalization_iterations/mdm2_merged_geneset_lps141_220505_p53_SInorm_tagnoNorm_10kb.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",12600)),dec=".") 

rownames <- myData[,4]

p53.matrix <- myData[,-c(1:6)]
rownames(p53.matrix) <- rownames


x <- 1:2100

y <- apply(p53.matrix[,c(1:2100)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(2101:4200)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(4201:6300)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(6301:8400)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(8401:10500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(p53.matrix[,c(10501:12600)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_p53Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#42BB6E', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='MDM2 regions - P53 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#42BB6E",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#00741F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#00741F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_p53Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_p53Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_p53Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#009052', lwd=2, type='l', ylim=c(0,2500), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#009052",0.2),border=F)
dev.off()



###############################################################################
## MDM2 binding on MDM2 regions

myData <- read.csv2('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/metagene_normalization_iterations/mdm2_merged_geneset_lps141_LPS141_220505_MDM2_SInorm_tagnoNorm_10kb.gz', 
                    sep='\t', header=F, skip=2, colClasses = c(rep("character",6), rep("numeric",12600)),dec=".") 

rownames <- myData[,4]

mdm2.matrix <- myData[,-c(1:6)]
rownames(mdm2.matrix) <- rownames


x <- 1:2100

y <- apply(mdm2.matrix[,c(1:2100)],2,mean)
NT_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(2101:4200)],2,mean)
NT_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(4201:6300)],2,mean)
T2h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(6301:8400)],2,mean)
T2h_2 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(8401:10500)],2,mean)
T24h_1 = smooth.spline(x, y, spar=0.6)

y <- apply(mdm2.matrix[,c(10501:12600)],2,mean)
T24h_2 = smooth.spline(x, y, spar=0.6)


library(scales)

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_mdm2Binding.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='palevioletred2', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='MDM2 regions - MDM2 binding')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("palevioletred2",0.2),border=F)


lines(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)


lines(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#B1210F', lwd=2, type='l')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#B1210F",0.2),border=F)
dev.off()


pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_mdm2Binding_NT.pdf', width = 5, height = 7)
plot(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='NT')
polygon(c(x, rev(x)), c(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)-apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd), 
                        rev(apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, mean)+apply(cbind(predict(NT_1)$y,predict(NT_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_mdm2Binding_2hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='2hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)-apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, mean)+apply(cbind(predict(T2h_1)$y,predict(T2h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()

pdf('/Volumes/seq_epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/plots/tagNorm/MDM2regions_mdm2Binding_24hr.pdf', width = 5, height = 7)
plot(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean), col='#F15A29', lwd=2, type='l', ylim=c(7,50), frame.plot=F, xlab='Region', ylab='Intensity', main='24hr')
polygon(c(x, rev(x)), c(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)-apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd), 
                        rev(apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, mean)+apply(cbind(predict(T24h_1)$y,predict(T24h_2)$y), 1, sd))),
        col = alpha("#F15A29",0.2),border=F)
dev.off()
