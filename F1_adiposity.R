# IMPORT data


x <- read.csv("../../../../../../Liposarcoma/FINAL_COHORT_SAMP_ID.csv", row.names = 1)

y <- read.delim("../../../../../../Liposarcoma/Computational analyses/Chadi_analysis/RNA-Seq_Final_Matrix/Cohort1/Liposarcoma_Cohort1_TPM_Matrix.txt", sep = "\t", header = TRUE, row.names = 1)

z <- read.delim("../../../../../../Liposarcoma/Computational analyses/Chadi_analysis/RNA-Seq_Final_Matrix/Cohort2/Liposarcoma_Cohort2_TPM_Matrix.txt", sep = "\t", header = TRUE, row.names = 1)

w <- read.delim("../../../../../../Liposarcoma/Computational analyses/Chadi_analysis/RNA-Seq_Cell_Lines/Cell_Lines_TPM_Matrix.txt", sep = "\t", header = TRUE, row.names = 1)

lps <- cbind(y,z)
cell_lines <- w[,c("LPS141_1","LPS853_2","T449_2","T778_1")]

IDs <- row.names(x)

exp <- lps[,colnames(lps) %in% IDs]

drops <- c("WD14.1", "N27", "N6")
exp <- exp[,!(names(exp) %in% drops)]

exp <- merge(exp, cell_lines, by = 'row.names')
row.names(exp) <- exp$Row.names
exp <- exp[,-1]

col_order <- c("N34", "N15", "N22", "N33", "N31", "N30", "N14", "N19", "N10", "N12", "WD36", "WD34", "WD24", "WD33", "WD12", "WD15", "WD23", "WD22", "WD4", "WD10", "WD14", "WD19", "WD13", "WD28", "WD30", "WD5", "DD30", "DD20", "DD10", "DD2", "DD1", "DD15", "DD35", "DD36", "DD31", "DD_SS", "LPS141_1","LPS853_2","T449_2","T778_1")

exp <- exp[,col_order]


#write.table(exp, "220804_exp_final_tumor_cohort_and_cell_lines.txt", sep = "\t")


### Adiposity scoring of LPS cohort - tumors and cell lines

library("GSVA")

Adipo <- scan("Adipo30.txt", what="")

GeneLists <- list(Adipo)

M <- as.matrix(exp)

G <- gsva(M, GeneLists)

row.names(G) = c("adipo_score")

G <- G[,col_order]

#write.table(G, "adiposcore_tumor_and_cell_line.txt", sep = "\t")


####FIGURE 1 plots

histology <- c(rep("A_N",10), rep("B_WD",16), rep("C_DD",10))

data <- read.delim("adiposcore_MDM2exp_tumor_and_cell_line.txt", header = TRUE, sep = "\t")

tumors <- data[which(data$histology == "DD" | data$histology == "WD"),]
tumors_NF <- data[which(data$histology == "NF"),]
tumors_WD <- data[which(data$histology == "WD"),] 
tumors_DD <- data[which(data$histology == "DD"),] 

save(tumors_NF,tumors_WD,tumors_DD, file = "lps_expression_by_histology.RData")

#pdf("WD_DD_tumors_adiposity_v_MDM2exp.pdf")
plot(tumors$adiposity, tumors$MDM2, col = ifelse(tumors$histology == "DD", "blue", "red"), cex = 1.5)
abline(lm(tumors_WD$MDM2 ~ tumors_WD$adiposity, col = "red"))
abline(lm(tumors_DD$MDM2 ~ tumors_DD$adiposity, col = "blue"))

#pdf("WD_tumors_adiposity_v_MDM2exp.pdf")
plot(tumors_WD$adiposity, tumors_WD$MDM2, col =  "red", abline(lm(tumors_WD$MDM2 ~ tumors_WD$adiposity)), xlim = c(-1,1), cex = 1.5)

#pdf("DD_tumors_adiposity_v_MDM2exp.pdf")
plot(tumors_DD$adiposity, tumors_DD$MDM2, col = "blue", abline(lm(tumors_DD$MDM2 ~ tumors_DD$adiposity)), xlim = c(-1,1), cex = 1.5)


cor(tumors$adiposity, tumors$MDM2)
cor(tumors_WD$adiposity, tumors_WD$MDM2)
cor(tumors_DD$adiposity, tumors_DD$MDM2)


#pdf("adiposity_score_bar_plot.pdf")
barplot(c(tumors_NF$adiposity, tumors_WD$adiposity, tumors_DD$adiposity))


##pdf("MDM2_boxplot_tumors_cell_lines.pdf")
##boxplot(tumors_NF$MDM2, tumors_WD$MDM2, tumors_DD$MDM2)


#pdf("MDM2_boxplot_tumors_cell_lines.pdf")
boxplot(MDM2~histology, data = data, col = "white")
stripchart(MDM2~histology,
           data = data, # Data
           method = "jitter", # Random noise
           pch = 19,          # Pch symbols
           col = c("#BEBEBE", "#D383B7", "#4698D3", "black"), # Color of the symbol
           vertical = TRUE,   # Vertical mode
           add = TRUE)

t.test(tumors_WD$MDM2, tumors_DD$MDM2, alternative = 'less')

#pdf("PPARG_CEBPA_expression_supplemental_tumors_with_histology.pdf")
barplot(c(exp[c("PPARG","CEBPA"), "WD12"], exp[c("PPARG","CEBPA"), "WD23"], exp[c("PPARG","CEBPA"), "WD4"], exp[c("PPARG","CEBPA"), "DD15"]))



WD_exp <- exp[,11:26]
WD_exp_t <- as.data.frame(t(WD_exp))
pdf("WD_tumors_MDM2vPPARG_exp.pdf")
plot(WD_exp_t$MDM2, WD_exp_t$PPARG, col =  "red", abline(lm(WD_exp_t$PPARG ~ WD_exp_t$MDM2)), cex = 1.5)
pdf("WD_tumors_MDM2vCEBPA_exp.pdf")
plot(WD_exp_t$MDM2, WD_exp_t$CEBPA, col =  "red", abline(lm(WD_exp_t$CEBPA ~ WD_exp_t$MDM2)), cex = 1.5)

cor(WD_exp_t$PPARG, WD_exp_t$MDM2, method = 'spearman')
cor(WD_exp_t$CEBPA, WD_exp_t$MDM2, method = 'spearman')


library(car)
pdf("histogram_MDM2_exp_NF_DD_WD.pdf")
densityPlot(log2(data$MDM2) ~ as.factor(data$histology)) 

t.test(tumors_DD$MDM2, tumors_WD$MDM2, alternative = "greater")$p.value
