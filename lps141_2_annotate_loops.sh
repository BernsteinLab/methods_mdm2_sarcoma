#!/bin/bash -l

###################################
## Generate EE EP PP tsv files


use .bedtools-2.29.0
use .r-3.6.0

sample=$1
score=$2

##################################
## E-P loops

cd /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/loop_annotation/LPS141_2/

mkdir EP_loops

bedtools pairtobed -a $sample"_E_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6,7,8 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq > totalLoops.tmp

#Get enhancer 
enhan_first_half=`bedtools pairtobed -a $sample"_E_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3 | \
bedtools intersect -v -wa -a stdin -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed | bedtools sort | uniq`
enhan_second_half=`bedtools pairtobed -a $sample"_E_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type either | awk -F"\t" '$8>'$score | cut -f 4,5,6 | \
bedtools intersect -v -wa -a stdin -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed | bedtools sort | uniq`

#Create enhancers file
printf "$enhan_first_half\n$enhan_second_half" | bedtools sort | uniq > enhancers_totalLoops.tmp

# Get promoters in loops with loop coordinates
prom_first_half=`bedtools pairtobed -a $sample"_E_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3 | \
bedtools intersect -wa -a stdin -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed | bedtools sort | uniq`
prom_second_half=`bedtools pairtobed -a $sample"_E_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type either | awk -F"\t" '$8>'$score | cut -f 4,5,6 | \
bedtools intersect -wa -a stdin -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed | bedtools sort | uniq`

#Create enhancers file
printf "$prom_first_half\n$prom_second_half" | bedtools sort | uniq > promoters_loopAnnot_totalLoops.tmp

# MDM2-bound promoters
bedtools intersect -wa -wb -a promoters_loopAnnot_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_mdm2.bed > promoter_totalLoops_mdm2Bound.tmp

# MDM2-bound enhancers
bedtools intersect -wa -wb -a enhancers_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_mdm2.bed > enhancers_totalLoops_mdm2Bound.tmp

# p53-bound promoters
bedtools intersect -wa -wb -a promoters_loopAnnot_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_p53.bed > promoter_totalLoops_p53Bound.tmp

# p53-bound enhancers
bedtools intersect -wa -wb -a enhancers_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_p53.bed > enhancers_totalLoops_p53Bound.tmp

# cjun-bound promoters
bedtools intersect -wa -wb -a promoters_loopAnnot_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_cjun.bed > promoter_totalLoops_cjunBound.tmp

# cjun-bound enhancers
bedtools intersect -wa -wb -a enhancers_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_cjun.bed > enhancers_totalLoops_cjunBound.tmp

# runx-bound promoters
bedtools intersect -wa -wb -a promoters_loopAnnot_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_runx.bed > promoter_totalLoops_runxBound.tmp

# runx-bound enhancers
bedtools intersect -wa -wb -a enhancers_totalLoops.tmp -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_runx.bed > enhancers_totalLoops_runxBound.tmp
#exit 0
##################################################
## Annotate promoters to genes by proximity
echo ------ Annotate Genes ------

mkdir promoter_annotation

# All promoters
annotatePeaks.pl promoters_loopAnnot_totalLoops.tmp hg38 > "promoter_annotation/EP_"$score"_promoter_annotation.txt"

#MDM2 bound
annotatePeaks.pl promoter_totalLoops_mdm2Bound.tmp hg38 > "promoter_annotation/EP_"$score"_mdm2Bound_promoter_annotation.txt"

#p53 bound
annotatePeaks.pl promoter_totalLoops_p53Bound.tmp hg38 > "promoter_annotation/EP_"$score"_p53Bound_promoter_annotation.txt"


echo ------ Create annotation file ------



##################################
## P-P loops

echo ------ Annotate P-P loops ------

mkdir PP_loops

bedtools pairtobed -a $sample"_P_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_peaks_all_noFilter.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6,7,8 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq > "PP_loops/PP_totalLoops.tmp"


prom_first_half=`bedtools pairtobed -a $sample"_P_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_peaks_all_noFilter.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3 | \
bedtools intersect -wa -a stdin -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed | bedtools sort | uniq`
prom_second_half=`bedtools pairtobed -a $sample"_P_P_loops.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_peaks_all_noFilter.bed -type either | awk -F"\t" '$8>'$score | cut -f 4,5,6 | \
bedtools intersect -wa -a stdin -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed | bedtools sort | uniq`

#Create enhancers file
printf "$prom_first_half\n$prom_second_half" | bedtools sort | uniq > "PP_loops/PP_promoters_loopAnnot_totalLoops.tmp"

# MDM2-bound promoters
bedtools intersect -wa -wb -a "PP_loops/PP_promoters_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_mdm2.bed > "PP_loops/PP_promoter_totalLoops_mdm2Bound.tmp"
# p53-bound promoters
bedtools intersect -wa -wb -a "PP_loops/PP_promoters_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_p53.bed > "PP_loops/PP_promoter_totalLoops_p53Bound.tmp"
# cJun-bound promoters
bedtools intersect -wa -wb -a "PP_loops/PP_promoters_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_cjun.bed > "PP_loops/PP_promoter_totalLoops_cjunBound.tmp"
# RUNX-bound promoters
bedtools intersect -wa -wb -a "PP_loops/PP_promoters_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_runx.bed > "PP_loops/PP_promoter_totalLoops_runxBound.tmp"

##################################################
## Annotate promoters to genes by proximity
echo ------ Annotate Genes ------

mkdir PP_loops/promoter_annotation

# All promoters
annotatePeaks.pl "PP_loops/PP_promoters_loopAnnot_totalLoops.tmp" hg38 > "PP_loops/promoter_annotation/PP_"$score"_promoter_annotation.txt"


##################################
## E-E loops

mkdir EE_loops

awk -F"\t" '$8>'$score $sample"_E_E_loops.bedpe" | cut -f 1,2,3,4,5,6,7,8 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq > "EE_loops/EE_totalLoops.tmp"

enh_first_half=`awk -F"\t" '$8>'$score $sample"_E_E_loops.bedpe" | cut -f 1,2,3 |  bedtools sort | uniq`
enh_second_half=`awk -F"\t" '$8>'$score $sample"_E_E_loops.bedpe" | cut -f 4,5,6 |  bedtools sort | uniq`

#Create enhancers file
printf "$enh_first_half\n$enh_second_half" | bedtools sort | uniq > "EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp"

# MDM2-bound enhancers
bedtools intersect -wa -wb -a "EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_mdm2.bed > "EE_loops/EE_enhancers_totalLoops_mdm2Bound.tmp"
# p53-bound enhancers
bedtools intersect -wa -wb -a "EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_p53.bed > "EE_loops/EE_enhancers_totalLoops_p53Bound.tmp"
# cJun-bound enhancers
bedtools intersect -wa -wb -a "EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_cjun.bed > "EE_loops/EE_enhancers_totalLoops_cjunBound.tmp"
# RUNX-bound enhancers
bedtools intersect -wa -wb -a "EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/lps141_runx.bed > "EE_loops/EE_enhancers_totalLoops_runxBound.tmp"

##################################################
## Annotate enhancers to genes by proximity
echo ------ Annotate Genes ------

mkdir EE_loops/enhancers_annotation

# All enhancers
annotatePeaks.pl "EE_loops/EE_enhancers_loopAnnot_totalLoops.tmp" hg38 > "EE_loops/enhancers_annotation/EE_"$score"_enhancers_annotation.txt"


cd ..

Rscript /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/lps141_loop_annotation_table.R $sample

cd $sample

rm *.tmp






