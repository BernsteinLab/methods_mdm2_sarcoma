~/deepTools/bin/computeMatrix scale-regions -S /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/NT/210223_LPS141_DMSO_cJun.bw /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps853/NT/210308_LPS853_DMSO_cJun.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/NT/210223_LPS141_DMSO_RUNX2.bw /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps853/NT/210308_LPS853_DMSO_RUNX2.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/NT/210223_LPS141_DMSO_p53.bw /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps853/NT/210308_LPS853_DMSO_p53.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/NT/210223_LPS141_DMSO_MDM2.bw /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps853/NT/210308_LPS853_DMSO_MDM2.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/NT/210127_LPS141_200910_201204_DMSO_201208_K27ac.bw /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps853/NT/210308_LPS853_DMSO_H3K27ac.bw \
-R  cjun_runx_enhancers_noAmplified_lps141.bed cjun_runx_promoters_noAmplified_lps141.bed p53_enhancers_noAmplified_lps141.bed p53_promoters_noAmplified_lps141.bed mdm2_enhancers_noAmplified_lps141.bed mdm2_promoters_noAmplified_lps141.bed \
-b 5000 -a 5000 --skipZeros --sortUsingSamples 7 8 5 6 --smartLabels --binSize 5 -o LPS141_853_enhancers_combined.gz;


~/deepTools/bin/plotHeatmap -m LPS141_853_enhancers_combined.gz -o LPS141_853_enhancers_combined.pdf \
--outFileSortedRegions LPS141_853_enhancers_combined.bed --zMax 20 --colorMap Blues Blues Purples Purples Greens Greens Reds Reds Greys Greys \
--startLabel . --endLabel . --whatToShow "heatmap and colorbar";
