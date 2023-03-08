## Heatmaps figure 3A

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


## Metagenes Figure 3K and L


##MDM2

#bed_file=$1

#bed_path=$2

out=${bed_file%.*}"_LPS141_220505_MDM2_SInorm"

~/deepTools/bin/computeMatrix scale-regions -S /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_untreated_rep1_MDM2_S1_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_untreated_rep2_MDM2_S2_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_2h_rep1_MDM2_S3_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_2h_rep2_MDM2_S4_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_24h_rep1_MDM2_S5_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_24h_rep2_MDM2_S6_SInorm.bw \
-R "$bed_path"$bed_file -b 10000 -a 10000 \
--skipZeros --smartLabels -o "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.gz";

~/deepTools/bin/plotProfile -m "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.gz" \
-o "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.pdf" \
--outFileSortedRegions "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.bed" \
--perGroup \
--colors black black orange orange red red \
--plotTitle "" --samplesLabel "NT" "NT" "2h" "2h" "24h" "24h" \
--refPointLabel "TSS" \
-T "MDM2 binding LPS141 TC SInorm" \
-z "";




## P53

#bed_file=$1

#bed_path=$2

out=${bed_file%.*}"_220505_p53_SInorm"

~/deepTools/bin/computeMatrix scale-regions -S /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_untreated_rep1_P53_S7_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_untreated_rep2_P53_S8_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_2h_rep1_P53_S9_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_2h_rep2_P53_S10_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_24h_rep1_P53_S11_SInorm.bw \
/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bigWig_files/lps141/220505/220505_LPS141_24h_rep2_P53_S12_SInorm.bw \
-R "$bed_path"$bed_file -b 10000 -a 10000 \
--skipZeros --smartLabels -o "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.gz";

~/deepTools/bin/plotProfile -m "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.gz" \
-o "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.pdf" \
--outFileSortedRegions "/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/scripts/metagenes/metagenes_HDM201/"$out"_10kb.bed" \
--perGroup \
--colors black black orange orange red red \
--plotTitle "" --samplesLabel "NT" "NT" "2h" "2h" "24h" "24h" \
--refPointLabel "TSS" \
-T "p53 binding LPS141 TC SInorm" \
-z "";

