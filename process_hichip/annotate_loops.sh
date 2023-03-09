#!/bin/bash -l
# First convert spaces to tabs in newBEDPE using sed -e 's/ /\t/g'
use .bedtools-2.29.0

file=$1
hichipperPath='/seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/hichipper/hichipper_down_2Samplesk27_out'
hichipperPath=$2

#cd /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/hiChip/hichipper/hichipper_out
mkdir ./$file

echo $file

bedtools pairtobed -a $hichipperPath"/new"$file".filt.intra.loop_counts.bedpe" \
-b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type either | cut -f 1,2,3,4,5,6,7,8 | bedtools pairtobed -a stdin \
-b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_peaks_all_noFilter.bed -type either | cut -f 1,2,3,4,5,6,7,8 | sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 | uniq > ./$file/$file"_E_P_loops.bedpe"

bedtools pairtobed -a $hichipperPath"/new"$file".filt.intra.loop_counts.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_peaks_all_noFilter.bed -type both | cut -f 1,2,3,4,5,6,7,8 | sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 | uniq > ./$file/$file"_E_E_loops.bedpe"

bedtools pairtobed -a $hichipperPath"/new"$file".filt.intra.loop_counts.bedpe" -b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/hg38_refGene_promoters_500bp.bed -type both | cut -f 1,2,3,4,5,6,7,8 | bedtools pairtobed -a stdin \
-b /seq/epiprod02/sgaldon/liposarcoma/chip_master_datasets/bed_files/H3k27ac_peaks_all_noFilter.bed -type either | cut -f 1,2,3,4,5,6,7,8 | sort -k1,1 -k2,2 -k3,3 -k4,4 -k5,5 -k6,6 | uniq > ./$file/$file"_P_P_loops.bedpe"

# # Prioritize P annotations over E annotations: Delete all redundant loops at Es
comm -23 ./$file/$file"_E_E_loops.bedpe" ./$file/$file"_P_P_loops.bedpe" > ./$file/$file"_E_E_loops.bedpe.tmp"
#bedtools pairtopair -a ./$file/$file"_E_E_loops.bedpe" -b ./$file/$file"_P_P_loops.bedpe" -type both
comm -23 ./$file/$file"_E_E_loops.bedpe.tmp" ./$file/$file"_E_P_loops.bedpe" > ./$file/$file"_E_E_loops.bedpe.tmp.tmp"
comm -23 ./$file/$file"_E_P_loops.bedpe" ./$file/$file"_P_P_loops.bedpe" > ./$file/$file"_E_P_loops.bedpe.tmp"

mv ./$file/$file"_E_E_loops.bedpe.tmp.tmp" ./$file/$file"_E_E_loops.bedpe"
rm ./$file/$file"_E_E_loops.bedpe.tmp"
mv ./$file/$file"_E_P_loops.bedpe.tmp" ./$file/$file"_E_P_loops.bedpe"