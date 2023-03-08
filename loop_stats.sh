#!/bin/bash -l

use .bedtools-2.29.0

file=$1
score=$2

cd /seq/epiprod02/sgaldon/loop_annotation

echo $file
echo
echo 
echo Using filter score of $score

EP=`wc -l < $file"/"$file"_E_P_loops.bedpe"`
PP=`wc -l < $file"/"$file"_P_P_loops.bedpe"`
EE=`wc -l < $file"/"$file"_E_E_loops.bedpe"`

EP=`awk -F"\t" '$8>'$score $file"/"$file"_E_P_loops.bedpe" | wc -l`
PP=`awk -F"\t" '$8>'$score $file"/"$file"_P_P_loops.bedpe" | wc -l`
EE=`awk -F"\t" '$8>'$score $file"/"$file"_E_E_loops.bedpe" | wc -l`

total=$((EP+PP+EE))

echo $total
#echo Number of E-P loops: $EP '('$((EP*100/total))'%)'
echo Number of E-P loops: $EP $(bc <<< "scale=3 ; $EP * 100 / $total")
#echo Number of E-E loops: $EE '('$((EE*100/total))'%)'
echo Number of E-E loops: $EE $(bc <<< "scale=3 ; $EE * 100 / $total")
#echo Number of P-P loops: $PP '('$((PP*100/total))'%)'
echo Number of P-P loops: $PP $(bc <<< "scale=3 ; $PP * 100 / $total")

exit 0
## Stats on TFs in the EP loops

mdm2_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
p53_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
cjun_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
runx_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`

echo MDM2 in E-P loops: $mdm2_EP '('$((mdm2_EP*100/EP))'%)'
echo p53 in E-P loops: $p53_EP '('$((p53_EP*100/EP))'%)'
echo cJun in E-P loops: $cjun_EP '('$((cjun_EP*100/EP))'%)'
echo RUNX in E-P loops: $runx_EP '('$((runx_EP*100/EP))'%)'


# MDM2-bound vs -unbound genes

mdm2_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $mdm2_EP
mdm2_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $mdm2_EE
mdm2_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $mdm2_PP
total=$((mdm2_EP+mdm2_PP+mdm2_EE))

echo Number of MDM2 E-P loops: $(bc <<< "scale=3 ; $mdm2_EP * 100 / $total")
echo Number of MDM2 E-E loops: $(bc <<< "scale=3 ; $mdm2_EE * 100 / $total")
echo Number of MDM2 P-P loops: $(bc <<< "scale=3 ; $mdm2_PP * 100 / $total")

## MDM2 unbound 

mdm2_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type neither | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $mdm2_EP
mdm2_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type neither | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $mdm2_EE
mdm2_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type neither | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $mdm2_PP
total=$((mdm2_EP+mdm2_PP+mdm2_EE))

echo Number of nonMDM2 E-P loops: $(bc <<< "scale=3 ; $mdm2_EP * 100 / $total")
echo Number of nonMDM2 E-E loops: $(bc <<< "scale=3 ; $mdm2_EE * 100 / $total")
echo Number of nonMDM2 P-P loops: $(bc <<< "scale=3 ; $mdm2_PP * 100 / $total")

# p53-bound vs -unbound genes

p53_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $p53_EP
p53_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $p53_EE
p53_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $p53_PP
total=$((p53_EP+p53_PP+p53_EE))

echo Number of p53 E-P loops: $(bc <<< "scale=3 ; $p53_EP * 100 / $total")
echo Number of p53 E-E loops: $(bc <<< "scale=3 ; $p53_EE * 100 / $total")
echo Number of p53 P-P loops: $(bc <<< "scale=3 ; $p53_PP * 100 / $total")

## p53 unbound 

p53_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type neither | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $p53_EP
p53_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_p53.bed -type neither | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $p53_EE
p53_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type neither | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
echo $p53_PP
total=$((p53_EP+p53_PP+p53_EE))

echo Number of nonp53 E-P loops: $(bc <<< "scale=3 ; $p53_EP * 100 / $total")
echo Number of nonp53 E-E loops: $(bc <<< "scale=3 ; $p53_EE * 100 / $total")
echo Number of nonp53 P-P loops: $(bc <<< "scale=3 ; $p53_PP * 100 / $total")


###########################################
## Percentage of mdm2 present at loops
echo MDM2 stats
mdm2_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 9,10,11 | uniq | bedtools intersect -wa -a TF_chips/lps141_mdm2.bed -b stdin | cut -f 1,2,3 | bedtools sort | uniq | wc -l`
mdm2_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 9,10,11 | uniq | bedtools intersect -wa -a TF_chips/lps141_mdm2.bed -b stdin | cut -f 1,2,3 | bedtools sort | uniq | wc -l`
mdm2_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 9,10,11 | uniq | bedtools intersect -wa -a TF_chips/lps141_mdm2.bed -b stdin | cut -f 1,2,3 | bedtools sort | uniq | wc -l`

echo $mdm2_EP
echo $mdm2_PP
echo $mdm2_EE
wc -l TF_chips/lps141_mdm2.bed

###########################################
## Percentage of p53 present at loops
echo p53 stats
p53_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 9,10,11 | uniq | bedtools intersect -wa -a TF_chips/lps141_p53.bed -b stdin | cut -f 1,2,3 | bedtools sort | uniq | wc -l`
p53_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 9,10,11 | uniq | bedtools intersect -wa -a TF_chips/lps141_p53.bed -b stdin | cut -f 1,2,3 | bedtools sort | uniq | wc -l`
p53_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 9,10,11 | uniq | bedtools intersect -wa -a TF_chips/lps141_p53.bed -b stdin | cut -f 1,2,3 | bedtools sort | uniq | wc -l`

echo $p53_EP
echo $p53_PP
echo $p53_EE
wc -l TF_chips/lps141_p53.bed

#############################################
## Percentage of MDM2 binding in both sides

mdm2_EP_total=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
mdm2_EE_total=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
mdm2_PP_total=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`

mdm2_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_mdm2.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
mdm2_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_mdm2.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
mdm2_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_mdm2.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_mdm2.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`


echo MDM2 in both ends of E-E loops: $(bc <<< "scale=3 ; $mdm2_EE * 100 / ($mdm2_EE_total)")
echo MDM2 in both ends of E-P loops: $(bc <<< "scale=3 ; $mdm2_EP * 100 / ($mdm2_EP_total)")
echo MDM2 in both ends of P-P loops: $(bc <<< "scale=3 ; $mdm2_PP * 100 / ($mdm2_PP_total)")

#############################################
## Percentage of p53 binding in both sides

p53_EP_total=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
p53_EE_total=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
p53_PP_total=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`

p53_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_p53.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
p53_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_p53.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
p53_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_p53.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_p53.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`


echo p53 in both ends of E-E loops: $(bc <<< "scale=3 ; $p53_EE * 100 / $p53_EE_total")
echo p53 in both ends of E-P loops: $(bc <<< "scale=3 ; $p53_EP * 100 / $p53_EP_total")
echo p53 in both ends of P-P loops: $(bc <<< "scale=3 ; $p53_PP * 100 / $p53_PP_total")


#############################################
## Percentage of cjun binding in both sides

cjun_EP_total=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
cjun_EE_total=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
cjun_PP_total=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`

cjun_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_cjun.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
cjun_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_cjun.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
cjun_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_cjun.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_cjun.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`


echo cjun in both ends of E-E loops: $(bc <<< "scale=3 ; $cjun_EE * 100 / $cjun_EE_total")
echo cjun in both ends of E-P loops: $(bc <<< "scale=3 ; $cjun_EP * 100 / $cjun_EP_total")
echo cjun in both ends of P-P loops: $(bc <<< "scale=3 ; $cjun_PP * 100 / $cjun_PP_total")


#############################################
## Percentage of runx binding in both sides

runx_EP_total=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
runx_EE_total=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
runx_PP_total=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`

runx_EP=`bedtools pairtobed -a $file"/"$file"_E_P_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_runx.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
runx_EE=`bedtools pairtobed -a $file"/"$file"_E_E_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_runx.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`
runx_PP=`bedtools pairtobed -a $file"/"$file"_P_P_loops.bedpe" -b TF_chips/lps141_runx.bed -type either | awk -F"\t" '$8>'$score | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | bedtools pairtobed -a stdin -b TF_chips/lps141_runx.bed -type both | cut -f 1,2,3,4,5,6 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n | uniq | wc -l`


echo runx in both ends of E-E loops: $(bc <<< "scale=3 ; $runx_EE * 100 / $runx_EE_total")
echo runx in both ends of E-P loops: $(bc <<< "scale=3 ; $runx_EP * 100 / $runx_EP_total")
echo runx in both ends of P-P loops: $(bc <<< "scale=3 ; $runx_PP * 100 / $runx_PP_total")

