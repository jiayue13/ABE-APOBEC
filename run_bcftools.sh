#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

res_dir="$HOME/myscratch/ProQ_STAR_results" # result dir
genomeFa="$HOME/reference/kpn/genome/kp.fa" # genome index dir
bed="$HOME/reference/kpn/genome/kp.bed" # genome bed dir (bedtools gtf2bed < genes.gtf > genes.bed)

# mpileup + filter** (Keep threshold logic: **any sample with ALT depth > 2** and **any sample with DP > 20**)

bcftools mpileup -f "$genomeFa" -R "$bed" -d 10000000 -I -a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR \
    ${res_dir}/ProQ_1_Aligned.out.FWD.bam \
    ${res_dir}/ProQ_2_Aligned.out.FWD.bam \
    ${res_dir}/ProQ_3_Aligned.out.FWD.bam\
    ${res_dir}/WT_1_Aligned.out.FWD.bam \
    ${res_dir}/WT_2_Aligned.out.FWD.bam \
    ${res_dir}/WT_3_Aligned.out.FWD.bam \
    | bcftools filter -i 'INFO/AD[1-]>2 & MAX(FORMAT/DP)>20' -O v - > /mnt/e/Kpn_data/seq_data/20251028/results/ProQ-rABE_FWD.vcf

bcftools mpileup -f "$genomeFa" -R "$bed" -d 10000000 -I -a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR \
    ${res_dir}/ProQ_1_Aligned.out.REV.bam \
    ${res_dir}/ProQ_2_Aligned.out.REV.bam \
    ${res_dir}/ProQ_3_Aligned.out.REV.bam\
    ${res_dir}/WT_1_Aligned.out.REV.bam \
    ${res_dir}/WT_2_Aligned.out.REV.bam \
    ${res_dir}/WT_3_Aligned.out.REV.bam \
    | bcftools filter -i 'INFO/AD[1-]>2 & MAX(FORMAT/DP)>20' -O v - > /mnt/e/Kpn_data/seq_data/20251028/results/ProQ-rABE_REV.vcf
	
echo "all samples REV finished!"
