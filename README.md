# ABE-APOBEC: RBP and Base Editor Analysis Pipeline

This repository supports the analysis of RNA-binding protein (RBP) effects on Base Editor activity.  
The `bactools` module processes strand-specific VCF files generated from RNA-seq or targeted sequencing data to quantify editing frequencies across samples.

---
## ✨ Quality Control
## ✒ bowtie2 alignment
## 🧬 bcftools mpileup
```wsl (Ubuntu 24.04.1 LTS)
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

res_dir="$HOME/myscratch/ProQ_STAR_results"
genomeFa="$HOME/reference/kpn/genome/kp.fa"
bed="$HOME/reference/kpn/genome/kp.bed"

# mpileup + 过滤（保持你原始阈值逻辑：任一样本 ALT 深度>2 且任一样本DP>20）
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

## 🧪 `vcf_process.py`, `vcf_process.R`: VCF Processing Tool

Processes strand-separated VCF files (e.g., from `mpileup` + strand filtering) into a unified, sample-annotated tabular format for downstream analysis. (The results are same when you use R or python to process the VCF files.)

### 📥 Input
- VCF files (optionally gzipped) named like:
  - `RBP_1_Aligned.out.FWD.vcf`
  - `WT_2_Aligned.out.REV.vcf.gz`
- Files should reside in a single input directory.
- Sample names are auto-detected from VCF headers unless explicitly provided.

### 📤 Output
- A single vcf file (tab-separated) with columns:
  - Genomic coordinates (`CHROM`, `POS`, `REF`, `ALT`)
  - Strand-specific editing counts per sample

---

### ▶️ Usage

```bash or powershell
python vcf_process.py [STRAND] \
  --input <INPUT_VCF_DIR> \
  --output <OUTPUT_DIR> \
  [--samples SAMPLE1 SAMPLE2 ...]

Rscript vcf_process.R [STRAND] \
  --input <INPUT_VCF_DIR> \
  --output <OUTPUT_DIR> \
  [--samples SAMPLE1 SAMPLE2 ...]
 
