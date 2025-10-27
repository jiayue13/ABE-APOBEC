# ABE-APOBEC: RBP and Base Editor Analysis Pipeline

This repository supports the analysis of RNA-binding protein (RBP) effects on Base Editor activity.  
The `bactools` module processes strand-specific VCF files generated from RNA-seq or targeted sequencing data to quantify editing frequencies across samples.

---
## ✨ Quality Control

## ✒ bowtie2 alignment
- I put all fastq files in a Folder.
  - `samtools version 1.19.2`
  - `bowtie2 version 2.5.4`
  - `gff2bed version 2.4.42`
```
# make genome fasta index
mkdir -p ~/reference/kpn/genome/bowtie2_index
bowtie2-build ~/reference/kpn/genome/kp.fa ~/reference/kpn/genome/bowtie2_index/kp

# make genome bed
gff2bed < kp.gff > kp.bed
```
```
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

indir=/mnt/e/Kpn_data/seq_data/ProQ-ABE-fastq
outdir=~/myscratch/ProQ_STAR_results
indexloc=~/reference/kpn/genome/bowtie2_index/kp

mkdir -p "$outdir"

# unzip
gunzip -dk ${indir}/*.fastq.gz

for r1 in "$indir"/*.R1.raw.fastq.gz; do
    sample=$(basename "$r1" .R1.raw.fastq.gz)
    r2="${r1/.R1.raw.fastq.gz/.R2.raw.fastq.gz}"
    [[ -f "$r2" ]] || { echo "[WARN] Missing R2 for $sample"; continue; }

    echo "[INFO] Processing sample: $sample"

    sam="${outdir}/${sample}.sam"
    bam="${outdir}/${sample}_Aligned.out.bam"
    sortbam="${outdir}/${sample}_Aligned.out.sorted.bam"

    # align
    bowtie2 -x "$indexloc" -p 16 -1 "$r1" -2 "$r2" -S "$sam"

    # convert, sort, index
    samtools view -@ 16 -bhS -o "$bam" "$sam"
    samtools sort -@ 16 -o "$sortbam" "$bam"
    samtools index -@ 16 "$sortbam"

    # split by strand
    samtools view -@ 16 -b -f 99  "$sortbam" > "${outdir}/${sample}_Aligned.out.R1F.bam"
    samtools view -@ 16 -b -f 147 "$sortbam" > "${outdir}/${sample}_Aligned.out.R2R.bam"
    samtools view -@ 16 -b -f 83  "$sortbam" > "${outdir}/${sample}_Aligned.out.R1R.bam"
    samtools view -@ 16 -b -f 163 "$sortbam" > "${outdir}/${sample}_Aligned.out.R2F.bam"

    samtools merge -@ 16 "${outdir}/${sample}_Aligned.out.FWD.bam" \
        "${outdir}/${sample}_Aligned.out.R1F.bam" \
        "${outdir}/${sample}_Aligned.out.R2R.bam"

    samtools merge -@ 16 "${outdir}/${sample}_Aligned.out.REV.bam" \
        "${outdir}/${sample}_Aligned.out.R1R.bam" \
        "${outdir}/${sample}_Aligned.out.R2F.bam"

    samtools index -@ 16 "${outdir}/${sample}_Aligned.out.FWD.bam"
    samtools index -@ 16 "${outdir}/${sample}_Aligned.out.REV.bam"

    echo "[DONE] $sample"
done

```
## 🧬 bcftools mpileup
- **Set workdir and bash**
  - `bcftools 1.19`
  - `Using htslib 1.19`
```
#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

res_dir="$HOME/myscratch/ProQ_STAR_results" # result dir
genomeFa="$HOME/reference/kpn/genome/kp.fa" # genome index dir
bed="$HOME/reference/kpn/genome/kp.bed" # genome bed dir (bedtools gtf2bed < genes.gtf > genes.bed)
```
- **mpileup + filter** (Keep threshold logic: **any sample with ALT depth > 2** and **any sample with DP > 20**)
```
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
```
- The colnames of vcf files count on the bam you submitted.
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
  - `Genomic coordinates (`CHROM`, `POS`, `REF`, `ALT`)`
  - `Strand-specific editing counts per sample`

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
 
