# ABE-APOBEC: RBP and Base Editor Analysis Pipeline

This repository supports the analysis of RNA-binding protein (RBP) effects on Base Editor activity.  
The `bcftools` module processes strand-specific VCF files generated from RNA-seq or targeted sequencing data to quantify editing frequencies across samples.

---
## ✨ Quality Control
- Use `fastp` to perform quality control on original fastq. 
  - `fastp version 0.23.4`
- Run `run_fastp.sh` to perform quality control.

## ✒ bowtie2 alignment
- I put all fastq files in a Folder. First make genome index and genome bed.
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
- Run `run_bowtie2.sh` to align fastq and split chains.
## 🧬 bcftools mpileup
- **Set workdir and bash**
  - `bcftools 1.19`
  - `Using htslib 1.19`
- Run `run_bcftools.sh` to perform variants calling.
- **mpileup + filter** (Keep threshold logic: **any sample with ALT depth > 2** and **any sample with DP > 20**)
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
 
