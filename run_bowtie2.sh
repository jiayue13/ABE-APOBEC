#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

indir=/mnt/e/Kpn_data/seq_data/ProQ-ABE-fastq
outdir=~/myscratch/ProQ_STAR_results
indexloc=~/reference/kpn/genome/bowtie2_index/kp

mkdir -p "$outdir"
#---------------------------------------------------------------#
# Decompression: decompress only if .gz exists; skip if already decompressed
gz_list=( "$indir"/*.fastq.gz )
if [ ${#gz_list[@]} -gt 0 ]; then
  echo "[INFO] Decompressing *.fastq.gz in $indir ..."
  for f in "${gz_list[@]}"; do
    gunzip "$f"
  done
else
  echo "[INFO] No .fastq.gz files. Assuming FASTQ already uncompressed."
fi

#---------------------------------------------------------------#
for r1 in "$indir"/*.R1.raw.fastq; do
    sample=$(basename "$r1" .R1.raw.fastq)
    r2="$indir/${sample}.R2.raw.fastq"
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
