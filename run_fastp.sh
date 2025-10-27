#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

indir=/mnt/e/Kpn_data/seq_data/ProQ-ABE-fastq
# sample list
SAMPLES=("ProQ_1" "ProQ_2" "ProQ_3" "WT_1" "WT_2" "WT_3")

# Fastp Quality Control
for SAMPLE in "${SAMPLES[@]}"; do
    echo "🧼 Running fastp for $SAMPLE..."
    fastp -i $indir/${SAMPLE}_R1.fastq.gz -I $FASTQ_DIR/${SAMPLE}_R2.fastq.gz \
          -o ${SAMPLE}.R1.raw.fastq.gz -O ${SAMPLE}.R1.raw.fastq.gz \
          -h ${SAMPLE}_fastp.html -j ${SAMPLE}_fastp.json \
          -q 20 -u 30 -n 5 -l 50
done
