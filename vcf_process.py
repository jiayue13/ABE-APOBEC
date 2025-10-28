#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Parse bcftools mpileup VCF, extract per-sample DP/AD and compute ratios.

Features:
- Preserves #CHROM column names; skips ## metadata
- Compatible with .vcf and .vcf.gz
- Automatic path normalization on Windows/WSL
- Use --input to specify the input directory; --output to specify the output directory
- Infers sample columns from headers (can also be specified with --samples)
- Prioritizes dynamic DP/AD location based on FORMAT; falls back to a fixed location if not found
"""

import os
import re
import gzip
import argparse
import platform
import pandas as pd
import numpy as np
from typing import List, Optional
import tempfile
import shutil

# ---------- Path Normalization (Windows -> WSL) ----------
def normalize_path(p: Optional[str]) -> Optional[str]:
    if p is None:
        return p
    p = p.strip().strip('"').strip("'")
    # Map "E:\path" or "e:/path" to "/mnt/e/path" under Linux/WSL
    if platform.system() == "Linux":
        m = re.match(r'^([a-zA-Z]):[\\/](.*)$', p)
        if m:
            drive = m.group(1).lower()
            rest = m.group(2).replace('\\', '/')
            return f"/mnt/{drive}/{rest}"
    return p.replace('\\', '/')

# ---------- Security Split ----------
def safe_split(s: Optional[str], sep: str, idx: int) -> Optional[str]:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return np.nan
    parts = str(s).split(sep)
    return parts[idx] if 0 <= idx < len(parts) else np.nan

def find_column_index(fmt: str, key: str) -> Optional[int]:
    """Find the 0-based position of the key in the FORMAT field"""
    if fmt is None or (isinstance(fmt, float) and pd.isna(fmt)):
        return None
    fields = str(fmt).split(':')
    for i, k in enumerate(fields):
        if k == key:
            return i
    return None

# ---------- Sample name cleaning (remove paths, remove .bam/.cram, and remove illegal characters -> _) ----------
def sanitize_sample(name: str) -> str:
    base = os.path.basename(name)
    base = re.sub(r"\.(bam|cram)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"[^\w.\-]+", "_", base)
    return base

# ---------- parameter ----------
def parse_args():
    ap = argparse.ArgumentParser(description="Fix strand and compute mutation ratios from VCF.")
    ap.add_argument("s", type=int, help="strand: 1=FWD/neg, 2=REV/pos")
    ap.add_argument("--input", required=True, help="Input VCF directory (including ProQ-rABE_FWD/REV.vcf[.gz])")
    ap.add_argument("--output", required=True, help="Output Directory")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="Optional: list of sample names; if not provided, they will be automatically inferred from the VCF header")
    return ap.parse_args()

# ---------- Main process ----------
def main():
    args = parse_args()

    # strand
    if args.s == 1:
        st, strand = "FWD", "neg"
    elif args.s == 2:
        st, strand = "REV", "pos"
    else:
        raise SystemExit(" s can be only 1 or 2.")

    input = normalize_path(args.input)
    output  = normalize_path(args.output)
    os.makedirs(output, exist_ok=True)

    # Input file (prefer .vcf, if not available, try .vcf.gz)
    vcf_path = os.path.join(input, f"ProQ-rABE_{st}.vcf")
    vcf_gz   = vcf_path + ".gz"
    if not os.path.exists(vcf_path):
        if os.path.exists(vcf_gz):
            vcf_path = vcf_gz
        else:
            raise FileNotFoundError(f"Input file not found: \n {vcf_path}\n or {vcf_gz}")

    print(f"[INFO] read file: {vcf_path}")

    # Read the header (#CHROM line) to get the column names and sample columns
    open_func = gzip.open if vcf_path.endswith(".gz") else open
    with open_func(vcf_path, "rt") as f:
        header_line = None
        for line in f:
            if line.startswith("#CHROM"):
                header_line = line.strip().lstrip("#").split("\t")
                break
    if header_line is None:
        raise ValueError("No #CHROM header row found in VCF.")

    # The sample columns start from the 10th column of the header.
    sample_cols_from_header = header_line[9:]

    # Read data: skip lines starting with # and use header_line as column names
    df = pd.read_csv(
        vcf_path,
        sep="\t",
        comment="#",
        names=header_line,
        dtype=str,
        na_filter=True,
        engine="python",
        compression="infer",
    )

    # Strand and ALT processing
    df["strand"] = strand
    df["ALT"] = df["ALT"].astype(str).str.slice(0, 1)

    # Drop rows with missing core columns (similar to R’s drop_na )
    core_cols = ["CHROM", "POS", "REF", "ALT", "FORMAT"]
    df = df.dropna(subset=[c for c in core_cols if c in df.columns]).copy()

    # RNAREF / RNAALT
    df["RNAREF"] = df["REF"].astype(str)
    df["RNAALT"] = df["ALT"].astype(str)
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if strand == "neg":
        df["RNAREF"] = df["RNAREF"].map(comp).fillna(df["RNAREF"])
        df["RNAALT"] = df["RNAALT"].map(comp).fillna(df["RNAALT"])

    # Sample name: User-provided is preferred, otherwise inferred from the header and cleaned
    if args.samples and len(args.samples) > 0:
        if len(args.samples) != len(sample_cols_from_header):
            print("[WARN] The number of --samples is inconsistent with the number of VCF sample columns. The VCF header will take precedence.")
            sList: List[str] = [sanitize_sample(x) for x in sample_cols_from_header]
        else:
            sList = args.samples
    else:
        sList = [sanitize_sample(x) for x in sample_cols_from_header]

    sample_cols = sample_cols_from_header  # 与 df 中样本列一一对应

    print("[INFO] Sample columns (raw):" + ", ".join(sample_cols))
    print("[INFO] Sample names (for column prefixes):" + ", ".join(sList))

    # Parse DP/AD for each sample by FORMAT (fallback to fixed position)
    def extract_sample_dp_ad(sample_series: pd.Series, format_series: pd.Series):
        dp_vals, ad_vals = [], []
        for fmt, cell in zip(format_series, sample_series):
            dp_idx = find_column_index(fmt, "DP")
            ad_idx = find_column_index(fmt, "AD")
            if dp_idx is None: dp_idx = 1   # Fallback: Paragraph 2  # noqa: E701
            if ad_idx is None: ad_idx = 5   # Fallback: Paragraph 6  # noqa: E701
            dp_str = safe_split(cell, ":", dp_idx)
            ad_str = safe_split(cell, ":", ad_idx)
            alt1  = safe_split(ad_str, ",", 1) if pd.notna(ad_str) else np.nan
            dp_vals.append(pd.to_numeric(dp_str, errors="coerce"))
            ad_vals.append(pd.to_numeric(alt1,  errors="coerce"))
        return pd.Series(dp_vals), pd.Series(ad_vals)

    for sample_col, sp in zip(sample_cols, sList):
        dp_s, ad_s = extract_sample_dp_ad(df[sample_col], df["FORMAT"])
        df[f"{sp}.DP"] = dp_s
        df[f"{sp}.AD"] = ad_s
        with np.errstate(divide="ignore", invalid="ignore"):
            df[f"{sp}.ratio"] = (df[f"{sp}.AD"].astype(float) /
                                 df[f"{sp}.DP"].astype(float)).replace([np.inf, -np.inf], np.nan)

    # Average coverage (DP average)
    dp_cols = [f"{sp}.DP" for sp in sList]
    df["avgDP"] = df[dp_cols].mean(axis=1, skipna=True)

# === Write out the file (Windows writes directly; Linux/WSL writes to /tmp first and then copies) ===


    out_path = os.path.join(output, f"mpileup_fixstrand_{strand}.vcf")
    out_path = os.path.normpath(out_path)
    print("[INFO] Target output path repr:", repr(out_path))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    if platform.system() == "Windows":
        # Windows Write the target path directly
        df.to_csv(out_path, sep="\t", index=False, header=True)
        print(f"[OK] Output File: {out_path}  (total {len(df)} lines）")
    else:
        # Linux/WSL：Write to /tmp first, then copy to /mnt/disk to avoid drvfs Errno 22
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False,
                                        dir="/tmp", prefix="vcf_", suffix=".vcf") as tf:
            tmp_out = tf.name
        df.to_csv(tmp_out, sep="\t", index=False, header=True)
        try:
            shutil.copyfile(tmp_out, out_path)
            print(f"[OK] Output File: {out_path}  (total {len(df)} lines）")
        except OSError as e:
            print("[WARN] Copying to the target directory failed:", repr(e))
            print(f"[HINT] Please copy manually: cp {tmp_out} {out_path}")

if __name__ == "__main__":
    main()




