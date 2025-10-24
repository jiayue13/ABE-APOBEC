#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
功能：
  解析 bcftools mpileup 生成的 VCF，按样本提取 DP、AD 并计算 ratio。
  支持正反链、指定输入和输出目录。

用法示例：
  python proq_fixstrand.py 1 \
    --workdir "e:/Kpn_data/seq_data/20251028/results" \
    --resdir "e:/Kpn_data/seq_data/20251028/output"

参数：
  s: 1 -> FWD/neg strand
     2 -> REV/pos strand
  --workdir: 输入文件所在目录（包含 ProQ-rABE_FWD.vcf / ProQ-rABE_REV.vcf）
  --resdir:  输出文件保存目录
  --samples: 样本名顺序，与 VCF 中一致
"""

import os
import argparse
import pandas as pd
import numpy as np

def parse_args():
    ap = argparse.ArgumentParser(description="Fix strand and compute mutation ratios from VCF.")
    ap.add_argument("s", type=int, help="strand code: 1=FWD/neg, 2=REV/pos")
    ap.add_argument("--workdir", required=True, help="输入 VCF 文件所在目录")
    ap.add_argument("--resdir", required=True, help="结果输出目录")
    ap.add_argument("--samples", nargs="*", default=["T.H.1","T.H.2","T.H.3","C.H.1","C.H.2","C.H.3"],
                    help="样本名列表，与 VCF 样本列顺序一致")
    return ap.parse_args()

def safe_split(s, sep, idx):
    if pd.isna(s): return np.nan  # noqa: E701
    parts = str(s).split(sep)
    return parts[idx] if 0 <= idx < len(parts) else np.nan

def main():
    args = parse_args()

    # strand setting
    if args.s == 1:
        st, strand = "FWD", "neg"
    elif args.s == 2:
        st, strand = "REV", "pos"
    else:
        raise SystemExit("参数 s 只能是 1 或 2。")

    # make sure output directory exists
    os.makedirs(args.resdir, exist_ok=True)

    # read vcf
    vcf_path = os.path.join(args.workdir, f"ProQ-rABE_{st}.vcf")
    print(f"[INFO] 读取文件: {vcf_path}")
    df = pd.read_csv(vcf_path, sep="\t", header=None, dtype=str, na_filter=True)

    # add strand column
    df["strand"] = strand

    # keep only first ALT base
    df.iloc[:, 4] = df.iloc[:, 4].astype(str).str.slice(0, 1)

    # drop rows with NA
    df = df.dropna(how="any").copy()

    # build RNAREF / RNAALT
    df["RNAREF"] = df.iloc[:, 3].astype(str)
    df["RNAALT"] = df.iloc[:, 4].astype(str)

    # complement on neg strand
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if strand == "neg":
        df["RNAREF"] = df["RNAREF"].map(comp).fillna(df["RNAREF"])
        df["RNAALT"] = df["RNAALT"].map(comp).fillna(df["RNAALT"])

    sList = args.samples
    sample_cols = list(df.columns[9:9 + len(sList)])

    # extract DP / AD and calculate ratios
    for col, sp in zip(sample_cols, sList):
        df[f"{sp}.DP"] = pd.to_numeric(df[col].apply(lambda x: safe_split(x, ":", 1)), errors="coerce")
        ad_raw = df[col].apply(lambda x: safe_split(x, ":", 5))
        df[f"{sp}.AD"] = pd.to_numeric(ad_raw.apply(lambda x: safe_split(x, ",", 1)), errors="coerce")

    for sp in sList:
        dp, ad = df[f"{sp}.DP"], df[f"{sp}.AD"]
        with np.errstate(divide="ignore", invalid="ignore"):
            df[f"{sp}.ratio"] = ad.astype(float) / dp.astype(float)

    df["avgDP"] = df[[f"{sp}.DP" for sp in sList]].mean(axis=1, skipna=True)

    # output
    out_path = os.path.join(args.resdir, f"mpileup_fixstrand_{strand}.vcf")
    df.to_csv(out_path, sep="\t", index=False, header=True)
    print(f"[OK] 输出文件: {out_path}")

if __name__ == "__main__":
    main()
