#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Parse bcftools mpileup VCF, extract per-sample DP/AD and compute ratios.

特点：
- 保留 #CHROM 行为列名；跳过 ## 元信息
- 兼容 .vcf 与 .vcf.gz
- Win/WSL 路径自动规范化
- 用 --input 指定输入目录；--output 指定输出目录
- 样本列从 header 推断（也可 --samples 指定）
- 优先根据 FORMAT 动态定位 DP/AD；找不到则回退到固定位置
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

# ---------- 路径规范化（Windows -> WSL） ----------
def normalize_path(p: Optional[str]) -> Optional[str]:
    if p is None:
        return p
    p = p.strip().strip('"').strip("'")
    # Linux/WSL 下把 "E:\path" 或 "e:/path" 映射到 "/mnt/e/path"
    if platform.system() == "Linux":
        m = re.match(r'^([a-zA-Z]):[\\/](.*)$', p)
        if m:
            drive = m.group(1).lower()
            rest = m.group(2).replace('\\', '/')
            return f"/mnt/{drive}/{rest}"
    return p.replace('\\', '/')

# ---------- 安全 split ----------
def safe_split(s: Optional[str], sep: str, idx: int) -> Optional[str]:
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return np.nan
    parts = str(s).split(sep)
    return parts[idx] if 0 <= idx < len(parts) else np.nan

def find_column_index(fmt: str, key: str) -> Optional[int]:
    """在 FORMAT 字段中找 key 的 0-based 位置"""
    if fmt is None or (isinstance(fmt, float) and pd.isna(fmt)):
        return None
    fields = str(fmt).split(':')
    for i, k in enumerate(fields):
        if k == key:
            return i
    return None

# ---------- 样本名清洗（去路径、去 .bam/.cram、非法字符 -> _） ----------
def sanitize_sample(name: str) -> str:
    base = os.path.basename(name)
    base = re.sub(r"\.(bam|cram)$", "", base, flags=re.IGNORECASE)
    base = re.sub(r"[^\w.\-]+", "_", base)
    return base

# ---------- 参数 ----------
def parse_args():
    ap = argparse.ArgumentParser(description="Fix strand and compute mutation ratios from VCF.")
    ap.add_argument("s", type=int, help="strand: 1=FWD/neg, 2=REV/pos")
    ap.add_argument("--input", required=True, help="Input VCF directory (including ProQ-rABE_FWD/REV.vcf[.gz])")
    ap.add_argument("--output", required=True, help="Output Directory")
    ap.add_argument("--samples", nargs="*", default=None,
                    help="Optional: list of sample names; if not provided, they will be automatically inferred from the VCF header")
    return ap.parse_args()

# ---------- 主流程 ----------
def main():
    args = parse_args()

    # strand
    if args.s == 1:
        st, strand = "FWD", "neg"
    elif args.s == 2:
        st, strand = "REV", "pos"
    else:
        raise SystemExit("参数 s 只能是 1 或 2。")

    input = normalize_path(args.input)
    output  = normalize_path(args.output)
    os.makedirs(output, exist_ok=True)

    # 输入文件（优先 .vcf，无则尝试 .vcf.gz）
    vcf_path = os.path.join(input, f"ProQ-rABE_{st}.vcf")
    vcf_gz   = vcf_path + ".gz"
    if not os.path.exists(vcf_path):
        if os.path.exists(vcf_gz):
            vcf_path = vcf_gz
        else:
            raise FileNotFoundError(f"未找到输入文件：\n  {vcf_path}\n或 {vcf_gz}")

    print(f"[INFO] 读取文件: {vcf_path}")

    # 读取 header（#CHROM 行）获取列名与样本列
    open_func = gzip.open if vcf_path.endswith(".gz") else open
    with open_func(vcf_path, "rt") as f:
        header_line = None
        for line in f:
            if line.startswith("#CHROM"):
                header_line = line.strip().lstrip("#").split("\t")
                break
    if header_line is None:
        raise ValueError("未在 VCF 中找到 #CHROM 表头行。")

    # 样本列按 header 的第 10 列起
    sample_cols_from_header = header_line[9:]

    # 读入数据：跳过以 # 开头的行，使用 header_line 作为列名
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

    # strand 与 ALT 处理
    df["strand"] = strand
    df["ALT"] = df["ALT"].astype(str).str.slice(0, 1)

    # 核心列缺失的行去掉（与 R 的 drop_na 类似）
    core_cols = ["CHROM", "POS", "REF", "ALT", "FORMAT"]
    df = df.dropna(subset=[c for c in core_cols if c in df.columns]).copy()

    # RNAREF / RNAALT（负链互补）
    df["RNAREF"] = df["REF"].astype(str)
    df["RNAALT"] = df["ALT"].astype(str)
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    if strand == "neg":
        df["RNAREF"] = df["RNAREF"].map(comp).fillna(df["RNAREF"])
        df["RNAALT"] = df["RNAALT"].map(comp).fillna(df["RNAALT"])

    # 样本名：优先用户提供，否则从 header 推断并清洗
    if args.samples and len(args.samples) > 0:
        if len(args.samples) != len(sample_cols_from_header):
            print("[WARN] --samples 数量与 VCF 样本列数量不一致，将以 VCF header 为准。")
            sList: List[str] = [sanitize_sample(x) for x in sample_cols_from_header]
        else:
            sList = args.samples
    else:
        sList = [sanitize_sample(x) for x in sample_cols_from_header]

    sample_cols = sample_cols_from_header  # 与 df 中样本列一一对应

    print("[INFO] 样本列（原始）：" + ", ".join(sample_cols))
    print("[INFO] 样本名（用于列前缀）：" + ", ".join(sList))

    # 按 FORMAT 解析每个样本的 DP/AD（回退到固定位置）
    def extract_sample_dp_ad(sample_series: pd.Series, format_series: pd.Series):
        dp_vals, ad_vals = [], []
        for fmt, cell in zip(format_series, sample_series):
            dp_idx = find_column_index(fmt, "DP")
            ad_idx = find_column_index(fmt, "AD")
            if dp_idx is None: dp_idx = 1   # 回退：第2段  # noqa: E701
            if ad_idx is None: ad_idx = 5   # 回退：第6段  # noqa: E701
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

    # 平均覆盖（DP 平均）
    dp_cols = [f"{sp}.DP" for sp in sList]
    df["avgDP"] = df[dp_cols].mean(axis=1, skipna=True)

# === 写出文件（Windows 直接写；Linux/WSL 先写 /tmp 再复制） ===


    out_path = os.path.join(output, f"mpileup_fixstrand_{strand}.vcf")
    out_path = os.path.normpath(out_path)
    print("[INFO] 目标输出路径 repr:", repr(out_path))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    if platform.system() == "Windows":
        # Windows 直接写目标路径
        df.to_csv(out_path, sep="\t", index=False, header=True)
        print(f"[OK] 输出文件: {out_path}  （共 {len(df)} 行）")
    else:
        # Linux/WSL：先写 /tmp，再复制到 /mnt/盘，规避 drvfs Errno 22
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False,
                                        dir="/tmp", prefix="vcf_", suffix=".vcf") as tf:
            tmp_out = tf.name
        df.to_csv(tmp_out, sep="\t", index=False, header=True)
        try:
            shutil.copyfile(tmp_out, out_path)
            print(f"[OK] 输出文件: {out_path}  （共 {len(df)} 行）")
        except OSError as e:
            print("[WARN] 复制到目标目录失败：", repr(e))
            print(f"[HINT] 请手动复制：cp {tmp_out} {out_path}")

if __name__ == "__main__":
    main()


