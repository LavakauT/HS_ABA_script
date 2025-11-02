
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fastp_qc_summary.py
-------------------
Summarize multiple fastp JSON summaries into one CSV and plots.

Usage:
  python fastp_qc_summary.py --input /path/to/jsons --output ./qc_out --pattern "*fastp.json"

Outputs:
  - qc_out/fastp_qc_summary.csv       (tidy table, before/after rows per sample)
  - qc_out/fastp_qc_wide.csv          (wide table, one row per sample with prefixed columns)
  - qc_out/fastp_qc_plots.pdf         (multi-page PDF with matplotlib charts)

Notes:
  - "mean_length" is the average of read1_mean_length and read2_mean_length.
  - If your libraries are single-end, the script will handle missing read2 fields.
"""

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", "-i", required=True, help="Directory containing fastp JSON files")
    ap.add_argument("--output", "-o", required=True, help="Output directory")
    ap.add_argument("--pattern", "-p", default="*fastp.json", help="Glob to match JSON files (default: *fastp.json)")
    ap.add_argument("--sample-from", choices=["stem", "parent"], default="stem",
                    help="Get sample name from file stem or parent dir name (default: stem)")
    return ap.parse_args()

def safe_get(d: Dict[str, Any], *keys, default=np.nan):
    cur = d
    for k in keys:
        if cur is None:
            return default
        cur = cur.get(k, None)
    return default if cur is None else cur

def extract_one(json_path: Path, sample_from: str = "stem") -> pd.DataFrame:
    with open(json_path, "r") as f:
        data = json.load(f)

    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})

    # Try to capture read lengths; handle single-end (missing read2 fields)
    r1_len_b = safe_get(before, "read1_mean_length", default=np.nan)
    r2_len_b = safe_get(before, "read2_mean_length", default=np.nan)
    r1_len_a = safe_get(after, "read1_mean_length", default=np.nan)
    r2_len_a = safe_get(after, "read2_mean_length", default=np.nan)

    def mean_len(r1, r2):
        if pd.isna(r1) and pd.isna(r2):
            return np.nan
        if pd.isna(r2):
            return float(r1)
        if pd.isna(r1):
            return float(r2)
        return (float(r1) + float(r2)) / 2.0

    sample = json_path.stem if sample_from == "stem" else json_path.parent.name

    rows = []
    rows.append({
        "sample": sample,
        "stage": "before",
        "total_reads": safe_get(before, "total_reads"),
        "total_bases": safe_get(before, "total_bases"),
        "q20_bases": safe_get(before, "q20_bases"),
        "q30_bases": safe_get(before, "q30_bases"),
        "q20_rate": safe_get(before, "q20_rate"),
        "q30_rate": safe_get(before, "q30_rate"),
        "read1_mean_length": r1_len_b,
        "read2_mean_length": r2_len_b,
        "mean_length": mean_len(r1_len_b, r2_len_b),
        "gc_content": safe_get(before, "gc_content"),
        "fastp_version": safe_get(summary, "fastp_version", default=""),
        "sequencing": safe_get(summary, "sequencing", default=""),
        "source_file": str(json_path),
    })
    rows.append({
        "sample": sample,
        "stage": "after",
        "total_reads": safe_get(after, "total_reads"),
        "total_bases": safe_get(after, "total_bases"),
        "q20_bases": safe_get(after, "q20_bases"),
        "q30_bases": safe_get(after, "q30_bases"),
        "q20_rate": safe_get(after, "q20_rate"),
        "q30_rate": safe_get(after, "q30_rate"),
        "read1_mean_length": r1_len_a,
        "read2_mean_length": r2_len_a,
        "mean_length": mean_len(r1_len_a, r2_len_a),
        "gc_content": safe_get(after, "gc_content"),
        "fastp_version": safe_get(summary, "fastp_version", default=""),
        "sequencing": safe_get(summary, "sequencing", default=""),
        "source_file": str(json_path),
    })
    return pd.DataFrame(rows)

def gather(input_dir: Path, pattern: str, sample_from: str) -> pd.DataFrame:
    files = sorted(input_dir.rglob(pattern))
    if not files:
        raise SystemExit(f"No files matched pattern '{pattern}' under {input_dir}")
    dfs = [extract_one(p, sample_from=sample_from) for p in files]
    return pd.concat(dfs, ignore_index=True)

def make_wide(df: pd.DataFrame) -> pd.DataFrame:
    # Pivot before/after to columns
    metrics = ["total_reads", "total_bases", "q20_bases", "q30_bases",
               "q20_rate", "q30_rate", "read1_mean_length", "read2_mean_length",
               "mean_length", "gc_content"]
    meta = ["fastp_version", "sequencing", "source_file"]
    pivot_df = df.pivot(index="sample", columns="stage", values=metrics)
    # Flatten multiindex columns
    pivot_df.columns = [f"{m}_{s}" for m, s in pivot_df.columns]
    # Get meta from 'before' rows
    meta_df = (df[df["stage"]=="before"]
               .set_index("sample")[meta]
               .drop_duplicates())
    out = pivot_df.join(meta_df, how="left")
    out.reset_index(inplace=True)
    return out

def plot_pdf(df: pd.DataFrame, out_pdf: Path):
    # df is tidy (before/after). We'll draw a separate page for each metric.
    metrics = [
        ("q30_rate", "Q30 rate"),
        ("q20_rate", "Q20 rate"),
        ("total_reads", "Total reads"),
        ("total_bases", "Total bases (bp)"),
        ("mean_length", "Mean read length (bp)"),
        ("gc_content", "GC content"),
    ]

    # Order samples by sample name for consistent plotting
    samples = sorted(df["sample"].unique())
    stage_order = ["before", "after"]

    with PdfPages(out_pdf) as pdf:
        for key, title in metrics:
            sub = df[["sample", "stage", key]].copy()
            # Ensure plotting order
            sub["sample"] = pd.Categorical(sub["sample"], categories=samples, ordered=True)
            sub["stage"] = pd.Categorical(sub["stage"], categories=stage_order, ordered=True)
            sub = sub.sort_values(["sample", "stage"])

            # Create a simple grouped bar plot with matplotlib only
            # one chart per metric (per tool rules).
            fig, ax = plt.subplots(figsize=(max(8, len(samples)*0.5), 5))
            # x positions for groups
            x = np.arange(len(samples))
            width = 0.35  # bar width for before/after

            # split values
            vals_before = sub[sub["stage"]=="before"][key].values
            vals_after  = sub[sub["stage"]=="after"][key].values

            ax.bar(x - width/2, vals_before, width, label="before")
            ax.bar(x + width/2, vals_after,  width, label="after")

            ax.set_title(title)
            ax.set_xlabel("Sample")
            ax.set_ylabel(key)
            ax.set_xticks(x)
            ax.set_xticklabels(samples, rotation=90)
            ax.legend()

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

def main():
    args = parse_args()
    in_dir = Path(args.input).expanduser().resolve()
    out_dir = Path(args.output).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    df = gather(in_dir, args.pattern, args.sample_from)
    # Save tidy
    tidy_csv = out_dir / "fastp_qc_summary.csv"
    df.to_csv(tidy_csv, index=False)

    # Save wide
    wide = make_wide(df)
    wide_csv = out_dir / "fastp_qc_wide.csv"
    wide.to_csv(wide_csv, index=False)

    # Plots
    out_pdf = out_dir / "fastp_qc_plots.pdf"
    plot_pdf(df, out_pdf)

    print(f"Wrote tidy CSV: {tidy_csv}")
    print(f"Wrote wide CSV: {wide_csv}")
    print(f"Wrote plots PDF: {out_pdf}")

if __name__ == "__main__":
    main()
