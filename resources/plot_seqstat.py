#!/usr/bin/env python3
"""
Script to model and visualize the relationship between sequencing reads and CpG counts
using an asymptotic growth function (arctangent). The script fits the model to downsampling data
and visualizes sequencing saturation across different samples.

Author: Edoardo Giuili
"""

import sys
import os
import warnings
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit, OptimizeWarning


# ===============================================================
# Mathematical Model
# ===============================================================

def asymptotic_growth(x, beta0, beta1):
    """Asymptotic growth model using the arctangent function."""
    return beta0 * np.arctan(beta1 * x)


def derivative_asymptotic_growth(x, beta0, beta1):
    """Derivative of the asymptotic growth model."""
    return beta0 * beta1 / (1 + (beta1 * x) ** 2)


def find_asymptote(params):
    """Return the asymptote value (y-limit as x → ∞)."""
    beta0, _ = params
    return beta0 * np.pi / 2


# ===============================================================
# Plotting Functions
# ===============================================================

def plot_error_data(x_data, y_data, reads, output_path, title, error):
    """Plot the data when curve fitting fails."""
    x_pred = np.append(x_data, [1.2, 1.4, 1.6, 1.8, 2.0])
    x_pred_reads = reads * x_pred

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(x_data, y_data, "o", color="blue")

    ax1.set_xticks(x_pred)
    ax1.set_title(f"{title}\n{error}", fontsize=12, color="red", loc="left")
    ax1.set_xlabel("Percentage of downsampling", fontsize=14)
    ax1.set_ylabel("Number of CpGs", fontsize=14)
    ax1.grid(True, linestyle="--", alpha=0.6)

    # Add secondary x-axis for read counts
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(x_pred)
    ax2.tick_params(axis="x", pad=7)
    ax2.set_xticklabels([f"{x:.1e}" for x in x_pred_reads], rotation=45)
    ax2.set_xlabel("Number of reads", fontsize=14)

    fig.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"⚠️  Curve fitting failed. Plot saved to {output_path}")


def plot_data(x_data, y_data, reads, asymptote, params, output_path, title):
    """Plot data and fitted curve for a single sample."""
    x_pred = np.append(x_data, [1.2, 1.4, 1.6, 1.8, 2.0])
    x_pred_reads = reads * x_pred
    y_pred = asymptotic_growth(x_pred, *params)
    y_diff = [0, 1] + [(y_pred[i] / asymptote) * 100 for i in range(2, len(y_pred))]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Raw data + fit
    ax1.plot(x_data, y_data, "o", color="blue")
    ax1.plot(x_pred, y_pred, "g-", label=f"fit: β₀={params[0]:.3f}, β₁={params[1]:.3f}")
    ax1.plot(x_pred, y_pred, "|", color="black", markersize=10)

    # Add percentage annotations
    for i in range(3, len(x_pred)):
        ax1.text(
            x_pred[i], y_pred[i] - 0.09 * y_pred[i],
            f"{y_diff[i]:.1f}%", color="black", fontsize=10, ha="center"
        )

    # Add asymptote line
    ax1.axhline(y=asymptote, color="grey", linestyle="--")
    ax1.text(0, asymptote - asymptote * 0.05,
             f"Asymptote = {asymptote:.2e}", color="grey", fontsize=12)

    # Formatting
    ax1.set_title(title, fontsize=16, loc="left")
    ax1.set_xlabel("Percentage of downsampling", fontsize=14)
    ax1.set_ylabel("Number of CpGs", fontsize=14)
    ax1.grid(True, linestyle="--", alpha=0.6)

    # Secondary axis for reads
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(x_pred)
    ax2.tick_params(axis="x", pad=7)
    ax2.set_xticklabels([f"{x:.1e}" for x in x_pred_reads], rotation=45)
    ax2.set_xlabel("Number of reads", fontsize=14)

    # Legend
    leg_patch = mpatches.Patch(
        label=r"% : Sequencing saturation ($\frac{\hat{y}}{\text{asymptote}} \times 100$)"
    )
    plt.legend(handles=[leg_patch], loc="lower right", handletextpad=-1.0, handlelength=0)

    fig.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"✅ Plot saved to {output_path}")


# ===============================================================
# Core Functions
# ===============================================================

def plot_reads_vs_cpgs(data, output_path, percentages):
    """Fit the asymptotic model and plot reads vs CpGs."""
    x_data = np.array([0] + percentages)
    y_data = np.array([0] + data["cpgs_counts"].tolist())

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            params, _ = curve_fit(asymptotic_growth, x_data, y_data, p0=[1, 1])
        asymptote = find_asymptote(params)
        fit_success = True
    except (RuntimeError, OptimizeWarning) as e:
        fit_success = False
        params, asymptote, error_msg = None, None, str(e)
        print(f"⚠️ Curve fitting error: {e}")

    title = data["sample"].iloc[0]
    reads = int(data["reads_counts"].iloc[-1])

    if fit_success:
        plot_data(x_data, y_data, reads, asymptote, params, output_path, title)
    else:
        plot_error_data(x_data, y_data, reads, output_path, title, error_msg)


def select_sample(cpgs_file, reads_file, percentages,out_dir):
    """Main routine to process each sample and generate plots."""
    col_cpgs = ["sample", "percentage", "min_counts", "cpgs_counts"]
    col_reads = ["sample", "percentage", "reads_counts"]

    cpg_df = pd.read_csv(cpgs_file, sep=",", header=None, names=col_cpgs)
    reads_df = pd.read_csv(reads_file, sep=",", header=None, names=col_reads)
    data = pd.merge(cpg_df, reads_df, on=["sample", "percentage"])

    for sample in data["sample"].unique():
        sample_data = data[data["sample"] == sample]
        for min_val in sample_data["min_counts"].unique():
            subset = sample_data[sample_data["min_counts"] == min_val].sort_values(by="percentage")
            #if out_dir is
            output_path = f"{out_dir}/{sample}_{min_val}x_plot.svg"

            if subset.duplicated().any():
                print(f"❌ Duplicated rows found in {sample}. "
                      "Check input files.")
                sys.exit(1)

            plot_reads_vs_cpgs(subset, output_path, percentages)


# ===============================================================
# Main
# ===============================================================

def main():
    parser = argparse.ArgumentParser(description="Plot Reads vs. CpGs from two CSV files.")
    parser.add_argument("--cpgs_file", required=True, help="Path to the CpGs CSV file.")
    parser.add_argument("--read_file", required=True, help="Path to the Reads CSV file.")
    parser.add_argument("--percentages", required=True,
                        help="Comma-separated list of downsampling percentages (e.g. 0.25,0.5,0.75,1).")
    parser.add_argument("--outdir", required=False, help="Output directory path to save the plots. Default is current directory.", default=".")

    args = parser.parse_args()
    percentages = [float(p) for p in args.percentages.split(",")]

    # Check if output directory exists, create if not
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)
        print(f"Created output directory: {args.outdir}")
    else:
        print(f"Using existing output directory: {args.outdir}")

    select_sample(args.cpgs_file, args.read_file, percentages, args.outdir)


if __name__ == "__main__":
    main()
