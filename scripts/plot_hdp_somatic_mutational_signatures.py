#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from natsort import natsorted
import pandas as pd
import plotnine as p9

NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATION = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read HDP somatic mutational signatures",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to plot HDP somatic mutational signatures",
    )
    args = args[1:]
    return parser.parse_args(args)


def plot_signature(df: Path, sig_name: str) -> p9.ggplot:
    df_subset = df[df["SIGNATURE"] == sig_name]
    plot = (
        p9.ggplot(df_subset, p9.aes(x="TRI", y="PROBABILITY", fill="SUB"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw(16)
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(
            values=(
                "#98D7EC",
                "#212121",
                "#FF003A",
                "#A6A6A6",
                "#83A603",
                "#F5ABCC"
                )
        )
        + p9.labs(
            title=f"\n{sig_name}\n",
            x="\nTrinucleotide\n",
            y="\nProportion of Mutations (%)\n"
        )
        + p9.theme(
            figure_size=(22, 12),
            text=p9.element_text(family="Helvetica"),
            plot_title=p9.element_text(ha="left", weight="bold"),
            legend_position="none",
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(
                family="monospace", angle=90, ha="center"
            ),
        )
    )
    return plot


def load_signatures(input_path: Path) -> pd.DataFrame:
    # Parse header
    with open(input_path) as f:
        header = f.readline().strip().split(",")
    sig_names = header[1:]
    sig_names = [sig_name.replace("X", "sToL") for sig_name in sig_names]

    # Collect mutational probabilities per signature
    probs_per_sig = {sig_name: [] for sig_name in sig_names}
    with open(input_path) as f:
        next(f)  # skip header
        for line in f:
            fields = line.rstrip().split(",")
            for cidx, prob in enumerate(fields[1:]):
                sig_name = sig_names[cidx]
                probs_per_sig[sig_name].append(float(prob))

    # Now convert into a long DataFrame using SBS96_LST mapping
    rows = []
    for sig, probs in probs_per_sig.items():
        for op_idx, prob in enumerate(probs):
            sbs96 = SBS96_CLASSIFICATION[op_idx]
            ubase, _, ref, _, alt, _, dbase = list(sbs96)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            rows.append({
                "SIGNATURE": sig,
                "SBS96": sbs96,
                "SUB": sub,
                "TRI": tri,
                "PROBABILITY": prob
            })
    return pd.DataFrame(rows)


def plot_somatic_mutational_signatures(input_path: Path, output_path: Path):
    # Load signatures
    sig_df = load_signatures(input_path)
    sig_names = natsorted(sig_df["SIGNATURE"].unique())

    # Plot each signature
    plots = []
    for sig_name in sig_names:
        p = plot_signature(sig_df, sig_name)
        plots.append(p)
    p9.save_as_pdf_pages(plots, output_path)


def main():
    options = parse_args(sys.argv)
    plot_somatic_mutational_signatures(options.input, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
