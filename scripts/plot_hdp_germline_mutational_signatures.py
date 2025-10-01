#!/usr/bin/env python

import argparse
from pathlib import Path
import sys
from typing import Dict, List

from natsort import natsorted
import pandas as pd
import plotnine as p9

SUBS = ["C>A", "C>G", "C>T", "T>A", "T>G"]


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="file to read germline mutational signatures",
    )
    parser.add_argument(
        "--sbs96-to-sbs52-lookup",
        type=Path,
        required=True,
        help="file to read SBS96 to SBS52 lookup table",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="file to plot germline mutational signatures",
    )
    args = args[1:]
    return parser.parse_args(args)


def load_sbs96_to_sbs52_lookup(sbs96_to_sbs52_lookup_table_path: Path) -> Dict[str, str]:
    sbs96_to_sbs52_lookup = {}
    with open(sbs96_to_sbs52_lookup_table_path) as f:
        for line in f:
            if line.startswith("SBS96"):
                continue
            fields = line.strip().split("\t")
            sbs96 = fields[0]
            sbs52 = fields[1]
            sbs96_to_sbs52_lookup[sbs96] = sbs52
    return sbs96_to_sbs52_lookup


def load_sbs52_classification(sbs96_to_sbs52_lookup: Dict[str, str]) -> List[str]:
    sub_to_sbs52_classification = {sub: [] for sub in SUBS}
    for sbs52 in list(set(sbs96_to_sbs52_lookup.values())):
        _ubase, _, ref, _, alt, _, _dbase = list(sbs52)
        sub = "{}>{}".format(ref, alt)
        sub_to_sbs52_classification[sub].append(sbs52)

    sbs52_classification = []
    for sub in SUBS:
        sbs52_classification.extend(natsorted(sub_to_sbs52_classification[sub]))
    return sbs52_classification


def plot_signature(df: Path, sig_name: str) -> p9.ggplot:
    df_subset = df[df["SIGNATURE"] == sig_name]
    plot = (
        p9.ggplot(df_subset, p9.aes(x="TRI", y="PROBABILITY", fill="SUB"))
        + p9.geom_bar(stat="identity")
        + p9.theme_bw(16)
        + p9.facet_grid(". ~ SUB", scales="free")
        + p9.scale_fill_manual(values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#F5ABCC"))
        + p9.labs(title=f"\n{sig_name}\n", x="\nTrinucleotide\n", y="\nProportion of Mutations (%)\n")
        + p9.theme(
            figure_size=(22, 12),
            text=p9.element_text(family="Helvetica"),
            plot_title=p9.element_text(ha="left", weight="bold"),
            legend_position="none",
            legend_title=p9.element_blank(),
            axis_text_x=p9.element_text(family="monospace", angle=90, ha="center"),
        )
    )
    return plot


def load_signatures(input_path: Path, sbs52_classification: List[str]) -> pd.DataFrame:
    # Parse header
    with open(input_path) as f:
        header = f.readline().strip().split(",")
    sig_names = header[1:]
    sig_names = [sig_name.replace("X", "gToL") for sig_name in sig_names]

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
            sbs52 = sbs52_classification[op_idx]
            ubase, _, ref, _, alt, _, dbase = list(sbs52)
            sub = f"{ref}>{alt}"
            tri = f"{ubase}{ref}{dbase}"
            rows.append({"SIGNATURE": sig, "SBS52": sbs52, "SUB": sub, "TRI": tri, "PROBABILITY": prob})
    return pd.DataFrame(rows)


def plot_germline_mutational_signatures(input_path: Path, sbs96_to_sbs52_lookup_path: Path, output_path: Path):
    # Load lookup table
    sbs96_to_sbs52_lookup = load_sbs96_to_sbs52_lookup(sbs96_to_sbs52_lookup_path)
    sbs52_classification = load_sbs52_classification(sbs96_to_sbs52_lookup)

    # Load signatures
    sig_df = load_signatures(input_path, sbs52_classification)

    # Plot each signature
    plots = []
    sig_names = natsorted(sig_df["SIGNATURE"].unique())
    for sig_name in sig_names:
        p = plot_signature(sig_df, sig_name)
        plots.append(p)
    p9.save_as_pdf_pages(plots, output_path)


def main():
    options = parse_args(sys.argv)
    plot_germline_mutational_signatures(options.input, options.sbs96_to_sbs52_lookup, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
