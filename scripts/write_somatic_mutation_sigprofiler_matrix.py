#!/usr/bin/env python

import argparse
import math
import os
from pathlib import Path
from typing import List
import sys

import natsort


NTS = ["A", "C", "G", "T"]
SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SBS96_CLASSIFICATIONS = [f"{nti}[{sub}]{ntj}" for sub in SUBS for nti in NTS for ntj in NTS]


def parse_args(args):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="file of file names to read somatic mutation counts from"
    )
    parser.add_argument("--samples", type=Path, required=True, help="list of samples to include into the HDP matrix")
    parser.add_argument("-o", "--output", type=Path, required=True, help="sigProfiler matrix")
    args = args[1:]
    return parser.parse_args(args)


def get_common_suffix(strings: List[str]) -> str:
    # Reverse each string
    rev_strings = [s[::-1] for s in strings]
    # Find common *prefix* of the reversed strings
    rev_prefix = os.path.commonprefix(rev_strings)
    # Reverse back to get the suffix
    return rev_prefix[::-1]


def write_sigprofiler_matrix(input_path: Path, sample_path: Path, output_path: Path):
    # Load file paths
    file_paths = [line.rstrip() for line in open(input_path).readlines()]

    # Get common suffix
    suffix = get_common_suffix(file_paths[0:5])

    # Load target sampels
    target_samples = set([line.strip() for line in open(sample_path).readlines()])

    # Aggregate counts
    query_samples = []
    sample_sbs96_counts = {}
    for file_path in file_paths:
        query_sample = file_path.replace(suffix, "")
        if query_sample not in target_samples:
            continue
        query_samples.append(query_sample)
        for line in open(file_path).readlines():
            if line.startswith("#"):
                continue
            if line.startswith("sub"):
                continue
            fields = line.strip().split()
            sub = fields[0]
            tri = fields[1]
            count = fields[3]
            upstream, _, downstream = list(tri)
            sbs96 = f"{upstream}[{sub}]{downstream}"
            sample_sbs96_counts[f"{query_sample}:{sbs96}"] = count

    # Write matrix
    target_samples = natsort.natsorted(list(target_samples))
    with open(output_path, "w") as outfile:
        outfile.write("{}\n".format("\t".join(["Mutation Types"] + query_samples)))
        for sbs96 in SBS96_CLASSIFICATIONS:
            row_values = [sbs96]
            for query_sample in query_samples:
                row_values.append(str(math.ceil(float(sample_sbs96_counts[f"{query_sample}:{sbs96}"]))))
            outfile.write("{}\n".format("\t".join(row_values)))


def main():
    options = parse_args(sys.argv)
    write_sigprofiler_matrix(options.input, options.samples, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
