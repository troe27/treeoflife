#!/usr/bin/env python

import argparse
import math
import os
from pathlib import Path
from typing import List
import sys

import natsort


def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="file of file names to read germline mutation counts from"
    )
    parser.add_argument("--samples", type=Path, required=True, help="list of samples to include into the HDP matrix")
    parser.add_argument("-o", "--output", type=Path, required=True, help="HDP matrix")
    args = args[1:]
    return parser.parse_args(args)


def get_common_suffix(strings: List[str]) -> str:
    # Reverse each string
    rev_strings = [s[::-1] for s in strings]
    # Find common *prefix* of the reversed strings
    rev_prefix = os.path.commonprefix(rev_strings)
    # Reverse back to get the suffix
    return rev_prefix[::-1]


def write_hdp_matrix(input_path: Path, sample_path: Path, output_path: Path):
    # Load file paths
    file_paths = [line.rstrip() for line in open(input_path).readlines()]

    # Get common suffix
    suffix = get_common_suffix(file_paths[0:5])

    # Load target sampels
    target_samples = set([line.strip() for line in open(sample_path).readlines()])

    # Aggregate counts
    sample_sbs52_counts = {}
    sbs52_classifications = set()
    for file_path in file_paths:
        query_sample = file_path.replace(suffix, "")
        if query_sample not in target_samples:
            continue
        for line in open(file_path).readlines():
            if line.startswith("SUB"):
                continue
            fields = line.rstrip().split("\t")
            sub = fields[0]
            tri = fields[1].split()[1]
            sbs52 = fields[2]
            count = fields[3]
            upstream, _, downstream = list(tri)
            sbs52 = "{},{}-{}".format(sub, upstream, downstream)
            sample_sbs52_counts["{}:{}".format(query_sample, sbs52)] = count
            sbs52_classifications.add(sbs52)
    sbs52_classifications = natsort.natsorted(list(sbs52_classifications))

    # Write matrix
    target_samples = natsort.natsorted(list(target_samples))
    with open(output_path, "w") as outfile:
        outfile.write("{}\n".format("\t".join(sbs52_classifications)))
        for sample in target_samples:
            row_values = [sample] + [
                str(math.ceil(float(sample_sbs52_counts["{}:{}".format(sample, sbs52)])))
                for sbs52 in sbs52_classifications
            ]
            outfile.write("{}\n".format(",".join(row_values)))


def main():
    options = parse_args(sys.argv)
    write_hdp_matrix(options.input, options.samples, options.output)
    sys.exit(0)


if __name__ == "__main__":
    main()
