#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import argparse


def parse_num(s):
    x = float(s)
    return int(x) if x.is_integer() else x


def fmt_num(x):
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x)


def read_rows(path):
    rows = []
    skipped_header = False

    with open(path, "r", encoding="utf-8") as f:
        for lineno, line in enumerate(f, start=1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue

            fields = line.split()
            if len(fields) < 2:
                continue

            gene_id = fields[0]
            try:
                num = parse_num(fields[1])
            except ValueError:
                if lineno == 1 and not skipped_header:
                    skipped_header = True
                    continue
                else:
                    raise ValueError(f"  {lineno}  2 :{fields[1]}")

            rows.append((gene_id, num))

    if not rows:
        raise ValueError(" .")

    return rows


def split_sizes(total, n_bins):
    base = total // n_bins
    rem = total % n_bins
    return [base + (1 if i < rem else 0) for i in range(n_bins)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("n_bins", type=int, help=" , ,  3,5,7")
    parser.add_argument("--input", default="H3K27ac_output_E5_gene_count.txt",
                        help="inputfile,  H3K27ac_output_E5_gene_count.txt")
    args = parser.parse_args()

    n_bins = args.n_bins
    if n_bins < 3:
        raise ValueError("n_bins   >= 3")
    if n_bins % 2 == 0:
        raise ValueError(" “ ”,n_bins  ,  3,5,7")

    rows = read_rows(args.input)

    rows.sort(key=lambda x: (x[1], x[0]))

    sizes = split_sizes(len(rows), n_bins)

    bins = []
    start = 0
    for size in sizes:
        end = start + size
        bins.append(rows[start:end])
        start = end

    low_bin = bins[0]
    mid_bin = bins[(n_bins - 1) // 2]
    high_bin = bins[-1]

    low_upper = low_bin[-1][1]      #  1 ,  <= low_upper
    mid_low = mid_bin[0][1]         #
    mid_high = mid_bin[-1][1]       #
    high_lower = high_bin[0][1]     #  N ,  >= high_lower

    print(f"{fmt_num(high_lower)} {fmt_num(mid_low)} {fmt_num(mid_high)} {fmt_num(low_upper)}")


if __name__ == "__main__":
    main()
