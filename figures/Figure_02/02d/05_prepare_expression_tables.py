#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
import sys

with open(sys.argv[1]) as f:
    samples = [x.strip("'") for x in f.readline().split()[3:]]
    print(','.join(sorted(samples)))
    for line in f:
        cols = line.split()
        name = (cols[0] + ":" + cols[1] + "-" + cols[2])
        line = {}
        for i, num in enumerate(cols[3:]):
            line[samples[i]] = num
        print(','.join([name] + [line[s] for s in sorted(samples)]))
