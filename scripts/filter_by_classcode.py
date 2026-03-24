#!/usr/bin/env python3
"""
Filter an annotated GTF (from gffcompare) to keep only transcripts with class codes i, u, x.
Usage: filter_by_classcode.py annotated.gtf output.gtf
"""
import sys
import re

def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: filter_by_classcode.py annotated.gtf output.gtf")
    annotated_gtf = sys.argv[1]
    out_file = sys.argv[2]

    keep_transcripts = set()
    # First pass: find transcript IDs with desired class codes
    with open(annotated_gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            attr = fields[8]
            # Extract class_code
            m = re.search(r'class_code "([^"]+)"', attr)
            if m and m.group(1) in ('i', 'u', 'x'):
                # Extract transcript_id
                t = re.search(r'transcript_id "([^"]+)"', attr)
                if t:
                    keep_transcripts.add(t.group(1))

    # Second pass: write all lines belonging to those transcripts
    with open(annotated_gtf) as fin, open(out_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            attr = fields[8]
            t = re.search(r'transcript_id "([^"]+)"', attr)
            if t and t.group(1) in keep_transcripts:
                fout.write(line)

if __name__ == "__main__":
    main()
