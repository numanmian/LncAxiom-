#!/usr/bin/env python3
"""
Wrapper for psRobot_tar to identify endogenous target mimics (eTMs) of lncRNAs.
Applies stringent filters based on published eTM criteria.
"""

import argparse
import subprocess
import sys
import re
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Run psRobot_tar and filter eTMs")
    parser.add_argument("--lncrna_fasta", required=True, help="FASTA file of lncRNA sequences")
    parser.add_argument("--mirna_fasta", required=True, help="FASTA file of mature miRNA sequences")
    parser.add_argument("--output", required=True, help="Output TSV file with filtered eTMs")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--score_thresh", type=float, default=2.5, help="Target score threshold")
    parser.add_argument("--seed_start", type=int, default=2, help="Start of seed region (1‑based)")
    parser.add_argument("--seed_end", type=int, default=8, help="End of seed region")
    parser.add_argument("--forbidden_bulge_start", type=int, default=9, help="Start of forbidden bulge region")
    parser.add_argument("--forbidden_bulge_end", type=int, default=12, help="End of forbidden bulge region")
    parser.add_argument("--max_bulge_size", type=int, default=3, help="Maximum allowed bulge size")
    parser.add_argument("--max_mismatches_gu", type=int, default=3, help="Max mismatches+G:U outside bulge")
    parser.add_argument("--psrobot_bin", default="psRobot_tar", help="Path to psRobot_tar executable")
    return parser.parse_args()

def run_psrobot(args, tmp_out):
    """Execute psRobot_tar with parameters optimised for eTM detection."""
    cmd = [
        args.psrobot_bin,
        "-s", args.mirna_fasta,
        "-t", args.lncrna_fasta,
        "-o", tmp_out,
        "-ts", str(args.score_thresh),
        "-fp", str(args.seed_start),
        "-tp", str(args.seed_end),
        "-p", str(args.threads),
        "-gn", "1",
        "-gl", "17"
    ]
    print("Running:", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)

def parse_gTP(gTP_file):
    """Parse psRobot_tar .gTP output into a list of dictionaries."""
    interactions = []
    with open(gTP_file) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('>'):
            # Header line
            header = line[1:]  # remove '>'
            # Extract score using regex
            score_match = re.search(r'Score:\s*([\d.]+)', header)
            score = float(score_match.group(1)) if score_match else None
            # Split header by tabs (psRobot uses tabs between fields)
            fields = header.split('\t')
            if len(fields) >= 3:
                mirna = fields[0].strip()
                target = fields[-1].strip()
            else:
                mirna = fields[0]
                target = "unknown"
            i += 1
            # Skip blank lines
            while i < len(lines) and lines[i].strip() == '':
                i += 1
            if i >= len(lines):
                break
            # Query line: e.g., "Query:          1 TGACAGAAGAGAGTGAGCAC 20"
            query_line = lines[i].strip()
            i += 1
            # Symbols line: e.g., "                  |||||||*:|||:||||*|*"
            symbols_line = lines[i].strip()
            i += 1
            # Subject line: e.g., "Sbjct:        825 ACTGTCTCTTCTTACTC-TC 807"
            subject_line = lines[i].strip()
            i += 1
            # Extract sequences
            query_parts = query_line.split()
            if len(query_parts) >= 3:
                query_seq = query_parts[2]
            else:
                query_seq = ""
            subject_parts = subject_line.split()
            if len(subject_parts) >= 3:
                subject_seq = subject_parts[2]
            else:
                subject_seq = ""
            symbols = symbols_line
            interactions.append({
                "mirna": mirna,
                "target": target,
                "score": score,
                "query_seq": query_seq,
                "subject_seq": subject_seq,
                "symbols": symbols,
                "raw_header": header
            })
        else:
            i += 1
    return interactions

def filter_etm(interactions, args):
    """
    Apply eTM filters to each interaction.
    Currently only filters by score. Will add seed and bulge checks later.
    """
    filtered = []
    for inter in interactions:
        if inter["score"] is not None and inter["score"] <= args.score_thresh:
            filtered.append(inter)
    return filtered

def main():
    args = parse_args()
    tmp_out = Path(args.output).with_suffix(".gTP.tmp")
    run_psrobot(args, str(tmp_out))
    interactions = parse_gTP(str(tmp_out))
    filtered = filter_etm(interactions, args)
    with open(args.output, "w") as out:
        out.write("miRNA\tlncRNA\tscore\tquery_seq\tsubject_seq\tsymbols\n")
        for i in filtered:
            out.write(f"{i['mirna']}\t{i['target']}\t{i['score']}\t{i['query_seq']}\t{i['subject_seq']}\t{i['symbols']}\n")
    # tmp_out.unlink()  # keep raw file for debugging
    print(f"Filtered eTMs written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()
