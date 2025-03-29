#!/usr/bin/env python3

import subprocess
import os
import sys
import tempfile
from collections import defaultdict

def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    hits = find_annealing(primer_file, assembly_file)
    good_hits = filter_blast(hits)
    sorted_hits = sorted(good_hits, key=lambda x: (x[1], int(x[8])))
    return sorted_hits

def step_two(
    sorted_hits: list[str],
    max_amplicon_size: int
) -> list[tuple[list[str]]]:
    hit_pairs = identify_paired_hits(sorted_hits, max_amplicon_size)
    return hit_pairs

def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    amplicons = get_amplicons(assembly_file, hit_pairs)
    return amplicons

def find_annealing(primer_file: str, assembly_file: str) -> list[str]:
    if not os.path.exists(primer_file):
        print(f"Error: Primer file {primer_file} does not exist", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(assembly_file):
        print(f"Error: Assembly file {assembly_file} does not exist", file=sys.stderr)
        sys.exit(1)

    blast_command = ["blastn"]
    blast_command += ["-query", primer_file]
    blast_command += ["-subject", assembly_file]
    blast_command += ["-task", "blastn-short"]
    blast_command += ["-word_size", "6"]
    blast_command += ["-penalty=-2"]
    blast_command += ["-outfmt", "6 std qlen"]
    blast_result, _ = run_external(blast_command)
    return blast_result

def run_external(command: list[str], stdin=None) -> tuple[str, str]:
    try:
        if stdin is None:
            result = subprocess.run(command, capture_output=True, text=True)
        else:
            result = subprocess.run(command, capture_output=True, text=True, input=stdin)
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(command)}", file=sys.stderr)
        print(f"Error message: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def filter_blast(blast_output: str) -> list[list[str]]:
    good_hits = []
    for line in blast_output.split("\n"):
        if len(line) == 0:
            continue
        cols = line.split()
        if cols[3] != cols[12]:
            continue
        good_hits.append(cols)
    return good_hits

def identify_paired_hits(
    good_hits: list[list[str]],
    max_amp_size: int
) -> list[tuple[list[str]]]:
    pairs = []
    for i in range(len(good_hits)-1):
        a_hit = good_hits[i]
        a_start, a_stop = [int(i) for i in a_hit[8:10]]
        a_dir = "fwd" if a_start < a_stop else "rev" #whether it is forward or reverse
        for j in range(i+1, len(good_hits)):
            b_hit = good_hits[j]
            if a_hit[1] != b_hit[1]: # if theey are different contigs
                continue

            b_start, b_stop = [int(i) for i in b_hit[8:10]]
            b_dir = "fwd" if b_start < b_stop else "rev"

            if a_dir == b_dir: #if they are in the same direction
                continue

            if a_start < b_start: #if the first sequence starts before the other sequence
                if not a_dir == "fwd": #if they are in the opposite direction
                    continue
                if not a_stop > b_start - max_amp_size: 
                    break 

                pairs.append((a_hit, b_hit))
                continue

            if not b_dir == "fwd": #if its in thw wrong direction
                continue
            if not b_stop > a_start - max_amp_size: #if its out of range
                continue

            pairs.append((a_hit, b_hit))
    return pairs

def get_amplicons(
    assembly: str,
    hit_pairs: list[tuple[list[str]]]
) -> str:
    amplicons = []
    bed = []
    for f_hit, r_hit in hit_pairs:
        contig = f_hit[1]
        start = int(f_hit[9])
        end = int(r_hit[9])-1

        bed.append(f"{contig}\t{start}\t{end}\n")

    bed_string = "".join(bed)

    #using stdin to run seqtk command
    seqtk_command = ["seqtk", "subseq", assembly]
    seqtk_command += ["-"]
    amplicons, stderr = run_external(seqtk_command, stdin=bed_string)
    
    return amplicons

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    step1 = step_one(primer_file, assembly_file)
    step2 = step_two(step1, max_amplicon_size)
    step3 = step_three(step2, assembly_file)
    return step3 