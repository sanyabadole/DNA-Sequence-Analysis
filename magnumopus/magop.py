#!/usr/bin/env python3

import argparse
import glob
import re
import subprocess
import os
import sys
from tempfile import TemporaryDirectory
from .scripts.sam import SAM, Read
from .scripts import ispcr
from .scripts import nw

def parse_arguments():
    parser = argparse.ArgumentParser(description="DNA sequence analysis tool for PCR primer analysis and sequence alignment")
    parser.add_argument("-a", "--assemblies", nargs='*', help="Path to assemblies", required=False)
    parser.add_argument("-p", "--primers", help="Path to primers", required=True)
    parser.add_argument("-r", "--reads", nargs='*', help="Path to reads", required=False)
    parser.add_argument("-s", "--reference", help="Path to reference", required=False)
    return parser.parse_args()

def pair_reads(read_files):
    read_pairs = {}
    for read in read_files:
        match = re.search(r"(.+)_\d.fastq", read)
        if not match:
            print(f"Warning: Could not parse read file name: {read}", file=sys.stderr)
            continue
        sample_id = match.group(1)
        if sample_id not in read_pairs:
            read_pairs[sample_id] = []
        read_pairs[sample_id].append(read)
    return read_pairs

def run_external(command: list[str], shell=False, stdin=None) -> tuple[str, str]:
    try:
        if stdin is None:
            result = subprocess.run(command, capture_output=True, text=True, shell=shell)
        else:
            result = subprocess.run(command, capture_output=True, text=True, input=stdin, shell=shell)
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {' '.join(command)}", file=sys.stderr)
        print(f"Error message: {e.stderr}", file=sys.stderr)
        sys.exit(1)

def map_reads_to_ref(ref: str, r1: str, r2: str) -> SAM:
    if not os.path.exists(ref):
        print(f"Error: Reference file {ref} does not exist", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(r1) or not os.path.exists(r2):
        print(f"Error: Read files {r1} or {r2} do not exist", file=sys.stderr)
        sys.exit(1)

    mapping_args = ["minimap2", "-ax", "sr", "-k", "10", "-B", "0"]
    mapping_args += [ref, r1, r2]
    map_cmd = " ".join([str(a) for a in mapping_args])
    map_cmd += "| samtools view -hF4"

    with TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "tmp.sam")
        map_cmd += f" > {outfile}"
        run_external(map_cmd, shell=True)
        sam = SAM.from_sam(outfile)
    return sam

def parse_fasta(filename):
    if not os.path.exists(filename):
        print(f"Error: Fasta file {filename} does not exist", file=sys.stderr)
        sys.exit(1)
        
    sequences = {}
    current_id = None
    current_seq = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = '\n'.join(current_seq)
                current_id = line
                current_seq = [line]
            else:
                current_seq.append(line)
    if current_id:
        sequences[current_id] = '\n'.join(current_seq)
    return sequences

def main():
    args = parse_arguments()
    amps_list = []
    
    if args.assemblies:
        assembly_files = []
        for item in args.assemblies:
            if '*' in item:
                assembly_files.extend(glob.glob(item))
            else:
                assembly_files.append(item)
                
        if not assembly_files:
            print("Error: No assembly files found", file=sys.stderr)
            sys.exit(1)
            
        for assembly in assembly_files:
            if not os.path.exists(assembly):
                print(f"Warning: Assembly file {assembly} does not exist", file=sys.stderr)
                continue
            amps = ispcr.ispcr(args.primers, assembly, 2000)
            amps_list.append(amps.split('\n')[0] + '\n' + amps.split('\n')[1] + '\n')
            
    if args.reads:
        read_files = []
        for item in args.reads:
            if '*' in item:
                read_files.extend(glob.glob(item))
            else:
                read_files.append(item)
                
        if not read_files:
            print("Error: No read files found", file=sys.stderr)
            sys.exit(1)
            
        read_pairs = pair_reads(read_files)
        
        if not read_pairs:
            print("Error: No valid read pairs found", file=sys.stderr)
            sys.exit(1)
            
        for pair in read_pairs.values():
            if len(pair) != 2:
                print(f"Warning: Skipping incomplete read pair: {pair}", file=sys.stderr)
                continue
                
            read1, read2 = pair
            sam1 = map_reads_to_ref(
                ref=args.reference,
                r1=read1,
                r2=read2
            )
            with open('temp.fasta', 'w') as file:
                file.write(sam1.best_consensus(fasta=True))
            seq = 'temp.fasta'
            amps_list.append(ispcr.ispcr(args.primers, seq, 2000))
            os.remove(seq)

        if args.reference:
            ref = parse_fasta(args.reference)
            for seq in ref.values():
                with open('temp.fasta', 'w') as file:
                    file.write(seq)
                fasta = 'temp.fasta'
                amps_list.append(ispcr.ispcr(args.primers, fasta, 2000))
                os.remove(fasta)

    if not amps_list:
        print("Error: No amplicons were generated", file=sys.stderr)
        sys.exit(1)

    anchor = amps_list[0].split('\n')[1]
    final_amplicons = []
    for amp in amps_list:
        amplicon_2 = amp.split('\n')[1]
    
        match = 1
        mismatch = -1
        gap = -1
        aln, score = nw.needleman_wunsch(anchor, amplicon_2, match, mismatch, gap)
        base_complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        amplicon_2_rc = ''.join([base_complements[nuc] for nuc in reversed(amplicon_2)])
        aln_rc, score_rc = nw.needleman_wunsch(anchor, amplicon_2_rc, match, mismatch, gap)
        if score_rc > score:
            score = score_rc
            final_amplicons.append(amp.split('\n')[0] + '\n' + amplicon_2_rc + '\n')
        else:
            final_amplicons.append(amp)
    
    for amp in final_amplicons:
        print(amp)

if __name__ == "__main__":
    main() 