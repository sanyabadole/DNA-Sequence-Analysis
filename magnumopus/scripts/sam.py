import re
from collections import defaultdict
from copy import deepcopy

class Read:
    def __init__(self, sam_line: str):
        (qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, *tags) = sam_line.strip().split("\t")
        self.qname: str = qname
        self.flag: int = int(flag)
        self.rname: str = rname
        self.pos: str = int(pos)
        self.mapq: str = int(mapq)
        self.cigar: str = cigar
        self.rnext: str = rnext
        self.pnext: str = int(pnext)
        self.tlen: str = int(tlen)
        self.seq: str = seq
        self.qual: str = qual
        self.tags: list[str] = tags
        self.is_mapped: bool = not bool(self.flag & 4) # 4 bit not in flag
        self.is_forward: bool = not bool(self.flag & 16)
        self.is_reverse: bool = bool(self.flag & 16)
        self.is_primary: bool = not (bool(self.flag & 256) or bool(self.flag & 2048)) # not secondary or supplemental

        #initializing empty variables to use for mapped reads
        self.cigar_bits: tuple[tuple[int, str]] = None
        self.mapped_len: int = None
        
        #if mapped then add to cigar bits and mapped length
        if self.is_mapped:
            self.cigar_bits = tuple([(int(n), cig) for n, cig in re.findall(r"(\d+)([A-Z])", self.cigar)])
            self.mapped_len = sum([n for n, cig in self.cigar_bits if cig in {"M", "D"}])

    def read_idx_at_pos(self, pos: int) -> list[None|int]:
        if not self.is_mapped:
            return []
        pos -= self.pos
        if pos < 0: 
            return []
        if pos >= self.mapped_len:
            return []
        
        cigar_bits = re.findall(r"(\d+)([A-Z])", self.cigar)
        #seq position based on cigar
        seq = 0
        mapped_count = 0
        for n, (size, cig_type) in enumerate(cigar_bits):
            size = int(size)
            if cig_type in {"S", "H", "I"}:
                seq += size
                continue
            if mapped_count + size >= pos+1:
                if cig_type == "M":
                    seq += (pos-mapped_count)
                    break
                if cig_type == "D":
                    return []
            else:
                if cig_type == "M":
                    mapped_count += size
                    seq += size
                if cig_type == "D":
                    mapped_count += size

        #if insertion then
        if mapped_count + size == pos+1:
            if n+1 != len(cigar_bits):
                size, cig_type = cigar_bits[n+1]
                size = int(size)
                if cig_type == "I":
                    return [i for i in range(seq, seq+size+1)]
        return [seq]

    def mapped_seq(self) -> str:
        if not self.is_mapped:
            return ""

        idx = 0 #what position we are in the sequence
        bases = [] 
        for n, cig in self.cigar_bits:
            if cig == "S":
                idx += n
            elif cig == "D":
                bases += ["-"]*n
            elif cig in {"M", "I"}:
                bases += [self.seq[i] for i in range(idx, idx+n)]
                idx += n
        return "".join(bases)

    def base_at_pos(self, pos: int) -> str:
        idx = self.read_idx_at_pos(pos)
        return "".join([self.seq[i] for i in idx])

    def qual_at_pos(self, pos: int) -> str:
        idx = self.read_idx_at_pos(pos)
        return "".join([self.seq[i] for i in idx])


class SAM:
    def __init__(self):
        self._reads: dict[dict[list[Read]]] = defaultdict(lambda: defaultdict(list))
    
    @property
    def reads(self):
        return deepcopy(self._reads)
    
    @property
    def best_ref(self):
        most_reads = 0
        best = ""
        for k,v in self.reads.items():
            if len(v) > most_reads:
                most_reads = len(v)
                best = k
        return best

    def add_read(self, read: Read):
        if read.is_mapped:
            for i in range(read.pos, read.pos + read.mapped_len):
                self._reads[read.rname][i].append(read)
        else:
            # don't store unmapped reads
            pass

    def reads_at_pos(self, seq_name: str, pos: int):
        return self._reads[seq_name][pos]

    def pileup_at_pos(self, seq_name: str, pos: int) -> tuple[list[str], list[str]]:
        s_pileup = []
        q_pileup = []
        for read in self._reads[seq_name][pos]:
            s_pileup.append(read.base_at_pos(pos))
            q_pileup.append(read.qual_at_pos(pos))
        return (s_pileup, q_pileup)

    def consensus_at_pos(self, seq_name: str, pos: int) -> str:
        if pos not in self.reads[seq_name]:
            return " "
        bases, quals = self.pileup_at_pos(seq_name, pos)
        base_count = defaultdict(int)
        for base, qual in zip(bases, quals):
            base_count[base] += 1
        sorted_base_count = sorted([(b, c) for b, c in base_count.items()], key=lambda x: x[1], reverse=True)
        if len(sorted_base_count) > 1:
            if sorted_base_count[0][1] <= sum([i[1] for i in sorted_base_count])/2:
                top_base = "N"
                return top_base
        top_base = sorted_base_count[0][0]
        return top_base

    def consensus_sequence(self, seq_name: str, start: int=None, end: int=None, fasta=False):
        if fasta:
            consensus = f">{seq_name}_consensus\n"
        else:
            consensus = ""
        if seq_name not in self.reads:
            return consensus+ "\n"
        if start is None or start < 1:
            start = 1
        if end is None or end > max(self._reads[seq_name].keys()):
            end = max(self._reads[seq_name].keys())
        if fasta:
            consensus = f">{seq_name}_consensus\n"
        consensus = "".join([consensus] + [self.consensus_at_pos(seq_name, i) for i in range(start, end + 1)] + ["\n"]) 
        return consensus
    
    def best_consensus(self, fasta=False):
        best = self.best_ref
        return self.consensus_sequence(best, fasta=fasta)

    @classmethod
    def from_sam(cls, sam_file: str):
        sam = cls()
        with open(sam_file) as fin:
            for line in fin:
                if line[0] != "@":
                    #reading the lines in the sam file
                    read = Read(line)
                    sam.add_read(read)
                    continue
        return sam 