# MagnumOpus - DNA Sequence Analysis Tool

MagnumOpus is a comprehensive bioinformatics tool designed for DNA sequence analysis, with a particular focus on PCR primer analysis and sequence alignment. This tool combines multiple bioinformatics techniques to provide a robust solution for DNA sequence analysis tasks.

## Features

### 1. PCR Primer Analysis
- Analyzes PCR primers against DNA assemblies
- Identifies potential amplicons within specified size limits
- Supports both single and paired-end read analysis
- Handles multiple assembly files simultaneously

### 2. Sequence Alignment
- Implements Needleman-Wunsch algorithm for global sequence alignment
- Supports both forward and reverse complement sequence alignment
- Configurable scoring parameters (match, mismatch, gap penalties)

### 3. Read Mapping
- Maps sequencing reads to reference sequences using minimap2
- Generates SAM format output
- Supports paired-end read mapping
- Includes consensus sequence generation

### 4. File Format Support
- FASTA file parsing and manipulation
- SAM/BAM file handling
- FASTQ read file processing
- Support for multiple input formats

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/magnumopus.git
cd magnumopus
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

3. Make sure you have the following external tools installed:
- minimap2 (for read mapping)
- samtools (for SAM/BAM file manipulation)
- blastn (for primer analysis)

## Usage

The main script `magnop.py` can be used with the following arguments:

```bash
python magop.py -p primers.fasta [-a assemblies/*.fasta] [-r reads/*.fastq] [-s reference.fasta]
```

### Arguments:
- `-p, --primers`: Path to primers file (required)
- `-a, --assemblies`: Path to assembly files (optional)
- `-r, --reads`: Path to read files (optional)
- `-s, --reference`: Path to reference sequence (optional)

### Example Usage

1. Analyze PCR primers against assemblies:
```bash
python magop.py -p primers.fasta -a assemblies/*.fasta
```

2. Map reads to reference and analyze:
```bash
python magop.py -p primers.fasta -r reads/*.fastq -s reference.fasta
```

3. Combined analysis:
```bash
python magop.py -p primers.fasta -a assemblies/*.fasta -r reads/*.fastq -s reference.fasta
```

## Project Structure

```
magnumopus/
├── magnumopus/
│   ├── __init__.py
│   ├── magop.py
│   └── scripts/
│       ├── __init__.py
│       ├── sam.py      # SAM file handling
│       ├── nw.py       # Needleman-Wunsch alignment
│       └── ispcr.py    # PCR primer analysis
├── data/
├── tests/
├── requirements.txt
├── setup.py
├── LICENSE
└── README.md
```

## Dependencies

### Python Packages
- numpy (>=1.21.0)
- pandas (>=1.3.0)
- biopython (>=1.79)
- setuptools (>=45.0.0)
- wheel (>=0.37.0)

### External Tools
- minimap2
- samtools
- blastn

## Output

The tool generates:
1. Aligned amplicon sequences
2. Consensus sequences from mapped reads
3. PCR primer analysis results
4. Sequence alignment scores and alignments

## Error Handling

The tool includes comprehensive error handling for:
- Missing input files
- Invalid file formats
- Failed external tool executions
- Incomplete read pairs
- Missing dependencies

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## Acknowledgments

- Uses minimap2 for read mapping
- Implements Needleman-Wunsch algorithm for sequence alignment
- Utilizes BLAST for primer analysis 
