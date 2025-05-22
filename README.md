# Infernommatic

This script was developed to address a problem at the SetuLab (University of Paulo), in which reads were still contaminated with adapters even after running multiple read trimming tools. This might be due to wetlab problems.  
Here we try to remove those fragmented adapters from many positions in the read sequence!   

## Installation

```bash
pip install git+https://github.com/FefoRossi/Infernommatic
```

## Basic usage

### Commands
```
usage: infernommatic [-h] -i INPUT_R1 -I INPUT_R2 [-adapt ADAPTERS_FILE] -o OUTPUT_R1 -O OUTPUT_R2 [-pct SEQUENCE_PCT]

Script for trimming reads with bad linked adapters

options:
  -h, --help            show this help message and exit
  -i INPUT_R1, --input_r1 INPUT_R1
                        R1 fastq file -- either gzipped or not
  -I INPUT_R2, --input_r2 INPUT_R2
                        R2 fastq file -- either gzipped or not
  -adapt ADAPTERS_FILE, --adapters_file ADAPTERS_FILE
                        Fasta file with all adapters to be processed
  -o OUTPUT_R1, --output_r1 OUTPUT_R1
                        R1 fastq output file path
  -O OUTPUT_R2, --output_r2 OUTPUT_R2
                        R2 fastq output file path
  -pct SEQUENCE_PCT, --sequence_pct SEQUENCE_PCT
                        percent_threshold': If an adaptor (or part) is within this percentage from an end, it's considered at the start/end of the sequence.

Ex. usage: infernommatic.py -i R1.fastq.gz -I R2.fastq.gz -o filtered_R1.fastq -O filtered_R2.fastq
```
### Basic run

```bash
infernommatic -i INPUT_R1 -I INPUT_R2 -o OUTPUT_R1 -O OUTPUT_R2
```

>By default, a list with all trimommatic adapters is available, for custom adapters, please provide by: `-adapt path/to/adapter.fasta`  

