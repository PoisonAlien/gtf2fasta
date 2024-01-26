# gtf2fasta
A minimal tool for sequence extraction for every transcript (only from exonic regions) in a [gtf](http://mblab.wustl.edu/GTF22.html) from a fasta file. 

It's written in [Julia](https://julialang.org/) and has no dependencies. See [here](https://julialang.org/downloads/) for Julia installation.

This tool requires [indexed fasta file](http://www.htslib.org/doc/faidx.html) for memory efficient sequence extraction.

## Usage

```bash
#Extracting fasta sequence from gtf file (sequences are written to stdout)

$ gtf2fasta.jl ens82.gtf hg19.fa | head
>ENST00000456328
GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCA
GACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCA
GGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGGAGTTTTCCTG
...

It will also generates a `ens82.gtf.transcript.dict.tsv` with transcript to gene mappings.
