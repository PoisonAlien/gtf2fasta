# gtfsanta
Tookit for [gtf](http://mblab.wustl.edu/GTF22.html) format conversion and sequence extraction. Still in development.

##Description

A set of functions for gtf file manipulation. It's written in [Julia](https://julialang.org/) and has no dependencies except for Julia itself.
Julia a fast dynamic programming language which offers speeds on par with C and Fortran. See [here](https://julialang.org/downloads/) for Julia installation.

This tool requires [indexed fasta file](http://www.htslib.org/doc/faidx.html) for memory efficient sequence extraction.

##Usage

```bash
#Extracting fasta sequence from gtf file

$ gtfsanta.jl gtf2fasta ens82.gtf hg19.fa | head
>ENST00000456328
GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCA
GACTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCA
GGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGGAGTTTTCCTG
...

#Extracting fasta sequence from bed file

$ gtfsanta.jl bed2fasta CEBPe_peaks.bed mm10.fa | head
>chr1:4769907-4770192
AGTGTCTGAACCTCTGAGGCAACTGCCTTTCTAATTTCAATGTGTCTAGCATTGTGCCTG
TGTGTGTGCCCTTGAGAACTGTCGATTGTAGAAACTCATGAATACTTAATATTATTAATC
AAACCTGGAATTTCCCAATGCTTGCATTTGGCCACCAGGGGGCAGTCTTGCCAATTTGAA
ATGCCGATTTCGTTTTTGGGAGACTCCAATTCTTCACCTTTGAGCATCTCAAACGTTACT
TCTTTTAGAACTTTTCACTTGAAAGTATTTCAAAGGTAACTCCTG
>chr1:4773351-4773679
AGTCTATCATGCATTTTCTTCTCATTAGTGAGACTCAAGAGTGCTTTACAGACTATATGT
AAATCCTAAGAGGGAAGGAAGTAGAATCCACAGGGTGGCAGTTACCTAAGGATCAGAAGA
...
```

##To do

- [x] gtf2fasta
- [x] bed2fasta
- [ ] gtf2bed6
- [ ] gtf2bed12
- [ ] bed2gtf
- [ ] gtf2gff
- [ ] txTable
