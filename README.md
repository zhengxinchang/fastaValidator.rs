# Introduction

fastaValitator is a fast validator for FASTA files. It is implemented by Rust and used in [GenBase](https://ngdc.cncb.ac.cn/genbase/) SARS-CoV-2 submission pipeline.

# System requirements

This program can run at any Linux distribution theoretically, but we only test it on Ubuntu(WSL) and CentOS 7.

# Validators

1. Unique Sequence ID.
2. Sequence_ID's may contain only the following characters - letters, digits, hyphens "-", underscores "_".
3. Max length of Sequence ID is 23.
4. Sequence length must be 50 - 30,000 bases, and <50% Ns, excluding terminal Ns which will be trimmed automatically, and the inside of the sequence cannot contain horizontal lines "-" (some of which are inserted after the software is run).
5. Only A,T,C,G and degenerative nucleotides are permitted.