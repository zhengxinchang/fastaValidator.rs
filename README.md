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
6. Percentage of unknown nucleotide (N/n) should be less than 50%.
7. Sequence should not be starts or ends with "N" or "n"


# Output

## source table for table2asn

<pre>
seqid   organism        genetic_code    moltype topology        strand
>ssss   Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
>ssss   Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
>eq4    Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
>Beijing-3107-2022      Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
>Beijing-AAA-2022       Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
>asdf   Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
>asdf   Severe acute respiratory syndrome coronavirus 2 1       genomic RNA     linear  single
</pre>

## errors in table format

<pre>
+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Error_type | Message                                                                                                                                                       |
+============+===============================================================================================================================================================+
| Nucleotide | Found invalid 'N' at start of sequence '>ssss'. It should not start with 'N' or 'n'                                                                           |
+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Nucleotide | Found invalid '>' at sequence(seqid:'>Beijing-AAA-2022'). This symbol is not allowed in the sequence. Please check whether the new-line character is missing. |
+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Sequence   | For SARS-CoV-2 submission, sequence length must be between 50 and 30000. SeqLength of '>Beijing-AAA-2022' is 39852                                            |
+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Defline    | Found duplicated sequence id: '>ssss'                                                                                                                         |
+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Defline    | Found duplicated sequence id: '>asdf'                                                                                                                         |
+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+
</pre>