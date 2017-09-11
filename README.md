### Introduction

Bioawk is an extension to [Brian Kernighan's awk][1], adding the support of
several common biological data formats, including optionally gzip'ed BED, GFF,
SAM, VCF, FASTA/Q and TAB-delimited formats with column names. It also adds a
few built-in functions and an command line option to use TAB as the
input/output delimiter. When the new functionality is not used, bioawk is
intended to behave exactly the same as the original BWK awk.

The original awk requires a YACC-compatible parser generator (e.g. Byacc or
Bison). Bioawk further depends on [zlib][zlib] so as to work with gzip'd files.

### New functionality

##### Command line option `-t`

Using this option is equivalent to

    bioawk -F'\t' -v OFS="\t"

##### Command line option `-c arg`

This option specifies the input format. When this option is in use, bioawk will
seamlessly add variables that name the fields, based on either the format or
the first line of the input, depending *arg*. This option also enables bioawk
to read gzip'd files. The argument *arg* may take the following values:

* `help`. List the supported formats and the naming variables.

* `hdr` or `header`. Name each column based on the first line in the input.
  Special characters in the first are converted to underscore. For example:

        grep -v ^## in.vcf | bioawk -tc hdr '{print $_CHROM,$POS}'

  prints the `CHROM` and `POS` columns of the input VCF file.

* `sam`, `vcf`, `bed` and `gff`. SAM, VCF, BED and GFF formats.

* `fastx`. This option regards a FASTA or FASTQ as a TAB delimited file with
  four columns: sequence name, sequence, quality and FASTA/Q comment, such that
  various fields can be retrieved with column names. See also example 4 in the
  following.

##### New built-in functions

See `awk.1`.

### Examples

1. List the supported formats:

        bioawk -c help

2. Extract unmapped reads without header:

        bioawk -c sam 'and($flag,4)' aln.sam.gz

3. Extract mapped reads with header:

        bioawk -Hc sam '!and($flag,4)'

4. Reverse complement FASTA:

        bioawk -c fastx '{print ">"$name;print revcomp($seq)}' seq.fa.gz

5. Create FASTA from SAM (uses revcomp if FLAG & 16)

        samtools view aln.bam | \
            bioawk -c sam '{s=$seq; if(and($flag, 16)) {s=revcomp($seq)} print ">"$qname"\n"s}'

6. Print the genotypes of sample `foo` and `bar` from a VCF:

        grep -v ^## in.vcf | bioawk -tc hdr '{print $foo,$bar}'


### Potential limitations

1. When option `-c` is in use, bioawk replaces the line reading module of awk.
   The new line reading function parses FASTA and FASTQ files and seamlessly
   reads gzip'ed files. However, the new code does not fully mimic the original
   code. It may fail in corner cases (though this has not happened yet). Thus
   when `-c` is not specified, awk falls back to the original line reading code
   and does not support gzip'ed input.

2. When `-c` is in use, several strings allocated in the new line reading
   module are not freed in the end. These will be reported by valgrind as
   "still reachable". To some extent, these are not memory leaks.


[1]: http://www.cs.princeton.edu/~bwk/btl.mirror/
[zlib]: http://zlib.net
