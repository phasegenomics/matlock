# Matlock 

 "I know the prosecutor in this case, and she is a very tough cookie." 
    --Matlock

Matlock is a set of tools designed to work with Hi-C data. Functions include bam filtering, conversion, and clustering.

## installing

```
git clone --recursive https://github.com/phasegenomics/matlock.git matlock ; cd matlock ; make
```

## usage
```
usage: matlock <command> [options]


commands:
 - bam2 - converts alignments to several useful hi-c formats.
   + usage: matlock bam2 [binmat|lachesis|juicer|counts] input output
   + details:
      The input file format is automatically determined [cram|bam|sam].
      The output is written to the fineame provided, no extention.


 - bamfilt - filter a hi-c bam.
   + usage: matlock bamfilt input.[cram|bam|sam] output.bam


 - motif - count motifs.
   + usage: matlock motif input.fasta ATGC TGCA ...
```
