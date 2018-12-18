#!/bin/sh -x

#conda install art
art_illumina -ss HS20 -p -f 30 -l 100 -m 200 -s 10 -nf 0 -na -i Sample.fa -o Sample
head -n 4000 Sample1.fq > Sample_1.fastq
head -n 4000 Sample2.fq > Sample_2.fastq
