# trans2genome
tophat-like method for transcriptome-aware spliced alignment using HISAT

The method performs alignment in two stages:
1. Using Bowtie2 or HISAT2 (Bowtie2 by default) to perform unspliced local alignment to transcriptome index
2. Using HISAT2 to align reads that were not mapped to the transcriptome in step #1

This method improves upon TopHat2 by optimizing transcriptome alignment, which reduces computational time
