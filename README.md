# T2G - conversion of transcriptomic alignments to genomic
T2G is a fast and efficient tool for conversion of transcriptomic alignments to genomic coordinates. 
T2G is developed to perform a variety of complementing and filtering tasks on transcriptomic alignments 
in order to reduce the inherent bias towards the provided annotation and improve both sensitivity and 
precision of the alignment results.

## Installation
In most cases the following T2G should can be isntalled from source with the following steps:
1. cmake -DCMAKE_BUILD_TYPE=Release .
2. make
3. make install

In order to use the transHISAT2.py several python packages may need to be installed on your system
for python3. The packages can be installed using the requirements.txt file located in the root directory.

## trans2genome alignment parser

#### Multimapper resolution
T2G can be used to process all possible multimappers to a given read. This is achieved by extracting
coordinates of all identical sequences of length n from the transcriptome during the index building step.
T2G can then use the information to retrieve a set of possible coordinates of reads in the alignment.
This step is supplementary to the commonly used read-aligners such as Bowtie2 or HISAT2 which are not
designed to search for all occurrences of a given read.

###### Multimapper reolution strategies
T2G provides several strategies for outputting multimappers
1. `-l` flag instructs T2G to output all multimappers regardless of their location or
 the number of occurences
2. `-a` argument allows the user to allocate multimappers by the likelihood of the read having come
 from a given transcript or locus. This option requires a user input of transcript abundances, such as
 a SALMON (Supported) or Kallisto (Experimental) report
3. `-k` argument allows specifying the ceil of the number of multimappers to report for each mapping
 (Experimental)
4. default -in default mode T2G will cycle through `-e` (percent reads) or `-n` (numner of reads) and use
 non-multimapping reads to compute starting transcript abundance estimates. Abundance estimates will then
 be used to assign multimappers based on likelihood. This computation is done largely in real-time, saving
 a lot of computation time by relying on a PID controller with a variable threshold to stabilize the values across
 multimapping loci. Note that in order to disable real-time abundance-based likelihood evaluation
 through the PID controller, a 2-pass over the data will be required and can be set with -e 100.
 
#### Misalignment detection
T2G can be used to screen for and reduce the impact of overfitting done by the mapping to a subset of a 
 genome.
 
`-s` flag instructs T2G to precompute a distribution of errors in the reads and use it to check whether
 a given read is an outlier or not. The error distribution is computed based on the edit distance if the
 `NM` tag is set in the record, or, alternatively, based on the `MD` tag.

`-t` argument can be set to detect outliers based on the number of standard deviations within the Poisson
distribution

`-r` argument can be set to detect outliers based on a strict threshold of raw edit distance as set by
 the user.
 
#### Default settings
 Besides the user-controlled parameters, T2G will automatically attempt at standardizing the alignments
 in a number of ways
 - T2G will automatically try to repair discordant alignments by checking whether two mates map back to
 the same locus with the expected fragment length
 - T2G will automatically collapse transcriptomic alignments which cover the same set of genomic coordinates.
 As such, it is not necessary for the base aligner to search for multimappers, and such tools can be instructed
 to search for and report only 1 alignment (`-k 1` in Bowtie2 and HISAT2, etc.).
 
#### Inputs
Currently T2G supports three types of inputs:
1. SAM
2. BAM
3. piped SAM - note that in this mode precomputation can only be performed using `-n`
 since the total number of reads in the sample is unknown. This changes the internal workings
 of the protocol with little effect on the final results and runtime.
 
#### Help Page: `trans2genome --help`
Trans2Genome Help Page

```
salmon2genome -o -x [-a -e -f -g -i -j -k -l -m -n -p -q -r -s -t -u ]
Arguments:
	-a/--abund	        use abundances precomputed. Only gene-level abundance estimates computed by salmon are supported at the moment
	-e/--perc	        percent of the reads to be evaluated when pre-loading data before parsing
	-f/--fraglen	        fragment length of the paired reads
	-g/--nosingle	        report all mates for which the second mate is unaligned as unaligned. If --un is specified, the reads will be deposited into respective fasta/fastq files. If not, SAM information will be set to indicate unaligned
	-i/--input	        input alignment SAM/BAM
	-j/--nodiscord	        report all dicscordant pairs as unaligned. If --un is specified, the reads will be deposited into respective fasta/fastq files. If not, SAM information will be set to indicate unaligned
	-k/--nmult	        The number of most likely multimappers to report
	-l/--all	        whether to output all multimappers and not assign them based on likelihood. This flag negates -k
	-m/--multi	        whether to search and evaluate multimappers
	-n/--nrp	        number of reads to precompute based on the streaming
	-o/--output	        output file path (BAM)
	-p/--threads	        number of threads (default = 1)
	-q/--uniq	        input alignment contains only 1 mapping per read (no secondary alignments present such as in bowtie k1 mode)
	-r/--outlier_raw	Edit distance threshold for misalignments. For paired alignments the threshold is applied to each read separately
	-s/--mis	        try to eliminate misaligned reads based on the error distribution
	-t/--outlier_stdv	Poisson threshold as the number of standard deviations. Everythin above the threshold will be discarded as a misalignment
	-u/--un	                output unaligned reads (singletons and pairs) to separate fastq.
	-x/--index	        path and basename of the index built by gtf_to_fasta
```
	
## Building the Index for T2G with trans2genome_indexer
A T2G index is comprised of a collection of files which describe several relationships between
 transcriptome and genome as well as within transcriptome and within genome. The index contains:
1. `<basename>.fasta` - FASTA-formatted sequences of all isoforms annotated in the provided annotation
2. `<basename>.genome_header` - SAM header of the genome for which the annotation is provided
3. `<basename>.glst` - list and genomic coordinates of all loci in the annotation
4. `<basename>.info` - auxilary data about the input data
5. `<basename>.multi` - optional file computed if `-m` flag is enabled. The file describes sets of unique genomic
coordinates which are share identical sequence within the transcriptome
6. `<basename>.tgmap` - relationship between transcripts and loci as observed from the annotation
7. `<basename>.tlst` - exon coordinates of all transcripts in the database

#### Help Page: `trans2genome --help`
```
   trans2genome_indexer -a -o -r [-k -m -u ]
   Arguments:
   	-a/--gff	    path to the annotation of the genome is GFF/GTF format
   	-k/--kmer	    kmer length to use for building the index
   	-m/--multi	    identify all multimappers present in the input transcriptome
   	-o/--output	    base name for the output files
   	-r/--ref	    path to the reference genome in the FASTA format
   	-u/--uniq	    get a separate output with uniq kmers per each transcript
```

## TransHISAT2.py - automated alignment protocol
A transHISAT2 protocol is also provided for the convenience of end-users. 
The protocol aims to combine all steps of a transcriptomic alignment by reference, 
emmulating in part Tophat2 while taking advantage of the faster and more efficient HISAT2 software 
for the denovo-discovery of spliced alignments.

#### Locus alignment
transHISAT2.py protocol makes it possible to force discovery of intron-retention events for the
 known transcriptome. Users may choose to create a locus-index, and use it with the `--locus` flag
 in `transHISAT2.py align` forcing alignment of reads that were unmapped to the transcriptome onto the
 introns of the annotated loci.
 
 `extractLocus.py` script is provided to aid users in the extraction of correct information
 
 ```
usage: extractLocus.py [-h] --input INPUT --faidx FAIDX --out OUT

Help Page

optional arguments:
  -h, --help     show this help message and exit
  --input INPUT  GFF or GTF annotation of the genome
  --faidx FAIDX  fasta index of the reference genome
  --out OUT      output file
```

#### Help Page: `transHISAT2.py align --help`
```
usage: transHISAT2.py align [-h] --m1 M1 --m2 M2 [--single SINGLE] [--fasta]
                            --db DB -o OUTPUT [--type {hisat,bowtie}]
                            [--tmp TMP] --genome-db GENOME_DB
                            [--threads THREADS] [-a] [-k K]
                            [--bowtie [BOWTIE [BOWTIE ...]]]
                            [--hisat [HISAT [HISAT ...]]] [--keep] [--mf]
                            [--locus] [--abunds ABUNDS] [--errcheck]
                            [--nounal] [-v] [--log] [--no-discord]
                            [--no-single] [--sleep SLEEP]

optional arguments:
  -h, --help            show this help message and exit
  --m1 M1               File with #1 reads of the pair
  --m2 M2               File with #2 reads of the pair
  --single SINGLE       File with singletons
  --fasta               Reads are in fasta format
  --db DB               path to the database directory
  -o OUTPUT, --output OUTPUT
                        output BAM file
  --type {hisat,bowtie}
                        which aligner to use
  --tmp TMP             directory for tmp files
  --genome-db GENOME_DB
                        path to the genome database for HISAT2
  --threads THREADS     number of threads to use
  -a, --all             run hisat/bowtie with the -a option
  -k K                  -k argument for bowtie/hisat first pass
  --bowtie [BOWTIE [BOWTIE ...]]
                        additional arguments to be passed over to bowtie
  --hisat [HISAT [HISAT ...]]
                        additional arguments to be passed over to hisat
  --keep                keep temporary files
  --mf                  use multimapper index to supplement alignments
  --locus               perform second bowtie run to align against loci -
                        parsimonious preMRNA mode
  --abunds ABUNDS       perform transcript-level and gene-level abundance
                        estimation using salmon
  --errcheck            perform error correction to remove misalignments from
                        the input alignment
  --nounal              do not include unaligned reads in the output BAM
  -v, --verbose         print all output of the tools to stderr. If not
                        specified - only the final report will be generated
  --log                 If enabled logs of the steps in the protocol will be
                        compiled into a single document and saved in
                        <output>.log
  --no-discord          realign discordant transcriptomic reads
  --no-single           realign single mates reported in transcriptomic
                        alignment
  --sleep SLEEP         instructs the protocol to sleep for the specified
                        number of seconds after each operation which writes
                        data to disk.
```

#### Help page: `transHISAT2.py build --help`

```
usage: transHISAT2.py build [-h] --gff GFF --ref REF -o OUTPUT
                            [--threads THREADS] [--type {hisat,bowtie}]
                            [--kmerlen KMERLEN] [--locus]

optional arguments:
  -h, --help            show this help message and exit
  --gff GFF             GFF or GTF annotation of the genome
  --ref REF             Fasta-formatted reference genome
  -o OUTPUT, --output OUTPUT
                        output database directory
  --threads THREADS     number of threads to use
  --type {hisat,bowtie}
                        which aligner should be used
  --kmerlen KMERLEN     kmer length to use in the search for multimappers
  --locus               build index for a second bowtie run to align against
                        loci - parsimonious preMRNA mode
```