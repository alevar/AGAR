import os
import sys
import subprocess


def main(args):
    assert os.path.exists(os.path.abspath(args.gff)), "gff annotation not found"
    assert os.path.exists(os.path.abspath(args.ref)), "reference file not found"
    # TODO: replace the samtools faidx and build_header function with an implementation within gtf_to_fasta
    # TODO: rename gtf_to_fasta to something else
    # TODO: create a new name for the tool
    if not os.path.exists(os.path.abspath(args.ref) + ".fai"):
        print("FASTA index for the reference genome not found. Building now.")
        subprocess.call(["samtools", "faidx", args.ref])

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # indexer
    print("Extracting transcriptome data and building an index")
    indexer_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     'trans2genome_indexer')  # get path to the gtf_to_fasta that was compiled with the package
    print(" ".join([indexer_path, "-k", str(args.kmerlen), "-a", args.gff, "-r", args.ref, "-o",
                    os.path.abspath(args.output) + "/db", "-m"]))
    subprocess.call([indexer_path, "-k", str(args.kmerlen), "-a", args.gff, "-r", args.ref, "-o",
                     os.path.abspath(args.output) + "/db", "-m"])

    if args.type == "bowtie":
        print("Building additional bowtie transcriptome index")
        subprocess.call(
            ["bowtie2-build", os.path.abspath(args.output) + "/db.fasta", os.path.abspath(args.output) + "/db"])
    elif args.type == "hisat":
        print("Building transcriptome database for HISAT2")
        subprocess.call(["hisat2-build", "-p", args.threads, os.path.abspath(args.output) + "/db.fasta",
                         os.path.abspath(args.output) + "/db"])
    else:
        print("unknown aligner specified")
    if args.locus:
        print("building a locus specific bowtie index for the second bowtie stage")
        # first extract the locus information
        subprocess.call["./extractLocus.py", args.gff, s.path.abspath(args.output) + "/db.locus.gff"]
        # second run the bedtools
        subprocess.call["bedtools", "maskfasta", "-fo", s.path.abspath(
            args.output) + "/db.locus.fasta", "-fi", args.ref, "-bed", s.path.abspath(args.output) + "/db.locus.gff"]
        # third build the index
        subprocess.call(["bowtie2-build", os.path.abspath(args.output) + "/db.locus.fasta",
                         os.path.abspath(args.output) + "/db.locus"])
