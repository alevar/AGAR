import os
import sys
import subprocess
from shutil import which


def main(args):
    assert os.path.exists(os.path.abspath(args.gff)), "gff annotation not found"
    assert os.path.exists(os.path.abspath(args.ref)), "reference file not found"
    # TODO: replace the samtools faidx with an implementation within gtf_to_fasta
    if not os.path.exists(os.path.abspath(args.ref) + ".fai"):
        print("FASTA index for the reference genome not found. Building now.")
        assert which("samtools") is not None, "samtools is not found - unable to proceed with building the index"
        subprocess.call(["samtools", "faidx", args.ref])

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    # indexer
    print("Extracting transcriptome data and building an index")
    indexer_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'agar_indexer')  # get path to the gtf_to_fasta that was compiled with the package
    print(" ".join([indexer_path, "-k", str(args.kmerlen), "-a", args.gff, "-r", args.ref, "-o",
                    os.path.abspath(args.output) + "/db", "-m"]))
    subprocess.call([indexer_path, "-k", str(args.kmerlen), "-a", args.gff, "-r", args.ref, "-o",
                     os.path.abspath(args.output) + "/db", "-m"])

    hisat2_rnasens_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                       'external/hisat2/hisat2-build')  # get path to the trans2genome that was compiled with the package

    if args.type == "bowtie":
        print("Building additional bowtie transcriptome index")
        subprocess.call(
            ["bowtie2-build", "--threads", str(args.threads), os.path.abspath(args.output) + "/db.fasta", os.path.abspath(args.output) + "/db"])
    elif args.type == "hisat":
        print("Building transcriptome database for HISAT2")
        subprocess.call([hisat2_rnasens_path, "-p", str(args.threads), os.path.abspath(args.output) + "/db.fasta",
                         os.path.abspath(args.output) + "/db"])
    else:
        print("unknown aligner specified")
        sys.exit()

    if args.locus:
        print("building a locus specific bowtie index for the second bowtie stage")
        # first extract the locus information
        extract_locus_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    'extractLocus.py')
        subprocess.call([extract_locus_path, args.gff, os.path.abspath(args.output) + "/db.locus.gff"])
        # print(" ".join([extract_locus_path,"--input",args.gff,"--out",os.path.abspath(args.output) + "/db.locus.gff","--faidx",args.ref+".fai"]))
        # second run the bedtools

        assert which("bedtools") is not None, "bedtools is not installed on the system - unable to complete building locus-specific index"

        subprocess.call(["bedtools", "maskfasta", "-fo", os.path.abspath(
            args.output) + "/db.locus.fasta", "-fi", args.ref, "-bed", os.path.abspath(args.output) + "/db.locus.gff"])
        # third build the index
        if args.type == "bowtie":
            print("building additional bowtie2 index for locus alignment")
            subprocess.call(["bowtie2-build", "--threads", str(args.threads), os.path.abspath(args.output) + "/db.locus.fasta", os.path.abspath(args.output) + "/db.locus"])
        else:
            print("building additional hisat2 index for locus alignment")
            subprocess.call([hisat2_rnasens_path, "-p". str(args.threads), os.path.abspath(args.output) + "/db.locus.fasta", os.path.abspath(args.output) + "/db.locus"])
