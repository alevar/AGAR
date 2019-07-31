#!/usr/bin/env python

# ./transHISAT2.py --1 ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/reads_1/sample_01_1.fasta
# --2 ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/reads_1/sample_01_2.fasta
# --gff ~/JHU/transcriptome/hisatTrans2Genome/errorData/me_mi/data/ann.gff --ref ~/genomicData/hg38/hg38_p8.fa
# -o ./test --tmp ./tmp

import os
import sys
import argparse
import subprocess

import build
import align


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')
    subparsers = parser.add_subparsers(help='sub-command help')

    # ==============================================================
    # ==============================================================
    # ==============================================================
    # ======================BUILD THE DATABASE======================
    # ==============================================================
    # ==============================================================
    # ==============================================================
    parser_build = subparsers.add_parser('build',
                                         help='build help')
    parser_build.add_argument("--gff",
                              type=str,
                              required=True,
                              help="GFF or GTF annotation of the genome")
    parser_build.add_argument("--ref",
                              type=str,
                              required=True,
                              help="Fasta-formatted reference genome")
    parser_build.add_argument("-o",
                              "--output",
                              type=str,
                              required=True,
                              help="output database directory")
    parser_build.add_argument("--threads",
                              type=str,
                              required=False,
                              default="1",
                              help="number of threads to use")
    parser_build.add_argument("--type",
                              type=str,
                              required=False,
                              default="hisat",
                              choices=["hisat", "bowtie"],
                              help="which aligner should be used")
    parser_build.add_argument("--kmerlen",
                              required=False,
                              type=int,
                              default=76,
                              help="kmer length to use in the search for multimappers")
    parser_build.add_argument("--locus",
                              action="store_true",
                              help="build index for a second bowtie run to align against loci - parsimonious preMRNA mode")

    parser_build.set_defaults(func=build.main)

    # ==============================================================
    # ==============================================================
    # ==============================================================
    # ======================PERFORM ALIGNMENT=======================
    # ==============================================================
    # ==============================================================
    # ==============================================================
    parser_align = subparsers.add_parser('align',
                                         help='alignment help')
    parser_align.add_argument("--m1",
                              type=str,
                              required=True,
                              help="File with #1 reads of the pair")
    parser_align.add_argument("--m2",
                              type=str,
                              required=True,
                              help="File with #2 reads of the pair")
    parser_align.add_argument("--single",
                              type=str,
                              required=False,
                              help="File with singletons")
    parser_align.add_argument("--fasta",
                              action="store_true",
                              help="Reads are in fasta format")
    parser_align.add_argument("--db",
                              type=str,
                              required=True,
                              help="path to the database directory")
    parser_align.add_argument("-o",
                              "--output",
                              type=str,
                              required=True,
                              help="output BAM file")
    parser_align.add_argument("--type",
                              type=str,
                              required=False,
                              default="bowtie",
                              choices=["hisat", "bowtie"],
                              help="which aligner to use")
    parser_align.add_argument("--tmp",
                              type=str,
                              required=False,
                              default="./tmp",
                              help="directory for tmp files")
    parser_align.add_argument("--genome-db",
                              type=str,
                              required=True,
                              help="path to the genome database for HISAT2")
    parser_align.add_argument("--threads",
                              type=str,
                              required=False,
                              default="1",
                              help="number of threads to use")
    parser_align.add_argument("-a",
                              "--all",
                              action="store_true",
                              help="run hisat/bowtie with the -a option")
    parser_align.add_argument("-k",
                              required=False,
                              type=int,
                              default=0,
                              help="-k argument for bowtie/hisat first pass")
    parser_align.add_argument("--bowtie",
                              required=False,
                              nargs='*',
                              help="additional arguments to be passed over to bowtie")
    parser_align.add_argument("--hisat",
                              required=False,
                              nargs='*',
                              help="additional arguments to be passed over to hisat")
    parser_align.add_argument("--keep",
                              action="store_true",
                              help="keep temporary files")
    parser_align.add_argument("--mf",
                              action="store_true",
                              help="use multimapper index to supplement alignments")
    parser_align.add_argument("--locus",
                              action="store_true",
                              help="perform second bowtie run to align against loci - parsimonious preMRNA mode")
    parser_align.add_argument("--abunds",
                              type=str,
                              required=False,
                              help="perform transcript-level and gene-level abundance estimation using salmon")
    parser_align.add_argument("--errcheck",
                              action="store_true",
                              help="perform error correction to remove misalignments from the input alignment")
    parser_align.add_argument("--nounal",
                              action="store_true",
                              help="do not include unaligned reads in the output BAM")

    parser_align.set_defaults(func=align.main)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
