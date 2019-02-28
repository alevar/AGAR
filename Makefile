CC=g++
CFLAGS=-Wall -std=c++11 -O3
TOPDIR := ./
SRCDIR := $(TOPDIR)src/

src := $(wildcard $(SRCDIR)*.cpp)
obj := $(src:.cpp=.o)

all: trans2genome gtf_to_fasta

trans2genome: src/arg_parse.o src/map2gff.o src/trans2genome.o src/codons.o src/GBase.o src/tokenize.o
	$(CC) $(CFLAGS) -o $@ $^ -lhts

gtf_to_fasta: src/arg_parse.o src/FastaTools.o src/gdna.o src/gff.o src/codons.o src/GBase.o src/GFaSeqGet.o src/GTFToFasta.o src/tokenize.o
	$(CC) $(CFLAGS) -o $@ $^ -lhts

clean:
	rm $(obj)
