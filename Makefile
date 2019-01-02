CC=g++
CXXFLAGS=-Wall -std=c++11 -O3
TOPDIR := ./
SRCDIR := $(TOPDIR)src/

src := $(wildcard $(SRCDIR)*.cpp)
obj := $(src:.cpp=.o)

all: trans2genome 

trans2genome: $(obj)
	$(CC) -o $@ /home/avaraby1/genomicTools/trans2genome/src/include/htslib/sam.h $^ /home/avaraby1/genomicTools/htslib-1.9/libhts.so

# trans2genome: trans2genome.o
# 	$(CC) $(CFLAGS) trans2genome.o -o trans2genome -lhts

clean:
	rm $(obj)
