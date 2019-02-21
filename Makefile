CC=g++
CFLAGS=-Wall -std=c++11 -O3
TOPDIR := ./
SRCDIR := $(TOPDIR)src/

src := $(wildcard $(SRCDIR)*.cpp)
obj := $(src:.cpp=.o)

all: trans2genome 

trans2genome: $(obj)
	$(CC) $(CFLAGS) -o $@ $^ -lhts

# trans2genome: trans2genome.o
# 	$(CC) $(CFLAGS) trans2genome.o -o trans2genome -lhts

clean:
	rm $(obj)
