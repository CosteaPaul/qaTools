# Using gcc for now
CC=g++
# Flags for the compiles
CFLAGS=-c -Wall -g
# Samtools path
SAMTOOLS=/home/pcostea/private/tools/samtools-0.1.8/

all: removeUnmapped qaCompute computeInsertSizeHistogram

removeUnmapped: removeUnmapped.o
	$(CC) removeUnmapped.o -o removeUnmapped -lz -L$(SAMTOOLS) -lbam

qaCompute: qaCompute.o
	$(CC) qaCompute.o -o qaCompute -lz -L$(SAMTOOLS) -lbam #-lefence

computeInsertSizeHistogram: computeInsertSizeHistogram.o
	$(CC) computeInsertSizeHistogram.o -o computeInsertSizeHistogram -lz -L$(SAMTOOLS) -lbam #-fopenmp

removeUnmapped.o: removeUnmapped.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) removeUnmapped.c

qaCompute.o: qaCompute.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) qaCompute.c

computeInsertSizeHistogram.o: computeInsertSizeHistogram.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) computeInsertSizeHistogram.c

clean:
	rm -rf *o removeUnmapped
	rm -rf *o qaCompute
	rm -rf *o computeInsertSizeHistogram