# Using gcc for now
CC=g++
# Flags for the compiles
CFLAGS=-c -Wall
# Samtools path
SAMTOOLS=/home/costea/samtools
# HTS lib path
HTSLIB=/home/costea/htslib
INCLUDE=include/

all: removeUnmapped qaCompute computeInsertSizeHistogram doBWAQualTrimming

fixPairedEnd: fixPairedEnd.o
	$(CC) fixPairedEnd.o -o fixPairedEnd -L$(SAMTOOLS) -L$(HTSLIB) -lbam -lhts -lz

removeUnmapped: removeUnmapped.o
	$(CC) removeUnmapped.o -o removeUnmapped -L$(SAMTOOLS) -L$(HTSLIB) -lbam -lhts -lz

qaCompute: qaCompute.o
	$(CC) qaCompute.o -o qaCompute -L$(SAMTOOLS) -L$(HTSLIB) -lbam -lhts -lz #-lefence

computeInsertSizeHistogram: computeInsertSizeHistogram.o
	$(CC) computeInsertSizeHistogram.o -o computeInsertSizeHistogram -lz -L$(SAMTOOLS) -L$(HTSLIB) -lbam -lhts #-fopenmp

doBWAQualTrimming: doBWAQualTrimming.o
	$(CC) doBWAQualTrimming.o -o doBWAQualTrimming -lz

fixPairedEnd.o:	fixPairedEnd.c
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) $(CFLAGS) fixPairedEnd.c

removeUnmapped.o: removeUnmapped.c
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) $(CFLAGS) removeUnmapped.c

qaCompute.o: qaCompute.c
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) $(CFLAGS) -o qaCompute.o qaCompute.c

computeInsertSizeHistogram.o: computeInsertSizeHistogram.c
	$(CC) -I$(SAMTOOLS) -I$(HTSLIB) $(CFLAGS) computeInsertSizeHistogram.c

doBWAQualTrimming.o: doBWAQualTrimming.c
	$(CC) -I$(INCLUDE) -I$(HTSLIB) $(CFLAGS) doBWAQualTrimming.c

clean:
	rm -rf *o removeUnmapped
	rm -rf *o qaCompute
	rm -rf *o computeInsertSizeHistogram
	rm -rf *o doBWAQualTrimming
	rm -rf *o fixPairedEnd
