# Using gcc for now
CC=g++
# Flags for the compiles
CFLAGS=-c -Wall
# Samtools path
SAMTOOLS=/g/bork8/costea/samtools/samtools-0.1.18
#SAMTOOLS=/g/bork8/costea/samtools/samtools-1.1
INCLUDE=include/

all: removeUnmapped qaCompute computeInsertSizeHistogram doBWAQualTrimming

fixPairedEnd: fixPairedEnd.o
	$(CC) fixPairedEnd.o -o fixPairedEnd -L$(SAMTOOLS) -lbam -lz

removeUnmapped: removeUnmapped.o
	$(CC) removeUnmapped.o -o removeUnmapped -L$(SAMTOOLS) -lbam -lz

qaCompute: qaCompute.o
	$(CC) qaCompute.o -o qaCompute -L$(SAMTOOLS) -lbam -lz #-lefence

computeInsertSizeHistogram: computeInsertSizeHistogram.o
	$(CC) computeInsertSizeHistogram.o -o computeInsertSizeHistogram -lz -L$(SAMTOOLS) -lbam #-fopenmp

doBWAQualTrimming: doBWAQualTrimming.o
	$(CC) doBWAQualTrimming.o -o doBWAQualTrimming -lz

fixPairedEnd.o:	fixPairedEnd.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) fixPairedEnd.c

removeUnmapped.o: removeUnmapped.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) removeUnmapped.c

qaCompute.o: qaCompute.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) -o qaCompute.o qaCompute.c

computeInsertSizeHistogram.o: computeInsertSizeHistogram.c
	$(CC) -I$(SAMTOOLS) $(CFLAGS) computeInsertSizeHistogram.c

doBWAQualTrimming.o: doBWAQualTrimming.c
	$(CC) -I$(INCLUDE) $(CFLAGS) doBWAQualTrimming.c

clean:
	rm -rf *o removeUnmapped
	rm -rf *o qaCompute
	rm -rf *o computeInsertSizeHistogram
	rm -rf *o doBWAQualTrimming
	rm -rf *o fixPairedEnd
