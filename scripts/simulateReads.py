#!/bin/python

import sys, itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import single_letter_alphabet
from random import randint

def read_fasta(fastaFile):
    handle = open(fastaFile, "rU")
    fastaRecords = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    return fastaRecords

def generate_intra_feature_reads(fastaRecords, numReads=5, minReadLength=20, maxReadLength=50, end="none"):
    fastqRecords = []
    numRecords = len(fastaRecords)
    for i in range(0, numReads):
        recordIdx = randint(0, numRecords - 1)
        srcRecord = fastaRecords[recordIdx]
        seqLength = len(srcRecord)
        if end == "left":
            subSeqStart = 0
            subSeqEnd = maxReadLength
        elif end == "right":
            subSeqEnd = seqLength
            subSeqStart = seqLength - maxReadLength
        else:
            subSeqStart = randint(0, seqLength)
            remainingSeq = seqLength - subSeqStart
            if remainingSeq < minReadLength:
                subSeqEnd = subSeqStart
                subSeqStart = randint(max(0, subSeqEnd - maxReadLength), 
                                      subSeqStart - minReadLength)

            else:
                subSeqEnd = randint(subSeqStart + minReadLength, 
                                    min(subSeqStart + maxReadLength, seqLength))

        fastqRecord = srcRecord[subSeqStart:subSeqEnd]
        fastqRecord.letter_annotations["phred_quality"] = \
            list(itertools.repeat(37, len(fastqRecord))) 
        fastqRecord.id = "%s[%s-%s]" % (fastqRecord.id, 
                                        str(subSeqStart), str(subSeqEnd))
        fastqRecords.append(fastqRecord)
    
    return fastqRecords

def generate_cross_feature_reads(fastaRecords, numReads=5, minReadLength=20, maxReadLength=50):
    fastqRecords = []
    numRecords = len(fastaRecords)
    for i in range(0, numReads):
        recordIdx = randint(0, numRecords - 2)
        srcRecord1 = fastaRecords[recordIdx]
        srcRecord2 = fastaRecords[recordIdx + 1]

        seqLength1 = len(srcRecord1)
        seqLength2 = len(srcRecord2)
        readLength1 = min(seqLength1, randint(minReadLength, maxReadLength - 10))
        readLength2 = min(seqLength2, randint(max(minReadLength - readLength1, 0),
                                              maxReadLength - readLength1))
        read1 = generate_intra_feature_reads([srcRecord1], 1, readLength1, readLength1, end="right")[0]
        read2 = generate_intra_feature_reads([srcRecord2], 1, readLength2, readLength2, end="left")[0]

        joinedRead = sum([read1, read2], Seq('', single_letter_alphabet))
        joinedRead.id = (':').join([read1.id, read2.id])
        joinedRead.name = (':').join([read1.name, read2.name])
        joinedRead.description = (':').join([read1.description, read2.description])
        joinedRead.letter_annotations["phred_quality"] = \
            list(itertools.repeat(37, len(joinedRead)))        

        fastqRecords.append(joinedRead)

    return fastqRecords

def write_fastq(fastqRecords, fastqOutFile):
    with open(fastqOutFile, 'wb') as f:
        [SeqIO.write(record, f, "fastq") for record in fastqRecords]

def main(argv):
    fastaFile = sys.argv[1]
    fastqOutFile = sys.argv[2]
    numReads = int(sys.argv[3])

    fastaRecords = read_fasta(fastaFile)
    fastqRecords = generate_intra_feature_reads(fastaRecords, numReads)
    fastqRecords = fastqRecords + generate_cross_feature_reads(fastaRecords, numReads)
    write_fastq(fastqRecords, fastqOutFile)

if __name__ == "__main__":
    main(sys.argv[1:])
