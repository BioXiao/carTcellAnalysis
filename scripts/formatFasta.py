#!/bin/python
import sys, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import single_letter_alphabet

def read_fasta(rawFastaFile):
    '''
    Simply read in raw (unformatted) sequence file as text.
    '''

    with open(rawFastaFile, 'rU') as f:
        fastaLines = f.readlines()
    return fastaLines

def get_parts(fastaLines):
    '''
    Parse unformatted sequences from gene construct to identify division points
    between individual 'parts'; store part names and corresponding lines of
    text.
    '''

    parts = []
    partLines = []
    for idx,line in enumerate(fastaLines):
        if len(line.split(' ')) > 1:
            tag = line.split(' ')[1]
            parts.append(re.sub("(\(|\)|\-|\n)", "", tag))
            partLines.append(idx)
    return (parts, partLines)

def build_parts_dict(fastaLines, parts, partLines):
    '''
    Extract sequences and labels from lists of part names and lines.
    '''

    partLines = [-1] + partLines
    partsDict = {}
    for partNum in range(len(partLines) - 1):
        part = parts[partNum]
        partRange = fastaLines[partLines[partNum]+1:partLines[partNum + 1]+1]
        partSeq = [ fastaLine.split(' ')[0] \
                    for fastaLine in partRange ]
        partsDict[part] = ('').join([ re.sub("\n", "", seq).lower() \
                                      for seq in partSeq ])
    return partsDict

def build_fasta_records(partsDict, gene, merge=True, transcript=False):
    '''
    Build FASTA records from dictionary of gene construct 'parts'.

    :param merge: If set to True, parts will be combined into a single FASTA
                  record; otherwise, a separate record will be written for each
                  part. Note: the True setting is most useful for building
                  test data.
    :param transcript: If set to True, the merged sequence is treated as a
                       processed (contiguous) transcript; if False, buffers of
                       'N' nucleotides will be inserted between parts to
                       simulate not-yet-spliced introns.
    '''

    fastaRecords = []
    # pad each individual sequence with 5 N's on each side
    if merge:
        if transcript:
            bufSeq = ''
            seqId = gene + '-1'
        else:
            bufSeq = 'NNNNN'
            seqId = '0'
        mergedSeq = ('').join([ bufSeq + partsDict[part] + bufSeq \
                                for part in partsDict ])
        record = SeqRecord(Seq(mergedSeq, single_letter_alphabet),
                           id=seqId,
                           description=("%s_full_gene" % gene))
        fastaRecords.append(record)
    else:
        for part in partsDict:
            record = SeqRecord(Seq(partsDict[part], single_letter_alphabet),
                               id=part,
                               description=("%s_gene_part:%s" % (gene, part)))
            fastaRecords.append(record)
    return fastaRecords

def write_fasta(fastaRecords, fastaOutFile):
    with open(fastaOutFile, 'wb') as f:
        [ SeqIO.write(record, f, "fasta") for record in fastaRecords ]

def main(argv):
    rawFastaFile = sys.argv[1]
    gene = sys.argv[2]
    fastaOutFile = sys.argv[3]
    if len(sys.argv) > 4:
        merge = sys.argv[4]
    else:
        merge = False
    if len(sys.argv) > 5:
        transcript = sys.argv[5]
    else:
        transcript = False

    fastaLines = read_fasta(rawFastaFile)
    parts,partLines = get_parts(fastaLines)
    partsDict = build_parts_dict(fastaLines, parts, partLines)
    fastaRecords = build_fasta_records(partsDict, gene, merge, transcript)
    write_fasta(fastaRecords, fastaOutFile)

if __name__ == "__main__":
    main(sys.argv[1:])
