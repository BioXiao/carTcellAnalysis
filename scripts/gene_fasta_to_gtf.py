#!/bin/python

'''
DEPRECATED!
Used to generate gene model (GTF/GFF file) from gene contstruct parts FASTA
for read quantification with htseq-count or featureCounts.
'''

import sys
from Bio import SeqIO


def read_fasta(fastaFile):
    handle = open(fastaFile, "rU")
    fastaRecords = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    return fastaRecords

def build_fasta_list(fastaRecords, gene, fastaList=[]):
    totalLength = 0
    spacer = 1 # use 5 for gapped reference
    for idx,record in enumerate(fastaRecords):
        recordStart = totalLength + spacer + 1
        recordStop = totalLength + len(record) + spacer
        # fastaList.append(('\t').join([gene + '-1', 'custom', 'exon',
        #                               str(recordStart), str(recordStop),
        #                               '.', '+', '.',
        #                               ('gene_id "%s"; '
        #                                'transcript_id "%s"; '
        #                                'exon_number "%s"; '
        #                                'exon_name "%s";\n') %
        #                                (gene, gene + '-1', idx + 1, record.name)]))
        # fastaList.append(('\t').join([gene + '-1', 'custom', 'CDS',
        #                               str(recordStart), str(recordStop),
        #                               '.', '+', '.',
        #                               ('gene_id "%s"; '
        #                                'transcript_id "%s"; '
        #                                'exon_number "%s"; '
        #                                'exon_name "%s";\n') %
        #                                (gene, gene + '-1', idx + 1, record.name)]))
        fastaList.append(('\t').join([gene + '-1', 'custom', 'transcript',
                                          str(recordStart), str(recordStop),
                                          '.', '+', '.',
                                          ('gene_id "%s"; '
                                           'transcript_id "%s";\n') %
                                           (gene, record.name)]))
        totalLength += len(record) + (spacer * 2)
    # fastaList.append(('\t').join([gene + '-1', 'custom', 'transcript',
    #                                   str(1), str(totalLength),
    #                                   '.', '+', '.',
    #                                   ('gene_id "%s"; '
    #                                    'transcript_id "%s";\n') %
    #                                    (gene, gene + '-1')]))
    return fastaList

def write_gtf(fastaList, gtfOutFile):
    with open(gtfOutFile, 'wb') as f:
        f.writelines(fastaList)

def main(argv):
    fastaFile = sys.argv[1]
    gene = sys.argv[2]
    gtfOutFile = sys.argv[3]

    fastaRecords = read_fasta(fastaFile)
    fastaList = build_fasta_list(fastaRecords, gene)

    write_gtf(fastaList, gtfOutFile)


if __name__ == "__main__":
    main(sys.argv[1:])
