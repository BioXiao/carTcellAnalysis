#!/bin/python

import sys, re
from Bio import SeqIO

def read_gtf(gtfFile):
    with open(gtfFile, 'rU') as f:
        gtfLines = f.readlines()
    return gtfLines

def build_gene_dict(gtfLines, geneDict={}):
    for idx,line in enumerate(gtfLines):
        lineParts = line.split('\t')

        geneIdRe = re.compile('(?<=(gene_id "))[A-Z]+')
        geneId = geneIdRe.search(line).group()
        if geneId in geneDict:
            geneDict[geneId]['start'].append(lineParts[3])
            geneDict[geneId]['stop'].append(lineParts[4])
        else:
            geneDict[geneId] = {'start': [lineParts[3]],
                                 'stop': [lineParts[4]]}
    return geneDict

def build_gene_refflat(geneDict, refFlatRecords=[]):
    for gene in geneDict:
        refFlatList = [gene, gene + '-1', gene, '+',
                       geneDict[gene]['start'][0],
                       geneDict[gene]['stop'][-1],
                       geneDict[gene]['start'][0],
                       geneDict[gene]['start'][0],
                       str(len(geneDict[gene]['start'])),
                       (',').join(geneDict[gene]['start']) + ',',
                       (',').join(geneDict[gene]['stop']) + ',',
                       '\n']
        refFlatRecords.append(('\t').join(refFlatList))
    return refFlatRecords

def write_refflat(refFlatRecords, refFlatOutFile):
    with open(refFlatOutFile, 'wb') as f:
        f.writelines(refFlatRecords)

def main(argv):
    gtfFile = sys.argv[1]
    refFlatOutFile = sys.argv[2]

    gtfLines = read_gtf(gtfFile)
    geneDict = build_gene_dict(gtfLines)
    refFlatRecords = build_gene_refflat(geneDict)
    write_refflat(refFlatRecords, refFlatOutFile)

if __name__ == "__main__":
    main(sys.argv[1:])


