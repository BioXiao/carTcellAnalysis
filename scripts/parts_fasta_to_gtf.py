#!/usr/bin/env python

import sys, re
from Bio import SeqIO


def read_fasta(fastaFile):
    handle = open(fastaFile, "rU")
    fastaRecords = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    return fastaRecords

def build_parts_dict(fastaRecords):
    partsDict = {}
    for idx,record in enumerate(fastaRecords):
        part_tag = re.search('(?<=gene_part=).*',
                             record.description).group()

        num_records = len(partsDict.get(part_tag, {}))
        part_record = '{}_r{}'.format(part_tag, num_records + 1)
        partsDict.setdefault(part_tag, {}).update({idx: part_record})

    return partsDict

def build_fasta_list(fastaRecords, spacer=0, contiguous=True):
    totalLength = 0
    fastaList = []
    partsDict = build_parts_dict(fastaRecords)
    part_lengths = map(lambda x: len(partsDict[x]), partsDict)
    if any([l > 1 for l in part_lengths]):
        multi = True
    else:
        multi = False

    for idx,record in enumerate(fastaRecords):
        if not contiguous:
            totalLength = 0

        recordStart = totalLength + spacer + 1
        recordStop = totalLength + len(record) - spacer
        gene_tag = re.search('(?<=gene=).*(?= )',
                             record.description).group()
        part_tag = re.search('(?<=gene_part=).*',
                             record.description).group()

        if multi:
            part_tag = partsDict[part_tag][idx]

        fastaList.append(('\t').join([record.id, 'custom', 'transcript',
                                          str(recordStart), str(recordStop),
                                          '.', '+', '.',
                                          ('gene_id "%s"; '
                                           'transcript_id "%s";\n') %
                                           (gene_tag, part_tag)]))
        totalLength += len(record)

    return fastaList

def write_gtf(fastaList, gtfOutFile):
    with open(gtfOutFile, 'wb') as f:
        f.writelines(fastaList)

def main(argv):
    fastaFile = sys.argv[1]
    gtfOutFile = sys.argv[2]
    if len(sys.argv) > 3:
        spacer = int(sys.argv[3])
    else:
        spacer = 0
    if len(sys.argv) > 4:
        contiguous = {'true': True, 'false': False}[sys.argv[4].lower()]
    else:
        contiguous = True


    fastaRecords = read_fasta(fastaFile)
    fastaList = build_fasta_list(fastaRecords, spacer, contiguous)

    write_gtf(fastaList, gtfOutFile)


if __name__ == "__main__":
    main(sys.argv[1:])
