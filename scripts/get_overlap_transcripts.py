#!/usr/bin/env python

import sys, re, yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import single_letter_alphabet

def read_gene_yaml(geneYaml):
    with open(geneYaml) as f:
        geneDict = yaml.load(f)
    return geneDict

def read_fasta(fastaFile):
    handle = open(fastaFile, "rU")
    fastaRecords = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    return fastaRecords

def get_fasta_subset(fastaRecords, geneDict):
    return [r for r in fastaRecords
            if re.search('(?<=gene=)\w*', r.description).group() in geneDict
            and re.search('^NM', r.name)]

def format_subset_records(fastaSubset, name, geneDict):
    for record in fastaSubset:
        record.id = record.name
        gene = re.search('(?<=gene=)\w*', record.description).group()
        record.description = "gene={} {}_gene_part={}".format(gene, name,
                                                              geneDict[gene])
    return fastaSubset


def write_fasta(fastaRecords, fastaOutFile):
    with open(fastaOutFile, 'wb') as f:
        writer = SeqIO.FastaIO.FastaWriter(f, wrap=70)
        writer.write_file(fastaRecords)

def main(argv):
    fastaFile = sys.argv[1]
    name = sys.argv[2]
    geneYaml = sys.argv[3]
    fastaOutFile = sys.argv[4]

    geneDict = read_gene_yaml(geneYaml)
    fastaRecords = read_fasta(fastaFile)
    fastaSubset = get_fasta_subset(fastaRecords, geneDict)
    fastaSubset = format_subset_records(fastaSubset, name, geneDict)

    write_fasta(fastaSubset, fastaOutFile)


if __name__ == "__main__":
    main(sys.argv[1:])
