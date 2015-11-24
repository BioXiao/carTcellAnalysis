#!/bin/python

'''
DEPRECATED!
Used to parse and format results from TopHat and featureCounts.
'''

import sys, re, csv

def build_junction_dict(gtfFile, juncDict={}):
    with open(gtfFile, 'rU') as f:
        gtfLines = f.readlines()
    gtfLines = gtfLines[:-1]
    for idx in range(len(gtfLines) - 1):
        feature1 = gtfLines[idx]
        feature2 = gtfLines[idx + 1]

        featureRe = re.compile('(?<=(exon_name ")).*(?=")')
        f1name = featureRe.search(feature1).group()
        f2name = featureRe.search(feature2).group()
        juncName = (':').join([f1name, f2name])
        juncDict[juncName] = [int(feature1.split('\t')[4]) - 1,
                              int(feature2.split('\t')[3]) - 1]
    return juncDict

def count_junction_hits(juncBedFile, juncDict, juncCountDict={}):
    with open(juncBedFile, 'rU') as f:
        juncBedLines = f.readlines()
    juncBedLines = juncBedLines[1:]
    for junc in juncDict:
        juncStart = juncDict[junc][0]
        juncEnd = juncDict[junc][1]
        juncCount = [line.split('\t')[4] for line in juncBedLines \
                     if int(line.split('\t')[1]) < juncStart \
                     and int(line.split('\t')[2]) > juncEnd]
        if len(juncCount):
            juncCountDict[junc] = juncCount[0]
        else:
            juncCountDict[junc] = '0'
    return juncCountDict

def build_feature_dict(gtfFile, featureDict={}):
    with open(gtfFile, 'rU') as f:
        gtfLines = f.readlines()
    gtfLines = gtfLines[:-1]
    for featureLine in gtfLines:
        featureRe = re.compile('(?<=(exon_name ")).*(?=")')
        fname = featureRe.search(featureLine).group()
        featureDict[fname] = [featureLine.split('\t')[3],
                              featureLine.split('\t')[4]]
    return featureDict

def count_feature_hits(countsFile, featureDict, featureCountDict={}):
    with open(countsFile, 'rU') as f:
        countLines = f.readlines()
    countLines = countLines[2:]
    for feature in featureDict:
        fStart = featureDict[feature][0]
        fEnd = featureDict[feature][1]
        fCount = [line.split('\t')[6] for line in countLines \
                  if line.split('\t')[2] == fStart \
                  and line.split('\t')[3] == fEnd]
        if len(fCount):
            featureCountDict[feature] = fCount[0].rstrip('\n')
        else:
            featureCountDict[feature] = '0'
    return featureCountDict

def write_results(outFile, featureCountDict, juncCountDict):
    w = csv.writer(open(outFile, 'wb'), delimiter='\t')
    for key, value in featureCountDict.items():
        w.writerow([key, value])

    for key, value in juncCountDict.items():
        w.writerow([key, value])

def main(argv):
    gtfFile = sys.argv[1]
    juncBedFile = sys.argv[2]
    countsFile = sys.argv[3]
    outFile = sys.argv[4]

    juncDict = build_junction_dict(gtfFile)
    juncCountDict = count_junction_hits(juncBedFile, juncDict)

    featureDict = build_feature_dict(gtfFile)
    featureCountDict = count_feature_hits(countsFile, featureDict)
    write_results(outFile, featureCountDict, juncCountDict)

if __name__ == "__main__":
    main(sys.argv[1:])
