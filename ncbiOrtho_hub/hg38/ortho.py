#! /cluster/home/jeltje/miniconda3/bin/python3.11

import os 
import csv
import sys
sys.path.append('/hive/groups/browser/pycbio/lib')
from pycbio.hgdata.bed import Bed, BedBlock, BedReader, intArraySplit



def makeGenesFromTx(infile):
    '''Extract gene symbols and create gene boundary blocks''' 
    with open(infile, 'r') as f:
        reader = csv.DictReader(f,  delimiter='\t')
        someRows = [next(reader)]
        prevID = someRows[0]['GeneID']
        geneObjDict = dict()
        for row in reader:
            if row['GeneID'] != prevID:
                 geneObj = makeGeneObject(someRows)
                 geneObjDict[prevID] = geneObj
                 prevID = row['GeneID']
                 someRows = []
            someRows.append(row)
        geneObj = makeGeneObject(someRows)
        geneObjDict[prevID] = geneObj
        return geneObjDict   
    

def makeGeneObject(rows):
    '''From a set of rows for one gene, get the best representatives, genome accession, location, and symbol'''
    returnObjects = dict()
    best = ['REVIEWED', 'VALIDATED', 'PROVISIONAL', 'INFERRED']
    ok = ['PREDICTED', 'MODEL']
    selected = []
    for status in best:
        selected = [row for row in rows if row['status'] in [best]]
    if len(selected) == 0:
        for status in ok:
            selected = [row for row in rows if row['status'] in [ok]]
    if len(selected) == 0:
        # all tx for this gene have status '-'
        selected = rows
    
    chroms = set(row['genomic_nucleotide_accession.version'] for row in selected)
    # there are duplicate IDs on chrX and Y TODO: make unique IDs
    for ch in chroms:
        subselect = [row for row in selected if row['genomic_nucleotide_accession.version'] == ch]
        # be aware we're only taking starts and stops of the best set, not of all tx
        genStart = min([int(row['start_position_on_the_genomic_accession']) for row in subselect])
        genStop = max([int(row['end_position_on_the_genomic_accession']) for row in subselect])
        geneObject = Gene(genStart, genStop, subselect[0])
	# can only do this if we fix duplicate IDs, hm
        #returnObjects[subselect['GeneID']] = geneObject
        # for now just return when we have one
        return geneObject

class Gene(object):
    def __init__(self, genStart, genStop, representativeRow):
        self.genStart = genStart
        self.genStop = genStop
        self.geneID = representativeRow['GeneID']
        self.symbol = representativeRow['Symbol']
        self.tx = representativeRow['RNA_nucleotide_accession.version']
        self.chrom = representativeRow['genomic_nucleotide_accession.version']
        self.name = representativeRow['Symbol']
        self.strand = representativeRow['orientation']
    def makeHtml(self, genome='mm39', species=mouse):
        '''Turn coordinate info into html'''
        # hgTracks?db=hg38&position=chr7:155799980-155812463&knownGene=pack
        baselink = '<a href="hgTracks?db=' + genome
        baselink += f'&position={self.chrom}:{self.genStart}-{self.genStop}'
        baselink += '&knownGene=pack">'
        baselink += f'{self.symbol}</a>'
        self.url = f'{species}:{baselink}<br>'
    def write(self, url, ortholog, outf):
        '''Create bed object and print'''
        bedObj = Bed(self.chrom, self.genStart, self.genStop, name=self.geneID, strand=self.strand,
                 score='0', numStdCols=9, extraCols=[self.name, url, ortholog])
        bedObj.write(outf)

# we want a field that gives the ortholog in mouse, with symbol and chromosome location



# Main

humGeneDict = makeGenesFromTx('humOrtho.genes')
musGeneDict = makeGenesFromTx('musOrtho.genes')
with open('humOrtho.bed', 'w') as h, open('musOrtho.bed', 'w') as m, open('humVsMouse.ortho', 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if fields[1] in humGeneDict and fields[4] in musGeneDict:
            humgene = humGeneDict[fields[1]]
            musgene = musGeneDict[fields[4]]
            musgene.makeHtml()
            humgene.makeHtml()
            humgene.write(musgene.url, musgene.symbol, h)
            musgene.write(humgene.url, humgene.symbol, m)
        



