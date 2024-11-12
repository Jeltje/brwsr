#! /cluster/home/jeltje/miniconda3/bin/python3.11

import os 
import csv
import sys
from collections import namedtuple
sys.path.append('/hive/groups/browser/pycbio/lib')
from pycbio.hgdata.bed import Bed, BedBlock, BedReader, intArraySplit



def makeGenesFromTx(infile, orthoDict):
    '''Extract gene symbols and create gene boundary blocks from multiple lines''' 
    with open(infile, 'r') as f:
        reader = csv.DictReader(f,  delimiter='\t')
        someRows = [next(reader)]
        prevID = someRows[0]['GeneID']
        geneObjDict = dict()
        for row in reader:
            # collect lines until we reach a  new gene
            if row['GeneID'] != prevID:
                 if prevID in orthoDict:
                     # make a gene object with the ortholog IDs added
                     geneObj = makeGeneObject(someRows, orthoDict[prevID])
                     geneObjDict[prevID] = geneObj
                 prevID = row['GeneID']
                 someRows = []
            someRows.append(row)
        # add the final gene
        if prevID in orthoDict:
            geneObj = makeGeneObject(someRows, orthoDict[prevID])
            geneObjDict[prevID] = geneObj
        return geneObjDict   
    

def makeGeneObject(rows, orthoList):
    '''From a set of rows for one gene, get the best representatives, genome accession, location, and symbol'''
    returnObjects = dict()
    best = ['REVIEWED', 'VALIDATED', 'PROVISIONAL', 'INFERRED']
    ok = ['PREDICTED', 'MODEL']
    selected = []
    for status in best:
        selected = [row for row in rows if row['status'] in [best]]
    if len(selected) == 0:
        # if we don't have a set of genes now, see if we have ok genes
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
        geneObject = Gene(genStart, genStop, subselect[0], orthoList)
    # can only do this if we fix duplicate IDs, hm
        #returnObjects[subselect['GeneID']] = geneObject
        # for now just return when we have one
        return geneObject

class Gene(object):
    '''Hold coordinate, ortholog, and species info per gene'''
    def __init__(self, genStart, genStop, representativeRow, orthoList):
        self.genStart = genStart
        self.genStop = genStop
        self.geneID = representativeRow['GeneID']
        self.orthoList = orthoList
        self.taxID = representativeRow['#tax_id']
        self.symbol = representativeRow['Symbol']
        self.tx = representativeRow['RNA_nucleotide_accession.version']
        self.chrom = representativeRow['genomic_nucleotide_accession.version']
        self.name = representativeRow['Symbol']
        self.strand = representativeRow['orientation']
    def makeHtml(self, speciesDict, geneDict):
        '''Turn coordinate info for ortholog genes into html'''
        # hgTracks?db=hg38&position=chr7:155799980-155812463&knownGene=pack
        self.url = ''
        baselink = '{species}:<a href="hgTracks?db={genome}&position={chrom}:{start}-{end}&knownGene=pack">{symbol}</a><br>'
        for geneName in self.orthoList:
             if geneName in geneDict:
                 gene = geneDict[geneName]
                 org = speciesDict[gene.taxID]
                 self.url += baselink.format(species=org.colloq, genome=org.ucsc, chrom=gene.chrom, start=gene.genStart, end=gene.genStop, symbol=gene.symbol)
    def write(self, outf):
        '''Create bed object and print'''
        bedObj = Bed(self.chrom, self.genStart, self.genStop, name=self.geneID, strand=self.strand,
                 score='0', numStdCols=9, extraCols=[self.name, self.url])
        bedObj.write(outf)

# get tax_id, ucsc name and common name for each species we want tracks for
# we're making all vs all tracks
speciesDict = dict()
species = namedtuple('species', ['colloq', 'ucsc']) 
with open('species.list', 'r') as f:
    for line in f:
        ncbi, colloq, ucsc = line.strip().split('\t')
        speciesDict[ncbi] = species(colloq, ucsc)

# gene IDs are unique, you don't need to keep the species. So it's a matter of keeping pairs
# we want to capture forward and reverse matches because they're not symmetrical
# meaning there are mouse-human matches that are not in human-mouse
orthoDict = dict()
with open('gene_orthologs', 'r') as f:
#with open('test_orthologs', 'r') as f:
    reader = csv.DictReader(f,  delimiter='\t')
    for row in reader:
        if row['#tax_id'] in speciesDict and row['Other_tax_id'] in speciesDict:
            if not row['GeneID'] in orthoDict:
                orthoDict[row['GeneID']] = set()
            orthoDict[row['GeneID']].add(row['Other_GeneID'])
            if not row['Other_GeneID'] in orthoDict:
                orthoDict[row['Other_GeneID']] = set()
            orthoDict[row['Other_GeneID']].add(row['GeneID'])


# now we know which genes we want, get their coordinate info
# add their ortholog set while we're at it
geneDict = makeGenesFromTx('gene2ucscAccession.txt', orthoDict)
 
# for each organism, print a bed file
for tax_id in speciesDict:
    organism = speciesDict[tax_id]
    organismGenes = {k:v for k, v in geneDict.items() if v.taxID == tax_id}
    outf = open(f'{organism.ucsc}.bed', 'w')
    for gene in organismGenes.values():
        gene.makeHtml(speciesDict, geneDict)
        print(gene.url)
        gene.write(outf)
    outf.close()
