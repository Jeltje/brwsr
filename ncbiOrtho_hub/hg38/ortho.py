#! /cluster/home/jeltje/miniconda3/bin/python3.11

import os 
import csv
import gzip
from multiprocessing import Pool
import sys
sys.path.append('/hive/groups/browser/pycbio/lib')
from pycbio.hgdata.bed import Bed, BedBlock, BedReader, intArraySplit

def process_gzipped_file(file_path):
  '''create csv dict'''
  
  with gzip.open(file_path, 'rt') as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    reader = list(reader)
  return reader


def orthoMatches(orthofile, org1='9606', org2='10090'):
    '''Using the NCBI taxon IDs, get ortholog combinations'''
    orthoDict = dict()
    genesOfInterest = set()
    reader = process_gzipped_file(orthofile)
    for row in reader:
        if row['#tax_id'] == org1:
            if row['Other_tax_id'] == org2:
                gid = row['GeneID']
                ogid = row['Other_GeneID']
                if gid in orthoDict:
                    orthoDict[gid].append(ogid)
                else:
                    orthoDict[gid] = ogid
                genesOfInterest.update([gid, ogid])
    return orthoDict, genesOfInterest

def accessions(accfile, genesOfInterest):
    '''Extract gene symbols and one best accession for all genes of interest'''
    geneObjs = list()
    reader = process_gzipped_file(accfile)
    prevId = None
    someRows = []
    for row in reader:
        if row['GeneID'] in genesOfInterest:
            if row['GeneID'] == prevId:
                geneObj = makeGeneObject(someRows)
                if geneObject:
                    geneObjs.append(geneObject)
                prevID = row['GeneID']
            someRows.append(row)
    geneObj = makeGeneObject(someRows)
    if geneObject:
        geneObjs.append(geneObject)
    return(geneObjs)
    

def makeGeneObject(rows):
    '''From a set of rows for one gene, get the best representative, genome accession, location, and symbol'''
    bestToWorst = ['REVIEWED', 'VALIDATED', 'PROVISIONAL', 'INFERRED', 'PREDICTED', 'MODEL']
    for status in bestToWorst:
        selected = [row for row in rows if row['status'] == status]
        if len(selected) > 0:
            break
    if len(selected) == 0:
        print('No good members found', rows[0]['GeneID'], file=sys.stderr)
        return False
    
    chroms = set(row['genomic_nucleotide_accession.version'] for row in selected)
    if len(chroms) > 1:
        print('Confusing chromosome location', selected[0]['geneID'], chroms, file=sys.stderr)
        return False
    # be aware we're only taking starts and stops of the best set, not of all tx
    genStart = min([int(row['start_position_on_the_genomic_accession']) for row in selected])
    genStop = max([int(row['end_position_on_the_genomic_accession']) for row in selected])
    geneObject = gene(genStart, genStop, selected[0])
    geneObject.write()
    return geneObject

class Gene(object):
    def __init__(self, genStart, genStop, representativeRow):
        self.genStart = genStart
        self.genStop = genStop
        self.geneID = representativeRow['GeneID']
        self.tx = representativeRow['RNA_nucleotide_accession.version']
        self.chrom = representativeRow['genomic_nucleotide_accession.version']
        self.name = representativeRow['Symbol']
        self.strand = representativeRow['orientation']
    def write(self):
        '''Create bed object and print'''
        bedObj = Bed(self.chrom, self.geneStart, self.geneStop, name=name, strand=self.strand,
                 score=score, numStdCols=6)
        bedOjb.write(sys.stdout)



# Main

orthofile = '/hive/users/markd/NCBI-SqlLite/2024-10-22/download/gene_orthologs.gz'
accfile = '/hive/users/markd/NCBI-SqlLite/2024-10-22/download/gene2accession.gz'
accfile = 'tryacc.gz'

orthoDict, genesOfInterest = orthoMatches(orthofile)
print('Step 2')
geneObjs = accessions(accfile, genesOfInterest)
