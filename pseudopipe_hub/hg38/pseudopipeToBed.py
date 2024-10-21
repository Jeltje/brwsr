#!/usr/bin/env python3
import sys
import csv
import os
import re
import argparse
import textwrap
sys.path.append('/hive/groups/browser/pycbio/lib')
from pycbio.hgdata.bed import Bed, BedBlock, intArraySplit

def main():
	parser = argparse.ArgumentParser(
	formatter_class=argparse.RawDescriptionHelpFormatter,
	description=textwrap.dedent('''\
Converts gtf to bed format.

        '''))
	parser.add_argument('--genes', type=str, help='ENST and hugo gene IDs')
	parser.add_argument('--gtf', type=str, required=True, default=None,  help="annotation gtf")
	parser.add_argument('--bed', type=str, required=True, default=None,  help="output bedfile")
	args = parser.parse_args()

	gdict = dict()
	with open(args.genes, 'r') as g:
		for line in g:
			tx, gene = line.strip().split(" ")
			gdict[tx] = gene
	g.close()
	gtf_to_bed(args.bed, args.gtf, gdict)


class Gene(object):
	def __init__(self, itemRgb, tx, chrom, start, end, strand, extrafields):
		self.itemRgb = itemRgb
		self.txID = tx
		self.chrom = chrom
		# replace this misguided CHR_ otherwise we can't run chromToUcsc
		if chrom.startswith('CHR_'):
			self.chrom = re.sub('CHR_', '', chrom)
		self.geneStart = start
		self.geneEnd = end
		self.strand = strand
		self.blockList = [BedBlock(start, end)]
		self.extraCols = (tx,) + extrafields
	def add(self, start, end):
		self.geneStart = min(start, self.geneStart)
		self.geneEnd = max(end, self.geneEnd)
		self.blockList.append(BedBlock(start, end))
	def write(self, FH):
		'''Create a Bed object, then write'''
		name = self.txID
		sorted_blocks = sorted(self.blockList, key=lambda x: x.start)
		self.merge(sorted_blocks)
		bedObj = Bed(self.chrom, self.geneStart, self.geneEnd, name=name, strand=self.strand, 
		blocks=self.blocks, itemRgb=self.itemRgb, numStdCols=12, extraCols = self.extraCols)
		bedObj.write(FH)
	def merge(self, sorted_blocks):
		merged_blocks = []
		for block in sorted_blocks:
		        if not merged_blocks:
		            merged_blocks.append(block)
		        else:
		            last_block = merged_blocks[-1]
		            if block.start <= last_block.end:
		                new_end = max(last_block.end, block.end)
		                merged_blocks[-1] = BedBlock(last_block.start, new_end)
		            else:
		                # No overlap, so just add the current block
		                merged_blocks.append(block)
		self.blocks = merged_blocks


def ids_from_gtf(descrField, gdict):
	'''Parse the last field of a gtf to extract transcript and gene IDs'''
	pairs = descrField.split("; ")
	data_dict = {pair.split(" ", 1)[0].strip(): pair.split(" ", 1)[1].strip('"') for pair in pairs}
	# gene_id, transcript_id, gene_type, transcript_type, gene_ens_id, transcript_ens_id, gene_ens_id, gene_parent_id, transcript_parent_id, protein_parent_id
	# note: data_dict['gene_type'] is always identical data_dict['transcript_type'], I checked.
	hugo = 'NA'
	if data_dict['transcript_parent_id'] in gdict:
		hugo = gdict[data_dict['transcript_parent_id']]
	itemRgb = '255,140,0' # dark orange (pseudogene)
	if data_dict['gene_type'] == 'unprocessed_pseudogene':
		itemRgb = '0,0,255' # blue
	elif data_dict['gene_type'] == 'processed_pseudogene':
		itemRgb = '85,107,47' # dark olive green 
	extrafields = tuple([hugo, data_dict['gene_type'], data_dict['gene_id'], data_dict['gene_ens_id'], data_dict['gene_parent_id'], 
		data_dict['transcript_parent_id'], data_dict['protein_parent_id']])
	return itemRgb, data_dict['transcript_id'], extrafields


def gtf_to_bed(outputfile, gtf, gdict):
	geneObj = None
	with open(outputfile, 'wt') as outfile:
		writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)

		prev_transcript, blockstarts, blocksizes, prev_gene, prev_chrom, prev_strand = [None, None, None, None, None, None]
		blockList = []
		# extract all exons from the gtf, keep exons grouped by transcript
		for line in open(gtf):  
			if line.startswith('#'):
				continue
			line = line.rstrip().split('\t')
			chrom, ty, start, end, strand = line[0], line[2], int(line[3]) - 1, int(line[4]), line[6]
			if ty in ['gene', 'transcript', 'UTR', 'start_codon']:
				continue
			if not chrom in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']:
				continue 
			chrom = 'chr'+ chrom

			itemRgb, tx, extrafields = ids_from_gtf(line[8], gdict)
			if geneObj is None:
				geneObj = Gene(itemRgb, tx, chrom, start, end, strand, extrafields)
			elif geneObj.txID == tx:
				if ty == 'exon':
					geneObj.add(start, end)
			# this line is a new gene, so print the previous
			else:
				geneObj.write(outfile)
				geneObj = Gene(itemRgb, tx, chrom, start, end, strand, extrafields)
				cdsObj = None

		# last entry
		geneObj.write(outfile)


if __name__ == "__main__":
    main()
