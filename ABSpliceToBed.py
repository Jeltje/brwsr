#!/usr/bin/env python3

import sys, os, re, argparse, textwrap
import csv

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

Convert ABSplice output to bed so that the highest score is presented (color coded)
and with a field that lists the top three tissues, and another field that lists all tissues and values
so that it can be converted to mouseover and a click-through table.
We also need ref and alt allele fields.
Todo: figure out if pos is 0 or 1 based
bedToBigBed -type=bed9+1 -tab -as=abSplice.as sorted.bed chrom.sizes output.bb

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('inputfile', type=str, help='DESCRIBE ME')
# optional flag
parser.add_argument('--outfile', type=str, default='tmp.bed',  help="Output filename (default tmp.bed)")
parser.add_argument('--strands', type=str, default='gene.strands',  help="Tab separated file with strand and gene ID (default gene.strands)")
parser.add_argument('-d', '--debug', help="Optional debugging output", action='store_true')


if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()

def itemRGB(score):
    '''Return color based on score'''
    rgb = '255,255,255' # black
    # blue, violet, pink, red
    for cutoff, color in [[0.5, '0,0,255'], [0.3, '170,0,204'], [0.1, '230,0,191'], [0.01, '255,0,0']]:
        if score >= cutoff:
            return color
    return rgb
# Main

# read the strands
strands = dict()
with open(args.strands, 'r') as infile:
    for line in infile:
        strand, gene = line.strip().split('\t')
        strands[gene] = strand


# Open the input file with the csv.DictReader
ct = 0
print(args.inputfile)
with open(args.inputfile, 'r', newline='', encoding='utf-8') as infile:
    # Create a csv.DictReader object with tab as the delimiter
    reader = csv.DictReader(infile, delimiter='\t')

    # Get the header from the input file
    header = reader.fieldnames

    # Find the indices of columns starting with 'AbSplice_DNA'
    indices_to_keep = [i for i, col in enumerate(header) if col.startswith('AbSplice_DNA')]

    # Open the output file with the csv.writer
    with open(args.outfile, 'w', newline='', encoding='utf-8') as outfile:
        # Create a csv.writer object with tab as the delimiter
        writer = csv.writer(outfile, delimiter='\t')

        # Write the updated header to the output file
        if args.debug:
            print('track name=ABsplice_test1 description="ABsplice_test1" itemRgb="On" visibility="dense"')
        else:
            print(f'remember to sort the output like so: sort -k1,1 -k2,2n {args.outfile} > sorted.bed')

        # Iterate over rows in the input file and write selected columns to the output file
        for row in reader:
            # omit in-file header rows (we're working on a concatenated file)
            if row['chrom'] == 'chrom':
                continue
            # Get the index and maximum value for each column in indices_to_keep
            #max_values = [(header[i].replace('DNA_', ''), float(row[header[i]])) for i in indices_to_keep]
            max_values = [(header[i].replace('DNA_', ''), float(row[header[i]]) if row[header[i]] else 0) for i in indices_to_keep]


            # Find the maximum value
            max_value = max(max_values, key=lambda x: x[1])[1]

            # Filter max_values to include only entries with the maximum value
            max_entries = [(column, value) for column, value in max_values if value == max_value]
            topValString=''
            name = f"{row['ref']}/{row['alt']}"

            # Write the filtered entries to the output file
            for column, value in max_entries:
                topValString += f'{column}: {value} <br> '
# the score must be 0-1000 so either we must convert or not use
# the strand is of interest, but must be pulled in from another source if we want it
# I think the mouseover should show the score,right. And color the intensity.
            startpos = int(row['pos'])
            if args.debug:
                if ct == 15:
                    sys.exit()
                ct += 1
                print(f"{row['chrom']}\t{startpos}\t{startpos+1}\t{name}\t0\t{strands[row['gene_id']]}\t{startpos}\t{startpos}\t{itemRGB(max_value)}\t{row['gene_id']}")
            writer.writerow([row['chrom'], startpos, startpos+1, name, 0, strands[row['gene_id']], startpos, startpos, itemRGB(max_value), row['gene_id']])

#topValString)

