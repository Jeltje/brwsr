#!/usr/bin/env python3

import sys, os, re, argparse, textwrap
import csv

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\

Convert ABSplice output to bed so that the highest score is presented as a color coded block with both alleles as name.
A mouseOver field lists the top score and all tissues that have this score, 
and another field presents all tissues and values in a HTML table format that will be shown
on click through.
Strand information is pulled in from a tab separated file extracted from gencode v38 (gene ID, strand)
bedToBigBed -type=bed9+3 -tab -as=abSplice.as sorted.bed chrom.sizes output.bb

        '''))
group = parser.add_argument_group('required arguments')
group.add_argument('inputfile', type=str, help='DESCRIBE ME')
# optional flag
parser.add_argument('--outfile', type=str, default='tmp.bed',  help="Output filename (default tmp.bed)")
parser.add_argument('--strands', type=str, default='gene.strands',  help="Tab separated file with strand and gene ID (default gene.strands)")
parser.add_argument('-d', '--debug', help="Optional debugging output", action='store_true')

# the score field must be 0-1000 so either we must convert or not use

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
    header = [column.replace('ABSplice_DNA_', '') for column in header]

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
            # Get the index and value for each column in indices_to_keep
            all_values = [(header[i], row[header[i]]) for i in indices_to_keep]
            # turn this information into a html table
            html_table = '<table>'
            for tvals in all_values:
                html_table += f'  <tr><td>{tvals[0]}</td><td>{tvals[1]}</td></tr>'
            html_table += '</table>'

            # Find the maximum value (first make sure all row entries are floats)
            all_values = [(x, float(y) if y else 0) for x, y in all_values]
            max_value = max(all_values, key=lambda x: x[1])[1]

            # Filter max_values to include only entries with the maximum value
            max_entries = [column for column, value in all_values if value == max_value]

            # mouseover information
            topValString = f'Max score {max_value} in <br>'
            topValString += '<br>'.join(max_entries)

            # this will be the item label
            name = f"{row['ref']}>{row['alt']}"

            # AB coordinates appear to be 1-based
            startpos = int(row['pos']) -1 
            if args.debug:
                if ct == 15:
                    sys.exit()
                ct += 1
                print(f"{row['chrom']}\t{startpos}\t{startpos+1}\t{name}\t0\t{strands[row['gene_id']]}\t{startpos}\t{startpos}\t{itemRGB(max_value)}\t{row['gene_id']}\t{topValString}\t{html_table}")
            writer.writerow([row['chrom'], startpos, startpos+1, name, 0, strands[row['gene_id']], startpos, startpos, itemRGB(max_value), row['gene_id'], topValString, html_table])


