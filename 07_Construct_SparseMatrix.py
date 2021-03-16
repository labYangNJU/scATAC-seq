## this script was downloaded from https://github.com/shendurelab/fly-atac with a little change
import sys
import pysam
import numpy 
import numpy as np
from scipy.sparse import coo_matrix
import pandas as pd


#'Usage: python Construct_SparseMatrix.py [Input Bam file] [Input Index table] [Window BED] [cell index Output] [region index output] [count output]')

inbam = sys.argv[1]
indextable = sys.argv[2]
windowbed = sys.argv[3]
cell_out = sys.argv[4]
region_out = sys.argv[5]
matrix_put = sys.argv[6]

print "readin cell ..."
descer = open(indextable,'r')
cells = [x.strip().split()[0] for x in descer.readlines() if '@' not in x]
descer.close()
print "convert cells to cell indexes ..."
cellsdic = {}
cellsdic_reverse = {}
for x,cell in enumerate(cells):
    cellsdic[cell] = x+1   # +1 for construct sparseMatrix
cellsdic_reverse = {v : k for k, v in cellsdic.items()}  #change key and value
outter = open(cell_out,'w')
print >> outter, 'Barcode\tCellIndex'
for num in cellsdic_reverse.keys():
	print >> outter, str(num) + '\t' + cellsdic_reverse[num]
outter.close()

print "readin window ..."
desrer = open(windowbed,'r')
regions = [line.strip() for line in desrer.readlines()]
desrer.close()
print "convert regions to region indexes ..."
regionsdic = {}
regionsdic_reverse = {}
anno = {}
for x,region in enumerate(regions):
	regionsdic[region] = x+1
regionsdic_reverse = {v : k for k, v in regionsdic.items()}	
for region in regionsdic_reverse.values():
	rec=region.split()[0:3]
	anno[regionsdic[region]] = rec[0] + "_" + rec[1] + "_" + rec[2]
outter = open(region_out,'w')
#print >> outter, 'Region\tRegionIndex'
for index in regionsdic_reverse.keys():
	print >> outter, str(index) + '\t' + regionsdic_reverse[index] + '\t' + anno[index]
outter.close()

print "read in bam"
bamfile = pysam.Samfile(inbam,'rb')
	
print "Counting window reads for each cell..."
outter = open(matrix_put,'w')
print >> outter, 'Region\tbarcode\tcount'
for region in regionsdic.keys():
	rec=region.split()[0:3]
	reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
	currcounts = {}
	for read in reads:
		readname = read.qname.split(':')[0]
		#readname = read.qname.split('.')[0]
		#readname = read.get_tags()
		#readname=readname[7][1]
		try:
			cellsdic[readname]  # make sure readname in cell barcode
			try:
				currcounts[readname] += 1   # currcounts save read count(in each cell) in current region
			except KeyError:
				currcounts[readname] = 1
		except KeyError:
			continue
	if len(currcounts)>0:
		curr_col = []
		curr_data = []
		curr_row = regionsdic[region] * len(currcounts)
		'''
		for key in currcounts.keys():
			curr_col.append(cellsdic[key])
			curr_data.append(currcounts[key])
		row.append(curr_row)
		col.extend(curr_col)
		data.extend(curr_data)
		'''
		for key in currcounts.keys():
			print >> outter,str(regionsdic[region]) + '\t' + str(cellsdic[key]) + '\t' + str(currcounts[key])
outter.close()		
