## this script was downloaded from https://github.com/shendurelab/fly-atac with a little change
#!/usr/bin/env python
import argparse
import subprocess
import sys
sys.path.append('~/bin/Levenshtein/')
import Levenshtein
import gzip
import io
import cStringIO
io_method = cStringIO.StringIO

parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
parser.add_argument('-1','--Read1', help='Input fastq file 1',dest='input1',required=True)
parser.add_argument('-2','--Read2', help='Input fastq file 2',dest='input2',required=True)
parser.add_argument('-3','--barcode', help='barcode file',dest='input3',required=True)
parser.add_argument('-O1','--output1', help='Output fastq file 1',dest='output1',required=True)
parser.add_argument('-O2','--output2', help='Output fastq file 2',dest='output2',required=True)
parser.add_argument('-O3','--output3', help='Output fastq file 3',dest='output3',required=True)
parser.add_argument('-L','--log', help='Output log file',dest='logfile',required=True)
parser.add_argument('-X','--nextseq',help='NextSeq run indicator',dest='nextseq',action="store_true")
parser.add_argument('-Z','--gzip',help='Gzip indicator',dest='gzip',action="store_true")
args = parser.parse_args()

f = open('737K-cratac-v1.txt', 'r')
a = ""
refbarc = {}
for line in f.readlines():
    line = line.strip()
    if not len(line):
        continue
    refbarc[line.split()[0]] = a
f.close()


def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G','N':'N'}
def reverse_complement(x):
	xrev = x[::-1]
	xrevcomp = ''.join([complements[z] for z in xrev])
	return xrevcomp


totreads = 0
exactmatch = 0
failed = 0
readcount = 0
prevouts = 0

outfasta1 = open(args.output1,'w')
outfasta2 = open(args.output2,'w')
outfasta3 = open(args.output3,'w')
if1 = gzip.open(args.input1,'rb')
infasta1 = io.BufferedReader(if1)
if2 = gzip.open(args.input2,'rb')
infasta2 = io.BufferedReader(if2)
if3 = gzip.open(args.input3,'rb')
infasta3 = io.BufferedReader(if3)
reads1 = []
reads2 = []
#barcode_dic = {}
#barcode_check = {}
'''
for line in infasta3:
	name3= line.strip('@').split(':')[0]
	barcode0 = infasta3.next().strip()
	plus3 = infasta3.next()
	qual3 = infasta3.next().strip()
	barcode_re = reverse_complement(barcode0)
	try:
		barcode_dic[barcode_re] += 1
	except KeyError:
		barcode_dic[barcode_re] = 1
if3.close()
barcode_len = len(barcode_dic)
print(barcode_len)
for keys in barcode_dic:
	outfasta3.writelines(keys + '\n' + str(barcode_dic[keys]) + '\n')
outfasta3.close()
'''
#for k in barcode_dic.keys():
#try:
#refbarc[k]
#barcode_check[k] = [k,"0"]
#except KeyError:
#barcode_ed = editcheck(k,refbarc)
#seqbarc = barcode_ed[0]
#seqed = barcode_ed[1]
#barcode_check[k] = [seqbarc,seqed]

#if3 = gzip.open(args.input3,'rb')
#infasta3 = io.BufferedReader(if3)
for line in infasta1:
	name1= line.strip('@').split(':')[0]
	read1 = infasta1.next().strip()
	plus1 = infasta1.next()
	qual1 = infasta1.next().strip()
	name2 = infasta2.readline()
	read2 = infasta2.readline().strip()
	plus2 = infasta2.readline()
	qual2 = infasta2.readline().strip()
	name3 = infasta3.readline()
	barcode0 = infasta3.readline().strip()
	plus3 = infasta3.readline()
	qual3 = infasta3.readline().strip()
	#barcode_re = reverse_complement(barcode0)
	barcode_re = barcode0
	try:
		refbarc[barcode_re]
		reads1.append('@' + barcode_re + ':' + str(totreads) + '#' + '0/1\n' + read1 + '\n+\n' + qual1 + '\n')
		reads2.append('@' + barcode_re + ':' + str(totreads) + '#' + '0/2\n' + read2 + '\n+\n' + qual2 + '\n')
		readcount += 1
		totreads += 1
		exactmatch += 1
	except KeyError:
		failed += 1
		totreads += 1

	if readcount == 250000:		
		outfasta1.writelines(reads1)
		outfasta2.writelines(reads2)
		readcount = 0
		reads1 = []
		reads2 = []
if readcount > 0:
	outfasta1.writelines(reads1)
	outfasta2.writelines(reads2)

if1.close()
if2.close()
if3.close()
outfasta1.close()
outfasta2.close()


logout = open(args.logfile,'w')
print >> logout,'total=' + str(totreads) + '\texact=' + str(round(float(exactmatch)*100/totreads,2)) + '%\tfail=' + str(round(float(failed)*100/totreads,2)) + '%'
logout.close()

if args.gzip:
	zipper = 'gzip -f ' + args.output1 + '; gzip -f ' + args.output2 + ';'
	submitter(zipper)

