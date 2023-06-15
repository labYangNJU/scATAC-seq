import sys
import pysam


#'Usage: python specific_barcode_selectBamFile.py [Input Bam file] [Input cluster_barcode list] [output_bamFile]')

inbam = sys.argv[1]
inClusterBarcode = sys.argv[2]
out1bam = sys.argv[3]
readsin = pysam.AlignmentFile(inbam,"rb")
readsout1 = pysam.Samfile(out1bam,"wb",template=readsin)

cluster_info = {}
descer = open(inClusterBarcode,'r')
for line in descer.readlines():
	cells = line.split()
	cluster_info[cells[0]] = ""
descer.close()



for read in readsin.fetch():
	readname = read.qname.split('.')[0]
	print readname
	try:
		cluster_info[readname]
		readsout1.write(read)
	except KeyError:
		continue
