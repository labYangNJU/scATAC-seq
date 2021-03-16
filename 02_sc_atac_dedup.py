## this script was downloaded from https://github.com/shendurelab/fly-atac with a little change
import pysam #v0.8.1
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]
outlog = sys.argv[3]

readsin = pysam.AlignmentFile(inbam,"rb")
readsout = pysam.AlignmentFile(outbam,"wb",template=readsin)
refs = readsin.references
unique_count = 0
total = 0
for refchrom in refs:
	if 'chrM' in refchrom or 'chrGL' in refchrom or 'chrNC' in refchrom or 'chrhs' in refchrom or 'random' in refchrom or 'chrU' in refchrom:
		continue
	readdic = {}
	print "Deduplicating " + refchrom + "..."
	for read in readsin.fetch(refchrom):
		readname = read.qname.split(':')[0]
		if 'CTF' in readname or 'AMBIG' in readname:
			continue
		if read.tlen < 0:
			fragstart = str(read.mpos - read.tlen)
			fragend = str(read.mpos)
		else:
			fragstart = str(read.pos)
			fragend = str(read.pos + read.tlen)
		try:
			readdic[readname + fragstart + fragend + unicode(read.is_read1)]
			total += 1
		except KeyError:
			readdic[readname + fragstart + fragend + unicode(read.is_read1)] = read
			unique_count += 1
			total += 1
			readsout.write(read)
logout = open(outlog,'w')
print >> logout, 'total=' + str(total) + '\tunique_count=' + str(unique_count)
logout.close()
