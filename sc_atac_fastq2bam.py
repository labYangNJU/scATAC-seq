import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='A program to convert NextSeq BCL files to fastq files for scATAC-seq analysis.')
parser.add_argument('-R1','--read1', help='R1 fastq file',dest='read1',required=True)
parser.add_argument('-R2','--read2', help='R2 fastq file',dest='read2',required=True)
parser.add_argument('-R3','--read3', help='barcode file',dest='read3',required=True)
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir',required=True)
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix',required=True)
parser.add_argument('-G','--genome',help='Bowtie genome',dest='genome',required=True)
args = parser.parse_args()

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

if args.outdir[-1] != '/':
	args.outdir = args.outdir + '/'

try:
	os.makedirs(args.outdir)
except OSError:
	print 'Outdir already exists...'

scriptdir = os.path.dirname(os.path.abspath(__file__))

print "Fixing barcodes..."
splitter = 'python ' + scriptdir + '/sc_atac_barcode_check_NoEdit.py -1 ' + args.read1 + ' -2 ' + args.read2 + ' -3 ' + args.read3 + ' -O1 ' + args.outdir + args.prefix + '_R1_fast.fastq -O2 ' + args.outdir + args.prefix + '_R2_fast.fastq -O3 ' + args.outdir + args.prefix + '_R3_fast.fastq' + ' -L ' + args.outdir + args.prefix + '.check_fast.log -X -Z'
print(splitter)
submitter(splitter)


print "Trimming adapters..."
trimmer = 'java -Xmx1G -jar ' + scriptdir + '/trimmomatic-0.32.jar PE -threads 8 ' + args.outdir + args.prefix + '_R1_fast.fastq.gz ' + args.outdir + args.prefix + '_R2_fast.fastq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz ILLUMINACLIP:' + scriptdir + '/NexteraPE-PE.fa:2:30:10:1:true TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:25 2> ' + args.outdir + args.prefix + '.split.trimmomatic.log'
submitter(trimmer)

print "Cleaning up..."
cleaner = 'rm ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz; rm ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz'
submitter(cleaner)


print "Mapping reads..."
mapper = "/Users/lab2116/Desktop/analysis/bowtie2-2.3.0-legacy/bowtie2 -p 11 -X 2000 -3 1 -x " + args.genome + " -1 " + args.outdir + args.prefix +".split.1.trimmed.paired.fastq.gz -2 " + args.outdir + args.prefix +".split.2.trimmed.paired.fastq.gz 2> " + args.outdir + args.prefix +".split.bowtie2.log | samtools view -bS - > " + args.outdir + args.prefix +".split.bam; samtools view -h -f 3 -F 12 -q 30 " + args.outdir + args.prefix +".split.bam | grep -v '[0-9]''\\t'chrM | grep -v '[0-9]''\\t'chrU | samtools view -Su - | samtools sort -@ 8 - -o " + args.outdir + args.prefix +".split.q30.sort.bam; samtools index " + args.outdir + args.prefix +".split.q30.sort.bam"
print(mapper)
submitter(mapper)


print "Deduplicating reads..."
dedup = "python " + scriptdir + "/sc_atac_true_dedup.py " + args.outdir + args.prefix +".split.q30.sort.bam " + args.outdir + args.prefix +".true.nodups.bam " + args.outdir + args.prefix +"_dedup.log"
submitter(dedup)

