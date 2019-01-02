import os
import sys
import subprocess

def buildHeader(fastaIDX,outputHeaderFP):
	headerFP=open(outputHeaderFP,"w+")
	headerFP.write("@HD\tVN:1.0\tSO:unsorted\n")
	with open(os.path.abspath(fastaIDX)) as inFP:
	    for line in inFP.readlines():
	        lineCols=line.split("\t")
	        newTAG="\t".join(["@SQ","SN:"+lineCols[0],"LN:"+lineCols[1]])+"\n"
	        headerFP.write(newTAG)
	headerFP.write('@PG\tID:transHISAT2\tPN:transHISAT2\tVN:1.0\tCL:"./transHISAT2.py"\n')
	headerFP.close()

def main(args):
	assert os.path.exists(os.path.abspath(args.gff)),"gff annotation not found"
	assert os.path.exists(os.path.abspath(args.ref)),"reference file not found"
	if not os.path.exists(os.path.abspath(args.ref)+".fai"):
		print("FASTA index for the reference genome not found. Building now.")
		subprocess.call(["samtools","faidx",args.ref])

	if not os.path.exists(args.output):
		os.mkdir(args.output)

	# gtf_to_fasta
	print("Extracting fasta from gtf")
	subprocess.call(["/home/avaraby1/genomicTools/tophat/bin/gtf_to_fasta",args.gff,args.ref,os.path.abspath(args.output)+"/db.fasta"])
	# buildGenomeHeader.py
	print("Building genome header file")
	buildHeader(os.path.abspath(args.ref)+".fai",os.path.abspath(args.output)+"/db.genome.header")

	if args.type=="bowtie":
		print("Building additional bowtie transcriptome index")
		subprocess.call(["bowtie2-build",os.path.abspath(args.output)+"/db.fasta",os.path.abspath(args.output)+"/db"])
	elif args.type=="hisat":
		print("Building transcriptome database for HISAT2")
		subprocess.call(["hisat2-build","-p",args.threads,os.path.abspath(args.output)+"/db.fasta",os.path.abspath(args.output)+"/db"])
	else:
		print("unknown aligner specified")
