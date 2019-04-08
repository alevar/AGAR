import os
import sys
import shutil
import subprocess

def main(args):
	for ifp in args.m1.split(","):
		assert os.path.exists(os.path.abspath(ifp)),"#1 reads not found"
	for ifp in args.m2.split(","):
		assert os.path.exists(os.path.abspath(ifp)),"#2 reads not found"
	if (not args.single==None):
		assert os.path.exists(os.path.abspath(args.single)),"singletons file does not exist"
	assert os.path.exists(os.path.abspath(args.db)),"database directory not found"

	if args.type=="hisat":
		for fp in ["1.ht2","2.ht2","3.ht2","4.ht2","5.ht2","6.ht2","7.ht2","8.ht2"]:
			assert os.path.exists(os.path.abspath(args.db)+"/db."+fp),args.db+"/db."+fp+" file not found"
	if args.type=="bowtie":
		for fp in ["1.bt2","2.bt2","3.bt2","4.bt2","rev.1.bt2","rev.2.bt2"]:
			assert os.path.exists(os.path.abspath(args.db)+"/db."+fp),args.db+"/db."+fp+" file not found"
	for fp in ["fasta.tlst","fasta","genome.header"]:
		assert os.path.exists(os.path.abspath(args.db)+"/db."+fp),args.db+"/db."+fp+" file not found"
	for fp in ["1.ht2","2.ht2","3.ht2","4.ht2","5.ht2","6.ht2","7.ht2","8.ht2"]:
		assert os.path.exists(args.genome_db+"."+fp),args.genome_db+"."+fp+" file not found"

	curTMP=args.tmp
	if not os.path.exists(args.tmp):
		os.mkdir(args.tmp)
	else:
		counter=0
		while True:
			curTMP=os.path.abspath(args.tmp)+"_"+str(counter)
			if not os.path.exists(curTMP):
				os.mkdir(curTMP)
				break
			else:
				counter+=1
				continue
	
	unalignedR1="" # keep track of the current unaligned reads governed by what stages of the run are active
	unalignedR2=""
	unalignedS =""
	alignmentStage=0 # indicates how many stages have been performed on the current run - important for the cleanup
	transcriptome_process=None; # process in which transcriptome alignment is performed. it opens a pipe to samtools
	genome_process=None; # process in which genome alignment is performed. it also opens a pipe to samtools
	#hisat2
	if args.type=="hisat":
		print("aligning with hisat2 against transcriptome")
		# perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments with the "-a" option enabled
		hisat2_cmd_trans_noA=["hisat2",
						      "--no-spliced-alignment",
						      "--end-to-end",
						      "--no-unal",
						      "-x",os.path.abspath(args.db)+"/db",
						      "-p",args.threads]
		if (args.fasta):
			hisat2_cmd_trans_noA.append("-f")
		hisat2_cmd_trans_noA.extend(("-1",args.m1,
						 			 "-2",args.m2))
		if (not args.single==None):
			hisat2_cmd_trans_noA.extend(("-U",args.single))
		hisat2_cmd_trans_noA.extend(("--un-conc",os.path.abspath(curTMP)+"/sample.trans.unconc_first.fq",
									 "--un",os.path.abspath(curTMP)+"/sample.trans.un_first.fq"))
		if (args.hisat):
			hisat2_cmd_trans_noA.extend(args.hisat)
		transcriptome_process=subprocess.Popen(hisat2_cmd_trans_noA,stdout=subprocess.PIPE)
		unalignedR1=os.path.abspath(curTMP)+"/sample.trans.unconc_first.1.fq"
		unalignedR2=os.path.abspath(curTMP)+"/sample.trans.unconc_first.2.fq"
		unalignedS =os.path.abspath(curTMP)+"/sample.trans.un_first.fq"
		alignmentStage=1

		if args.all:
			print("Aligning with hisat2 against transcriptome with -a")
			hisat2_cmd_trans_A=["hisat2",
							    "--no-spliced-alignment",
							    "--end-to-end",
							    "--no-unal",
							    "-a",
							    "-x",os.path.abspath(args.db)+"/db",
							    "-p",args.threads]
			if (args.fasta):
				hisat2_cmd_trans_A.append("-f")
			hisat2_cmd_trans_A.extend(("-1",os.path.abspath(curTMP)+"/sample.trans.unconc_first.1.fq",
									   "-2",os.path.abspath(curTMP)+"/sample.trans.unconc_first.2.fq",
									   "-U",os.path.abspath(curTMP)+"/sample.trans.un_first.fq",
									   "-S",os.path.abspath(curTMP)+"/sample.trans_second.sam",
									   "--un-conc",os.path.abspath(curTMP)+"/sample.trans.unconc.fq",
									   "--un",os.path.abspath(curTMP)+"/sample.trans.un.fq"))
			if (args.hisat):
				hisat2_cmd_trans_A.extend(args.hisat)
			subprocess.call(hisat2_cmd_trans_A)
			unalignedR1=os.path.abspath(curTMP)+"/sample.trans.unconc.1.fq"
			unalignedR2=os.path.abspath(curTMP)+"/sample.trans.unconc.2.fq"
			unalignedS =os.path.abspath(curTMP)+"/sample.trans.un.fq"
			alignmentStage=2

	elif args.type=="bowtie":
		print("aligning with bowtie2 against transcriptome")
		# perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments with the "-a" option enabled
		bowtie2_cmd_trans_noA=["bowtie2",
						      "--end-to-end",
						      "--no-unal",
						      "-x",os.path.abspath(args.db)+"/db",
						      "-p",args.threads]
		if (args.k):
			bowtie2_cmd_trans_noA.extend(("-k",str(args.k)))

		if (args.fasta):
			bowtie2_cmd_trans_noA.append("-f")
		bowtie2_cmd_trans_noA.extend(("-1",args.m1,
						 			 "-2",args.m2))
		if (not args.single==None):
			bowtie2_cmd_trans_noA.extend(("-U",args.single))
		bowtie2_cmd_trans_noA.extend(("--un-conc",os.path.abspath(curTMP)+"/sample.trans.unconc_first.fq",
									 "--un",os.path.abspath(curTMP)+"/sample.trans.un_first.fq"))
		if (args.bowtie):
			bowtie2_cmd_trans_noA.extend(args.bowtie)
		transcriptome_process=subprocess.Popen(bowtie2_cmd_trans_noA,stdout=subprocess.PIPE)
		unalignedR1=os.path.abspath(curTMP)+"/sample.trans.unconc_first.1.fq"
		unalignedR2=os.path.abspath(curTMP)+"/sample.trans.unconc_first.2.fq"
		unalignedS =os.path.abspath(curTMP)+"/sample.trans.un_first.fq"
		alignmentStage=1

		if args.all:
			print("Aligning with bowtie2 against transcriptome with -a")
			bowtie2_cmd_trans_A=["bowtie2",
							    "--end-to-end",
							    "--no-unal",
							    "-a",
							    "-x",os.path.abspath(args.db)+"/db",
							    "-p",args.threads]
			if (args.fasta):
				bowtie2_cmd_trans_A.append("-f")
			bowtie2_cmd_trans_A.extend(("-1",os.path.abspath(curTMP)+"/sample.trans.unconc_first.1.fq",
									   "-2",os.path.abspath(curTMP)+"/sample.trans.unconc_first.2.fq",
									   "-U",os.path.abspath(curTMP)+"/sample.trans.un_first.fq",
									   "-S",os.path.abspath(curTMP)+"/sample.trans_second.sam",
									   "--un-conc",os.path.abspath(curTMP)+"/sample.trans.unconc.fq",
									   "--un",os.path.abspath(curTMP)+"/sample.trans.un.fq"))
			if (args.bowtie):
				bowtie2_cmd_trans_noA.extend(args.bowtie)
			subprocess.call(bowtie2_cmd_trans_A)
			unalignedR1=os.path.abspath(curTMP)+"/sample.trans.unconc.1.fq"
			unalignedR2=os.path.abspath(curTMP)+"/sample.trans.unconc.2.fq"
			unalignedS =os.path.abspath(curTMP)+"/sample.trans.un.fq"
			alignmentStage=2

	subprocess.Popen(["samtools","view","-h","--output-fmt=BAM","-@",args.threads,"-o",os.path.abspath(curTMP)+"/sample.trans_first.bam"],stdin=transcriptome_process.stdout)
	transcriptome_process.wait()
	transcriptome_process.stdout.close()

	print("aligning with hisat2 against the genome")
	hisat2_cmd_genome=["hisat2",
					   "--very-sensitive",
                       "--dta",
					    "--no-unal",
					   "-x",os.path.abspath(args.genome_db),
					   "-p",str(int(args.threads)-1)]
	if (args.fasta):
		hisat2_cmd_genome.append("-f")
	hisat2_cmd_genome.extend(("-1",unalignedR1,
							  "-2",unalignedR2,
							  "-U",unalignedS))
	if (args.hisat):
		hisat2_cmd_genome.extend(args.hisat)
	genome_process=subprocess.Popen(hisat2_cmd_genome,stdout=subprocess.PIPE)
	subprocess.Popen(["samtools","view","-h","--output-fmt=BAM","-@",args.threads,"-o",os.path.abspath(curTMP)+"/sample.genome.bam"],stdin=genome_process.stdout)
	
	#trans2genome
	print("converting coordinates to genomic")
	trans2genome_cmd=["/home/avaraby1/genomicTools/trans2genome/trans2genome",
					 "-t",os.path.abspath(args.db)+"/db.fasta.tlst",
					 "-g",os.path.abspath(args.db)+"/db.fasta.glst",
					 "-m",os.path.abspath(args.db)+"/db.fasta.multi",
					 "-i",os.path.abspath(curTMP)+"/sample.trans_first.bam",
					 "-s",os.path.abspath(args.db)+"/db.genome.header",
					 "-o",os.path.abspath(curTMP)+"/sample.trans2genome_first.bam"]
	if args.mf:
		trans2genome_cmd.append("-f")
	subprocess.call(trans2genome_cmd)
	if os.path.exists(os.path.abspath(curTMP)+"/sample.trans_first.sam") and not args.keep:
		os.remove(os.path.abspath(curTMP)+"/sample.trans_first.sam")

	if os.path.exists(os.path.abspath(curTMP)+"/sample.trans_second.bam"):
		trans2genome_a_cmd=["/home/avaraby1/genomicTools/trans2genome/trans2genome",
						 "-t",os.path.abspath(args.db)+"/db.fasta.tlst",
						 "-g",os.path.abspath(args.db)+"/db.fasta.glst",
						 "-m",os.path.abspath(args.db)+"/db.fasta.multi",
						 "-i",os.path.abspath(curTMP)+"/sample.trans_second.bam",
						 "-s",os.path.abspath(args.db)+"/db.genome.header",
						 "-o",os.path.abspath(curTMP)+"/sample.trans2genome_second.bam"]
		if args.mf:
			trans2genome_a_cmd.append("-f")
		subprocess.call(trans2genome_a_cmd)


	genome_process.wait() # allows trans2genome to run at the same time as hisat2
	genome_process.stdout.close()

	if not args.keep:
		if os.path.exists(os.path.abspath(curTMP)+"/sample.trans.unconc_first.1.fq"):
			os.remove(os.path.abspath(curTMP)+"/sample.trans.unconc_first.1.fq")
		if os.path.exists(os.path.abspath(curTMP)+"/sample.trans.unconc_first.2.fq"):
			os.remove(os.path.abspath(curTMP)+"/sample.trans.unconc_first.2.fq")
		if os.path.exists(os.path.abspath(curTMP)+"/sample.trans.un_first.fq"):
			os.remove(os.path.abspath(curTMP)+"/sample.trans.un_first.fq")

	if alignmentStage==2:
		subprocess.call(["samtools","view","-h","--output-fmt=BAM","-@",args.threads,"-o",os.path.abspath(curTMP)+"/sample.trans_second.bam",os.path.abspath(curTMP)+"/sample.trans_second.sam"])
		if not args.keep:
			os.remove(os.path.abspath(curTMP)+"/sample.trans_second.sam")

	print("merging all sub-alignments")
	merge_cmd=["samtools","merge",
			   "-@",args.threads,
			   args.output,
			   os.path.abspath(curTMP)+"/sample.trans2genome_first.bam"]
	if os.path.exists(os.path.abspath(curTMP)+"/sample.trans2genome_second.bam"):
		merge_cmd.append(os.path.abspath(curTMP)+"/sample.trans2genome_second.bam")

	merge_cmd.append(os.path.abspath(curTMP)+"/sample.genome.bam")

	if len(merge_cmd)>6:
		subprocess.call(merge_cmd)
	else:
		if os.path.exists(os.path.abspath(curTMP)+"/sample.genome.bam") and not args.keep:
			os.remove(os.path.abspath(curTMP)+"/sample.genome.bam")
		os.rename(os.path.abspath(curTMP)+"/sample.trans2genome_first.bam",args.output)
		if not args.keep:
			os.remove(os.path.abspath(curTMP)+"/sample.trans2genome_first.bam")

	if not args.keep:
		shutil.rmtree(os.path.abspath(curTMP))
