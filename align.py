import os
import sys
import shutil
import subprocess

def main(args):
	assert os.path.exists(os.path.abspath(args.m1)),"#1 reads not found"
	assert os.path.exists(os.path.abspath(args.m2)),"#2 reads not found"
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
	
	#hisat2
	if args.type=="hisat":
		print("aligning with hisat2 against transcriptome")
		# perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments with the "-a" option enabled
		hisat2_cmd_trans_noA=["hisat2",
						      "--no-spliced-alignment",
						      "--very-sensitive",
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
		hisat2_cmd_trans_noA.extend(("-S",os.path.abspath(curTMP)+"/sample.trans_first.sam",
									 "--un-conc",os.path.abspath(curTMP)+"/sample.trans.unconc_first.fq",
									 "--un",os.path.abspath(curTMP)+"/sample.trans.un_first.fq"))
		subprocess.call(hisat2_cmd_trans_noA)

		print("Aligning with hisat2 against transcriptome with -a")
		hisat2_cmd_trans_A=["hisat2",
						    "--no-spliced-alignment",
						    "--very-sensitive",
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
		subprocess.call(hisat2_cmd_trans_A)

	elif args.type=="bowtie":
		print("aligning with bowtie2 against transcriptome")
		# perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments with the "-a" option enabled
		hisat2_cmd_trans_noA=["bowtie2",
						      "--very-sensitive",
						      "--end-to-end",
						      "--no-unal",
                              "-a",
						      "-x",os.path.abspath(args.db)+"/db",
						      "-p",args.threads]
		if (args.fasta):
			hisat2_cmd_trans_noA.append("-f")
		hisat2_cmd_trans_noA.extend(("-1",args.m1,
						 			 "-2",args.m2))
		if (not args.single==None):
			hisat2_cmd_trans_noA.extend(("-U",args.single))
		hisat2_cmd_trans_noA.extend(("-S",os.path.abspath(curTMP)+"/sample.trans_first.sam",
									 "--un-conc",os.path.abspath(curTMP)+"/sample.trans.unconc_first.fq",
									 "--un",os.path.abspath(curTMP)+"/sample.trans.un_first.fq"))
		subprocess.call(hisat2_cmd_trans_noA)

		print("Aligning with bowtie2 against transcriptome with -a")
		hisat2_cmd_trans_A=["bowtie2",
						    "--very-sensitive",
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
		subprocess.call(hisat2_cmd_trans_A)


	print("aligning with hisat2 against the genome")
	hisat2_cmd_genome=["hisat2",
					   "--very-sensitive",
					   "--no-unal",
					   "-x",os.path.abspath(args.genome_db),
					   "-p",args.threads]
	if (args.fasta):
		hisat2_cmd_genome.append("-f")
	hisat2_cmd_genome.extend(("-1",os.path.abspath(curTMP)+"/sample.trans.unconc.1.fq",
							  "-2",os.path.abspath(curTMP)+"/sample.trans.unconc.2.fq",
							  "-U",os.path.abspath(curTMP)+"/sample.trans.un.fq",
							  "-S",os.path.abspath(curTMP)+"/sample.genome.sam"))
	subprocess.call(hisat2_cmd_genome)

	#samtools view -h -b 
	print("converting alignments to BAM and sorting by read name")
	subprocess.call(["samtools","sort","-n","--output-fmt=BAM","-@",args.threads,"-o",os.path.abspath(curTMP)+"/sample.trans_first.bam",os.path.abspath(curTMP)+"/sample.trans_first.sam"])
	subprocess.call(["samtools","sort","-n","--output-fmt=BAM","-@",args.threads,"-o",os.path.abspath(curTMP)+"/sample.trans_second.bam",os.path.abspath(curTMP)+"/sample.trans_second.sam"])
	# these can be outsourced to a different thread for removal
	# in order to save some time
	os.remove(os.path.abspath(curTMP)+"/sample.trans_first.sam")
	os.remove(os.path.abspath(curTMP)+"/sample.trans_second.sam")
	
	#trans2genome
	print("converting coordinates to genomic")
	subprocess.call(["trans2genome",
					 "-g",os.path.abspath(args.db)+"/db.fasta.tlst",
					 "-i",os.path.abspath(curTMP)+"/sample.trans_first.bam",
					 "-s",os.path.abspath(args.db)+"/db.genome.header",
					 "-o",os.path.abspath(curTMP)+"/sample.trans2genome_first.bam"])
	# if os.path.exists(os.path.abspath(curTMP)+"/sample.trans_second.bam"):
	# 	subprocess.call(["trans2genome",
	# 					 "-g",os.path.abspath(args.db)+"/db.fasta.tlst",
	# 					 "-i",os.path.abspath(curTMP)+"/sample.trans_second.bam",
	# 					 "-s",os.path.abspath(args.db)+"/db.genome.header",
	# 					 "-o",os.path.abspath(curTMP)+"/sample.trans2genome_second.bam"])

	print("merging all sub-alignments")
	merge_cmd=["samtools","merge",
			   "-@",args.threads,
			   args.output,
			   os.path.abspath(curTMP)+"/sample.trans2genome_first.bam"]
	if os.path.exists(os.path.abspath(curTMP)+"/sample.trans2genome_second.bam"):
		merge_cmd.append(os.path.abspath(curTMP)+"/sample.trans2genome_second.bam")

	if os.path.exists(os.path.abspath(curTMP)+"/sample.genome.sam"):
		subprocess.call(["samtools","sort","-n","--output-fmt=BAM","-@",args.threads,"-o",os.path.abspath(curTMP)+"/sample.genome.bam",os.path.abspath(curTMP)+"/sample.genome.sam"])
		os.remove(os.path.abspath(curTMP)+"/sample.genome.sam")
		merge_cmd.append(os.path.abspath(curTMP)+"/sample.genome.bam")

	if len(merge_cmd)>6:
		subprocess.call(merge_cmd)
	else:
		os.rename(os.path.abspath(curTMP)+"/sample.trans2genome_first.bam",args.output)

	shutil.rmtree(os.path.abspath(curTMP))
