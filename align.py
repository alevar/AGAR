import os
import sys
import shutil
import subprocess


# for the salmon mode - here are the commands
# ~/soft/salmon-latest_linux_x86_64/bin/salmon index -t db.fasta -i ./db
# ~/soft/salmon-latest_linux_x86_64/bin/salmon quant -l A -i ./data/salmonHisatIDX_p8/db -1 ./data/SRR1071717_r1.fq -2 ./data/SRR1071717_r2.fq --writeMappings ./outBladder_salmon/SRR1071717.sam --writeUnmappedNames -o ./outBladder_salmon/SRR1071717 -p 24 --validateMappings
# samtools sort -n --output-fmt=BAM -@ 3 -o ./outBladder_salmon/SRR1071717.sorted.bam ./outBladder_salmon/SRR1071717.bam

# /ccb/salz4-4/avaraby/tools_under_development/trans2genome_gtftofasta/salmon2genome -t ./salmonHisatIDX_p8/db.fasta.tlst -i ./outBladder_salmon/SRR1071717_salmon_stdout.sorted.bam -s ./salmonHisatIDX_p8/db.genome.header -o ./outBladder_salmon/SRR1071717.trans2genome.bam -n 154368 -a ./outBladder_salmon/SRR1071717/quant.sf

# cut -d' ' -f1 SRR1071717/aux_info/unmapped_names.txt | seqtk subseq ../data/SRR1071717_r2.fq - > SRR1071717_salmon_unmapped_r2.fq
# cut -d' ' -f1 SRR1071717/aux_info/unmapped_names.txt | seqtk subseq ../data/SRR1071717_r2.fq - > SRR1071717_salmon_unmapped_r2.fq

# hisat2 -x /home/gpertea/ncbi/dbGaP-10908/hg38/genome_tran -1 ./SRR1071717_salmon_unmapped_r1.fq -2 SRR1071717_salmon_unmapped_r2.fq --very-sensitive --no-unal -k 30 -p 24 -S ./SRR1071717_salmon_unmapped_hisat.sam

# samtools view -h --output-fmt=BAM -@ 24 ./SRR1071717_salmon_unmapped_hisat.sam -o ./SRR1071717_salmon_unmapped_hisat.bam

# samtools merge -@ 24 SRR1071717_salmon_hisat.bam SRR1071717_salmon_unmapped_hisat.bam SRR1071717.trans2genome.bam

# samtools sort --output-fmt=BAM -@ 24 -o ./SRR1071717_salmon_hisat.sorted.bam ./SRR1071717_salmon_hisat.bam

# stringtie ./SRR1071717_salmon_hisat.sorted.bam -p 24 -m 150 -G ../data/hg38_p8.biotype_flt.cls.gff3 -o SRR1071717_stringtie.gt

def main(args):
    for ifp in args.m1.split(","):
        assert os.path.exists(os.path.abspath(ifp)), "#1 reads not found"
    for ifp in args.m2.split(","):
        assert os.path.exists(os.path.abspath(ifp)), "#2 reads not found"
    if args.single is not None:
        assert os.path.exists(os.path.abspath(args.single)), "singletons file does not exist"
    assert os.path.exists(os.path.abspath(args.db)), "database directory not found"

    if args.type == "hisat":
        for fp in ["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"]:
            assert os.path.exists(os.path.abspath(args.db) + "/db." + fp), args.db + "/db." + fp + " file not found"
    if args.type == "bowtie":
        for fp in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
            assert os.path.exists(os.path.abspath(args.db) + "/db." + fp), args.db + "/db." + fp + " file not found"
    if args.locus:
        for fp in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
            assert os.path.exists(
                os.path.abspath(args.db) + "/db.locus." + fp), args.db + "/db." + fp + " file not found"
    # assert (bool(args.type=="hisat") != bool(args.bowtie)), "can not define both type hisat and bowtie2"
    for fp in ["tlst", "fasta", "genome_header", "info", "glst"]:
        assert os.path.exists(os.path.abspath(args.db) + "/db." + fp), args.db + "/db." + fp + " file not found"
    if args.mf:
        assert os.path.exists(os.path.abspath(args.db) + "/db.multi"), args.db + "/db.multi file not found"
    for fp in ["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"]:
        assert os.path.exists(args.genome_db + "." + fp), args.genome_db + "." + fp + " file not found"

    cur_tmp = args.tmp
    if not os.path.exists(args.tmp):
        os.mkdir(args.tmp)
    else:
        counter = 0
        while True:
            cur_tmp = os.path.abspath(args.tmp) + "_" + str(counter)
            if not os.path.exists(cur_tmp):
                os.mkdir(cur_tmp)
                break
            else:
                counter += 1
                continue

    unaligned_r1 = ""  # keep track of the current unaligned reads governed by what stages of the run are active
    unaligned_r2 = ""
    unaligned_s = ""
    transcriptome_process = None  # process in which transcriptome alignment is performed. it opens a pipe to samtools
    genome_process = None  # process in which genome alignment is performed. it also opens a pipe to samtools
    locus_process = None  # process in which alignment of non-transcriptomic reads is performed against loci (introns)
    if args.type == "hisat":
        print("aligning with hisat2 against transcriptome")
        # perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments
        # with the "-a" option enabled
        hisat2_cmd_trans_noA = ["hisat2",
                                "--no-spliced-alignment",
                                "--end-to-end",
                                "--no-unal",
                                "-x", os.path.abspath(args.db) + "/db",
                                "-p", args.threads]
        if args.fasta:
            hisat2_cmd_trans_noA.append("-f")
        hisat2_cmd_trans_noA.extend(("-1", args.m1,
                                     "-2", args.m2))
        if not args.single is None:
            hisat2_cmd_trans_noA.extend(("-U", args.single))
        hisat2_cmd_trans_noA.extend(("--un-conc", os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.fq",
                                     "--un", os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"))
        if args.hisat:
            hisat2_cmd_trans_noA.extend(args.hisat)
        transcriptome_process = subprocess.Popen(hisat2_cmd_trans_noA, stdout=subprocess.PIPE)
        unaligned_r1 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq"
        unaligned_r2 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq"
        unaligned_s = os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"

    elif args.type == "bowtie":
        print("aligning with bowtie2 against transcriptome")
        # perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments
        # with the "-a" option enabled
        bowtie2_cmd_trans_no_a = ["bowtie2",
                                  "--end-to-end",
                                  "--no-unal",
                                  "-x", os.path.abspath(args.db) + "/db",
                                  "-p", args.threads]
        if args.k:
            bowtie2_cmd_trans_no_a.extend(("-k", str(args.k)))
        else:
            bowtie2_cmd_trans_no_a.extend(("-k", "1"))

        if args.fasta:
            bowtie2_cmd_trans_no_a.append("-f")
        bowtie2_cmd_trans_no_a.extend(("-1", args.m1,
                                       "-2", args.m2))
        if not args.single is None:
            bowtie2_cmd_trans_no_a.extend(("-U", args.single))
        bowtie2_cmd_trans_no_a.extend(("--un-conc", os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.fq",
                                       "--un", os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"))
        if args.bowtie:
            bowtie2_cmd_trans_no_a.extend(args.bowtie)
        transcriptome_process = subprocess.Popen(bowtie2_cmd_trans_no_a, stdout=subprocess.PIPE)
        unaligned_r1 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq"
        unaligned_r2 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq"
        unaligned_s = os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"

    subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
                      os.path.abspath(cur_tmp) + "/sample.trans_first.bam"], stdin=transcriptome_process.stdout)

    transcriptome_process.wait()
    transcriptome_process.stdout.close()

    # trans2genome
    print("converting coordinates to genomic")
    trans2genome_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     'trans2genome')  # get path to the trans2genome that was compiled with the package
    trans2genome_cmd = [trans2genome_path,
                        "-x", os.path.abspath(args.db) + "/db",
                        "-i", os.path.abspath(cur_tmp) + "/sample.trans_first.bam",
                        "-o", os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam",
                        "-q", "-p", "1"]
    if args.mf:
        trans2genome_cmd.append("-m")
    if args.abunds is not None:
        trans2genome_cmd.extend(["-a", args.abunds])
    trans2genome_process = subprocess.Popen(trans2genome_cmd)

    # locus-level alignment
    bowtie2_cmd_locus = None
    if args.locus and args.type == "bowtie":
        print("performing the locus lookup using bowtie")
        bowtie2_cmd_locus = ["bowtie2",
                             "--end-to-end",
                             "--no-unal",
                             "--very-sensitive",
                             "-k", "5",
                             "-x", os.path.abspath(args.db) + "/db.locus",
                             "-p", args.threads]

        if args.fasta:
            bowtie2_cmd_locus.append("-f")
        bowtie2_cmd_locus.extend(("-1", unaligned_r1,
                                  "-2", unaligned_r2,
                                  "-U", unaligned_s,
                                  "--un-conc", os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.fq",
                                  "--un", os.path.abspath(cur_tmp) + "/sample.trans.un_first.locus.fq"))
        if args.bowtie:
            bowtie2_cmd_locus.extend(args.bowtie)

        unaligned_r1 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.1.fq"
        unaligned_r2 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.2.fq"
        unaligned_s = os.path.abspath(cur_tmp) + "/sample.trans.un_first.locus.fq"

        locus_process = subprocess.Popen(bowtie2_cmd_locus, stdout=subprocess.PIPE)
        subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
                          os.path.abspath(cur_tmp) + "/sample.trans_second.locus.bam"], stdin=locus_process.stdout)

    if args.locus:
        locus_process.wait()
        locus_process.stdout.close()

    print("aligning with hisat2 against the genome")
    hisat2_cmd_genome = ["hisat2",
                         "--very-sensitive",
                         "--no-unal",
                         "-k", "30",
                         "-x", os.path.abspath(args.genome_db),
                         "-p", str(int(args.threads) - 1)]
    if args.fasta:
        hisat2_cmd_genome.append("-f")
    hisat2_cmd_genome.extend(("-1", unaligned_r1,
                              "-2", unaligned_r2,
                              "-U", unaligned_s))
    if args.hisat:
        hisat2_cmd_genome.extend(args.hisat)
    genome_process = subprocess.Popen(hisat2_cmd_genome, stdout=subprocess.PIPE)
    subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
                      os.path.abspath(cur_tmp) + "/sample.genome.bam"], stdin=genome_process.stdout)

    genome_process.wait()  # allows trans2genome to run at the same time as hisat2
    genome_process.stdout.close()

    trans2genome_process.wait()

    if not args.keep:
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq"):
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq")
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq"):
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq")
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"):
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq")

    if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans_first.sam") and not args.keep:
        os.remove(os.path.abspath(cur_tmp) + "/sample.trans_first.sam")

    print("merging all sub-alignments")
    merge_cmd = ["samtools", "merge",
                 "-@", args.threads,
                 args.output,
                 os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam"]
    if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans_second.locus.bam") and args.locus:
        merge_cmd.append(os.path.abspath(cur_tmp) + "/sample.trans_second.locus.bam")

    if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans2genome_second.bam"):
        merge_cmd.append(os.path.abspath(cur_tmp) + "/sample.trans2genome_second.bam")

    merge_cmd.append(os.path.abspath(cur_tmp) + "/sample.genome.bam")

    if len(merge_cmd) > 6:
        subprocess.call(merge_cmd)
    else:
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.genome.bam") and not args.keep:
            os.remove(os.path.abspath(cur_tmp) + "/sample.genome.bam")
        os.rename(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam", args.output)
        if not args.keep:
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam")

    if not args.keep:
        shutil.rmtree(os.path.abspath(cur_tmp))
