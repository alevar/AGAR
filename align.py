import os
import sys
import time
import shutil
import datetime
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


def parse_logs():
    res = "======================\n\n===== T2G REPORT =====\n======================\n\n"

    return res


def main(args):
    start_total = time.time()
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

    # FIRST INITIALIZE FILE HANDLERS TO STAGE OUTPUTS
    stage1_transcriptome_fh = open(os.path.abspath(cur_tmp) + "/stage1_transcriptome.tmp", "w+")
    stage2_translate_fh = open(os.path.abspath(cur_tmp) + "/stage2_translate.tmp", "w+")
    stage3_locus_fh = open(os.path.abspath(cur_tmp) + "/stage3_locus.tmp", "w+")
    stage4_genome_fh = open(os.path.abspath(cur_tmp) + "/stage4_locus.tmp", "w+")
    stage5_merge_fh = open(os.path.abspath(cur_tmp) + "/stage5_locus.tmp", "w+")
    final_fh = None
    if args.output.split(".")[-1] in ["bam", "sam", "cram"]:
        final_fname = ".".join(args.output.split(".")[:-1]) + ".stats"
        final_fh = open(final_fname, "w+")
    else:
        final_fh = open(args.output + ".stats", "w+")

    transcriptome_cmd = None
    if args.type == "hisat":
        print("aligning with hisat2 against transcriptome")
        # perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments
        # with the "-a" option enabled
        transcriptome_cmd = ['/home/varabyou/soft/hisat2/hisat2',
                             # get path to the trans2genome that was compiled with the package,
                             "--no-spliced-alignment",
                             "--end-to-end",
                             "--no-unal",
                             "-x", os.path.abspath(args.db) + "/db",
                             "-p", args.threads]

        if args.hisat:
            transcriptome_cmd.extend(args.hisat)

    elif args.type == "bowtie":
        print("aligning with bowtie2 against transcriptome")
        # perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments
        # with the "-a" option enabled
        transcriptome_cmd = {"bowtie2",
                             "--end-to-end",
                             "--no-unal",
                             "-x", os.path.abspath(args.db) + "/db",
                             "-p", args.threads}
        if args.bowtie:
            transcriptome_cmd.extend(args.bowtie)

    if args.k:
        transcriptome_cmd.extend(("-k", str(args.k)))
    else:
        transcriptome_cmd.extend(("-k", "1"))
    if args.fasta:
        transcriptome_cmd.append("-f")
    transcriptome_cmd.extend(("-1", args.m1,
                              "-2", args.m2))
    if not args.single is None:
        transcriptome_cmd.extend(("-U", args.single))
    transcriptome_cmd.extend(("--un-conc", os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.fq",
                              "--un", os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"))

    start_transcriptome = time.time()
    transcriptome_process = subprocess.Popen(transcriptome_cmd, stdout=subprocess.PIPE, stderr=stage1_transcriptome_fh)
    unaligned_r1 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq"
    unaligned_r2 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq"
    unaligned_s = os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"
    trans2genome_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     'cmake-build-release/trans2genome')  # get path to the trans2genome that was compiled with the package
    trans2genome_cmd = [trans2genome_path,
                        "-x", os.path.abspath(args.db) + "/db",
                        "-o", os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam",
                        "-q", "-p", str(args.threads)]
    if args.mf:
        trans2genome_cmd.append("-m")
    if args.abunds is not None:
        trans2genome_cmd.extend(["-a", args.abunds])
    if args.all is not None:
        trans2genome_cmd.extend(["-l"])
    if args.errcheck:
        trans2genome_cmd.extend(["-s"])

    translate_process = subprocess.Popen(trans2genome_cmd, stdin=transcriptome_process.stdout,
                                         stderr=stage2_translate_fh)

    transcriptome_process.wait()
    transcriptome_process.stdout.close()
    # transcriptome_process.stderr.close()

    stage1_transcriptome_fh.close()

    # translate_process.stderr.close()
    translate_process.wait()

    stage2_translate_fh.close()

    end_transcriptome = time.time()
    print("Transcriptome alignment time: "+str(datetime.timedelta(seconds=int(end_transcriptome-start_transcriptome))))

    if args.errcheck:
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam.unal_r1.fastq"):
            unaligned_r1 = unaligned_r1 + "," + os.path.abspath(
                cur_tmp) + "/sample.trans2genome_first.bam.unal_r1.fastq"
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam.unal_r1.fastq"):
            unaligned_r2 = unaligned_r2 + "," + os.path.abspath(
                cur_tmp) + "/sample.trans2genome_first.bam.unal_r2.fastq"
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam.unal_r1.fastq"):
            unaligned_s = unaligned_s + "," + os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam.unal_s.fastq"

        # locus-level alignment
    locus_cmd = None
    if args.locus:
        if args.type == "bowtie":
            print("performing the locus lookup using bowtie2")
            locus_cmd = ["bowtie2",
                         "--end-to-end",
                         "--no-unal",
                         "-k", "5",
                         "-x", os.path.abspath(args.db) + "/db.locus",
                         "-p", args.threads]

            if args.bowtie:
                locus_cmd.extend(args.bowtie)
        elif args.type == "hisat":
            print("performing the locus lookup using hisat2")
            locus_cmd = [os.path.abspath("/home/varabyou/soft/hisat2/hisat2"),
                         "--rna-sensitive",
                         "--end-to-end",
                         "--no-unal",
                         "-x", os.path.abspath(args.db) + "/db.locus",
                         "-p", args.threads]

            if args.hisat:
                locus_cmd.extend(args.hisat)

        if args.fasta:
            locus_cmd.append("-f")
        locus_cmd.extend(("-S", os.path.abspath(cur_tmp) + "/sample.trans_second.locus.sam",
                          "-1", unaligned_r1,
                          "-2", unaligned_r2,
                          "-U", unaligned_s,
                          "--un-conc", os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.fq",
                          "--un", os.path.abspath(cur_tmp) + "/sample.trans.un_first.locus.fq"))

        unaligned_r1 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.1.fq"
        unaligned_r2 = os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.2.fq"
        unaligned_s = os.path.abspath(cur_tmp) + "/sample.trans.un_first.locus.fq"

        start_locus = time.time()
        locus_process = subprocess.Popen(locus_cmd, stderr=stage3_locus_fh)
        # subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
        #                   os.path.abspath(cur_tmp) + "/sample.trans_second.locus.bam"], stdin=locus_process.stdout)

        locus_process.wait()
        # locus_process.stderr.close()

        stage3_locus_fh.close()

        end_locus = time.time()
        print("Locus alignment time: "+str(datetime.timedelta(seconds=int(end_locus-start_locus))))

    print("aligning with hisat2 against the genome")
    hisat2_cmd_genome = ["/home/varabyou/soft/hisat2/hisat2",
                         "--rna-sensitive",
                         "-x", os.path.abspath(args.genome_db),
                         "-p", args.threads]
    if args.nounal:
        hisat2_cmd_genome.append("--no-unal")
    if args.fasta:
        hisat2_cmd_genome.append("-f")
    hisat2_cmd_genome.extend(("-S", os.path.abspath(cur_tmp) + "/sample.genome.sam",
                              "-1", unaligned_r1,
                              "-2", unaligned_r2,
                              "-U", unaligned_s))
    if args.hisat:
        hisat2_cmd_genome.extend(args.hisat)

    start_genome = time.time()
    genome_process = subprocess.Popen(hisat2_cmd_genome, stderr=stage4_genome_fh)
    # convert_process = subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
    #                                     os.path.abspath(cur_tmp) + "/sample.genome.bam"], stdin=genome_process.stdout)

    genome_process.wait()  # allows trans2genome to run at the same time as hisat2
    # genome_process.stderr.close()

    stage4_genome_fh.close()

    # convert_process.wait()

    end_genome = time.time()
    print("Genome alignment time: "+str(datetime.timedelta(seconds=int(end_genome-start_genome))))

    if not args.keep:
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq"):
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.1.fq")
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq"):
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.2.fq")
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"):
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq")

    print("merging all sub-alignments")
    merge_cmd = ["samtools", "merge", "-f",
                 "-@", args.threads,
                 args.output,
                 os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam"]
    if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans_second.locus.sam") and args.locus:
        merge_cmd.append(os.path.abspath(cur_tmp) + "/sample.trans_second.locus.sam")

    if os.path.exists(os.path.abspath(cur_tmp) + "/sample.trans2genome_second.sam"):
        merge_cmd.append(os.path.abspath(cur_tmp) + "/sample.trans2genome_second.sam")

    merge_cmd.append(os.path.abspath(cur_tmp) + "/sample.genome.sam")

    if len(merge_cmd) > 6:
        start_merge = time.time()
        subprocess.call(merge_cmd, stderr=stage5_merge_fh)
        stage5_merge_fh.close()
        end_merge = time.time()
        print("Merge time: "+str(datetime.timedelta(seconds=int(end_merge-start_merge))))
    else:
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.genome.sam") and not args.keep:
            os.remove(os.path.abspath(cur_tmp) + "/sample.genome.sam")
        os.rename(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam", args.output)
        if not args.keep:
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam")

    # lastly process the outputs of the stages and create final report
    subprocess.call(["cat",
                     os.path.abspath(cur_tmp) + "/stage1_transcriptome.tmp",
                     os.path.abspath(cur_tmp) + "/stage2_translate.tmp",
                     os.path.abspath(cur_tmp) + "/stage3_locus.tmp",
                     os.path.abspath(cur_tmp) + "/stage4_locus.tmp",
                     os.path.abspath(cur_tmp) + "/stage5_locus.tmp"], stdout=final_fh)
    final_fh.close()

    print(parse_logs(), file=sys.stderr)

    if not args.keep:
        shutil.rmtree(os.path.abspath(cur_tmp))

    end_total = time.time()
    print("Total time elapsed: "+str(datetime.timedelta(seconds=int(end_total-start_total))))

    # TODO: need to remove all output from different internal tools and report only general alignment (and realignment report)
    # TODO: create a unified alignment rate report for the total number of alignments across both transcriptome and genome searches
    # TODO: automatically detect the directory where all tools are installed and verify that everything is installed
