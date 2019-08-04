import os
import re
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

def parse_hisat(stage_fname, log_fname, log_fh):
    hisat_container = {"total_reads": 0,
                       "total_pair": 0,
                       "conc0": 0,
                       "conc1": 0,
                       "concN": 0,
                       "disc": 0,
                       "conc_disc0": 0,
                       "al0": 0,
                       "al1": 0,
                       "alN": 0,
                       "total_unpair": 0,
                       "unpaired0": 0,
                       "unpaired1": 0,
                       "unpairedN": 0}

    if not os.path.exists(stage_fname):
        return hisat_container

    reg1 = re.compile(' *\d+ [A-z]+')
    reg2 = re.compile(' *\d+ \(\d+\.\d+\%\) [A-z]+')

    num_tab_stack = []

    # Stages - Level 0
    reads = False
    # Stages - Level 1
    paired = False
    unpaired = False
    # Stages - Level 2
    concord = False
    concord0 = False
    concord_disc_0 = False
    # Stages - Level 3
    mates = False

    with open(stage_fname, "r") as inFP:
        for line in inFP.readlines():
            if log_fname is not None:
                log_fh.write(line)

            if re.match(reg1, line) is not None or re.match(reg2, line) is not None:
                line = line.rstrip()
                lineCols = line.split(" ")
                num_tab = lineCols.count("")
                if num_tab == 0:  # total number of reads in the sample
                    hisat_container["total_reads"] = int(lineCols[num_tab])
                    reads = True
                    num_tab_stack.append(num_tab)
                elif num_tab == 2:
                    if "were paired; of these:" in line:
                        hisat_container["total_pair"] = int(lineCols[num_tab])
                        paired = True
                        unpaired = False
                    elif "were unpaired; of these:" in line:
                        hisat_container["total_unpair"] = int(lineCols[num_tab])
                        unpaired = True
                        paired = False
                elif num_tab == 4 and paired and not unpaired:
                    if "aligned concordantly 0 times" in line and not concord:
                        concord = True
                        hisat_container["conc0"] = int(lineCols[num_tab])
                    elif "aligned concordantly exactly 1 time" in line:
                        hisat_container["conc1"] = int(lineCols[num_tab])
                    elif "aligned concordantly >1 times" in line:
                        hisat_container["concN"] = int(lineCols[num_tab])
                    elif "pairs aligned concordantly 0 times; of these:" in line and concord:
                        concord0 = True
                        concord_disc_0 = False
                    elif "pairs aligned 0 times concordantly or discordantly; of these:" in line:
                        hisat_container["conc_disc0"] = int(lineCols[num_tab])
                        concord_disc_0 = True
                        concord0 = False
                    else:
                        print("ERROR parsing hisat2 summary", line)
                elif num_tab == 4 and unpaired and not paired:
                    if "aligned 0 times" in line:
                        hisat_container["unpaired0"] = int(lineCols[num_tab])
                    elif "aligned exactly 1 time" in line:
                        hisat_container["unpaired1"] = int(lineCols[num_tab])
                    elif "aligned >1 times" in line:
                        hisat_container["unpairedN"] = int(lineCols[num_tab])
                elif num_tab == 6 and paired and concord0:
                    if "aligned discordantly 1 time" in line:
                        hisat_container["disc"] = int(lineCols[num_tab])
                elif num_tab == 6 and paired and concord_disc_0:
                    if "mates make up the pairs; of these:" in line:
                        mates = True
                elif num_tab == 8 and paired and concord_disc_0 and mates:
                    if "aligned 0 times" in line:
                        hisat_container["al0"] = int(lineCols[num_tab])
                    elif "aligned exactly 1 time" in line:
                        hisat_container["al1"] = int(lineCols[num_tab])
                    elif "aligned >1 times" in line:
                        hisat_container["alN"] = int(lineCols[num_tab])
                    else:
                        print("ERROR parsing hisat2 summary #5", num_tab, line)
                else:
                    print("ERROR parsing hisat2 summary #6", num_tab, line)

    return hisat_container

def parse_hisat_new(stage_fname, log_fname, log_fh):
    hisat_container = {"total_pair": 0,
                       "conc_disc0": 0,
                       "conc1": 0,
                       "concN": 0,
                       "disc": 0,
                       "total_unpair": 0,
                       "al0": 0,
                       "al1": 0,
                       "alN": 0}

    if not os.path.exists(stage_fname):
        return hisat_container

    with open(stage_fname, "r") as stage_fh:
        for line in stage_fh.readlines():
            if log_fname is not None:
                log_fh.write(line)

            line = line.strip()
            reg = re.compile('Total pairs: \d+')
            if re.match(reg, line) is not None:
                hisat_container["total_pair"] = int(line.split(": ")[1])

            reg = re.compile('Aligned concordantly or discordantly 0 time: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["conc_disc0"] = int(line.split(": ")[1].split(" ")[0])

            reg = re.compile('Aligned concordantly 1 time: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["conc1"] = int(line.split(": ")[1].split(" ")[0])

            reg = re.compile('Aligned concordantly >1 times: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["concN"] = int(line.split(": ")[1].split(" ")[0])

            reg = re.compile('Aligned discordantly 1 time: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["disc"] = int(line.split(": ")[1].split(" ")[0])

            reg = re.compile('Total unpaired reads: \d+')
            if re.match(reg, line) is not None:
                hisat_container["total_unpair"] = int(line.split(": ")[1])

            reg = re.compile('Aligned 0 time: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["al0"] = int(line.split(": ")[1].split(" ")[0])

            reg = re.compile('Aligned 1 time: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["al1"] = int(line.split(": ")[1].split(" ")[0])

            reg = re.compile('Aligned >1 times: \d+ .*')
            if re.match(reg, line) is not None:
                hisat_container["alN"] = int(line.split(": ")[1].split(" ")[0])

    return hisat_container


def parse_t2g(stage_fname, log_fname, log_fh):
    assert os.path.exists(stage_fname), "stage2: translation log is not found"
    t2g_container = dict()
    with open(stage_fname, "r") as stage_fh:
        for line in stage_fh.readlines():
            if log_fname is not None:
                log_fh.write(line)

            if "single reads were precomputed" in line and line[:6] == "@STATS":
                t2g_container["single_precomp"] = int(line.split(":")[1].split(" ")[1])
            elif "discordant paired reads were precomputed" in line and line[:6] == "@STATS":
                t2g_container["disc_precomp"] = int(line.split(":")[1].split(" ")[1])
            elif "concordant paired reads were precomputed" in line and line[:6] == "@STATS":
                t2g_container["paired_precomp"] = int(line.split(":")[1].split(" ")[1])
            elif "minimum fragment length threshold" in line and line[:6] == "@STATS":
                t2g_container["min_fraglen_thresh"] = int(line.split(":")[1].split(" ")[1])
            elif "maximum fragment length threshold" in line and line[:6] == "@STATS":
                t2g_container["max_fraglen_thresh"] = int(line.split(":")[1].split(" ")[1])
            elif "singleton reads were evaluated" in line and line[:6] == "@STATS":
                t2g_container["single_al"] = int(line.split(":")[1].split(" ")[1])
            elif "discordant reads were evaluated" in line and line[:6] == "@STATS":
                t2g_container["disc_al"] = int(line.split(":")[1].split(" ")[1])
            elif "paired reads were evaluated" in line and line[:6] == "@STATS":
                t2g_container["pair_al"] = int(line.split(":")[1].split(" ")[1])
            elif "singletons discarded due to high edit distance" in line and line[:6] == "@STATS":
                t2g_container["single_un_err"] = int(line.split(":")[1].split(" ")[1])
            elif "discordant pairs discarded due to high edit distance" in line and line[:6] == "@STATS":
                t2g_container["disc_un_err"] = int(line.split(":")[1].split(" ")[1])
            elif "concordant pairs discarded due to high edit distance" in line and line[:6] == "@STATS":
                t2g_container["paired_un_err"] = int(line.split(":")[1].split(" ")[1])
            elif "paired reads were unaligned" in line and line[:6] == "@STATS":
                t2g_container["paired_un"] = int(line.split(":")[1].split(" ")[1])
            elif "singletons were unaligned" in line and line[:6] == "@STATS":
                t2g_container["single_un"] = int(line.split(":")[1].split(" ")[1])
            elif "singleton edit distance threshold used" in line and line[:6] == "@STATS":
                try:
                    t2g_container["single_edit_thresh"] = int(line.split(":")[1].split(" ")[1])
                except ValueError:
                    t2g_container["single_edit_thresh"] = None
            elif "concordant paired edit distance threshold used" in line and line[:6] == "@STATS":
                try:
                    t2g_container["pair_edit_thresh"] = int(line.split(":")[1].split(" ")[1])
                except ValueError:
                    t2g_container["pair_edit_thresh"] = None
            elif "multimapping singleton reads detected" in line and line[:6] == "@STATS":
                t2g_container["single_multi"] = int(line.split(":")[1].split(" ")[1])
            elif "multimapping concordantly paired reads detected" in line and line[:6] == "@STATS":
                t2g_container["paired_multi"] = int(line.split(":")[1].split(" ")[1])
            elif "mean singleton multimapping rate" in line and line[:6] == "@STATS":
                try:
                    t2g_container["single_rate"] = float(line.split(":")[1].split(" ")[1])
                except ValueError:
                    t2g_container["single_rate"] = None
            elif "concordant pair multimapping rate" in line and line[:6] == "@STATS":
                try:
                    t2g_container["paired_rate"] = float(line.split(":")[1].split(" ")[1])
                except ValueError:
                    t2g_container["paired_rate"] = None
            else:
                continue

    return t2g_container


def parse_logs(tmpDir_tmp, log_fname=None):
    log_fh = None
    if log_fname is not None:
        log_fh = open(log_fname, "w+")

    # begin by analyzing the alignment stats from transcriptomic alignment
    tmpDir = os.path.abspath(tmpDir_tmp)

    stage1_fname = tmpDir + "/stage1_transcriptome.tmp"
    assert os.path.exists(stage1_fname), "stage1: transcriptome alignment log is not found"
    stage1_res = parse_hisat(stage1_fname, log_fname, log_fh)

    stage2_fname = tmpDir + "/stage2_translate.tmp"
    assert os.path.exists(stage2_fname), "stage2: alignment translation log is not found"
    stage2_res = parse_t2g(stage2_fname, log_fname, log_fh)

    # STAGE3 - LOCUS ALIGNMENT - IF AVAILABLE
    stage3_fname = tmpDir + "/stage3_locus.tmp"
    stage3_res = dict()
    stage3_res = parse_hisat(stage3_fname, log_fname, log_fh)

    # STAGE4 - GENOME ALIGNMENT
    stage4_fname = tmpDir + "/stage4_genome.tmp"
    assert os.path.exists(stage4_fname), "stage4: genome alignment log is not found"
    stage4_res = parse_hisat(stage4_fname, log_fname, log_fh)

    # STAGE5 - MERGING OF THE ALIGNMENT FILES
    if log_fname is not None:
        stage5_fname = tmpDir + "/stage5_merge.tmp"
        assert os.path.exists(stage5_fname), "stage5: merge alignment log is not found"
        with open(stage5_fname, "r") as stage5_fh:
            for line in stage5_fh.readlines():
                log_fh.write(line)

    # COMPUTE FINAL STATS

    report = "======================\n===== T2G REPORT =====\n======================\n\n"
    total_paired = stage1_res["total_pair"]
    report += "Total pairs: %d\n" % total_paired

    trans_conc_N = stage2_res["paired_multi"]
    trans_conc_1 = stage2_res["pair_al"] - trans_conc_N
    total_conc_1 = trans_conc_1 + stage4_res["conc1"]
    total_conc_N = trans_conc_N + stage4_res["concN"]
    total_conc = total_conc_1 + total_conc_N

    total_conc_0 = total_paired - total_conc
    perc_conc_0 = (total_conc_0 / total_paired) * 100
    report += "\tAligned concordantly 0 time: %d (%.2f%%)\n" % (total_conc_0, perc_conc_0)

    trans_disc = stage2_res["disc_al"]
    genome_disc = stage4_res["disc"]
    total_disc = trans_disc + genome_disc
    perc_disc = (total_disc / total_conc_0) * 100
    report += "\t\tAligned discordantly: %d (%.2f%%)\n" % (total_disc, perc_disc)

    perc_disc_trans = (trans_disc / total_disc) * 100
    report += "\t\t\tAligned to the transcriptome: %d (%.2f%%)\n" % (trans_disc, perc_disc_trans)

    perc_disc_genome = (genome_disc / total_disc) * 100
    report += "\t\t\tAligned to the genome: %d (%.2f%%)\n" % (genome_disc, perc_disc_genome)

    total_conc_disc_0 = total_conc_0 - total_disc
    perc_conc_disc_0 = (total_conc_disc_0 / total_conc_0) * 100
    report += "\t\tAligned concordantly or discordantly 0 times: %d (%.2f%%)\n" % (total_conc_disc_0, perc_conc_disc_0)

    total_conc_disc_0_mates = total_conc_disc_0 * 2
    report += "\t\t\tNumber of mates that make up these pairs: %d\n" % total_conc_disc_0_mates

    trans_single = stage2_res["single_al"]
    genome_single = stage4_res["al1"] + stage4_res["alN"] + stage4_res["unpaired1"] + stage4_res["unpairedN"]
    total_single_genome = trans_single + genome_single
    perc_single = (total_single_genome / total_conc_disc_0_mates) * 100
    report += "\t\t\t\tAligned: %d (%.2f%%)\n" % (total_single_genome, perc_single)

    perc_trans_single = (trans_single / total_single_genome) * 100
    report += "\t\t\t\t\tAligned to the transcriptome: %d (%.2f%%)\n" % (trans_single, perc_trans_single)

    perc_genome_single = (genome_single / total_single_genome) * 100
    report += "\t\t\t\t\tAligned to the genome: %d (%.2f%%)\n" % (genome_single, perc_genome_single)

    total_unal = stage4_res["al0"] + stage4_res["unpaired0"]
    perc_unal = (total_unal / total_conc_disc_0_mates) * 100
    report += "\t\t\t\tAligned 0 times: %d (%.2f%%)\n" % (total_unal, perc_unal)

    perc_conc_1 = (total_conc_1 / total_paired) * 100
    report += "\tAligned concordantly 1 time: %d (%.2f%%)\n" % (total_conc_1, perc_conc_1)

    perc_conc_1_trans = (trans_conc_1 / total_conc_1) * 100
    report += "\t\tAligned to the transcriptome: %d (%.2f%%)\n" % (trans_conc_1, perc_conc_1_trans)

    perc_conc_1_genome = (stage4_res["conc1"] / total_conc_1) * 100
    report += "\t\tAligned to the genome: %d (%.2f%%)\n" % (stage4_res["conc1"], perc_conc_1_genome)

    perc_conc_N = (total_conc_N / total_paired) * 100
    report += "\tAligned concordantly >1 times: %d (%.2f%%)\n" % (total_conc_N, perc_conc_N)

    perc_conc_N_trans = (trans_conc_N / total_conc_N) * 100
    report += "\t\tAligned to the transcriptome: %d (%.2f%%)\n" % (trans_conc_N, perc_conc_N_trans)

    perc_conc_N_genome = (stage4_res["concN"] / total_conc_N) * 100
    report += "\t\tAligned to the genome: %d (%.2f%%)\n" % (stage4_res["concN"], perc_conc_N_genome)

    # lastly need to get the final alignment rate
    total_reads = total_paired*2 + stage1_res["total_unpair"]
    total_aligned = total_conc*2 + total_disc + total_single_genome
    al_rate = (total_aligned / total_reads) * 100
    report += "%.2f%% overall alignment rate\n" % al_rate

    print(report, file=sys.stderr)
    log_fh.write(report)

    if log_fname is not None:
        log_fh.close()

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

    if args.output.split(".")[-1] in ["bam", "sam", "cram"]:
        final_fname = ".".join(args.output.split(".")[:-1]) + ".logs"
    else:
        final_fname = args.output + ".stats"

    transcriptome_cmd = None
    if args.type == "hisat":
        print("aligning with hisat2 against transcriptome")
        # perform a two-pass alignment one for less redundant transcripts and one for more redundant alignments
        # with the "-a" option enabled
        transcriptome_cmd = ['hisat2',
                             # get path to the trans2genome that was compiled with the package,
                             "--no-spliced-alignment",
                             "--end-to-end",
                             "--rna-sensitive",
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
    # transcriptome_cmd.extend(("--un-conc", os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.fq",
    #                           "--un", os.path.abspath(cur_tmp) + "/sample.trans.un_first.fq"))

    start_transcriptome = time.time()
    stage1_transcriptome_fh = open(os.path.abspath(cur_tmp) + "/stage1_transcriptome.tmp", "w+")
    transcriptome_process = subprocess.Popen(transcriptome_cmd, stdout=subprocess.PIPE, stderr=stage1_transcriptome_fh)
    trans2genome_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     'trans2genome')  # get path to the trans2genome that was compiled with the package
    trans2genome_cmd = [trans2genome_path,
                        "-x", os.path.abspath(args.db) + "/db",
                        "-o", os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam",
                        "-u", os.path.abspath(cur_tmp) + "/sample.trans_first",
                        "-q", "-p", str(args.threads)]
    if args.mf:
        trans2genome_cmd.append("-m")
    if args.abunds is not None:
        trans2genome_cmd.extend(["-a", args.abunds])
    if args.all is not None:
        trans2genome_cmd.extend(["-l"])
    if args.errcheck:
        trans2genome_cmd.extend(["-s"])
    if args.no_discord:
        trans2genome_cmd.extend(["-j"])
    if args.no_single:
        trans2genome_cmd.extend(["-g"])

    stage2_translate_fh = open(os.path.abspath(cur_tmp) + "/stage2_translate.tmp", "w+")
    translate_process = subprocess.Popen(trans2genome_cmd, stdin=transcriptome_process.stdout,
                                         stderr=stage2_translate_fh)

    # make sure the sam header has the exact information to run the tool

    unaligned_r1 = os.path.abspath(cur_tmp) + "/sample.trans_first.unal_r1.fastq"
    unaligned_r2 = os.path.abspath(cur_tmp) + "/sample.trans_first.unal_r2.fastq"
    unaligned_s = os.path.abspath(cur_tmp) + "/sample.trans_first.unal_s.fastq"

    transcriptome_process.wait()
    transcriptome_process.stdout.close()
    # transcriptome_process.stderr.close()

    stage1_transcriptome_fh.close()
    if args.verbose:
        with open(os.path.abspath(cur_tmp) + "/stage1_transcriptome.tmp", "r") as inFP:
            for line in inFP.readlines():
                print(line.rstrip("\n"))

    # translate_process.stderr.close()
    translate_process.wait()

    stage2_translate_fh.close()
    if args.verbose:
        with open(os.path.abspath(cur_tmp) + "/stage2_translate.tmp", "r") as inFP:
            for line in inFP.readlines():
                print(line.rstrip("\n"))

    end_transcriptome = time.time()
    print("Transcriptome alignment time: "+str(datetime.timedelta(seconds=int(end_transcriptome-start_transcriptome))))

    # locus-level alignment
    locus_cmd = None
    if args.locus:
        stage3_locus_fh = open(os.path.abspath(cur_tmp) + "/stage3_locus.tmp", "w+")
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
            locus_cmd = [os.path.abspath("hisat2"),
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

        unaligned_r1 = unaligned_r1+","+os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.1.fq"
        unaligned_r2 = unaligned_r2+","+os.path.abspath(cur_tmp) + "/sample.trans.unconc_first.locus.2.fq"
        unaligned_s = unaligned_s+","+os.path.abspath(cur_tmp) + "/sample.trans.un_first.locus.fq"

        start_locus = time.time()
        locus_process = subprocess.Popen(locus_cmd, stderr=stage3_locus_fh)
        # subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
        #                   os.path.abspath(cur_tmp) + "/sample.trans_second.locus.bam"], stdin=locus_process.stdout)

        locus_process.wait()
        # locus_process.stderr.close()

        stage3_locus_fh.close()
        if args.verbose:
            with open(os.path.abspath(cur_tmp) + "/stage3_locus.tmp", "r") as inFP:
                for line in inFP.readlines():
                    print(line.rstrip("\n"))

        end_locus = time.time()
        print("Locus alignment time: "+str(datetime.timedelta(seconds=int(end_locus-start_locus))))

    print("aligning with hisat2 against the genome")
    hisat2_cmd_genome = ["hisat2",
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
    stage4_genome_fh = open(os.path.abspath(cur_tmp) + "/stage4_genome.tmp", "w+")
    genome_process = subprocess.Popen(hisat2_cmd_genome, stderr=stage4_genome_fh)
    # convert_process = subprocess.Popen(["samtools", "view", "-h", "--output-fmt=BAM", "-@", args.threads, "-o",
    #                                     os.path.abspath(cur_tmp) + "/sample.genome.bam"], stdin=genome_process.stdout)

    genome_process.wait()  # allows trans2genome to run at the same time as hisat2
    # genome_process.stderr.close()

    stage4_genome_fh.close()
    if args.verbose:
        with open(os.path.abspath(cur_tmp) + "/stage4_genome.tmp", "r") as inFP:
            for line in inFP.readlines():
                print(line.rstrip("\n"))

    # convert_process.wait()

    end_genome = time.time()
    print("Genome alignment time: "+str(datetime.timedelta(seconds=int(end_genome-start_genome))))

    if not args.keep:
        for item in unaligned_r1.split(","):
            if os.path.exists(item):
                os.remove(item)
        for item in unaligned_r2.split(","):
            if os.path.exists(item):
                os.remove(item)
        for item in unaligned_s.split(","):
            if os.path.exists(item):
                os.remove(item)

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

    stage5_merge_fh = open(os.path.abspath(cur_tmp) + "/stage5_merge.tmp", "w+")
    if len(merge_cmd) > 6:
        start_merge = time.time()
        subprocess.call(merge_cmd, stderr=stage5_merge_fh)
        stage5_merge_fh.close()
        if args.verbose:
            with open(os.path.abspath(cur_tmp) + "/stage5_merge.tmp", "r") as inFP:
                for line in inFP.readlines():
                    print(line.rstrip("\n"))
        end_merge = time.time()
        print("Merge time: "+str(datetime.timedelta(seconds=int(end_merge-start_merge))))
    else:
        if os.path.exists(os.path.abspath(cur_tmp) + "/sample.genome.sam") and not args.keep:
            os.remove(os.path.abspath(cur_tmp) + "/sample.genome.sam")
        os.rename(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam", args.output)
        if not args.keep:
            os.remove(os.path.abspath(cur_tmp) + "/sample.trans2genome_first.bam")

    if args.log:
        parse_logs(cur_tmp, final_fname)
    else:
        parse_logs(cur_tmp)

    if not args.keep:
        shutil.rmtree(os.path.abspath(cur_tmp))

    end_total = time.time()
    print("Total time elapsed: "+str(datetime.timedelta(seconds=int(end_total-start_total))))

    # TODO: need to remove all output from different internal tools and report only general alignment (and realignment report)
    # TODO: create a unified alignment rate report for the total number of alignments across both transcriptome and genome searches
    # TODO: automatically detect the directory where all tools are installed and verify that everything is installed
