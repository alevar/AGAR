#!/usr/bin/env python

# compute equivalence classes based on a gff annotation
import argparse
import sys
import os


def evaluate_interval(iv, oivs):  # takes current interval and a sorted list of original intervals
    res = []
    for o in oivs:
        if iv[1] > o[0] and o[1] > iv[0]:
            res.append(o[2])
    return res

# TODO: need to add coordinates of each equivalence class
# TODO: ideally should also include splice jucntions, but good enough for now
def process_locus(eps, oivs):
    interval = []
    res = []
    for i, e in enumerate(eps):  # iterate over and collect intervals
        if len(interval) == 2:  # full interval found, can be evaluated and reset
            segs = evaluate_interval(interval, oivs)
            if not len(segs) == 0:
                res.append(interval + segs)
                tmp = interval[1]
                if not e[0] == 0:
                    interval = [tmp, i]
                else:
                    interval = [tmp]
            else:
                tmp = interval[1]
                if not e[0] == 0:
                    interval = [tmp, i]
                else:
                    interval = [tmp]
        elif not e[0] == 0:
            interval.append(i)
        else:
            continue
    if len(interval) == 2:
        segs = evaluate_interval(interval, oivs)
        if not len(segs) == 0:
            res.append(interval + segs)

    return res

# TODO: need to write to a file
def report(segs):
    if len(segs) > 0:
        print(segs)


def gff2ec(args):
    endpoints = []

    original_intervals = set()  # these contain the transcript boundaries (without separation into exons and introns).
    # This is necessary in order to preserve splicejunction information whenver possible
    chrom = ""
    strand = "+"
    ls = 0  # locus start
    le = 0  # locus end
    geneID = ""
    ts = 0
    te = 0
    transID = ""
    transPar = ""

    assert os.path.exists(os.path.abspath(args.gff)), "unable to open the GFF file"

    with open(os.path.abspath(args.gff), "r") as gff:
        for line in gff.readlines():
            if line[0] == "#":
                continue
            else:
                lineCols = line.split("\t")
                if lineCols[2] == "gene":
                    final = process_locus(endpoints,
                                          original_intervals)  # process the intervals computed for the previous locus
                    final = [[x[0] + ls, x[1] + ls, tuple(x[2:])] for x in final]
                    report(final)
                    # get some basic information about the locus and initialize things
                    chrom = lineCols[0]
                    strand = lineCols[6]
                    ls = int(lineCols[3])
                    le = int(lineCols[4])
                    geneID = lineCols[8].strip().split("ID=")[1].split(";")[0]
                    print(geneID)

                    endpoints = [[0, set()] for x in range(ls, le + 1, 1)]
                    original_intervals = set()

                elif lineCols[2] == "transcript":
                    transID = lineCols[8].strip().split("ID=")[1].split(";")[0]
                    transPar = lineCols[8].strip().split("Parent=")[1].split(";")[0]
                    ts = int(lineCols[3])
                    te = int(lineCols[4])

                    # add to original intervals
                #                 original_intervals.add((ts-1,te-1,))

                elif lineCols[2] == "exon":
                    parent = lineCols[8].strip().split("Parent=")[1].split(";")[0]
                    start = int(lineCols[3])
                    end = int(lineCols[4])
                    endpoints[start - ls][0] += 1
                    endpoints[start - ls][1].add(parent)
                    endpoints[end - ls][0] += 1
                    endpoints[end - ls][1].add(parent)

                    # add to original intervals
                    original_intervals.add((start - ls, end - ls, parent))
                else:
                    continue

        final = process_locus(endpoints, original_intervals)
        final = [[x[0] + ls, x[1] + ls, tuple(x[2:])] for x in final]
        report(final)


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument("--gff",
                        type=str,
                        required=True,
                        help="GFF or GTF annotation of the genome")
    parser.add_argument("-o",
                        "--output",
                        type=str,
                        required=True,
                        help="output file base name for the equivalence class definitions")
    parser.set_defaults(func=gff2ec)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
