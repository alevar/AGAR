#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys
import os

def extractLocus(args):
    gff3Cols=["seqid","source","type","start","end","score","strand","phase","attributes"]

    chrDict = dict()
    with open(args.faidx) as inFP:
        for line in inFP.readlines():
            lineCols = line.split("\t")
            chrDict.setdefault(lineCols[0],0)
            chrDict[lineCols[0]]=int(lineCols[1])
            
    df=pd.read_csv(args.input,sep="\t",names=gff3Cols,comment="#")
    df["parent"]=df.attributes.str.split("gene_id \"",expand=True)[1].str.split("\";",expand=True)[0]
    df=df[~df["parent"].isnull()].reset_index(drop=True)
    df.drop_duplicates(["seqid","start","end","strand"],inplace=True)
    df["id"]=df.attributes.str.split("transcript_id \"",expand=True)[1].str.split("\";",expand=True)[0]
    df.drop_duplicates("id",inplace=True)
    df.reset_index(drop=True,inplace=True)
    df.drop(["id","parent"],axis=1)

    df["start"]=df.start.astype(int)
    df["end"]=df.end.astype(int)
    # add start and end for each strand
    for chrID in set(df["seqid"]):
        df=df.append({"seqid":chrID,"source":"new","type":"locus","start":0,"end":1,"score":".","strand":"+","phase":".","attributes":"none=0"}, ignore_index=True)
        df=df.append({"seqid":chrID,"source":"new","type":"locus","start":0,"end":1,"score":".","strand":"-","phase":".","attributes":"none=0"}, ignore_index=True)
        df=df.append({"seqid":chrID,"source":"new","type":"locus","start":chrDict[chrID]-1,"end":chrDict[chrID]-2,"score":".","strand":"+","phase":".","attributes":"none=0"}, ignore_index=True)
        df=df.append({"seqid":chrID,"source":"new","type":"locus","start":chrDict[chrID]-1,"end":chrDict[chrID]-2,"score":".","strand":"-","phase":".","attributes":"none=0"}, ignore_index=True)

    df.drop(["parent","id"],axis=1,inplace=True)

    df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)

    df["next_start"]=df["start"].shift(-1)
    df["next_seqid"]=df["seqid"].shift(-1)
    df["next_strand"]=df["strand"].shift(-1)
    df["next_end"] = df["end"].shift(-1)
    df=df.dropna(axis=0)
    df["next_start"]=df.next_start.astype(int)
    df["next_end"]=df.next_end.astype(int)
    
    before = 1
    after  = 0

    while(before>after):
        df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)
        
        df["next_start"]=df["start"].shift(-1)
        df["next_seqid"]=df["seqid"].shift(-1)
        df["next_strand"]=df["strand"].shift(-1)
        df["next_end"] = df["end"].shift(-1)
        df=df.dropna(axis=0)
        df["next_start"]=df.next_start.astype(int)
        df["next_end"]=df.next_end.astype(int)
        
        df.reset_index(drop=True,inplace=True)
        df["prev_start"]=df["start"].shift(1)
        df["prev_end"]=df["end"].shift(1)
        df["prev_seqid"]=df["seqid"].shift(1)
        df["prev_strand"]=df["strand"].shift(1)

        # now need to extend them as far as possible
        df["end"]=np.where((df["seqid"]==df["next_seqid"]) & \
                            (df["strand"]==df["next_strand"]) & \
                            (df["end"]>df["next_start"]),df["next_end"],df["end"])
        df["prev_start"]=np.where(df["prev_start"].isnull(),df["start"],df["prev_start"])
        df["prev_end"]=np.where(df["prev_end"].isnull(),df["end"],df["prev_end"])
        df["prev_seqid"]=np.where(df["prev_seqid"].isnull(),df["seqid"],df["prev_seqid"])
        df["prev_strand"]=np.where(df["prev_strand"].isnull(),df["strand"],df["prev_strand"])
        df["prev_start"]=df["prev_start"].astype(int)
        df["prev_end"]=df["prev_end"].astype(int)
        before=len(df)
        df_first = df.iloc[0]
        df=df[~((df["seqid"]==df["prev_seqid"]) & \
                  (df["strand"]==df["prev_strand"]) & \
                  (df["start"]>=df["prev_start"]) & \
                  (df["end"]<=df["prev_end"]))]
        df=df.append(df_first)
        after=len(df)

    # finally create new introns and write everything to the new locus gff for masking
    df.sort_values(by=["seqid","strand","start","end"],ascending=True,inplace=True)
    df.reset_index(drop=True,inplace=True)
    df["next_start"]=df["start"].shift(-1)
    df["next_seqid"]=df["seqid"].shift(-1)
    df["next_strand"]=df["strand"].shift(-1)
    df["start"]=df["end"]
    df["end"]=np.where((df["seqid"]==df["next_seqid"]) & \
                            (df["strand"]==df["next_strand"]),df["next_start"],df["end"])

    last_chr = df.iloc[len(df)-1].seqid
    df.set_value(len(df)-1, 'end', chrDict[last_chr]-1)

    df["start"]=df.start.astype(int)
    df["end"]=df.end.astype(int)
    df=df[~(df["start"]==df["end"])].reset_index(drop=True)
    
    df[gff3Cols].to_csv(args.out,sep='\t',index=False,header=False)


def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument("--input",
                        type=str,
                        required=True,
                        help="GFF or GTF annotation of the genome")
    parser.add_argument("--faidx",
                        type=str,
                        required=True,
                        help="fasta index of the reference genome")
    parser.add_argument("--out",
                        type=str,
                        required=True,
                        help="output file")

    parser.set_defaults(func=extractLocus)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])