#!/usr/bin/python3

import argparse
import yaml
import csv
import re
#import numpy

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("yaml",help="orig yaml",type=str)
    parser.add_argument("ranks",help="grinder ranks",type=str)
    parser.add_argument("bartender_barcode",help="bartender barcode file",type=str)
    parser.add_argument("bartender_cluster",help="bartender cluster file",type=str)
    parser.add_argument("actual_barcode_pattern",help="the real barcode pattern because bartender is written in a way that obscures honest benchmarking",type=str)
    args = parser.parse_args()

    orig_seq_to_id = dict()
    with open(args.yaml,"r") as f:
        for i in list(yaml.load_all(f))[0]:
            for key in i.keys():
                orig_seq_to_id[key] = i[key]

    grinder_ranks = dict()
    with open(args.ranks,"r") as f:
        library_ranks = csv.reader(f,delimiter="\t")
        next(library_ranks)
        for i in library_ranks:
            grinder_ranks[i[1]] = i[2]

    clusterid_to_orig_seqs = dict()
    orig_seqs_to_clusterid = dict()
    with open(args.bartender_barcode,"r") as f:
        bartender_barcode = csv.DictReader(f, delimiter=",")
        for i in bartender_barcode:
            real_barcode = args.actual_barcode_pattern
            for j in i["Unique.reads"]:
                real_barcode = re.sub("N",j,real_barcode,count=1)
            try:
                clusterid_to_orig_seqs[i["Cluster.ID"]].push(real_barcode)
            except:
                clusterid_to_orig_seqs[i["Cluster.ID"]] = [real_barcode]
            try:
                orig_seqs_to_clusterid[real_barcode].push(i["Cluster.ID"])
            except:
                orig_seqs_to_clusterid[real_barcode] = [i["Cluster.ID"]]

    clusterid_to_counts = dict()
    with open(args.bartender_cluster,"r") as f:
        bartender_cluster = csv.DictReader(f, delimiter=",")
        for i in bartender_cluster:
            clusterid_to_counts[i["Cluster.ID"]] = i["time_point_1"]

    for barcode, barcode_id in orig_seq_to_id.items():
        try:
            grinder_ranks[str(barcode_id)]
        except:
            grinder_ranks[str(barcode_id)] = 0
            orig_seqs_to_clusterid[barcode] = "NA"

    print(orig_seqs_to_clusterid)


    clusterid_to_counts["NA"] = 0

    for barcode, barcode_id in orig_seq_to_id.items():
        print()
        print(barcode)
        print(barcode.upper())
        print(barcode_id)
        print(grinder_ranks[str(barcode_id)])
        print(orig_seqs_to_clusterid[barcode.upper()])
        #print(clusterid_to_counts[orig_seqs_to_clusterid[barcode]])

