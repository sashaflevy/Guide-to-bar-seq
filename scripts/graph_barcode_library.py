#!/usr/bin/python3

# using
import argparse
import yaml
import re
import numpy
import jellyfish
#import os
#import random
#import datetime
#import csv


complement_dictionary = {
    'A':'T', 'a':'t',
    'C':'G', 'c':'g',
    'T':'A', 't':'a',
    'G':'C', 'g':'c'
    }


# recycling code from that SoBaSeq optimization business
# this returns a random base
# TODO should reimplement as being a list of letters to return and relative 
# frequencies, just so it plays nice for testing informative nature of biased
# base mixing. Should take that in as a base_mix dictionary.
def randBase():
    return str(["a","c","t","g"][int(numpy.floor(numpy.random.uniform(0,4,1)))])


# This takes a pattern, and replaces every "N" in there with one of 4 bases.
# TODO implement passing the base_mix dictionary through to the randBase().
def randBarcode(pattern):
    for i in range(re.subn("N","N",pattern)[1]):
        pattern = re.sub("N",randBase(),pattern,count=1)
    return(pattern)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("barcode_file",help="file",type=str)
    parser.add_argument("--yaml",help="is yaml?",action="store_true")
    parser.add_argument("outputbase",help="out",type=str)
    parser.add_argument("--store-graph",help="store a graph?",action="store_true")
    parser.add_argument("--print-as-matrix",help="do an upper triangle?",
        action="store_true")
    args = parser.parse_args()

    all_barcodes = []
    with open(args.barcode_file,"r") as f:
        if args.yaml:
            barcode_yaml = yaml.load(f)
            for i in barcode_yaml:
                all_barcodes.append(list(i.keys())[0])
        else:
            for i in f:
                if re.match("[ATCG]",i):
                    all_barcodes.append(i.strip())

    results_matrix = []
    graph = []
    for i in range(len(all_barcodes)-1):
        results_row = []
        for j in range(len(all_barcodes)-1):
            if args.print_as_matrix:
                if j >= i :
                    dist = jellyfish.levenshtein_distance(all_barcodes[i],all_barcodes[j])
                    results_row.append(dist)
                    graph.append([all_barcodes[i],all_barcodes[j],dist])
                else:
                    results_row.append(" ")
            else:
                if j != i :
                    dist = jellyfish.levenshtein_distance(all_barcodes[i],all_barcodes[j])
                    results_row.append(dist)
                    if j >= i :
                        graph.append([all_barcodes[i],all_barcodes[j],dist])
        results_matrix.append(results_row)

    with open(args.outputbase+".tsv","w") as f:
        for i in results_matrix:
            f.write( " ".join(str(j) for j in i)+"\n" );

    if args.store_graph:
        with open(args.outputbase+".graph","w") as f:
            for i in graph:
                f.write( " ".join(str(j) for j in i)+"\n" );
