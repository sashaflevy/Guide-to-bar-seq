#!/usr/bin/python3

# using
import argparse
import yaml
# not sure
import os
import sys
import maps
import hashlib 
import datetime
import isodate
import re 
import shutil
from Bio import Seq, SeqRecord
import numpy
import itertools
import re


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


# Probably won't need this one
## takes a string, mutates those bases based on a lambda
#def mutator(x,mutations=0):
#    numberToMutate = numpy.random.poisson(lam=mutations)
#    positionsToMutate = numpy.floor(numpy.random.uniform(0,len(x),numberToMutate))
#    outputString = list(x)
#    if len(outputString) == 0:
#        return ""
#    for eachMutation in positionsToMutate:
#        outputString[int(eachMutation)] = randBase()
#    return "".join(outputString)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This is going to make a "+
            "degenerate barcode library file, so each barcode and the variant "+
            "name that it maps to.")
    parser.add_argument("pattern",help="Pattern of the barcode, specified in "+
            "IUPAC codes. So, this can be the full construct.",
        type=str)
    parser.add_argument("outbase",help="The base filename to write out to. "+
            "This should describe the parameters, basically. It'll then have "+
            "a CSV and a YAML generated onto that basename.",
        type=str)
    parser.add_argument("--number-variants",help="The number of different "+
            "variants to try to barcode.",
        default=100, type=float)
    parser.add_argument("--number-barcoded-clones",help="Essentially, the "+
            "number of barcodes in the library. This is probably going to be "+
            "larger than the number of variants",
        default=1000, type=float)
    args = parser.parse_args()

    barcode_to_variant = list()

    for i in range(args.number_barcoded_clones):
        barcode_to_variant.append(
            {   randBarcode(args.pattern):
                int(numpy.floor(numpy.random.uniform(0,args.number_variants,1)))
                }
            )

    with open(args.outbase+".csv","w") as f:
        for i in barcode_to_variant:
            f.write( list(i.keys())[0] +","+ str(list(i.values())[0]) +"\n" )

    with open(args.outbase+".yaml","w") as f:
        f.write(yaml.dump(barcode_to_variant))

    with open(args.outbase+".fasta","w") as f:
        for i in barcode_to_variant:
            f.write( "> "+ str(list(i.values())[0]) +"\n"+ list(i.keys())[0] +"\n")


