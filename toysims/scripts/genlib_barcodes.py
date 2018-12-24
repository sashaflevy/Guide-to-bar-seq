#!/usr/bin/python3

import argparse
import os
import sys
import yaml
import maps
import hashlib 
import datetime
import isodate
import re 
import shutil
import argparse
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
def randBase():
    return str(["a","c","t","g"][int(numpy.floor(numpy.random.uniform(0,4,1)))])

# takes a string, mutates those bases based on a lambda
def mutator(x,mutations=0):
    numberToMutate = numpy.random.poisson(lam=mutations)
    positionsToMutate = numpy.floor(numpy.random.uniform(0,len(x),numberToMutate))
    outputString = list(x)
    if len(outputString) == 0:
        return ""
    for eachMutation in positionsToMutate:
        outputString[int(eachMutation)] = randBase()
    return "".join(outputString)

def randBarcode(pattern):
    for i in range(re.subn("N","N",pattern)[1]):
        pattern = re.sub("N",randBase(),pattern,count=1)
    return(pattern)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="this thing makes libraries of barcodes, with or without simulated reads per")
    parser.add_argument("pattern",help="pattern of the barcode, with N as random bases",type=str)
    parser.add_argument("outfile",help="where to write out",type=str)
    parser.add_argument("--number-variants",default=100,type=float)
    parser.add_argument("--number-barcoded-clones",default=1000,type=float)
#    parser.add_argument("--counts-exponent",default=1,type=float)
#    parser.add_argument("--counts-mean",default=10,type=float)
#    parser.add_argument("--counts-noise",default=1,type=float)
#    parser.add_argument("--per-base-error-rate",default=0.01,type=float)
    args = parser.parse_args()

    barcode_to_variant = list()

    for i in range(args.number_barcoded_clones):
        barcode_to_variant.append(
            {   randBarcode(args.pattern):
                int(numpy.floor(numpy.random.uniform(0,args.number_variants,1)))
                }
            )

    with open(args.outfile,"w") as f:
        for i in barcode_to_variant:
            f.write( list(i.keys())[0] +"\t"+ str(list(i.values())[0]) +"\n" )


#        for i, afile in enumerate(this_entry["f"]):
#            if not os.path.isfile(afile):
#                raise Exception("there ain't no file "+afile)
#
#        this_entry["d"] = str(this_entry["d"])
#        if not re.match("^20",this_entry["d"]):
#            this_entry["d"] = "20"+this_entry["d"]
#
#        this_entry["d"] = isodate.parse_date(this_entry["d"])
#
#        entry_dir = config["notebook_directory"]+"/"+this_hash_hex
#        print(entry_dir)
#        os.makedirs(entry_dir,exist_ok=True)
#
#        for i, afile in enumerate(this_entry["f"]):
#            try:
#                shutil.copyfile(
#                    this_entry["f"][i],
#                    entry_dir+"/"+this_hash_hex+"_"+this_entry["f"][i])
#            except:
#                raise Exception("copying "+afile+" failed!")
#            this_entry["n"] = re.sub(afile,this_hash_hex+"_"+afile,this_entry["n"])
#            this_entry["f"][i] = this_hash_hex+"_"+this_entry["f"][i]
#
#        with open(entry_dir+"/"+this_hash_hex+"_notes.txt","w") as notefile:
#            notefile.write(yaml.dump(this_entry))
#            
#
#    os.makedirs(config["trash"],exist_ok=True)
#    for this_entry in these_entries:
#        for i, afile in enumerate(this_entry["f"]):
#            this_real_filename = re.sub(r"^[0-9a-z]{10}_","",this_entry["f"][i])
#            try:
#                shutil.move(
#                    this_real_filename,
#                    config["trash"]+"/"+this_real_filename
#                    )
#            except:
#                pass
#
#    shutil.move(
#        config["intake_file"],
#        config["trash"]+"/"+this_hash_hex+"_notes.txt"
#        )
#
#    shutil.copyfile(
#        config["template_file"],
#        config["intake_file"]
#        )
#
#
#    with open(vars(args)["input-fasta"],"r") as input_fasta, open(vars(args)["output-basename"]+"_read1.fastq","w") as r1, open(vars(args)["output-basename"]+"_read2.fastq","w") as r2:
#
#            record_name = input_fasta.readline().strip()
#            record_seq  = input_fasta.readline().strip()
#
#            if record_name is "" or record_seq is "":
#                break
#
#            start_indicies = numpy.random.normal(
#                float(args.avg_insert),
#                float(args.insert_spread),
#                int(args.depth_per)
#                )
#
#            for i in start_indicies:
#                i = numpy.maximum(i,0)
#                i = numpy.minimum(i,len(record_seq)-args.read_length)
#                i = int(i)
#                record_name = re.sub(r">", "", record_name)
#                #
#                r1.write("@"+record_name+"\n")
#                r1.write(mutator(record_seq[(len(record_seq)-i):(len(record_seq)-i+args.read_length)],float(args.errors_per_base))+"\n")
#                r1.write("+"+"\n")
#                r1.write("E"*args.read_length+"\n")
#                #
#                r2.write("@"+record_name+"\n")
#                r2.write(mutator("".join(comp_dict.get(base,base) for base in reversed(record_seq[ (len(record_seq)-args.read_length):(len(record_seq)) ])),float(args.errors_per_base))+"\n")
#                r2.write("+"+"\n")
#                r2.write("E"*args.read_length+"\n")
#
