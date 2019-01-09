#!/usr/bin/python3

import argparse
import random
import re

if __name__ == '__main__':

    parser = argparse.ArgumentParser("")
    parser.add_argument("fasta",type=str)
    parser.add_argument("depth",type=int)
    parser.add_argument("output",type=str)
    args = parser.parse_args()

    that_lib = list()
    with open(args.fasta,"r") as f:
        for i in f:
            that_lib.append([ re.sub("> ","",str(i.strip())), next(f).strip() ])

    with open(args.output,"w") as f:
        for i in range(args.depth):
            j = random.randint(0,len(that_lib)-1)
            f.write("> "+that_lib[j][0]+"\n"+that_lib[j][1]+"\n")

