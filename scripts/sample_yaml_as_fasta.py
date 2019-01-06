#!/usr/bin/python3

import argparse
import yaml
import random

if __name__ == '__main__':

    parser = argparse.ArgumentParser("")
    parser.add_argument("yaml",type=str)
    parser.add_argument("depth",type=int)
    parser.add_argument("output",type=str)
    parser.add_argument("abundances",type=str)
    args = parser.parse_args()

    that_lib = list()
    with open(args.yaml,"r") as f:
        for i in list(yaml.load_all(f))[0]:
            that_lib.append([ list(i.keys())[0], i[list(i.keys())[0]]])

    with open(args.output,"w") as f, open(args.abundances,"w") as a:
        for i in range(args.depth):
            j = random.randint(0,len(that_lib)-1)
            f.write(">"+str(that_lib[j][1])+"\n"+that_lib[j][0]+"\n")
            a.write(str(that_lib[j][1])+"\n")

