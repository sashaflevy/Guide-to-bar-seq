#!/usr/bin/python3

import argparse
import re

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("abundances",help="",type=str)
    parser.add_argument("starcoded",help="",type=str)
    parser.add_argument("output",help="",type=str)
    args = parser.parse_args()

    clone_to_codes = dict()
    code_to_clone = dict()
    codes_to_input_counts = dict()
    codes_to_output_counts = dict()

    with open(args.abundances,"r") as f:
        for i in f:
            (count, rest) = re.split(r'\s+',i.strip())
            (clone, code) = re.split(r'_',rest)
            try:
                clone_to_codes[clone].append(code)
            except:
                clone_to_codes[clone] = [code]
            try:
                code_to_clone[code].append(clone)
            except:
                code_to_clone[code] = [clone]
            try:
                codes_to_input_counts[code] += count
            except:
                codes_to_input_counts[code] = count

    with open(args.starcoded,"r") as f:
        for i in f:
            (code,count,all_codes,idz) = re.split(r'\s+',i.strip())
            try:
                codes_to_output_counts[code] += count
            except:
                codes_to_output_counts[code] = count

    # Don't need this first one, it already handles htis okay
    #codes_input_only = set(map(lambda x: x.upper(), codes_to_input_counts.keys())).difference(set(codes_to_output_counts.keys()))
    codes_output_only = set(codes_to_output_counts.keys()).difference(set(map(lambda x: x.upper(), codes_to_input_counts.keys())))


    with open(args.output,"w") as f:
        f.write("clone_id\tcode\tinput_counts\toutput_counts\n")
        for clone, codes in clone_to_codes.items():
            for each_code in codes:
                f.write(clone+'\t')
                f.write(each_code+'\t')
                try:
                    f.write(codes_to_input_counts[each_code]+'\t')
                except:
                    f.write("0"+'\t')
                try:
                    f.write(codes_to_output_counts[each_code.upper()]+'')
                except:
                    f.write("0"+'')
                f.write("\n")
        for each_code in codes_output_only:
            f.write("NA"+'\t')
            f.write(each_code+'\t')
            f.write("0"+'\t')
            f.write(codes_to_output_counts[each_code.upper()]+'')
            f.write("\n")

