#!/usr/bin/python3

# using
import argparse
import yaml
import re
import numpy


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
            "degenerate barcode library file, so each barcode and the lineage "+
            "name that it maps to.")
    parser.add_argument("pattern",help="Pattern of the barcode, specified in "+
            "IUPAC codes. So, this can be the full construct.",
        type=str)
    parser.add_argument("outbase",help="The base filename to write out to. "+
            "This should describe the parameters, basically. It'll then have "+
            "a FASTA and a YAML generated onto that basename.",
        type=str)
    parser.add_argument("--number-lineages",help="The number of different "+
            "lineages to try to barcode.",
        default=100, type=int)
    parser.add_argument("--number-barcoded-clones",help="Essentially, the "+
            "number of barcodes in the library. This is probably going to be "+
            "larger than the number of lineages",
        default=1000, type=int)
    parser.add_argument("--replicate",help="A replicate ID to prevent "+
            "filename collisions.",
        default="", type=str)
    args = parser.parse_args()

    baseNameOut = args.outbase+"_"+args.pattern+"_"+str(args.number_lineages)+"lineages_"+str(args.number_barcoded_clones)+"clones_replicate"+args.replicate


    barcode_to_lineage = list()

    for i in range(args.number_barcoded_clones):
        barcode_to_lineage.append(
            {   randBarcode(args.pattern):
                int(numpy.floor(numpy.random.uniform(0,args.number_lineages,1)))
                }
            )

    with open(baseNameOut+".yaml","w") as f:
        f.write(yaml.dump(barcode_to_lineage))

    with open(baseNameOut+".fasta","w") as f:
        for i in barcode_to_lineage:
            f.write( "> "+ str(list(i.values())[0]) +"\n"+ list(i.keys())[0] +"\n")


