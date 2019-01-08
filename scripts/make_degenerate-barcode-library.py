#!/usr/bin/python3

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


iupac_dictionary = {
    "N": ["A","T","C","G"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"]
    }


# recycling code from that SoBaSeq optimization business
# this returns a random base
# TODO should reimplement as being a list of letters to return and relative 
# frequencies, just so it plays nice for testing informative nature of biased
# base mixing. Should take that in as a base_mix dictionary.
def randBase(possible=iupac_dictionary["N"],
        mixing_dict={ "A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}
        ):
    # Initializes two arrays, then subsets out only the nucleotides used
    ( letters, proportions) = ( [], [] )
    for i in possible:
        letters.append(i)
        proportions.append(mixing_dict[i])
    # Normalized to 1
    proportions = proportions / numpy.sum(proportions)
    # Picks a letter using these proportions
    return str(numpy.random.choice(letters,size=1,p=proportions)[0])


# This takes a pattern, and replaces every "IUPAC code" in there with one of
# the appropriate bases, using the mixing dictionary of proportions
def randBarcode(pattern, 
        mixing_dict={ "A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}
        ) :
    # For each code in the dictionary, go looking for it
    for each_code in iupac_dictionary.keys():
        # For each time we find it
        for i in range(re.subn(each_code,each_code,pattern)[1]):
            # Sub it out with a random base of the appropriate type!
            pattern = re.sub(each_code,
                randBase(iupac_dictionary[each_code],mixing_dict),
                pattern,count=1
                )
    return(pattern)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This is going to make a "+
            "degenerate barcode library, then associate that with lineages by "+
            "1 to 1 mapping or sampling barcodes for each lineage clone."
        )
    parser.add_argument("pattern",help="Pattern of the barcode, specified in "+
            "IUPAC codes. So, this can be the full construct if you want.",
        type=str)
    parser.add_argument("--mix",help="The nucelotide mix of A,T,C,G in the "+
            "barcodes. The script normalizes these to unity. For clarity's "+
            "sake, there's no option to alter the mix of individual "+
            "positions, except for just using other IUPAC codes in the "+
            "pattern.",
        type=str,default="0.25,0.25,0.25,0.25")
    parser.add_argument("--number-lineages",help="The number of different "+
            "lineages to try to barcode, so clones or genotypes barcoded with "+
            "a clone.",
        default=100, type=int)
    parser.add_argument("--barcodes-per-lineage",help="Barcodes per each "+
            "lineage Essentially, the number of barcodes in the library is "+
            "the product of this and the number of lineages.",
        default=3, type=float)
    parser.add_argument("--fixed-barcodes-per-clone",help="This turns on the "+
            "mode where a clone gets a fixed number of barcodes (above). "+
            "Otherwise, it gets a poisson distribution of barcodes with mean "+
            "as the number of barcodes-per-lineage.",
        type=bool)
    parser.add_argument("outbase",help="The base filename to write out to. "+
            "This should describe the parameters, basically. It'll then have "+
            "a FASTA generated onto that basename. The ID is the lineage ID, "+
            "then it's the barcode.",
        type=str)
    parser.add_argument("--replicate",help="A replicate ID to prevent "+
            "filename collisions, this is just tacked onto the end.",
        default="", type=str)
    args = parser.parse_args()

    baseNameOut = (
        args.outbase+"_"+args.pattern+"_"+
            re.sub(",","-",args.mix)+"_"+str(args.number_lineages)+"lineages_"+
            str(args.barcodes_per_lineage)+"barcodesper_fixedBarcodesPer"+
            str(args.fixed_barcodes_per_clone)+"_replicate"+args.replicate
        )

    print("Printing this out to : "+baseNameOut)

    barcode_to_lineage = list()

    mix_num = numpy.array( [float(i) for i in re.split(",",args.mix)] )
    mixing_dict = dict(zip( ["A","T","C","G"], mix_num ))

    for i in range(args.number_lineages):
        if args.fixed_barcodes_per_clone:
            number_of_barcodes = int(args.barcodes_per_lineage)
        else:
            number_of_barcodes = int(numpy.floor(numpy.random.poisson(args.barcodes_per_lineage,1)))
        for j in range(number_of_barcodes):
            this_barcode = randBarcode(args.pattern,mixing_dict=mixing_dict)
            barcode_to_lineage.append( [ this_barcode, i ] )

    with open(baseNameOut+".fasta","w") as f:
        for i in barcode_to_lineage:
            f.write( "> "+ str(i[1]) +"\n"+ str(i[0]) +"\n" )


