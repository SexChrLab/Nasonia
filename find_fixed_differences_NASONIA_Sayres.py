"""
what: find fixed differences in Ngir vs Nvit using nasion.STAR.raw.vcf
    
input: nasion.STAR.raw.vcf

output: nasion.STAR.NgirFixedandDifferentSites.vcf

"""
#! /bin/python
# import modules
import csv # imports module to make and work with csv files
import sys # import the "sys" module (sys is short for "system")
from sys import argv # sys.argv is simply a list of all of the parameters on the command line when Python was run
import numpy as np

print(" ****** The program will keep only sites that are fixed and different **********")
rawVCF= argv[1] # gff orignal input file from NCBI
onlySitesWithFixedDifferences = argv[2] # output file of our new lines that we want to keep, save as a gtf
out_list=[] # empty list where the output lines will be added to 

def inPutAndFilter(): # function to take the input, gff GCF_000002325.3_Nvit_2.1_genomic.gff, and reformat lines so that they are in cufflinks format
    """inPutAndFilter is a function that reads in the gene_exp.diff file.
    Keeps lines that haven an fpkm value that is greater than the user defined value
    """
    with open(rawVCF, "r") as rawVCFInput: # open and read gene_exp.diff, our argv[1])
        for line in rawVCFInput: # for each item in the gff file/list do the following
            line = line.strip()
            splt = line.split("\t") # split the items as tab delimated
            homREF= "0/0"
            homALT= "1/1"
            keep_info = splt[0:]
            genotypes = splt[9:]
            header = splt[0]
            header0 =header.split("=")
           # for i in line:
 #               print(line[0:180])
            for k in genotypes:
               # print(k[0:180])
                sample1 = genotypes[0].split(":")
                sample2 = genotypes[1].split(":")
                sample3 = genotypes[2].split(":")
                sample4 = genotypes[3].split(":")
                sample5 = genotypes[4].split(":")
                sample6 = genotypes[5].split(":")
                if (sample1[0] == sample2[0] == sample3[0]) and (sample4[0] == sample5[0] == sample6[0]) and (sample1[0] != sample4[0]) and sample1[0] == homREF and sample4[0] == homALT:
#                    print(line)
                    out_list.append(list(keep_info))
                    break
            for i in header:
                if i[0]== "#":
                    out_list.append(list(keep_info))
                    break

        return out_list

def outlist(out):
    with open(onlySitesWithFixedDifferences, "w") as outfile:
        for i in range(len(out)-1):
            i2 = "\t".join(out[i])+"\n"
            outfile.write(i2)
        i3 = "\t".join(out[-1])
        outfile.write(i3)


def main():
    if len(sys.argv)!= 3:
        print("\n",
              "********** This program will **************","\n","\n",
              "Please follow the directions below to have this program run properly:","\n","\n",
              "This program will read in 1 input files IN THIS ORDER! all inputs are required:", "\n","\n"
              " argv[0] is the python script","\n"
                  "\t","fixedDifferences.py","\n"
              " argv[1] will be a raw unfiltered vcf file containing information for a samples in the project","\n"
                  "\t", "Will contain a HEADER", "\n"
                  "\t", "Example input: nasonia.STAR.raw.vcf","\n"
              " argv[2] will be the output file","\n"
                 "\t", "Example output: nasion.STAR.NgirFixedandDifferentSites.vcf","\n")
        print("Usage example: $ python fixedDifferences.py nasonia.STAR.raw.vcf nasion.STAR.NgirFixedandDifferentSites.vcf")
    else:
        out=inPutAndFilter()
        outlist(out)


if __name__ == "__main__": # progam will not begin running until user has invoked the main function directly 
    main()

