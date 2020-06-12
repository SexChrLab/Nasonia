# In this script, we are converting the output of ASE from tsv to BED file format

from optparse import OptionParser

# Parsing arguments from command line
parser = OptionParser(usage='python convert_ase_to_bed.py <options>')
parser.add_option('--ase_filename', dest='ase_filename', action='store')
parser.add_option('--output_filename', dest='output_filename', action='store')
(options, args) = parser.parse_args()

outfile = open(options.output_filename, "w")
with open(options.ase_filename, "r") as f:
    for line in f:
        if not line.startswith("contig"):
            items = line.rstrip().split('\t')
            items.insert(1, str(int(items[1])-1))
            print ("\t".join(items), file=outfile)

