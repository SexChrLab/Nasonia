# In this script, we are converting the GTF to BED file format

from optparse import OptionParser

# Parsing arguments from command line
parser = OptionParser(usage='python convert_gtf_to_bed.py <options>')
parser.add_option('--gtf_filename', dest='gtf_filename', action='store')
parser.add_option('--output_filename', dest='output_filename', action='store')
(options, args) = parser.parse_args()

outfile = open(options.output_filename, "w")
with open(options.gtf_filename, "r") as f:
    for line in f:
        if not line.startswith("#"):
            items = line.rstrip().split('\t')
            gene_information = items[8].split(';')
            for i in gene_information:
                if i.startswith(" gene_name"):
                    gene_name = i.split(' ')
                    # print (gene_name)
                    new_line = [items[0], items[3], items[4], items[1], items[2], items[5], items[6], items[7], gene_name[2]]
                    # print (new_line)
                    if new_line[3] == 'ENSEMBL':
                        print ("\t".join(new_line), file=outfile)

