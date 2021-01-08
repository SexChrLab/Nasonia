#!/bin/sh
#$ -cwd #tells OGS to run the script in your current working directory
#$ -N download_NasoniaFastq #name of the job 
#$ -m abe # email with the script Aborts, Begins, and Ends. Use -m e to only be altered when the job ends. 
#$ -M olney.kimberly@mayo.edu # email address
#$ -l h_vmem=8G #All OGS resource requests start with an l -- this example is requesting 8 GB of RAM for this run.
#$ -q 1-hour # submit is job for 1 hour # -q 4-days may also be used
#$ -o /research/labs/neurology/fyer/m239830/stdout # path for ells OGS where you want your standard output 
#$ -e /research/labs/neurology/fyer/m239830/stdout # same as above for error output  
#$ -notify

# In the currect directory
cd /research/labs/neurology/fyer/m239830/Nasonia
conda activate nasonia_environment_210108

# We can now efetch and format the search result as so-called runinfo data type. The runinfo type is in CSV, comma separated value, format:
esearch -db sra -query PRJNA260391 | efetch -format runinfo > runinfo.csv

# The resulting runinfo.csv file is quite rich in metadata and contains a series of attributes for the run. 

#To filter for lines that match SRR and to store the first ten run ids in a file we could write:

cat runinfo.csv | cut -f 1 -d ',' | grep SRR > runids.txt

# download the fastq files using parallel
cat runids.txt | parallel fastq-dump --split-files {}

# gzip fastq files
gzip *.fastq
mkdir RNA
mv *fastq* RNA
rm runids.txt
rm runinfo.csv

# We can now efetch and format the search result as so-called runinfo data type. The runinfo type is in CSV, comma separated value, format:
esearch -db sra -query PRJNA299670 | efetch -format runinfo > runinfo.csv

# The resulting runinfo.csv file is quite rich in metadata and contains a series of attributes for the run. 

#To filter for lines that match SRR and to store the first ten run ids in a file we could write:

cat runinfo.csv | cut -f 1 -d ',' | grep SRR > runids.txt

# download the fastq files using parallel
cat runids.txt | parallel fastq-dump --split-files {}

# gzip fastq files
gzip *.fastq
mkdir RNA
mv *fastq* RNA
