import re
import os
import numpy as np
import pandas as pd

#The working directory to store all the sequencing data, will be created if not exists
os.environ['WORKDIR'] = '/Users/helen/Desktop/bioinformatic/script_rna_pipeline'
#The directory for the reference genome
##os.environ['GENOMEDIR'] = 'RNA_seq/genomes/Porphyra/'
print(os.getcwd())
#examine the RNA_seq data
##ls -lth $WORKDIR/*/*.fq

##bash rna_analysis.sh -gtf -fa -w -CPU

###################################STEP1: quality control###########################
# raw_reads_count = open('raw_reads_summary.txt', 'r')
# for raw in raw_reads_count:
# 	raw_detail = raw.split(" ")
# 	library = raw_detail[0]
# 	print (library)