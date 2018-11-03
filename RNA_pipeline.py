## this is used to summary the quality check step during trimming, error correction and alignment
## aim to get all quality control data into the final table for easier review

import re
import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict

##getting raw reads information [ library / # of left reads / # of right reads ]
raw_reads_count = open('raw_reads_summary.txt', 'r')
raw_reads_dict = dict()
for raw in raw_reads_count:
	raw = raw.replace('\n', '')
	raw_detail = raw.split(" ")
	library = raw_detail[0]
	## left reads number, dict element 0
	raw_reads_dict.setdefault (library,[]).append (raw_detail[1])
	## right reads number, dict element 1
	raw_reads_dict.setdefault (library,[]).append (raw_detail[2])
#print(raw_reads_dict)

trimming = open('trimming_summary.txt', 'r')
trim_dict = dict()
for line in trimming:
	if line.startswith ("remove adapters by bbduk for"):
		library = re.search("remove adapters by bbduk for (.*)", line)
		#print (library.group(1))
	## number of total reads (left + right)
	elif line.startswith ("Input:"):
		Input = line.split()
		#print (Input[1])
	elif line.startswith ("Result"):
		trimmed = line.split()
		trimmed_reads = trimmed[1]
		trimmed_reads_percentage = trimmed[1]+trimmed[3]
		trimmed_bases = trimmed [4]
		## number of reads after trimming, dict element 2
		raw_reads_dict.setdefault(library.group(1),[]).append(trimmed_reads)
		## number of reads and the percentage of reads kept after trimming, dict element 3
		raw_reads_dict.setdefault(library.group(1),[]).append(trimmed_reads_percentage)
		## number of bases kept after trimming, dict element 4
		raw_reads_dict.setdefault(library.group(1),[]).append(trimmed_bases)
#print (raw_reads_dict)

rcorrector = open('rcorrector_summary.txt', 'r')
for line1 in rcorrector:
	if line1.startswith("perform error correction by rcorrector for"):
		library = re.search("perform error correction by rcorrector for (.*)", line1)
	elif "Corrected" in line1:
		corrected = line1.split()
		percentage = round (float(corrected[1])/float(raw_reads_dict[library.group(1)][4])*100, 2)
		corrected_base = corrected[1] + "(" + str(percentage) +"%" + ")"
		## number of bases corrected and the percentage of the corrected bases, dict element 5
		raw_reads_dict.setdefault(library.group(1),[]).append(corrected_base)

root = os.path.abspath("./STAR")
for root, dirs, files in os.walk(root):
	for name in files:
		if name.endswith (".stats"):
			library2 = name.split(".")[0]
			full_path = os.path.join(root, name)
			#print (full_path)
			stats_file = open(full_path,'r')
			for line2 in stats_file:
				if re.search("raw total sequences:",line2):
					reads_aligned = line2.split(":")[1]
					## number of reads mapped to the reference
					reads_aligned = reads_aligned.replace('\n', '')
					## percentage of reads that mapped to the reference
					reads_aligned_percentage = round (float(reads_aligned)/float(raw_reads_dict[library2][2])*100, 2)
					reads_aligned_reads_percentage = reads_aligned + "(" + str(reads_aligned_percentage) +"%" + ")"
					## number of reads aligned to the reference and the percentage of the mapped reads, dict element 6
					raw_reads_dict.setdefault(library2,[]).append(reads_aligned_reads_percentage) 

print (raw_reads_dict)

sys.stdout = open ('mapped_quality_check.txt', 'w')
print ('{:<8} {:<15} {:<20} {:<25} {:<15} {:<15}'.format('library', 'left_raw_reads','right_raw_reads','reads_after_trimming','bases corrected','reads_mapped')) 
for k, v in raw_reads_dict.items():
    left_raw_reads,right_raw_reads,reads_after_trimming,reads_percentage_after_trimming, bases_after_trimming, bases_corrected,reads_mapped = v
    print ('{:<8} {:<15} {:<20} {:<25} {:<15} {:<15}'.format(k, left_raw_reads,right_raw_reads, reads_percentage_after_trimming, bases_corrected,reads_mapped))
