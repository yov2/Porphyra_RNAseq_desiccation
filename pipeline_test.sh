#!/bin/bash
while [[ $# -gt 0 ]]; do
	key = "$1"
	case $key in
		-g|--genome)
		GENOME="$2"
		shift
		;;
		-w|--workdir)
		WORKDIR="$2"
		shift
		;;
		-CPU|--cpu)
		CPU="$2"
		shift
		;;
		-h|--help)
		echo "Usage: ./raw_fastq.sh -g < GENOME> -w <WORKDIR> -CPU <CPU>"
		exit
		;;
		#unknow args
		*)
		echo "unknown option: $key, exiting."
		echo "usage: ./raw_fastq.sh -g < GENOME> -w <WORKDIR> -CPU <CPU>"
		exit
		;;
	esac
	shift
done

## check $WORKDIR
if [[ ! -d $WORKDIR ]]; then
	echo "Could not find working directory: $WORKDIR, exiting. Please make sure the working directory exists"
	exit 1
else
	##convert to absolute path
	WORKDIR = $(readlink -e $WROKDIR)
	echo "WORKDIR = $WORKDIR"
fi

##check $GENOME
if [[ ! -d $GENOME ]]; then
	echo "Could not find reference genome direcotry: $GENOME, exiting. Please make sure the reference genome directory exists"
	exit 1
else 
	## convert to absolute path
	GENOME = $(readlink -e $GENOME)
	echo "GENOME = $GENOME"

	#check if annotation gff/gtf file and sequence .fa file exist in the genome reference directory
	GENOME_GTF = "$GENOME/*.gtf"
	GENOME_GFF = "$GENOME/*.gff"
	GENOME_FA = "$GENOME/*.fa"
	if [[ ! -f $GENOME_GTF ]]; then
		echo "$GENOME_GTF not found in the genome directory, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_FA ]]; then
		echo "$GENOME_FA not found in the genome directory, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_GFF]]; then
		echo "$GENOME_GFF not found in the genome directory, exiting"
		exit 1
	fi
	STAR_INDEX="$GENOME/star/genome_index"
fi

##CHECK if STAR genome index is exist
if [[ ! -d $STAR_INDEX ]]; then
	echo "STAR index does not exist, building STAR index"
	mkdir $STAR_INDEX
	STAR \
		--runThreadN 24 \
		--runMode genomeGenerate \
		--genomeDir $STAR_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbGTFfile $GENOME_GTF \
		--sjdbOverhang 99 
fi

mkdir -p fastq
mkdir -p fastqc_raw
mkdir -p fastqc_trimming
mkdir -p trimming
mkdir -p rcorrector
mkdir -p STAR

fastq="$WORKDIR/fastq"
fastqc_raw="$WORKDIR/fastqc_raw"
raw_reads="$WORKDIR/raw_reads"
bbmap="$WORKDIR/bbmap"
trimming="$WORKDIR/trimming"
PREFIX="Sample_"
BBDUK_PATH="$bbmap/bbmap/bbduk.sh"
ADAPTER_PATH="$bbmap/bbmap/resources/adapters.fa"
rcorrector="$WORKDIR/rcorrector"
fastqc_trimming="$WORKDIR/fastqc_trimming"
STAR="$WORKDIR/STAR"

##CHECK if bbmap is installed in the right place (the wget didn't work in server premise)
if [[ ! -d $bbmap ]]; then
	echo "need to install bbduk from the source"
	mkdir bbmap
	cd bbmap
	URL_bbmap="https://sourceforge.net/projects/bbmap/files/BBMap_37.93.tar.gz"
	wget $URL_bbmap 
	echo "extract bbmap ... "
	tar xzvf "BBMap_*.tar.gz"
	cd $WORKDIR
fi

echo "checking the format of input files and file names"
if [[ ! -d $raw_reads ]]; then
	echo "Input files should be in the file named as raw_reads under the working direcotry"
	echo "NO input files found, exiting"
	exit 1
else
	cd $raw_reads
	for fq in $(ls); do
		if [[ ! $fq  == $PREFIX* ]]; then 
			echo "raw reads for each library should be in sub_file inside file raw_reads"
			echo "each sub_file should named as Samplel+ibraryname, connected with '_', S needs to be upper case, for example Sample_D1"
			echo "$fq doesn't qualify for naming rules, exiting"
			exit 1
		else 
			###########STEP 1: combine all fastq reads from the same individual ############
			basename=$(echo $fq | cut -f2 -d '_')
			echo "Start to combine raw fastq reads from the same library for $basename"
			cd $raw_reads/$PREFIX$basename
			fq1="_1.fastq.gz"
			fq2="_2.fastq.gz"
			fq1=$basename$fq1
			fq2=$basename$fq2
			mkdir $fastq/$PREFIX$basename
			cat *R1*.fastq.gz > $fastq/$PREFIX$basename/$fq1
			cat *R2*.fastq.gz > $fastq/$PREFIX$basename/$fq2
			
			##summary of # of raw reads in each library, and input summary results into reads_track.txt
			No_left=$(gzip -dc $fastq/$PREFIX$basename/$fq1 | grep -c '@')
			echo "the number of left raw reads for $basename is $No_left" 
			No_right=$(gzip -dc $fastq/$PREFIX$basename/$fq2 | grep -c '@')
			echo "the number of right raw reads for $basename is $No_right" 
			## check the number of left reads is equal to the number of right reads in each library
			if [ ! $No_left == $No_right ]; then 
				echo "WARNNING: The # of left reads is not equal to the # of right reads in $basename"
				echo "PLEASE check raw files"
				exit 1
			fi
			## output reads number into a text file for quality check purpose
			echo "$basename $No_left $No_right" >> $WORKDIR/raw_reads_summary.txt

			########################STEP 2: FASTQC of raw reads ######################
			echo "fastqc check of the raw reads for $basename"
			mkdir $fastqc_raw/$PREFIX$basename
			fastqc -f fastq \
				-t 24 \
				$fastq/$PREFIX$basename/$fq1 \
				$fastq/$PREFIX$basename/$fq2 \
				--outdir="$fastqc_raw/$PREFIX$basename/"

			#########################STEP 3: remove adapters using bbduk ############
			echo "remove adapters by bbduk for $basename" >> trimming_summary.txt
			suffix="BBDUK_"
			mkdir $trimming/$PREFIX$basename/
			source $BBDUK_PATH \
				in1=$fastq/$PREFIX$basename/$fq1 \
				in2=$fastq/$PREFIX$basename/$fq2 \
				out1=$trimming/$PREFIX$basename/$suffix$fq1 \
				out2=$trimming/$PREFIX$basename/$suffix$fq2 \
				ref=$ADAPTER_PATH \
				ktrim=r k=23 mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=2 \
				1>&2 2>> $trimming/trimming_summary.txt

			#####################STEP 4: perform error correction using rcorrect##########
			echo "perform error correction by rcorrector for $basename" >> $rcorrector/rcorrector_summary.txt
			mkdir $rcorrector/$PREFIX$basename/
			run_rcorrector.pl \
				-1 $trimming/$PREFIX$basename/$suffix$fq1 \
				-2 $trimming/$PREFIX$basename/$suffix$fq2 \
				-od  $rcorrector/$PREFIX$basename/ \
				-t 24 \
				1>&2 2>> $rcorrector/rcorrector_summary.txt

			## the squence names were changed by rcorrecotor, not accepted by Trinity
			gunzip $rcorrector/$PREFIX$basename/*.cor.fq.gz
			sfq1="_1.cor.fq"
			sfq2="_2.cor.fq"
			fq1=$suffix$basename$sfq1
			fq2=$suffix$basename$sfq2
			rename="rcorrector_renamed_"
			cut -d ' ' -f1 $rcorrector/$PREFIX$basename/$fq1 | sed '1~4 s/$/\/1/g' > $rcorrector/$PREFIX$basename/$rename$basename$sfq1
			cut -d ' ' -f1 $rcorrector/$PREFIX$basename/$fq2 | sed '1~4 s/$/\/2/g' > $rcorrector/$PREFIX$basename/$rename$basename$sfq2

			#####################STEP 5: FASTQC for clean reads #########################
			echo "perform FASTQC quality check for adapters removed reads for $basename"
			mkdir $fastqc_trimming/$PREFIX$basename
			fastqc -f fastq -t 24 \
				$rcorrector/$PREFIX$basename/$rename$basename$sfq1 \
				$rcorrector/$PREFIX$basename/$rename$basename$sfq2 \
				--outdir="$fastqc_trimming/$PREFIX$basename/"

			#####################STEP6: reads alignment by STAR #########################
			echo "Aligning $basename to the reference genome"
			mkdir $STAR/$PREFIX$basename
			STAR \
				--runThreadN 24 \
				--genomeDir $STAR_INDEX \
				--sjdbGTFfile $GENOME_GTF \
				--outSAMtype BAM Unsorted SortedByCoordinate\
				--outFileNamePrefix $STAR/$PREFIX$basename/$basename \
				--readFilesIn $rcorrector/$PREFIX$basename/$rename$basename$sfq1 $rcorrector/$PREFIX$basename/$rename$basename$sfq2 \
				--outReadsUnmapped Fastx \
				--quantMode GeneCounts
		fi
	done
fi
