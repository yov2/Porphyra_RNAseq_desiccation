#!/bin/bash
# This script start the pipeline from raw fastq files assuming them are located in ? relative to WORKDIR

## input raw reads should be fair-end fastq reads from ILLUMINA sequencing
## need to make sure that bbduk is installed
########step 0: check command line args and if fastq files exist ###############

###input the data and genome directory
WORKDIR = "RNA_seq/desiccation/raw_reads/"
GENOME = "$HOME/genome/porphyra/"

##Parse args
while [[$# -gt 0]]; do
	key = "$1"
	case $key in 
		-g|--genome)
		GENOME = "$2"
		shift
		;;
		-w|--workdir)
		WORKDIR = "$2"
		shift
		;;
		-CPU|-cpu)
		CPU = "$2"
		shift
		;;
		-h|--help)
		echo "Usage: ./raw_fastq.sh -g < GENOME> -w <WORKDIR> -CPU <CPU>"
		exit
		;;
		#unknow args
		*)
		echo "unknown option: $key, exiting."
		echo "usage: ./raw_fastq.sh -g < GENOME> -w <WORKDIR>"
		exit
		;;
	esac
	shift
done

## if it is running on Premise, the max of CPU is 24 per nodes
N_CPUS=$CPU

## check $WORKDIR
if [[ ! -d $WORKDIR ]]; then
	echo "Could not find working directory: $WORKDIR, exiting. Please make sure the working directory exists"
	exit 1
else
	##convert to absolute paths
	GENOME = $(readlink -e $GENOME)
	WORKDIR = $(readlink -e $WROKDIR)
	echo "GENOME = $GENOME"
	echo "WORKDIR = $WORKDIR"
fi

##check $GENOME
if [[ ! -d $GENOME ]]; then
	echo "Could not find reference genome: $GENOME, exiting. Please make sure the working directory exists"
	exit 1
else 
	#check if annotation gff/gtf file and sequence .fa file exist
	GENOME_GTF = "$GENOME/*.gtf"
	GENOME_GFF = "$GENOME/*.gff"
	GENOME_FA = "$GENOME/*.fa"
	if [[ ! -f $GENOME_GTF ]]; then
		echo "$GENOME_GTF not found, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_FA ]]; then
		echo "$GENOME_FA not found, exiting"
		exit 1
	fi
	if [[ ! -f $GENOME_GFF]]; then
		echo "$GENOME_GFF not found, exiting"
		exit 1
	fi
	STAR_INDEX="$GENOME/star/genome_index"
fi

## index genome reference if STAR index not exists
if [ ! -d $STAR_INDEX]; then
	echo "STAR index not exist, build STAR index"
	STAR \
		--runThreadN $CPU \
		--runMode genomeGenerate \
		--genomeDir $STAR_INDEX \
		--genomeFastaFiles $GENOME_FA \
		--sjdbGTFfile $GENOME_GTF \
		--sjdbOverhang 99
fi

##check if raw fastq files is complete fair-end sequences?
## create dirs if not exists
HOME = "/mnt/lustre/klein/yov2/pipeline_test/"
mkdir -p fastq
mkdir -p fastqc_raw
mkdir -p fastqc_trimming
fastq="$HOME/fastq"
fastqc_raw="$HOME/fastqc_raw"

###########STEP 1: combine all fastq reads from the same individual ############
echo "Start to combine all raw fastq reads from the same library"
cd raw_reads
for fq in $(ls); do
	basename=$(echo $fq | cut -f2 -d '_')
	echo "combining raw reads for $basename"
	cd Sample_$basename
	fq1="_1.fastq.gz"
	fq2="_2.fastq.gz"
	fq1=$basename$fq1
	fq2=$basename$fq2
	PREFIX="Sample_"
	mkdir $fastq/$PREFIX$basename
	cat *R1*.fastq.gz > $fastq/$PREFIX$basename/$fq1
	cat *R2*.fastq.gz *R3*.fastq.gz > $fastq/$PREFIX$basename/$fq2
	##summary of # of raw reads in each library
	No_left=$(gzip -dc $fastq/$PREFIX$basename/$fq1 | grep -c '@')
	echo "the number of left reads for $basename is $No_left" 
	No_right=$(gzip -dc $fastq/$PREFIX$basename/$fq2 | grep -c '@')
	echo "the number of left reads for $basename is $No_right" 
	## check raw reads fastqc
	mkdir $fastqc_raw/$basename
	fastqc -f fastq -t 24 \
		$fastq/$PREFIX$basename/$fq1 \
		$fastq/$PREFIX$basename/$fq2 \
		--outdir="$fastqc_raw/$basename/"
	cd ..
done

##########STEP 2: remove adapters using bbduk ############
echo "remove adapters by bbduk"
bbmap="$HOME/bbmap/"
## the wget didn't work in server premise
if [ ! -d $bbmap ]; then
	echo "need to install bbduk from the source"
	mkdir bbmap
	cd bbmap
	URL_BASE="https://sourceforge.net/projects/bbmap/files/BBMap_37.93.tar.gz"
	wget "$URL_BASE"
	tar xzvf "BBMap_*.tar.gz"
else
	cd fastq
	for fq in $(ls); do
		basename=$(echo $fq | cut -f2 -d '_')
		echo "perform adapter trimming for $basename" 
		fq1="_1.fastq.gz"
		fq2="_2.fastq.gz"
		fq1=$basename$fq1
		fq2=$basename$fq2
		suffix="trimming."
		prefix="Sample_"
		SCRIPT_PATH="$bbmap/bbmap/bbduk.sh"
		ADAPTER_PATH="$bbmap/bbmap/resources/adapters.fa"
		mkdir $bbmap/bbmap/$prefix$basename/
		source $SCRIPT_PATH \
			in1=$prefix$basename/$fq1 \
			in2=$prefix$basename/$fq2 \
			out1=../bbmap/$prefix$basename/$suffix$fq1 \
			out2=../bbmap/$prefix$basename/$suffix$fq2 \
			ref=$ADAPTER_PATH \
			ktrim=r k=23 mink=11 hdist=1 tbo=f tpe=f qtrim=r trimq=2 
	done
	cd ..

#######STEP 3: perform error correction using rcorrect##########
echo "perform error correction by rcorrector"
mkdir rcorrector
RCORRECT="$HOME/rcorrector"
## combine all these basename/fq1/fq2 stuff outside the loop
cd $HOME/bbmap
for fq in $(ls); do
		basename=$(echo $fq | cut -f2 -d '_')
		echo "perform error correction for $basename"
		fq1="_1.fastq.gz"
		fq2="_2.fastq.gz"
		fq1=$basename$fq1
		fq2=$basename$fq2
		suffix="trimming."
		prefix="Sample_"
		SCRIPT_PATH="$bbmap/bbmap/bbduk.sh"
		ADAPTER_PATH="$bbmap/bbmap/resources/adapters.fa"
run_rcorrector.pl 
	-1 D${fn}/bbduk_D${fn}_R1.fastq.gz \
	-2 D${fn}/bbduk_D${fn}_R2.fastq.gz \
	-od  D${fn}/rcorrector_bbduk_D${fn} \
	-t 24












