
#!/bin/bash #

#DATA ORGANIZATION#

## Directories were created which will be used for the rest of the workflow.

## First we should make sure that we are in our home directory.
pwd

## Create ngs_course folder within home directory
mkdir ngs_course

## Create dnaseq directory for DNA-seq analysis within the ngs_course folder
mkdir ngs_course/dnaseq

## Change the current directory to dnaseq
cd ngs_course/dnaseq

## Create four directories(data, meta, results, log) inside dnaseq
mkdir data meta results log

## Verify that you have created the directories
ls -lF

## Change the current directory and go into data
cd ~/ngs_course/dnaseq/data

## Create untrimmed_fastq subdirectory for untrimmed reads
mkdir untrimmed_fastq

## Create trimmed_fastq subdirectory for trimmed reads
mkdir trimmed_fastq

## Download the raw_fastq and BED data to be used in the workflow
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

## Move the raw_fastq files to the untrimmed_fastq directory. Use * sign before fastq.qz to move multiple files in a single command
mv *fastq.qz ~/ngs_course/dnaseq/data/untrimmed_fastq

## Move the bed file to the data directory
mv annotation.bed ~/ngs_course/dnaseq/data

## Download the reference genome which is required to perform alignment later in the workflow
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

## Move the reference genome into data directory
mv hg19.fa.gz ~/ngs_course/dnaseq/data/

## Go to home directory
cd ~/

## Download miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

## Give execute permissions to the script (make it executable)
chmod +x ./ Miniconda3-latest-Linux-x86_64.sh

## Run the script with bash command to install miniconda
bash Miniconda3-latest-Linux-x86_64.sh

## Run the source command with the .bashrc file as an argument for the changes to be applied
source ~/.bashrc

## Add ('defaults' already in 'channels' list, moving to the top)
conda config --add channels defaults

## Add the channel "bioconda" to the top of th channel list
conda config --add channels bioconda

## Add a new value to channels so conda looks for packages in this location
conda config --add channels conda-forge

## Install conda packages to be used in the downstream analyses
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
conda install vcflib


##QUALITY_ASSESSMENT##

## Change to untrimmed_fastq directory
cd ~/ngs_course/dnaseq/data/untrimmed_fastq

## List the file details in a given directory
ls -lart

## fastq.qz is a compressed format so uncompress R1 and R2 fastq.qz files
zcat NGS0001.R1.fastq.qz > NGS0001.R1.fastq
zcat NGS0001.R2.fastq.qz > NGS0001.R2.fastq
fastqc NGS0001.R1.fastq NGS0001.R2.fastq

## Create fastqc_trimmed_reads directory for the results
mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads


## Move the fastqc files into the newly created directory
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/

##TRIMMOMATIC##

## Change directories to the untrimmed fastq data
cd ~/ngs_course/dnaseq/data/untrimmed_fastq

## Run trimmomatic command line for the paired-end fastq files
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq /home/ubuntu/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq \
  -baseout /home/ubuntu/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R \
 ILLUMINACLIP:/home/ubuntu//miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50


##QUALITY ASSESSMENT OF PAIRED TRIMMED SEQUENCING DATA##

## Perform basic quality assessment of paired trimmed sequencing data by fastqc
fastqc ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P NGS0001.R1.fastq_trimmed_R_2P


##ALIGNMENT##

## Create reference folder for the reference genome and its index files
mkdir -p ~/ngs_course/dnaseq/data/reference

## Move the reference genome into the reference folder
mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/

## Run bwa to generate the index files(take ~45 mins)
bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz

## List the contents of the reference folder
ls ~/ngs_course/dnaseq/data/reference

## Due to the limited space on the terminal(40Gigabytes), remove the untrimmed_fastq files
rm -r ~/ngs_course/dnaseq/data/untrimmed_fastq

## Create aligned_data folder 
mkdir ~/ngs_course/dnaseq/data/aligned_data

## Run bwa mem with the read group information for the alignment
bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_1P ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001.R1.fastq_trimmed_R_2P > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam

## Change directories into the aligned_data folder
cd ~/ngs_course/dnaseq/data/aligned_data

## Convert the SAM file into BAM file format, sort it and generate an index using Samtools
samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam

## List the contents of the aligned_data to check the generated files
cd ~/ngs_course/dnaseq/data/aligned_data

## Remove SAM files (as it is too large and not needed anymore)
cd ~/ngs_course/dnaseq/data/aligned_data
rm NGS0001.sam

##POST_ALIGNMENT_QUALITY_CONTROL_&_FILTERING##

## Use Picard to perform duplicate marking
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam

## Quality filter the duplicate marked BAM reads according to mapping quality and bitwise flags using samtools
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam

## Use Samtools and Picard to generate standard alignment statistics( flagstats, idxstats)
samtools flagstat NGS0001_sorted.bam
samtools idxstats NGS0001_sorted.bam

##VARIANT_CALLING_&_FILTERING##

## Uncompress the gizip reference file
zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa

## Index the reference with samtools faidx
samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa

## Call variants with Freebayes
freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/NGS0001.vcf

## Compress the generated variant call file (VCF)
bgzip ~/ngs_course/dnaseq/results/NGS0001.vcf

## Index the VCF with tabix
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

## Filter the VCF to remove the bad/poor calls
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf

## Change from current directory to aligned_data
cd ~/ngs_course/dnaseq/data/aligned_data

## Filter the VCF file for the regions in the annotation.bed file
bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf -b ../annotation.bed > ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

## Compress the generated filtered VCF dataset
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf

## Index the filtered vcf.gz
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz

## Remove trimmed_fastq file to create more space in the terminal
rm -r ~/ngs_course/dnaseq/data/trimmed_fastq

## Go to home directory
cd ~/

##VARIANT_ANNOTATION_&_PRIORITIZATION##

## After uploading the annovar.latest.tar.gz file, unpack it and set annovar up
tar -zxvf annovar.latest.tar.gz

## Go to annovar directory
cd annovar

## Download annovar databases it uses for the annotation
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/

## Convert VCF to Annovar input format
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/NGS0001_filtered_R.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_R.avinput

## Run Annovar table function to create csv format file, download it via FileZilla
./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_R.avinput humandb/ -buildver hg19 -out ~/ngs_course/dnaseq/results/NGS0001 -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

