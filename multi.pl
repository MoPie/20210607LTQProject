#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

# take the list of commands as input and run through qsub
############# qsub parameter ################
$PM =
'#!/bin/sh
#$ -cwd 
#$ -o ./
#$ -M yixinguo.19@intl.zju.edu.cn
#$ -V
#$ -S /bin/bash
#$ -l h_data=6G,h_rt=1:00:00
#$ -e ./
# -pe shared 4
';

############# BSMAP parameter ########################
#$refG ='/u/home/w/wanluliu/reference/mm9.fa';
#MM9 microgenome '/u/home/w/wanluliu/reference/mm9_repMicroGenome/repRegion_microGenome.fa';
#single_end: using 8 threads [-p 8], and allow up to 1 multiple hits [-w 1], map to all four possible strands [-n 1], allow 2 mismatches [-v 2], 
$BSMAP_PM = '-p 4 -w 1 -n 1 -v 2';

#################### open folder contains all fastq files#################
my $USAGE = "\nUSAGE: multiple_BSMAP.pl 
-seqdir directory_with_fastq_files
";
my $options = {};
GetOptions($options, "-seqdir=s"); #, "-out=s" 
die $USAGE unless defined ($options->{seqdir});

############################# Grobal Variables #############################
my $seqdir = $options->{seqdir};

opendir (DIR, $seqdir) or die "Couldn't open $seqdir: $!\n";

my $command;
while (defined(my $seqfile = readdir(DIR))) {
unless ($seqfile =~ /(\S+)_lambdaDNA_srt.bam$/) {next;}
#unless ($seqfile =~ /(\S+)\.fastq$/) {next;}
$lib_name = $1;
# $r=$2;
# $chr=$2;
#bsmap -a /u/home/j/jxzhai/MCDB_folder/ALLDATA/BS-seq/Raw/BC5_ACCTCA.fastq.gz -d /u/home/j/jxzhai/MCDB_folder/ALLDATA/Genomic/TAIR10/GenomicSEQ/Ath_ChrAll.fa -o BC5.bam -p 8 -w 1 -n 0 -v 2
$command =
"
refDir='/u/home/w/wanluliu/wanluliu/refData/lambda.fa'
python /u/home/w/wanluliu/tools/bsmap/bsmap-2.74/methratio.py -d ${refDir} --out=${lib_name}\_lambdaDNA_methratio.txt -r -u -z $seqfile -s /u/home/w/wanluliu/tools/samtools-0.1.19

";
print "$command\n";
open(SH, ">$lib_name\_tmp_BSMAPstep_lambdaDNA.sh");
print SH "$PM";
print SH "$command\n";
# system "qsub $lib_name$chr\_tmp_BSMAPstep3.sh\n";
close SH;
#system "rm $lib_name$chr\_tmp_BSMAPstep1.sh";

}
#python /u/home/w/wanluliu/tools/bsmap/bsmap-2.74/methratio.py -d ~/project-mcdb/reference/hg19_ref/$chr.fa --out=$lib_name$chr\_methratio\.txt -r -u -z $lib_name$chr\.bam -s /u/home/w/wanluliu/tools/samtools/samtools-0.1.19
#python /u/home/w/wanluliu/tools/bsmap/bsmap-2.74/methratio.py -d ~/project-mcdb/reference/hg19_ref/hg19.fa --out=$lib_name\_methratio\.txt -r -u -z $lib_name\.bam -s /u/home/w/wanluliu/tools/samtools/samtools-0.1.19
#bsmap -a $seqdir$seqfile -d ~/project-mcdb/reference/hg19_ref/hg19.fa -o $lib_name\.bam $BSMAP_PM
#bsmap -a $seqdir$seqfile -d /u/home/w/wanluliu/project-mcdb/reference/lambda.fa -o $lib_name\_lambda.bam $BSMAP_PM
#perl ~/tools/script/extractType.pl --i $lib_name\.txt --type CG  --o ./CGmethratio/
#/u/home/w/wanluliu/tools/bamtools-master/bin/bamtools split -in $seqfile -reference
#python /u/home/w/wanluliu/tools/bsmap/bsmap-2.74/methratio.py -d ~/project-mcdb/reference/hg19_ref/$chr.fa --out=$lib_name$chr\_methratio\.txt -r -u -z $lib_name$chr\.bam -s /u/home/w/wanluliu/tools/samtools/samtools-0.1.19