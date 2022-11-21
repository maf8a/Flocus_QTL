#!/bin/bash

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE


bin_RMasker='/home/user/RepeatMasker/'
export PATH=$bin_RMasker:$PATH

bin_RMasker_util='/home/user/RepeatMasker/util/'
export PATH=$bin_RMasker_util:$PATH

bin_RModeler='/home/user/RepeatModeler-2.0.1/'
export PATH=$bin_RModeler:$PATH


Genome='XINB3'
cd
mkdir -p RepeatModeler_output_${Genome}
cd RepeatModeler_output_${Genome}

# BuildDatabase command
# create NCBI Blast database of the .fasta genome file
echo 'BuildDatabase for' $Genome
BuildDatabase -name ${Genome}_database /home/user/XINB3/FI-${Genome}_3.0_ref_genome.fasta -engine ncbi

# RepeatModeler command
# identify de-novo repeat families and create species-specific repeat library
echo 'RepeatModeler for' $Genome
nohup RepeatModeler -database ${Genome}_database -engine ncbi -pa 5 -LTRStruct >& ${Genome}_run.out & 


Genome='IINB1'
cd
mkdir -p RepeatModeler_output_${Genome}
cd RepeatModeler_output_${Genome}

# BuildDatabase command
# create NCBI Blast database of the .fasta genome file
echo 'BuildDatabase for' $Genome
BuildDatabase -name ${Genome}_database /home/user/IINB1/DE-${Genome}_draft_genome.fasta -engine ncbi

# RepeatModeler command
# identify de-novo repeat families and create species-specific repeat library
echo 'RepeatModeler for' $Genome
nohup RepeatModeler -database ${Genome}_database -engine ncbi -pa 5 -LTRStruct >& ${Genome}_run.out &


DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE