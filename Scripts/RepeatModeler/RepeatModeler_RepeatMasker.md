# RepeatModeler_RepeatMasker
***
**Purpose:** Install and run RepeatModeler followed by RepeatMasker on XINB3 and IINB1 genomes to identify repeat elements
*follow instructions on page: http://www.repeatmasker.org/RepeatModeler/*
**Written by:** Virginie Ricci and Maridel Fredericksen
**Updated:** 17.01.20


### Install prerequisites
***

- perl5.16.3
- RepeatMasker-open-4.1.0
- trf409
- rmblast-2.10.0
- RepeatScout-1.0.6
- RepeatModeler-2.0.1




##### *also install programs needed for LTR structural search pipeline*
***

- LtrHarvest (genometools-1.5.9)
 - cdhit-4.8.1
- ncbi-blast-2.10.0+-src
- hmmer-3.3
- LTR_retriever-2.8
- mafft-7.453-with-extensions
- NINJA-0.95-cluster_only
 

### Run RepeatModeler

***
To run the following script, use these commands:
```bash
chmod +x ./RepeatModeler.sh
 ./RepeatModeler.sh
 ```
#### Script: *RepeatModeler.sh*
```bash
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
```


### RepeatMasker
***

##### *combine repeat libraries:*
- note that RepBase library already has some annotated repeat families for *Daphnia pulex*
- also RepeatModeler may not detect older ancestral repeats that could be present in *Daphnia magna* genome

first extract the *Daphnia* sequences from RepeatMasker Database 
```bash
cd RepeatMasker/
perl ./util/queryRepeatDatabase.pl -species "daphnia" > RepeatMaskerLib_daphnia.fasta
```
now combine with my consensi.fa.classified files
```bash
cat /home/user/RepeatMasker/RepeatMaskerLib_daphnia.fasta /home/user/RepeatModeler_output_XINB3/RM_20190.WedJan291042012020/consensi.fa.classified > combined_repeat_libs_XINB3.fasta
cat /home/user/RepeatMasker/RepeatMaskerLib_daphnia.fasta /home/user/RepeatModeler_output_IINB1/RM_20244.WedJan291042022020/consensi.fa.classified > combined_repeat_libs_IINB1.fasta
```

##### *create output folders:*
#
```bash
cd
mkdir RepeatMasker_output_XINB3
mkdir RepeatMasker_output_IINB1
```

##### *run RepeatMasker with the following specifications:*

 -pa 5: run on 5 cores
 -s : slow search (Slow search; 0-5% more sensitive, 2-3 times slower than default)
 -lib: specify library I just created
 -dir: repeatmasker output folder that I created before 
```bash
RepeatMasker -pa 5 -s -lib combined_repeat_libs_XINB3.fasta -dir /home/user/RepeatMasker_output_XINB3/ -e ncbi  /home/user/XINB3/FI-XINB3_3.0_ref_genome.fasta 
RepeatMasker -pa 5 -s -lib combined_repeat_libs_IINB1.fasta -dir /home/user/RepeatMasker_output_IINB1/ -e ncbi  /home/user/IINB1/DE-IINB1_draft_genome.fasta
```

##### *interpret output:*
- use file genome.fasta.out, which is a table of where each element was found

cut it down to just the scaffold of interest
```bash
grep 000011F FI-XINB3_3.0_ref_genome.fasta.out > XINB3_11F.fasta.out
```

cut that down to just F locus
```bash
awk '{if ($6 >= 2329948) {print} }' XINB3_11F.fasta.out > XINB3_11F_Flocusstart.fasta.out
awk '{if ($7 <= 2358729) {print} }' XINB3_11F_Flocusstart.fasta.out > XINB3_11F_Flocus.fasta.out
```
repeat for IINB1
```bash
grep -w 129 DE-IINB1_draft_genome.fasta.out > IINB1_129.fasta.out
awk '{if ($6 >= 715834) {print} }' IINB1_129.fasta.out > IINB1_129_Flocusstart.fasta.out
awk '{if ($7 <= 762881) {print} }' IINB1_129_Flocusstart.fasta.out > IINB1_129_Flocus.fasta.out
```