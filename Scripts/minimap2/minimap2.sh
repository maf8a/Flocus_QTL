#!/bin/bash

##############################################################
#This bash script uses minimap2 (version 2.17-r941) to map and align
# de novo-assembled transcripts from each consensus transcriptome
# (from Daphnia magna QTL parent clones Xinb3 and Iinb1) to its
# respecitve whole genome.
#Written by Maridel Fredericksen
#Updated 20.07.2022
##############################################################


##############################################################
#For Xinb3
##############################################################

#align consensus transcriptome to Xinb3 genome
minimap2 -ax splice:hq -L -uf FI-XINB3_3.0_ref_genome.fasta stressflea_Xinb3_merged_transcriptome.fasta.gz > minimap2_spadeVSXinb3.sam

#convert sam to bam
samtools view -b -S -o minimap2_spadeVSXinb3.bam minimap2_spadeVSXinb3.sam

#sort bam file
samtools sort minimap2_spadeVSXinb3.bam -o minimap2_spadeVSXinb3.sorted.bam

#index sorted bam file
samtools index minimap2_spadeVSXinb3.sorted.bam

#apply filters (MapQ score 60, remove secondary alignments: bit flag 0x100)
samtools view -bq 60 -F 0x100 minimap2_spadeVSXinb3.sorted.bam > minimap2_spadeVSXinb3.sorted.filtered60.bam
#index new bam
samtools index minimap2_spadeVSXinb3.sorted.filtered60.bam

#only F locus:
samtools view -b minimap2_spadeVSXinb3.sorted.filtered60.bam "000011F:2329948-2358729"> minimap2_spadeVSXinb3.sorted.filtered60.Flocus.final.bam

#only F locus as fasta:
samtools fasta minimap2_spadeVSXinb3.sorted.filtered60.Flocus.final.bam > minimap2_spadeVSXinb3.sorted.filtered60.Flocus.final.fasta


##############################################################
#For Iinb1
##############################################################

#align consensus transcriptome to Iinb1 genome
minimap2 -ax splice:hq -L -uf DE-IINB1_draft_genome.fasta iinb1_transcriptome_assembly_31012020.fasta.gz > new_minimap2_spadeVSIinb1.sam

#convert sam to bam
samtools view -b -S -o new_minimap2_spadeVSIinb1.bam new_minimap2_spadeVSIinb1.sam

#sort bam file
samtools sort new_minimap2_spadeVSIinb1.bam -o new_minimap2_spadeVSIinb1.sorted.bam

#index sorted bam file
samtools index new_minimap2_spadeVSIinb1.sorted.bam

#apply filters (MapQ score 60, remove secondary alignments: bit flag 0x100)
samtools view -bq 60 -F 0x100 new_minimap2_spadeVSIinb1.sorted.bam > new_minimap2_spadeVSIinb1.sorted.filtered60.bam
#index new bam
samtools index new_minimap2_spadeVSIinb1.sorted.filtered60.bam

#only F locus:
samtools view -b new_minimap2_spadeVSIinb1.sorted.filtered60.bam "129:715834-762881"> new_minimap2_spadeVSIinb1.sorted.filtered60.Flocus.final.bam

#only F locus as fasta (default option removes supplementary alignments)
samtools fasta new_minimap2_spadeVSIinb1.sorted.filtered60.Flocus.final.bam > new_minimap2_spadeVSIinb1.sorted.filtered60.Flocus.final.fasta