#Oral microbiome profiling in smokers with and without head and neck cancer reveals variations between health and disease.
#The datasets generated for the current study have been deposited in the NCBI BioProject database under project number **PRJNA668810** 

#-- 16S rRNA (v4 region) data Analysis 

#Step1: Remove primers

for i in *_R1_001.fastq.gz;
do
  SAMPLE=$(echo ${i} | sed "s/_R1_\001\.fastq.gz//")
  echo ${SAMPLE}_R1_001.fastq.gz ${SAMPLE}_R2_001.fastq.gz
  cutadapt -a GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC -A GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
    --discard-untrimmed -o /home/gomeza/sharm646/Raw/Primer_Removed/${SAMPLE}_R1_001_primerRemoved.fastq \
    -p /home/gomeza/sharm646/Raw/Primer_Removed/${SAMPLE}_R2_001_primerRemoved.fastq  ${SAMPLE}_R1_001.fastq.gz  ${SAMPLE}_R2_001.fastq.gz

done

#Step2----Remove Low quality reads
cd Primer_Removed
module load fastx_toolkit
fastq_quality_filter -h
for i in *fastq; do fastq_quality_filter -q 30 -p 50 -i "$i" -o "$i".QC.gz -z; done;
mkdir QC_30
mv *gz QC_30

#Step3---- Make pairs again
sh repair.sh
~/Softwares/bbmap/repair.sh in1=IS013_S198_R1_001_primerRemoved.fastq.QC.gz in2=IS013_S198_R2_001_primerRemoved.fastq.
QC.gz out1=IS013_S198_R1_001_filtered.fastq.gz out2=IS013_S198_R2_001_filtered.fastq.gz outsingle=IS013_S198_Singleton
s.gz

#Step4--- Check the quality of your sequernces
cat *R1* > forward.fastq
cat *R2* > reverse.fastq
mkdir quality_check
moudle load fastqc
fastqc forward.fastq reverse.fastq -o quality_check

#Step5----- ASV picking and taxonomic assignment using Qiime2
#-- Prepare Manifest file
sh Sharma-qiime2.sh

#Step6--- Statistical Analysis- Scripts to generate figures
Main figure- follow Sharma-MainFigure.R
Supplementary figures - follow Sharma_SuppFigure.R

#---- End
