#!/bin/bash
#SBATCH --account=naiss2024-22-540
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=assig7_bioinfo
#SBATCH --error=assig7_%J.err
#SBATCH --output=assig7_%J.out
#SBATCH --mail-type=all
#SBATCH --mail-user=paula.camargo.romera@ki.se

module load sra-tools
module load seqkit
module load fastqc
module load flash2
module load bowtie2

#download accession numbers
sqlite3 -batch -noheader -csv /proj/applied_bioinformatics/common_data/sample_collab.db "select run_accession from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='x_pauca';" > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt

#download sequencing data
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif xargs -a /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt -n 1 -I {} fastq-dump -X 10 --readids --gzip --outdir /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/ --disable-multithreading --split-3 {}

#raw FASTQ seqkit
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit stats --threads 1 /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/* > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_stats.txt

#check for de-replicated
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit rmdup -n -d -i --threads 1 /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/* > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_dereplicated.txt
#output: [INFO] 0 duplicated records removed

#check for trimmed
#whole adapters
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit --threads 1 grep -s -i -p 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA' /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/* > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_adaptread1.txt
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit --threads 1 grep -s -i -p 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT' /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/* > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_adaptread2.txt
#short adapters
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit --threads 1 grep -s -i -p 'AGATCGGAAGAGCACA' /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/* > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_adaptread1s.txt
singularity exec /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit --threads 1 grep -s -i -p 'AGATCGGAAGAGCGTC' /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/data/sra_fastq/* > /proj/applied_bioinformatics/users/x_pauca/Assig7/MedBioinfo_assig7/analyses/x_pauca_adaptread2s.txt
#no output for any, so the adapters seems to be already trimmed

#quality control fastQC
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif xargs -a /mnt/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt -n 1 -I {} fastqc --outdir /mnt/MedBioinfo_assig7/analyses/fastqc --threads 2 --noextract /mnt/MedBioinfo_assig7/data/sra_fastq/{}_1.fastq.gz /mnt/MedBioinfo_assig7/data/sra_fastq/{}_2.fastq.gz
#Per base sequence quality has high values all the time in the different samples
#The adapter content is flat and 0% in the samples, which means the adaptors were trimmed
#Per sequence quality scores indicates good overall read quality in all samples
#Overrepresented Sequences appears normally bad, but the adaptors has been check and dont appear there

#merge reads, flash2
#the meta.def for the container was a bit tricky, I had to add the path of flash
#also problems with file output from tee with the mount container, so I executed the command inside the container
#I had to increase the max overlap of flash to 150 because of a warning, and now it is all combined
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif xargs -a /mnt/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt -I{} /usr/local/bin/FLASH-1.2.11/flash --threads=2 -z -M 150 --output-directory=/mnt/MedBioinfo_assig7/data/merged_pairs --output-prefix={}.flash /mnt/MedBioinfo_assig7/data/sra_fastq/{}_1.fastq.gz /mnt/MedBioinfo_assig7/data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a /mnt/MedBioinfo_assig7/analyses/x_pauca_flash.log

#checking the histograms, in general, the histogram suggests that the DNA library insert sizes are not uniform but rather display variability

#compare how many base pairs you had in your initial unmerged reads, versus how many you have left after merging: have you lost information, or was it redundant information?
# Count initial base pairs
initial_bp=$(singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit stats /mnt/MedBioinfo_assig7/data/sra_fastq/*_1.fastq.gz /mnt/MedBioinfo_assig7/data/sra_fastq/*_2.fastq.gz | awk 'NR>1 {sum += $6} END {print sum}')
#1376
# Count merged base pairs
merged_bp=$(singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif seqkit stats /mnt/MedBioinfo_assig7/data/merged_pairs/*.flash.extendedFrags.fastq.gz | awk 'NR>1 {sum += $6} END {print sum}')
#688
# Calculate difference
difference=$((initial_bp - merged_bp))
echo "Difference: $difference" #1376-688=688
#The result is the half of the initial, therefore, there is a loss of information


#read mapping PhiX, bowtie2
#create index
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif bowtie2-build -f /mnt/MedBioinfo_assig7/data/reference_seqs/PhiX_NC_001422.fna /mnt/MedBioinfo_assig7/data/bowtie2_DBs/PhiX_bowtie2_DB

#align
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif xargs -a /mnt/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt -I{} bowtie2 -x /mnt/MedBioinfo_assig7/data/bowtie2_DBs/PhiX_bowtie2_DB -U /mnt/MedBioinfo_assig7/data/merged_pairs/{}.flash.extendedFrags.fastq.gz -S /mnt/MedBioinfo_assig7/analyses/bowtie/x_pauca_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee /mnt/MedBioinfo_assig7/analyses/bowtie/x_pauca_bowtie_merged2PhiX.log
#I did not have any hit alignment, so it seems not to contain PhiX

#align with SARS-CoV-2
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif xargs -a /mnt/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt -I{} bowtie2 -x /mnt/MedBioinfo_assig7/data/bowtie2_DBs/SARS-CoV-2_bowtie2_DB -U /mnt/MedBioinfo_assig7/data/merged_pairs/{}.flash.extendedFrags.fastq.gz -S /mnt/MedBioinfo_assig7/analyses/bowtie/x_pauca_merged2SARS-CoV-2.sam --threads 8 --no-unal 2>&1 | tee /mnt/MedBioinfo_assig7/analyses/bowtie/x_pauca_bowtie_merged2SARS-CoV-2.log
#Again, no hit aligments with SARS-CoV-2

#new bacterial: Fusobacterium_nucleatum_AE009951
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif xargs -a /mnt/MedBioinfo_assig7/analyses/x_pauca_run_accessions.txt -I{} bowtie2 -x /mnt/MedBioinfo_assig7/data/bowtie2_DBs/Fusobacterium_nucleatum_bowtie2_DB -U /mnt/MedBioinfo_assig7/data/merged_pairs/{}.flash.extendedFrags.fastq.gz -S /mnt/MedBioinfo_assig7/analyses/bowtie/x_pauca_merged2Fuso.sam --threads 8 --no-unal 2>&1 | tee /mnt/MedBioinfo_assig7/analyses/bowtie/x_pauca_bowtie_merged2Fuso.log
#ALIGN!!!
#10 reads; of these:
#  10 (100.00%) were unpaired; of these:
#    8 (80.00%) aligned 0 times
#    0 (0.00%) aligned exactly 1 time
#    2 (20.00%) aligned >1 times
#20.00% overall alignment rate

#combine quality contorl, multiqc
singularity exec -B /proj/applied_bioinformatics/users/x_pauca/Assig7:/mnt /proj/applied_bioinformatics/users/x_pauca/Assig7/meta.sif multiqc --force --title "x_pauca sample sub-set" /mnt/MedBioinfo_assig7/data/merged_pairs/ /mnt/MedBioinfo_assig7/analyses/fastqc/ /mnt/MedBioinfo_assig7/analyses/x_pauca_flash.log /mnt/MedBioinfo_assig7/analyses/bowtie/


