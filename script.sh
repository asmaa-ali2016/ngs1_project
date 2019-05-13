###########################################
#  - NGS1 Course - Assignment             #
#  - Bash Script                          #
#  - Copyright: Asmaa Ali                 #
#  - Nile University                      #
###########################################

#!/bin/bash
# file: assignment.sh


cd Assignment
wget https://wsend.net/8012c94d7b1acec3e3be0c1305b78b60/SRR8797509_X.zip #Paired-end data
#download the shuffled data https://uploadfiles.io/aydom

mv shuffled_SRR8797509_2.part_001.tar.xz SRR8797509_1.part_002.tar.xz
tar xvf shuffled_SRR8797509_1.part_002.tar.xz

mv shuffled_SRR8797509_2.part_002.tar.xz SRR8797509_2.part_002.tar.xz
tar xvf shuffled_SRR8797509_2.part_002.tar.xz

gzip UnShuffled/SRR8797509_1.part_001.fastq 
gzip UnShuffled/SRR8797509_2.part_001.fastq

gzip Shuffled/SRR8797509_1.part_002.fastq 
gzip Shuffled/SRR8797509_2.part_002.fastq 


for i in 1 2
do  
seqkit split UnShuffled/SRR8797509_${i}.part_001.fastq.gz -p 5 
seqkit split Shuffled/SRR8797509_${i}.part_002.fastq.gz -p 5
done


fastqc -t 1 -f fastq -noextract Shuffled/SRR8797509_1.part_002.fastq.gz.split/SRR8797509_1.part_002.part_001.fastq.gz 
fastqc -t 1 -f fastq -noextract UnShuffled/SRR8797509_1.part_001.fastq.gz.split/SRR8797509_1.part_001.part_001.fastq.gz 


mkdir shuffled_trimmed && mkdir unshuffled_trimmed
for i in {1..5}
do
f1="~/NU_Bioinformatics_Diploma/NGS2/Assignment/UnShuffled/SRR8797509_1.part_001.fastq.gz.split/SRR8797509_1.part_001.part_00${i}.fastq.gz"
f2="~/NU_Bioinformatics_Diploma/NGS2/Assignment/UnShuffled/SRR8797509_2.part_001.fastq.gz.split/SRR8797509_2.part_001.part_00${i}.fastq.gz" 
newf1="~/NU_Bioinformatics_Diploma/NGS2/Assignment/unshuffled_trimmed/SRR8797509_1.part_001.part_00${i}.pe.trim.fq.gz" 
newf2="~/NU_Bioinformatics_Diploma/NGS2/Assignment/unshuffled_trimmed/SRR8797509_2.part_001.part_00${i}.pe.trim.fq.gz" 
newf1U="~/NU_Bioinformatics_Diploma/NGS2/Assignment/unshuffled_trimmed/SRR8797509_1.part_001.part_00${i}.se.trim.fq.gz" 
newf2U="~/NU_Bioinformatics_Diploma/NGS2/Assignment/unshuffled_trimmed/SRR8797509_2.part_001.part_00${i}.se.trim.fq.gz" 
adap="~/anaconda3/share/trimmomatic-0.39-1/adapters" 
 
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:15 MINLEN:36
done

for i in {1..5}
do
f1="~/NU_Bioinformatics_Diploma/NGS2/Assignment/Shuffled/SRR8797509_1.part_002.fastq.gz.split/SRR8797509_1.part_002.part_00${i}.fastq.gz"
f2="~/NU_Bioinformatics_Diploma/NGS2/Assignment/Shuffled/SRR8797509_2.part_002.fastq.gz.split/SRR8797509_2.part_002.part_00${i}.fastq.gz" 
newf1="~/NU_Bioinformatics_Diploma/NGS2/Assignment/shuffled_trimmed/SRR8797509_1.part_002.part_00${i}.pe.trim.fq.gz" 
newf2="~/NU_Bioinformatics_Diploma/NGS2/Assignment/shuffled_trimmed/SRR8797509_2.part_002.part_00${i}.pe.trim.fq.gz" 
newf1U="~/NU_Bioinformatics_Diploma/NGS2/Assignment/shuffled_trimmed/SRR8797509_1.part_002.part_00${i}.se.trim.fq.gz" 
newf2U="~/NU_Bioinformatics_Diploma/NGS2/Assignment/shuffled_trimmed/SRR8797509_2.part_002.part_00${i}.se.trim.fq.gz" 
adap="~/anaconda3/share/trimmomatic-0.39-1/adapters" 
 
trimmomatic PE -threads 1 -phred33 -trimlog trimLogFile -summary statsSummaryFile  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
ILLUMINACLIP:$adap/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:2:30 MINLEN:36
done


mkdir bwa_align && cd bwa_align
mkdir bwa_index && cd bwa_index
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.pc_transcripts.fa.gz 
gunzip gencode.v29.pc_transcripts.fa.gz 
bwa index -a bwtsw gencode.v29.pc_transcripts.fa
cd ..

for i in {1..5} 
do 
R1="~/NU_Bioinformatics_Diploma/NGS2/Assignment/unshuffled_trimmed/SRR8797509_1.part_001.part_00${i}.pe.trim.fq.gz" 
R2="~/NU_Bioinformatics_Diploma/NGS2/Assignment/unshuffled_trimmed/SRR8797509_2.part_001.part_00${i}.pe.trim.fq.gz" 

/usr/bin/time -v bwa mem bwa_index/gencode.v29.pc_transcripts.fa $R1 $R2 > unshuffled_trimmed.pe_00${i}.sam 
done 

cd ..
mkdir hisat_align && cd hisat_align 
mkdir hisat_index && cd hisat_index
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf 

hisat2_extract_splice_sites.py chr22_with_ERCC92.gtf > splicesites.tsv  
hisat2_extract_exons.py chr22_with_ERCC92.gtf > exons.tsv  
hisat2-build -p 1 --ss splicesites.tsv --exon exons.tsv chr22_with_ERCC92.fa chr22_with_ERCC92 

cd ..
for i in {1..5} 
do 
R1="~/NU_Bioinformatics_Diploma/NGS2/Assignment/shuffled_trimmed/SRR8797509_1.part_002.part_00${i}.pe.trim.fq.gz" 
R2="~/NU_Bioinformatics_Diploma/NGS2/Assignment/shuffled_trimmed/SRR8797509_2.part_002.part_00${i}.pe.trim.fq.gz" 

hisat2 -p 1 -x hisat_index/chr22_with_ERCC92 --dta --rna-strandness RF -1 $R1 -2 $R2 -S shuffled_trimmed.pe_00${i}.sam 
done


cd ..
for i in {1..5}
do
samtools flagstat bwa_align/unshuffled_trimmed.pe_00${i}.sam > bwa_align/unshuffled_trimmed.pe_00${i}.sam.stats
samtools view -bS bwa_align/unshuffled_trimmed.pe_00${i}.sam > bwa_align/unshuffled_trimmed.pe_00${i}.bam
samtools sort bwa_align/unshuffled_trimmed.pe_00${i}.bam -o bwa_align/unshuffled_trimmed.pe_00${i}.sorted.bam
done

for i in {1..5}
do
samtools flagstat hisat_align/shuffled_trimmed.pe_00${i}.sam > hisat_align/shuffled_trimmed.pe_00${i}.sam.stats
samtools view -bS hisat_align/shuffled_trimmed.pe_00${i}.sam > hisat_align/shuffled_trimmed.pe_00${i}.bam
samtools sort hisat_align/shuffled_trimmed.pe_00${i}.bam -o hisat_align/shuffled_trimmed.pe_00${i}.sorted.bam
done


mkdir gtfs
for i in {1..5}
do
stringtie hisat_align/shuffled_trimmed.pe_00${i}.sorted.bam --rf -l ref_free -o gtfs/ref_free_00${i}.gtf
cat gtfs/ref_free_00${i}.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l
stringtie hisat_align/shuffled_trimmed.pe_00${i}.sorted.bam --rf -l ref_sup -G hisat_align/hisat_index/chr22_with_ERCC92.gtf -o gtfs/ref_sup_00${i}.gtf 
cat gtfs/ref_sup_00${i}.gtf | grep -v "^@" | awk '$3=="transcript"' | wc -l
done


cd gtfs
for i in {1..5}
do
gffcompare -r ref_sup_00${i}.gtf ref_free_00${i}.gtf
done


mkdir diff_exp && cd diff_exp
conda install -y bioconductor-deseq
featureCounts -a ../hisat_align/hisat_index/chr22_with_ERCC92.gtf -g gene_name -o counts.txt  ../hisat_align/shuffled_trimmed.pe*.sorted.bam  ../bwa_align/unshuffled_trimmed.pe*.sorted.bam

cat counts.txt | cut -f 1,7-12 > simple_counts.txt
cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
cat results_deseq1.tsv | awk ' $8 < 0.05 { print $0 }' > filtered_results_deseq1.tsv
cat filtered_results_deseq1.tsv | Rscript draw-heatmap.r > hisat_output.pdf
