# Project details

## 1- Data Download

Download ~5M Human RNA-Seq fragements from the SRA.
Link to the suggested data: https://www.ncbi.nlm.nih.gov/sra/?term=SRR8797509

## 2- Prepare the data

- Construct 5 Samples, each sample contains two parts.
    - **Part 1:** 1M reads from the main reads file. (split the main file into 5 parts)
    - **Part 2:** Shuffle the main reads file, and take random 1M Reads. (split the shuffled file into 5 parts)
    - **Example:** Sample 1 will be divided into S1_1 & S1_2,  unshuffled and shuffled, respectively.

## 3- FASTQ Quality Control

- For **Sample 1** only, use FASTQC to report the difference between S1_1 and S1_2

## 4- Trimming

- For **All Samples** , Apply:
    - Mild Trimming for SX_1. {shuffled}
    - Aggressive Trimming for SX_2. {unshuffled}

## 5- Alignment

- Align all the samples (1:5) using **BWA** and **Hisat** against the human reference file.
	- **BWA** for SX_1 and **HISAT** for SX_2
- Export some useful statistics report for **each sample indvidually**.

## 6- Assembly

- Apply reference-based trasncriptome assembly using stringTie.
    -  **Step1.** For the 5 samples **unshuffled**.
    -  **Step2.** For the 5 samples **shuffled**.

## 7- Using GTF-Compare to Compare the Generated Annotation Files to a Reference Annotation.

## 8- Apply Differential Expression.
