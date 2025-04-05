# Transcriptomics
A pipeline for de novo transcriptome assembly and differential gene expression analysis using FastQC, Trimmomatic, Trinity, BUSCO, CD-HIT-EST, TransDecoder, and GFold.

# De Novo Transcriptome Assembly and Differential Expression Pipeline

This repository contains a reproducible pipeline for transcriptome assembly and differential gene expression analysis using widely adopted tools including FastQC, Trimmomatic, Trinity, BUSCO, CD-HIT-EST, TransDecoder, and GFold.

## 🔧 Tools and Steps

### 1. 🔍 Quality Check with FastQC
FastQC provides an overview of read quality, including per-base quality scores, GC content, sequence duplication levels, and potential adapter contamination. 

fastqc read*.fq.gz

### 2. ✂️ Adapter and Quality Trimming with Trimmomatic
Trimmomatic removes adapter sequences, trims low-quality bases from the start and end of reads, and uses a sliding window approach to ensure the remaining reads maintain a certain quality threshold. It also discards reads that fall below a minimum length.

java -jar trimmomatic-0.39.jar SE -phred33 reads.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50

### 3. 🧬 Transcriptome Assembly using Trinity
Trinity assembles short RNA-seq reads into full-length transcripts without the need for a reference genome. It uses a three-step process (Inchworm, Chrysalis, Butterfly) to reconstruct transcript sequences and handle alternative splicing.

Ensure Salmon is in your PATH:
export PATH=$PATH:/Path/To/salmon/bin/

Run Trinity:
Trinity --seqType fq --max_memory 10G --KMER-SIZE 25 --output Trinity_output --left S1_R1.fq.gz --right S2_R2.fq.gz --no_bowtie --SS_lib_type FR --CPU 4

### 4. ✅ Assembly Quality Check using BUSCO
Trinity assembles short RNA-seq reads into full-length transcripts without the need for a reference genome. It uses a three-step process (Inchworm, Chrysalis, Butterfly) to reconstruct transcript sequences and handle alternative splicing.

busco -i Trinity.fasta -l viridiplantae -f -o BP_Busco -m transcriptome

### 5. 🧹 Redundancy Removal using CD-HIT-EST
CD-HIT-EST clusters highly similar transcript sequences (e.g., >95% identity), retaining only representative ones. This helps reduce redundancy in the transcriptome assembly and produces a non-redundant set of unigenes.

cdhit-est -i Trinity.fasta -o unigenes -c 0.95 -n 8

### 6. 🧬 ORF Prediction with TransDecoder
TransDecoder identifies candidate coding regions within transcript sequences. It first finds long ORFs and then optionally evaluates their coding potential using homology evidence or sequence composition models.

TransDecoder.LongOrfs -t unigenes.fasta

### 7. 📊 Differential Expression Analysis using GFold
GFold (Generalized Fold Change) provides reliable ranking of gene expression differences, especially for datasets with few or no replicates. It generates read count files from aligned reads and calculates fold changes between conditions or samples.

Generate read counts:
gfold count -ann hg19Ref.gtf -tag sample1.sam -o sample1.read_cnt  
gfold count -ann hg19Ref.gtf -tag sample2.sam -o sample2.read_cnt  

Perform differential expression analysis:
gfold diff -s1 sample1 -s2 sample2 -suf .read_cnt -o sample2VSsample1.diff

## 📁 Output Structure
- `Trinity_output/` – Assembled transcripts  
- `BP_Busco/` – BUSCO completeness report  
- `unigenes.fasta` – Non-redundant transcripts  
- `sample*.read_cnt` – Read count files  
- `sample2VSsample1.diff` – Differential expression output  

## 📌 Requirements
- Java (for Trimmomatic)  
- Python & Biopython (for GFold)  
- Tools: FastQC, Trinity, Salmon, CD-HIT, BUSCO, TransDecoder, GFold  

## 📬 Citation
If you use this pipeline in your research, please consider citing the respective tools used.

## 🧠 License
This repository is released under the MIT License.
