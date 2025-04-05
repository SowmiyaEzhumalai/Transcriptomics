# Transcriptomics
A pipeline for de novo transcriptome assembly, differential gene expression analysis and its visualization using FastQC, Trimmomatic, Trinity, BUSCO, CD-HIT-EST, TransDecoder, GFold and RStudio.

# De Novo Transcriptome Assembly and Differential Expression Pipeline

This repository contains a reproducible pipeline for transcriptome assembly and differential gene expression analysis using widely adopted tools, including FastQC, Trimmomatic, Trinity, BUSCO, CD-HIT-EST, TransDecoder, and GFold.

## 🔧 Workflow:

### 1. 🔍 Quality Check with FastQC
FastQC provides an overview of read quality, including per-base quality scores, GC content, sequence duplication levels, and potential adapter contamination. 

_fastqc read*.fq.gz_


### 2. ✂️ Adapter and Quality Trimming with Trimmomatic
Trimmomatic removes adapter sequences, trims low-quality bases from the start and end of reads, and uses a sliding window approach to ensure the remaining reads maintain a certain quality threshold. It also discards reads that fall below a minimum length.

_java -jar trimmomatic-0.39.jar SE -phred33 reads.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50_


### 3. 🧬 Transcriptome Assembly using Trinity
Trinity assembles short RNA-seq reads into full-length transcripts without the need for a reference genome. It uses a three-step process (Inchworm, Chrysalis, Butterfly) to reconstruct transcript sequences and handle alternative splicing.

#### Ensure Salmon is in your PATH:
_export PATH=$PATH:/Path/To/salmon/bin/_

#### Run Trinity:
_Trinity --seqType fq --max_memory 10G --KMER-SIZE 25 --output Trinity_output --left S1_R1.fq.gz --right S2_R2.fq.gz --no_bowtie --SS_lib_type FR --CPU 4_


### 4. ✅ Assembly Quality Check using BUSCO
BUSCO evaluates the completeness of the assembled transcriptome by checking for the presence of expected single-copy orthologous genes that are highly conserved across specific lineages.

_busco -i Trinity.fasta -l viridiplantae -f -o BP_Busco -m transcriptome_


### 5. 🧹 Redundancy Removal using CD-HIT-EST
CD-HIT-EST clusters highly similar transcript sequences (e.g., >95% identity), retaining only representative ones. This helps reduce redundancy in the transcriptome assembly and produces a non-redundant set of unigenes.

_cdhit-est -i Trinity.fasta -o unigenes -c 0.95 -n 8_


### 6. 🧬 ORF Prediction with TransDecoder
TransDecoder identifies candidate coding regions within transcript sequences. It first finds long ORFs and then optionally evaluates their coding potential using homology evidence or sequence composition models.

_TransDecoder.LongOrfs -t unigenes.fasta_


### 7.  📊 Differential Expression Analysis using GFold
GFold (Generalized Fold Change) provides reliable ranking of gene expression differences, especially for datasets with few or no replicates. It generates read count files from aligned reads and calculates fold changes between conditions or samples.

#### Generate read counts:
gfold count -ann hg19Ref.gtf -tag sample1.sam -o sample1.read_cnt  
gfold count -ann hg19Ref.gtf -tag sample2.sam -o sample2.read_cnt  

#### Perform differential expression analysis:
gfold diff -s1 sample1 -s2 sample2 -suf .read_cnt -o sample2VSsample1.diff 


### 8. Functional Annotation using OmicsBox
The assembled transcripts can be annotated using functional annotation software such as OmicsBox in which the predicted unigenes was taken to BLASTX search against the specific lineage database followed by GO mapping and GO annotation. 


### 9.Pathway Identification using KAAS server
Pathway Mapping can be performed for whole nucleotide sequences with a   corresponding representative set of organisms. The downstream analysis can be narrowed down to any selected metabolic pathways or secondary metabolite biosynthesis pathways according to the research objective.


### 10. Visualization of Differential Expression using R
The expression of genes or transcripts involved in the selected pathway have been taken to generate visualization. It has been done by 'pheatmap', an R package. 'pheatmap' is a visualization tool for creating complex and customizable heatmaps, often used in transcriptomics or gene expression studies.

#### Load necessary library
_library(heatmap)_
#### Read RPKM data from CSV
_rpkm_data <- read.table("/path/to/the/rpkm.csv", sep=',', header=TRUE, row.names = 1, col.names = c("col_1", "col_2"))_
#### Convert the data into matrix
_data<-as.matrix(rpkm_data)_
#### Plot heatmap
_pheatmap(data, scale = "none", show_rownames = F, color = colorRampPalette(c("navy", "green","pink" , "white","orange","yellow", "firebrick3" ))(3000), cellwidth=100, cellheight=5, fontsize_row=4, border_color = "black")_


## 📌 Requirements
- Java (for Trimmomatic)  
- Python & Biopython (for GFold)
- RStudio
- Tools: FastQC, Trinity, Salmon, CD-HIT, BUSCO, TransDecoder, GFold  

## 📬 Citation
1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120.
2. Brown, J., Pirrung, M., & McCue, L. A. (2017). FQC Dashboard: integrates FastQC results into a web-based, interactive, and extensible FASTQ quality control tool. Bioinformatics, 33(19), 3137-3139.
3. Feng, J., Meyer, C. A., Wang, Q., Liu, J. S., Shirley Liu, X., & Zhang, Y. (2012). GFOLD: a generalized fold change for ranking differentially expressed genes from RNA-seq data. Bioinformatics, 28(21), 2782-2788.
4. Haas, B. J., Papanicolaou, A., Yassour, M., Grabherr, M., Blood, P. D., Bowden, J., ... & Regev, A. (2013). De novo transcript sequence reconstruction from RNA-seq using the Trinity platform for reference generation and analysis. Nature protocols, 8(8), 1494-1512.
5. Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature biotechnology, 37(8), 907-915.
6. Li, W., & Godzik, A. (2006). Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658-1659.
7. Simão, F. A., Waterhouse, R. M., Ioannidis, P., Kriventseva, E. V., & Zdobnov, E. M. (2015). BUSCO: assessing genome assembly and annotation completeness with single-copy orthologs. Bioinformatics, 31(19), 3210-3212.

