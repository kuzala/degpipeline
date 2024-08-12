

# RNA-seq Analysis Pipeline

This repository contains a Python script designed to automate the processing and analysis of bulk RNA-seq data. The pipeline begins with raw FASTQ files and performs quality control, read trimming, genome alignment, and differential gene expression analysis using DESeq2 tasks consecutively. Its modular and flexible design handles both paired-end and single-end RNA-seq data.

## Table of Contents
- [Purpose](#purpose)
- [Prerequisites](#prerequisites)
- [Script Overview](#script-overview)
- [Usage](#usage)
  - [Clone the Repository](#clone-the-repository)
  - [Naming Input Files](#naming-input-files)
  - [Directory, Parameter and Sample Condition Setup](#directory-parameter-and-sample-condition-setup)
  - [Run the Script](#run-the-script)
- [Output Files](#output-files)

## Purpose

This pipeline is designed to perform analysis of bulk RNA-seq data to find differentially expressed genes, covering the following steps:

1. **Quality Control**: Using FASTQC to assess the quality of raw reads.
2. **Read Trimming**: Using Trimmomatic to remove adapters and low-quality regions.
3. **Read Alignment**: Aligning reads to a reference genome using HiSAT2.
4. **SAM to BAM Conversion**: Converting SAM files to BAM, sorting, and indexing them using Samtools.
5. **Gene Count Matrix Generation**: Counting reads mapped to each gene using featureCounts.
6. **DESeq2 Analysis**: Performing differential gene expression analysis with DESeq2.
7. **Visualization**: Visualizing differential gene expression analysis results.

[Back to the top](#rna-seq-analysis-pipeline)

## Prerequisites

Before running this script, ensure that you have the following tools and libraries installed:

- Python (version 3.x)
- R (version 4.x)
- FASTQC
- Trimmomatic
- HiSAT2
- samtools
- subread (Includes featureCounts)
- DESeq2 (R package)
- ggplot2 (R package)
- readr (R package)
- Pandas (Python library, install it with `pip install pandas` if needed)
- rpy2 (Python library, install it with `pip install rpy2` if needed)

[Back to the top](#rna-seq-analysis-pipeline)

## Script Overview

### 1. Directory and Parameter Setup
The script starts by defining paths to input and output directories, to required files such as adapter files for Trimmomatic and gene annotation file, and to the provided R script. It also defines the Trimmomatic parameters. Users should modify these parts based on their system configuration and specific needs. Additionally, it also creates necessary output directory and subdirectories if they do not already exist.

### 2. Gathering FASTQ Files and Identifying Paired-End and Single-End Data
The script collects all FASTQ files in the specified input directory. It assumes the files follow a specific naming convention to identify whether the data is paired-end or single-end. If the data is paired-end, it matches the correct files with the correct samples. Naming convention will be explained in detail in the [Usage](#naming-input-files) section.

### 3. Quality Control with FASTQC
FASTQC is run on each FASTQ file to generate quality control reports, and results are saved in the designated output directory.

### 4. Trimming Adapters and Low-Quality Reads
The script uses Trimmomatic to trim adapters and low-quality regions. It automatically selects the appropriate mode based on whether the data is paired-end or single-end, and saves the trimmed reads in the output directory (paired and unpaired for paired-end data).

### 5. Aligning Reads to the Reference Genome
Trimmed reads (paired ones for paired-end data) are aligned to a reference genome using HiSAT2 by automatically selecting the appropriate mode based on the data type. It outputs a SAM file for each sample, and these SAM files are stored in the alignment directory.

### 6. Converting SAM Files to Sorted BAM Files
SAM files are converted to BAM format, sorted, and indexed using samtools. The script deletes intermediate files to save space.

### 7. Generating Gene Count Matrix
A gene count matrix is generated using sorted BAM files by featureCounts, which counts the number of reads mapping to each gene per each sample. Samples are alphabetically ordered while producing the matrix. The resulting matrix is saved in the output directory.

### 8. Preparing Data for DESeq2 Analysis
The gene count matrix is pre-processed to make it compatible with DESeq2. This step involves dropping unnecessary columns, properly renaming the sample columns, and saving the processed matrix in a new file.

### 9. Creating a Sample Sheet for DESeq2
A metadata file called a sample sheet is created. It lists each sample and its condition (e.g., treated vs. control). This part needs the user to modify the script based on each experiment. Points to pay attention to will be explained in detail in the [Usage](#sample-condition-definition) section. This sheet is required for DESeq2 analysis.

### 10. Running DESeq2 for Differential Gene Expression Analysis and Visualization
The script executes the provided R script to perform differential gene expression analysis based on the specified conditions in the sample sheet. The resulting standard DESeq2 table and generated figures (PCA plot, MA plot, and Volcano plot) are saved in the featureCounts’ output directory. The script automatically sets this directory as R’s working directory. Therefore, users typically don’t need to modify the given R script except for their specific needs.

[Back to the top](#rna-seq-analysis-pipeline) [Back to the Script Overview](#script-overview)

## Usage

### Clone the Repository
First, users need to clone the repository:
```bash
git clone https://github.com/kuzala/degpipeline.git
cd degpipeline
```
[Back to the top](#rna-seq-analysis-pipeline)

### Naming Input Files
The script needs original input FASTQ files to follow a specific naming convention given in the table below:
| Sequencing Type | Correct Naming       | Incorrect Naming       |
|:---------------:|----------------------|------------------------|
| **Single-End**  |                      |                        |
|                 | samplename1.fastq    | sample_name1.fastq      |
|                 | samplename2.fastq    | sample_name2.fastq     |
| **Paired-End**  |                      |                        |
| Forward         | samplename1_1.fastq    | sample_name1_1.fastq    |
| Reverse         | samplename1_2.fastq    | sample_name1_2.fastq    |
| Forward         | samplename2_1.fastq    | sample_name2_1.fastq    |
| Reverse         | samplename2_2.fastq    | sample_name2_2.fastq    |
| Forward         | samplename3_1.fastq    | sample_name_3_1.fastq    |
| Reverse         | samplename3_2.fastq    | sample_name_3_2.fastq    |
| Forward         | samplename4_1.fastq    | sample_name4_forward.fastq    |
| Reverse         | samplename4_2.fastq    | sample_name4_reverse.fastq    |


- Paired-end reads **MUST** have `_1.fastq` and `_2.fastq` in their filenames after sample names. Users should keep that in mind and stick with the examples given in the table. The script relies on this convention to identify paired-end reads.

- **DO NOT** use underscore (`_`) characters in sample names as the script uses underscores to parse filenames. Users need to be absolutely sure that filenames do not contain any additional underscores except those used to separate the sample name and the read number as just explained.

[Back to the top](#rna-seq-analysis-pipeline) [Back to the Usage](#usage)

### Directory, Parameter and Sample Condition Setup
Before running the script, users need to adjust several paths and parameters according to their system configuration and specific needs.

#### Paths to Modify:
- `input_dir`: Path to the directory containing input FASTQ files.
- `output_dir`: Path to the directory where users want to save the output files. The directory will be created automatically.
- `trimmomatic_path`: Full path to the Trimmomatic `.jar` file.
- `adapter_file`: Full path to the Trimmomatic adapter file.
- `hisat2_index`: Path to the HiSAT2 index files for the reference genome.
- `gtf_file`: Full path to the gene annotation file (GTF format).
- `deseq2_path`: Full path to the provided DESeq2 R script.

Within the output directory, the following subdirectories will be automatically created:
- `fastqc`: For FastQC quality control reports.
- `trimmed`: For trimmed FASTQ files.
- `aligned`: For aligned SAM and sorted BAM files after conversion.
- `counted`: For the gene count matrix, sample sheet, and DESeq2 files.

#### Parameters to Modify:
Users can modify the following Trimmomatic parameters according to their needs or leave them as default:
- `leading`: Minimum quality required to keep a base at the start (default: 3).
- `trailing`: Minimum quality required to keep a base at the end (default: 3).
- `slidingwindow`: Sliding window size and required quality (default: "4:15").
- `minlen`: Minimum length of reads to keep after trimming (default: 36).

[Back to the top](#rna-seq-analysis-pipeline) [Back to the Usage](#usage) [Back to the Directory, Parameter and Sample Condition Setup](#directory-parameter-and-sample-condition-setup)

### Sample Condition Definition
The script is designed to create a CSV file containing sample names and their conditions (e.g., treated vs. control) which is required for the DESeq2 analysis. FeatureCounts generates the gene count matrix by using alphabetically ordered sample names. Therefore, conditions must be given in a specific manner for correct sample assignment.

**Example:**

Let the samples and their conditions be as follows:
- SampleC is control
- SampleE is treated
- SampleA is control
- SampleD is treated
- SampleB is treated

Since FeatureCounts generates the count matrix by using alphabetically ordered sample names, the final pre-processed gene count matrix will look like the following (numbers are arbitrary raw counts):


| Geneid | SampleA | SampleB | SampleC | SampleD | SampleE |
|:------:|---------|---------|---------|---------|---------|
| Gene_1 | 14      | 10      | 22      | 18      | 9       |
| Gene_2 | 2       | 7       | 6       | 12      | 5       |

Therefore, the related part of the script should be modified as follows:

```python
conditions = [
    'control',  # Condition of SampleA
    'treated',  # Condition of SampleB
    'control',  # Condition of SampleC
    'treated',  # Condition of SampleD
    'treated',  # Condition of SampleE
]
```
[Back to the top](#rna-seq-analysis-pipeline) [Back to the Usage](#usage) [Back to the Directory, Parameter and Sample Condition Setup](#directory-parameter-and-sample-condition-setup)

### Run the Script

After modifying the script, users can simply run it using Python:

```bash
python rna_seq_pipeline.py
```


This will:
1. **Gather all FASTQ Files**: The script will select all the FASTQ files in the input directory.
2. **Perform Quality Control**: FASTQC is run on all FASTQ files to assess the quality.
3. **Trim Reads**: Trimmomatic trims adapters and low-quality sequences from the reads.
4. **Align Reads**: HISAT2 aligns the reads to the reference genome and produces a SAM file for each sample.
5. **Generate BAM Files**: Samtools converts and sorts the aligned SAM files into BAM files.
6. **Count Genes**: FeatureCounts creates a gene count matrix from the aligned reads.
7. **Prepare DESeq2 Input**: featureCounts output is pre-processed to a DESeq2-ready count matrix and a sample sheet is generated.
8.	**Run DESeq2**: The script calls the provided R script to perform differential expression analysis using DESeq2 and saves the resulting figures.

[Back to the top](#rna-seq-analysis-pipeline) [Back to the Usage](#usage)

## Output Files
Upon successful execution, the following output files will be available:
- **FastQC Reports**: Quality control reports for each initial FASTQ file.
- **Trimmed FASTQ Files**: Adapter-trimmed and quality-filtered reads.
- **SAM and BAM Files**: Alignment results, including sorted BAM files and an index file for each BAM file.
- **Gene Count Matrix**: A text file (`gene_count_matrix.txt`) containing the number of reads mapping to each gene.
- **Processed Gene Count Matrix**: A refined count matrix (`deseq2ready_gene_count_matrix.txt`) ready for differential expression analysis with DESeq2.
- **Sample Sheet**: A CSV file (sample_sheet.csv) listing each sample along with its associated condition.
- **Differential Expression Analysis Results**: The results of DESeq2 analysis, including differential expression table and figures.

[Back to the top](#rna-seq-analysis-pipeline)