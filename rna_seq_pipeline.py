import os
import subprocess
import glob
import csv
import pandas as pd
import rpy2.robjects as robjects

# --- Definion of paths and parameters ---

input_dir = 'Input/directory/path' # Path of the input directory where FASTQ files are located
output_dir = 'Output/folder/path' # Path of the output directory
trimmomatic_path = 'Trimmomatic/.jar file/path' # Complete path to the trimmomatic .jar file including file name
adapter_file = 'Trimmomatic/adapter file/path' # Complete path to the trimmomatic adapter file including file name 
hisat2_index = 'HiSAT2/index file/path' # Path of the HiSAT2 index file
gtf_file = 'Gene annotation file/path' # Complete path to the gene annotation (gene model) file including file name 
deseq2_path = 'DESeq2/R Script/path' #Complete path to the provided R script including file name

#***** Trimmomatic parameters *****
leading = 3
trailing = 3
slidingwindow = "4:15"
minlen = 36

# --- Creation of necessary output directories if they don't exist ---

os.makedirs(output_dir, exist_ok=True)
fastqc_dir = os.path.join(output_dir, "fastqc")
trimmed_dir = os.path.join(output_dir, "trimmed")
aligned_dir = os.path.join(output_dir, "aligned")
counted_dir = os.path.join(output_dir, "counted")
os.makedirs(fastqc_dir, exist_ok=True)
os.makedirs(trimmed_dir, exist_ok=True)
os.makedirs(aligned_dir, exist_ok=True)
os.makedirs(counted_dir, exist_ok=True)

# --- Gathering all FASTQ files alltogether ---

input_fastq = glob.glob(os.path.join(fastq_dir, "*.fastq*"))
input_fq = glob.glob(os.path.join(fastq_dir, "*.fq*"))
fastq_files = input_fastq + input_fq

# --- Identifying paired-end and single-end samples ---
## This step requires a particular file name convention.
## DO NOT use underscore "_" character except specifying forward end reverse end pairs.
## Forward and reverse end pairs should be named as:
## "samplename_1.extension" and "samplename_2.extension" for forward and reverse end pairs.
## "samplename" part of the paired end reads MUST be the same for end pairs

samples = {}
for fq in fastq_files:
    base_name = os.path.basename(fq)
    sample_name = base_name.split("_")[0]
    if sample_name not in samples:
        samples[sample_name] = {"paired": False, "files": []}
    if "_1.fastq" in base_name or "_2.fastq" in base_name:
        samples[sample_name]["paired"] = True
    samples[sample_name]["files"].append(fq)
    
# --- Step 1: Running FASTQC on all FASTQ files ---

for fq in fastq_files:
    subprocess.run(["fastqc", "-o", fastqc_dir, fq])

print("FASTQC for all files completed!")

# --- Step 2: Trimming adapters and low-quality regions using Trimmomatic ---
    
for sample_name, sample_info in samples.items():
    #***** For paired-end data *****
    if sample_info["paired"]:
        # Picking fastq file pairs of the sample
        fq1 = [f for f in sample_info["files"] if "_1.fastq" in f][0] 
        fq2 = [f for f in sample_info["files"] if "_2.fastq" in f][0]
        # Defining output directories for pairs
        output_paired_1 = os.path.join(trimmed_dir, os.path.basename(fq1).replace(".fastq", "_paired.fastq"))
        output_unpaired_1 = os.path.join(trimmed_dir, os.path.basename(fq1).replace(".fastq", "_unpaired.fastq"))
        output_paired_2 = os.path.join(trimmed_dir, os.path.basename(fq2).replace(".fastq", "_paired.fastq"))
        output_unpaired_2 = os.path.join(trimmed_dir, os.path.basename(fq2).replace(".fastq", "_unpaired.fastq"))
        # Trimmomatic command for paired-end reads
        trimmomatic_command = [
            "java", "-jar", trimmomatic_path, "PE", fq1, fq2, output_paired_1, output_unpaired_1,
            output_paired_2, output_unpaired_2, f"ILLUMINACLIP:{adapter_file}:2:30:10",
            f"LEADING:{leading}", f"TRAILING:{trailing}", f"SLIDINGWINDOW:{slidingwindow}", f"MINLEN:{minlen}"]
    #***** For single-end data *****
    else:
        # Picking fastq file of the sample
        fq1 = sample_info["files"][0]
        # Defining output directory for trimmed fastq file
        output_trimmed = os.path.join(trimmed_dir, os.path.basename(fq1).replace(".fastq", "_trimmed.fastq"))
        # Trimmomatic command for single-end read
        trimmomatic_command = [
            "java", "-jar", trimmomatic_path, "SE", fq1, output_trimmed,
            f"ILLUMINACLIP:{adapter_file}:2:30:10",
            f"LEADING:{leading}", f"TRAILING:{trailing}", f"SLIDINGWINDOW:{slidingwindow}", f"MINLEN:{minlen}"]
    # Running Trimmomatic based on chosen command
    subprocess.run(trimmomatic_command)

print("Trimming Done!")

# --- Step 3: Aligning trimmed reads to reference genome using HiSAT2 ---

for sample_name, sample_info in samples.items():
    #***** For paired-end data *****
    if sample_info["paired"]:
        # Picking trimmed fastq file pairs of the sample
        fq1 = os.path.join(trimmed_dir, f"{sample_name}_1_paired.fastq.gz")
        fq2 = os.path.join(trimmed_dir, f"{sample_name}_2_paired.fastq.gz")
        # Defining output directory for resulting SAM file of the sample
        output_sam = os.path.join(aligned_dir, f"{sample_name}.sam")
        # HiSAT2 command for paired-end reads
        hisat2_command = ["hisat2", "-x", hisat2_index, "-1", fq1, "-2", fq2, "-S", output_sam, "--dta"]
    #***** For single-end data *****
    else:
        fq1 = os.path.join(trimmed_dir, f"{sample_name}_trimmed.fastq")
        # Defining output directory for resulting SAM file of the sample
        output_sam = os.path.join(aligned_dir, f"{sample_name}.sam")
        # HiSAT2 command for single-end read
        hisat2_command = ["hisat2", "-x", hisat2_index, "-U", fq1, "-S", output_sam, "--dta"]
    # Running HiSAT2 based on chosen command
    subprocess.run(hisat2_command)

print("Mapping Done!")

# --- Step 4: Converting SAM files to Sorted BAM files using Samtools ---

sam_files = glob.glob(os.path.join(aligned_dir, "*.sam"))
for sam in sam_files:
    bam = sam.replace(".sam", ".bam")
    sorted_bam = bam.replace(".bam", "_sorted.bam")
    subprocess.run(["samtools", "view", "-bS", "-o", bam, sam]) #SAM file converted to BAM files
    subprocess.run(["samtools", "sort", "-o", sorted_bam, bam]) #BAM file sorted according to genomic location
    subprocess.run(["samtools", "index", sorted_bam]) #Sorted BAM file indexed
    os.remove(bam) # Intermediate BAM file (not sorted) removed from disk to save spaace

print("All SAM files converted to Sorted BAM Files!")

# --- Step 5: Generate gene count matrix using featureCounts ---

sorted_bam_files = glob.glob(os.path.join(aligned_dir, "*_sorted.bam"))
count_matrix_file = os.path.join(counted_dir, "gene_count_matrix.txt")
## IMPORTANT NOTICE: This section assumes that all the FASTQ files are coming from a single sequencing experiment.
## So make sure that all the FASTQ files in the input directory is either paired-en or single end

## Checking if any of the samples are paired-end by looking at the 'paired' key in the samples dictionary
## This will be used to decide whether featureCounts run in single-end or paired-end mode
paired_end = any(info["paired"] for info in samples.values())
# General featureCounts command
featurecounts_command = [
    "featureCounts", "-T", "4",  
    "-t", "exon", "-g", "gene_id", "-a", gtf_file, "-o", count_matrix_file]
# Adding the paired-end flag based on the paired_end decision
if paired_end:
    featurecounts_command.append("-p")  # Paired-end mode

# Adding Sorted BAM files to the command
## IMPORTANT NOTICE: This is where samples are sorted alphabetically while generating gene count matrix. 
## This should be keep in mind as this information will be needed while generating the sample sheet
featurecounts_command += sorted(sorted_bam_files)
# Running featureCounts based on chosen command
subprocess.run(featurecounts_command)

print("featureCounts executed and count matrix is generated!")

# --- Step 6: Pre-processing the resulting gene count matrix for DESeq2 anaysis  ---
## Since featureCounts generates a gene count matrix that has bulilt-in 6 columns first being the Gene ID's
## Other five annotation columns needs to be dropped as DESeq2 expected.
## So final count matrix which is ready for analyzing in DESeq2 will be consist of:
## Gene ID's in the first column and alphabetically sorted samples with gene counts

bam_base_names = [os.path.basename(bam).replace("_sorted.bam", "") for bam in sorted(sorted_bam_files)]
df = pd.read_csv(count_matrix_file, sep="\t")
df.columns = list(df.columns[:6]) + bam_base_names
df = df.drop(df.columns[1:6], axis=1)
deseq2_ready_file = os.path.join(counted_dir, "deseq2ready_gene_count_matrix.txt")
df.to_csv(deseq2_ready_file, sep="\t", index=False)

print(f"DESeq2-ready count matrix saved as {deseq2_ready_file}.")

# --- Step 7: Creating a sample sheet for DESeq2 ---
## IMPORTANT NOTICE: This section assigns conditions to samples like "control vs treated" or "healthy vs diseased" etc.
## So when writing contions for each sample, NEVER forget that samples are alphabetically ordered and type accordingly.
conditions = [
    'treated',  # Example condition for first sample
    'treated',  # Example condition for second sample
    'treated',  # Example condition for third sample
    'treated',  # Example condition for fourth sample
    'control',  # Example condition for fifth sample
    'control',  # Example condition for sixth sample
    # Add conditions if needed
]
# Ensuring the number of conditions matches with the number of samples and if it doesn't, it raising an error
assert len(conditions) == len(bam_base_names), "Number of conditions must match the number of samples."

sample_sheet_data = list(zip(bam_base_names, conditions))
sample_sheet_file = os.path.join(counted_dir, "sample_sheet.csv")
with open(sample_sheet_file, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["sample", "condition"])
    writer.writerows(sample_sheet_data)

print(f"Sample sheet saved as {sample_sheet_file}.")

# --- Step 8: Running DESeq2 for differential gene expression analysis and visualization
# Setting R's working directory to the relevant directory
working_directory = counted_dir
robjects.r(f'setwd("{working_directory}")')
subprocess.run(["Rscript", deseq2_path])

print("DESeq2 is completed and resulting figures saved!")
print("Differential Gene Expression analysis of bulk RNA-seq pipeline is completed!")
