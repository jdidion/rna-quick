# Overview

RNA Quick is a pipeline for rapid analysis of gene and transcript quantification and differential expression. Given 1) FASTQ files for one or more samples partitioned between one or more experimental conditions, and 2) a metadata file describing the samples, the following steps are performed:

1. If necessary, download the reference genome (default: GRCh37) and annotation catalogue (default: GENCODE v19); existing files can be used instead if available.
2. Build a transcriptome index using Kallisto.
3. Perform initial quantification. By default, only the first 1M reads from each sample are used.
4. Isoform filtering: based on the suggestions in Soneson et al. (2016), the annotation catalogue is filtered to exclude isoforms with relative abundance less than a threshold value (default: 0.1) in all conditions. For example, if gene X has 3 isoforms with relative abundance (0.7,0.25,0.05) and (0.92,0.07,0.01) in the two conditions being tested, isoform 3 would be excluded because its relative abundance is < 0.1 in both conditions.
5. Generate a new Kallisto index using the filtered annotation catalogue.
6. Perform transcript quantification of all samples using Kallisto.
7. Perform differential transcript expression testing using Sleuth.
8. Aggregate transcript counts into gene counts using tximport.
9. Perform differential gene expression analysis using DESeq2.

# Installation

RNA Quick is implemented as a NextFlow pipeline and uses Docker containers to run all software packages. Follow the NextFlow installation instructions here: http://www.nextflow.io/.

# Running RNA Quick

## 1. Prepare the metadata file

The metadata file is a tab-delimited file (default name: libraries.txt) describing the sequencing data to be processed. It has the following columns:

1. UniqueID: A unique identifier for the sequencing library.
2. SourceID: Identifier for the source material from which the sequencing library was derived; this is to identify technical replicates derived from the same source material.
3. LibraryType (optional): Describes the sequencing method used (single or paired). If all libraries were generated using the same method, this column can be omitted and library type information can be specified in the configuration.
4. NumFiles (optional): Number of FASTQ files (or file pairs) associated with this library. If missing, all libraries are assumed to have a single file/pair.
5. FileNames (optional): Rather than provide NumFiles, you can list the full path names of all of the FASTQ files for each library. This is a comma-delimited string. For paired end files, each pair must appear one after the other, e.g. pairA_1.fq,pairA_2.fq,pairB_1.fq,pairB_2.fq.
6. All remaining columns are treated as phenotypes. These can either be experimental conditions or covariates to be included in the expression model. By default, the column named "Condition" will be used as the experimental condition to test for differential expression; if this behavior is not desired, or if there is more than one model to be tested, the models should be defined in the configuration.

If optional columns are included, then each row must have a valid value.

## 2. Execution using defaults

From the directory containing your metadata file, run:

```
nextflow run jdidion/rna-quick
```

## 3. Execution using custom configuration

In the same directory that contains your metadata file, create a copy of `nextflow.config`. Edit the file to change any parameter values. Then run:

```
nextflow run jdidion/rna-quick
```

Make sure to publish this file along with your metadata file, in order to enable others to perform exact replication of your analysis.

# Project directory structure

RNA Quick will perform all operations within the working directory (default: current directory).