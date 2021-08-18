# umi_variant_calling : Output <!-- omit in toc -->

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished.
All paths are relative to the top-level results directory.

## Pipeline overview <!-- omit in toc -->

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [UMI processing STAGE ONE](#umi-processing-stage-one)
  - [align group consensus filter](#align-group-consensus-filter)
    - [fastq to bam](#fastq-to-bam)
      - [FGBIO FastqToBam](#fgbio-fastqtobam)
    - [mark illumina adapters](#mark-illumina-adapters)
      - [PICARD MarkIlluminaAdapters](#picard-markilluminaadapters)
    - [unaligned marked bam to fastq](#unaligned-marked-bam-to-fastq)
      - [PICARD SamToFastq](#picard-samtofastq)
    - [align fastq](#align-fastq)
      - [BWA align](#bwa-align)
    - [merge aligned reads](#merge-aligned-reads)
      - [PICARD MergeBamAlignment](#picard-mergebamalignment)
    - [map bam to consensus reads](#map-bam-to-consensus-reads)
      - [FGBIO GroupReadsByUmi](#fgbio-groupreadsbyumi)
    - [call molecular consensus reads](#call-molecular-consensus-reads)
      - [FGBIO CallMolecularConsensusReads](#fgbio-callmolecularconsensusreads)
    - [filter consensus reads](#filter-consensus-reads)
      - [FGBIO FilterConsensusReads](#fgbio-filterconsensusreads)
- [UMI processesing STAGE TWO](#umi-processesing-stage-two)
  - [align markduplicates](#align-markduplicates)
    - [unaligned filtered bam to fastq](#unaligned-filtered-bam-to-fastq)
      - [PICARD SamToFastq](#picard-samtofastq-1)
    - [align fastq](#align-fastq-1)
      - [BWA align](#bwa-align-1)
    - [sort sam](#sort-sam)
      - [PICARD SortSam](#picard-sortsam)
      - [PICARD SortSam from FilterConsensusReads](#picard-sortsam-from-filterconsensusreads)
    - [merge bam](#merge-bam)
      - [PICARD MergeBam](#picard-mergebam)
  - [Mark Duplicates](#mark-duplicates)
    - [GATK MarkDuplicates](#gatk-markduplicates)
  - [TSV files](#tsv-files)
- [QC and reporting](#qc-and-reporting)
  - [QC](#qc)
    - [FastQC](#fastqc)
    - [bamQC](#bamqc)
    - [GATK MarkDuplicates reports](#gatk-markduplicates-reports)
    - [samtools stats](#samtools-stats)
  - [Reporting](#reporting)
    - [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

## UMI processing STAGE ONE

`umi-variant-caller` pre-processes raw `FASTQ` files or `unmapped BAM` files, based on [IDT](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/user-guide-manual/demultiplexing-data-containing-unique-molecular-indexes-(umis).pdf?sfvrsn=99953207_26).

### align group consensus filter

#### fastq to bam

##### FGBIO FastqToBam
[fgbio](http://fulcrumgenomics.github.io/fgbio/) fgbio is a command line toolkit for working with genomic and particularly next generation sequencing data.
[FastqToBam](http://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html)
Generates an unmapped BAM (or SAM or CRAM) file from fastq files. Takes in one or more fastq files (optionally gzipped), each representing a different sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read structures to allocate bases in those reads to template reads, sample indices, unique molecular indices, or to designate bases to be skipped over.

#### mark illumina adapters

##### PICARD MarkIlluminaAdapters
[MarkIlluminaAdapters](https://gatk.broadinstitute.org/hc/en-us/articles/360037434291-MarkIlluminaAdapters-Picard-)
Reads the file and rewrites it with new adapter-trimming tags. 

#### unaligned marked bam to fastq

##### PICARD SamToFastq
[SamToFastq](https://gatk.broadinstitute.org/hc/en-us/articles/360041416572-SamToFastq-Picard-)
Convert BAM files to FASTQ filesPicard’s SamToFastq tool extracts read sequences and base quality scores from the input SAM/BAM file and writes them into the output file in FASTQ format. Sequencing adapter sequences present in the reads in the BAM file can be trimmed before writing to the FASTQ files.

#### align fastq

##### BWA align
[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome.
Such files are intermediate and not kept in the final files delivered to users.

#### merge aligned reads

##### PICARD MergeBamAlignment
[MergeBamAlignment](https://gatk.broadinstitute.org/hc/en-us/articles/360037225832-MergeBamAlignment-Picard-)
Picard’s MergeBamAlignment produces a new BAM file that includes all aligned and unaligned reads, and also contains additional read attributes from the uBAM that may otherwise be lost during the alignment process. Note that MergeBamAlignment expects to find a sequence dictionary in the same directory as REFERENCE_SEQUENCE and expects it to have the same base name as the reference FASTA, but with the extension “.dict”. 

#### map bam to consensus reads

##### FGBIO GroupReadsByUmi
[GroupReadsByUmi](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)
Identify which reads come from the same source molecule by using fgbio’s GroupReadsByUmi tool, which assigns a unique source molecule ID to each applicable read, stores the ID in the MI tag, and outputs a BAM file that is sorted by the MI tag and ready for consensus calling. The source molecule is identified using a combination of UMI sequence and mapping positions from reads 1 and 2.

#### call molecular consensus reads

##### FGBIO CallMolecularConsensusReads
[CallMolecularConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html)
This step generates unmapped consensus reads from the output of GroupReadsByUmi.

#### filter consensus reads

##### FGBIO FilterConsensusReads
[FilterConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html)
Filter consensus reads generated by CallMolecularConsensusReads using fgbio’s FilterConsensusReads. There are two kinds of filtering: 1) masking or filtering individual bases in reads, and 2) filtering reads (i.e., not writing them to the output file). Base-level masking/filtering is only applied if per-base tags are present (see the documentation for CallMolecularConsensusReads for tag descriptions). Read-level filtering is always applied.When filtering reads, secondary alignments and supplementary records may be removed independently if they fail one or more filters. If either R1 or R2 primary alignments fail a filter, then all records for the template will be filtered out.

## UMI processesing STAGE TWO

### align markduplicates

#### unaligned filtered bam to fastq

##### PICARD SamToFastq

Return unmapped filtered marked bam to fastq

#### align fastq

##### BWA align

align to reference again

#### sort sam



##### PICARD SortSam

sort bam from alinment by queryname

##### PICARD SortSam from FilterConsensusReads

sort filtered bam by query name

#### merge bam

##### PICARD MergeBam

merge the two bams.

### Mark Duplicates

#### GATK MarkDuplicates
By default, `umi_variant_caller` will use [GATK MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360042912511-MarkDuplicatesSpark), `Spark` implementation of [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360042477492-MarkDuplicates-Picard), which locates and tags duplicate reads in a `BAM` or `SAM` file, where duplicate reads are defined as originating from a single fragment of DNA.

Specify `--no_gatk_spark` to use `GATK MarkDuplicates` instead.

This directory is the location for the `BAM` files delivered to users.
Besides the `duplicates-marked BAM` files, the recalibration tables (`*.recal.table`) are also stored, and can be used to create `recalibrated BAM` files.

For all samples:

**Output directory: `results/Preprocessing/[SAMPLE]/DuplicatesMarked`**

- `[SAMPLE].md.bam` and `[SAMPLE].md.bai`
  - `BAM` file and index

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).


### TSV files

The `TSV` files are auto-generated and can be used by `Sarek` for further processing and/or variant calling.

For further reading and documentation see the [`--input`](usage.md#--input) section in the usage documentation.

For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `duplicates_marked_no_table.tsv`, `duplicates_marked.tsv` and `recalibrated.tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps.


For all samples:

**Output directory: `results/Preprocessing/TSV`**

- `mapped.tsv`, `mapped_no_duplicates_marked.tsv` and `recalibrated.tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps.
- `mapped_[SAMPLE].tsv`, `mapped_no_duplicates_marked_[SAMPLE].tsv` and `recalibrated_[SAMPLE].tsv`
  - `TSV` files to start `Sarek` from `prepare_recalibration`, `recalibrate` or `variantcalling` steps for a specific sample.


## QC and reporting

### QC

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads.
It provides information about the quality score distribution across your reads, per base sequence content (`%A/T/G/C`), adapter contamination and overrepresented sequences.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/fastqc`**

- `sample_R1_XXX_fastqc.html` and `sample_R2_XXX_fastqc.html`
  - `FastQC` report containing quality metrics for your untrimmed raw `FASTQ` files
- `sample_R1_XXX_fastqc.zip` and `sample_R2_XXX_fastqc.zip`
  - Zip archive containing the FastQC report, tab-delimited data file and plot images

> **NB:** The `FastQC` plots displayed in the `MultiQC` report shows _untrimmed_ reads.
> They may contain adapter sequence and potentially regions with low quality.

For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

#### bamQC

[Qualimap bamqc](http://qualimap.bioinfo.cipf.es/) reports information for the evaluation of the quality of the provided alignment data.
In short, the basic statistics of the alignment (number of reads, coverage, GC-content, etc.) are summarized and a number of useful graphs are produced.

Plot will show:

- Stats by non-reference allele frequency, depth distribution, stats by quality and per-sample counts, singleton stats, etc.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/bamQC`**

- `VariantCaller_[SAMPLE].bcf.tools.stats.out`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [Qualimap bamqc manual](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#id7)

#### GATK MarkDuplicates reports

More information in the [GATK MarkDuplicates section](#gatk-markduplicates)

Duplicates can arise during sample preparation _e.g._ library construction using PCR.
Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.
These duplication artifacts are referred to as optical duplicates.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/MarkDuplicates`**

- `[SAMPLE].bam.metrics`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [MarkDuplicates manual](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/picard_sam_markduplicates_MarkDuplicates.php).

#### samtools stats

[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from `BAM` files and outputs in a text format.

Plots will show:

- Alignment metrics.

For all samples:

**Output directory: `results/Reports/[SAMPLE]/SamToolsStats`**

- `[SAMPLE].bam.samtools.stats.out`
  - Raw statistics used by `MultiQC`

For further reading and documentation see the [`samtools` manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTION

### Reporting

#### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project.
Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the `MultiQC` output for future traceability.

**Output files:**

- `multiqc/`  
  - `multiqc_report.html`
    - Standalone HTML file that can be viewed in your web browser
  - `multiqc_data/`
    - Directory containing parsed statistics from the different tools used in the pipeline
  - `multiqc_plots/`
    - Directory containing static images from the report in various formats

For more information about how to use `MultiQC` reports, see [https://multiqc.info](https://multiqc.info).

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  - Documentation for interpretation of results in HTML format: `results_description.html`.
