This pipeline was initially designed to investigate potential Delta and Omicron co-infections after anomalous sequences were detected during routine genomic surveillance. More recent developments have expanded its capability to identify potential co-infections involving any SARS-CoV-2 variants including newer variants and cases involving two lineages within the same clade.

<br />
<br />

![Katmon](https://github.com/lanadelrea/Katmon/blob/main/assets/Katmon.jpg)

<br />
<br />

## Prerequisites
[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [Docker](https://docs.docker.com/engine/install/ubuntu/)

## How to use Katmon
1) Download this github repository using `git clone`.
```
git clone https://github.com/lanadelrea/Katmon
```

2) Navigate to the directory where `Katmon` is downloaded. Run the pipeline by indicating the path to the input and output folders. The input directory should contain FASTA, FASTQ, BAM, and BAM index files.

```
nextflow run Katmon --indir <input dir> --outdir <output dir> --sequence <illumina OR ont> -profile <docker OR conda>
```

For example:
```
nextflow run Katmon --indir /mnt/d/Samples --outdir /mnt/d/Results --sequence illumina -profile docker
```

Currently available parameters are as follows:
```
    Required arguments:
                 
        --indir              Input directory containing FASTA, FASTQ, BAM, and BAM index files
        --outdir             Output directory for results
        --sequence           Can be 'illumina' or 'ont'

    Optional arguments:
        --bammix_thresh      Set the bammix threshold for the proportion of the major allele
                                        Default: 0.8
        -profile             Can be 'docker' or 'conda'
        -resume              To resume the pipeline
        -w                   The NextFlow work directory 
                             Delete this directory once the process is finished
                                        Default: ${workDir} 
```
```
    Use `nextflow run Katmon --help` to view the help message
```

## Results
Output directory will have 8 folders containing results from each step used to test for SARS-CoV-2 variants co-infection. The summary report is in the final folder named `08-Report` with file name as `summary-report.html`. You can view sample results [here](https://github.com/lanadelrea/simKatmon/tree/main/katmon-results).