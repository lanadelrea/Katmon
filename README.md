This pipeline was initially designed to investigate possible Delta and Omicron co-infections after we encountered anomalous sequences during routine genomic surveillance. More recent developments have expanded its capability in detecting potential co-infections of any SARS-CoV-2 variants – including co-infections of newer variants, and cases involving two lineages under the same clade.

## Prerequisites
[Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation) and [Docker](https://docs.docker.com/engine/install/ubuntu/)

## How to use Katmon
1) Download this github repository using `git clone`.
```
git clone https://github.com/lanadelrea/Katmon
```

2) Navigate to the directory where `Katmon` is downloaded. Run the pipeline by indicating the path to the input and output folders. The input directory should contain FASTA, FASTQ, BAM, and BAM index files.

```
nextflow run Katmon --in_dir <input dir> --out_dir <output dir>
```

For example:
```
nextflow run Katmon --in_dir /Samples --out_dir /Results
```

The currently available parameters are as follows:
```
    Required arguments:
                 
        --in_dir             Input directory containing FASTA, FASTQ, BAM, and BAM index files
        --out_dir            Output directory for results
                  
    Optional arguments:
        --bammix_thresh      Set the bammix threshold for the proportion of the major allele
                                Default: 0.8
        -profile             Can be docker or conda
        -resume              To resume the pipeline
        -w                   The NextFlow work directory 
                             Delete this directory once the process is finished
                                Default: ${workDir} 
        --help               To view this help message
```

## Results
Output directory will have 8 folders containing results from each step used to test for SARS-CoV-2 variants coinfection. The summary report is in the final folder named `08-Report` with file name as `summary-report.html`. You can view sample results [here](https://github.com/lanadelrea/simKatmon/katmon-results).
