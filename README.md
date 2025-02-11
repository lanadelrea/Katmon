## :clipboard: Prerequisites
1) This pipeline is written in nextflow. Install nextflow by following the steps at [Nextflow Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2) The processes run on Docker containers. Install docker engine by following steps at [Docker Installation Guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/)

## :computer: How to use Katmon
1) Download this github repository using `git clone`.
```
git clone https://github.com/lanadelrea/Katmon
```

2) Navigate to the directory where `Katmon` is downloaded. Run the pipeline by indicating the path to the input and output folders. The input directory should contain fasta, fastq, bam, and bam index files. Remove the slash "/" after the input and output directory. You can run the pipeline using docker or conda by indicating the `-profile`. 

```
nextflow run Katmon --in_dir <input dir> --out_dir <output dir> -profile <docker or conda>
```

For example:
```
nextflow run Katmon --in_dir /Samples --out_dir /Results -profile docker
```

You can also view the help documentations by:
```
nextflow run Katmon --help
```

## :open_file_folder: Results
Output directory will have 7 folders containing results from each step used to test for SARS-CoV-2 variant coinfection. The summary report is in the final folder named `07-Report` with the file name as `summary-report.html`.
