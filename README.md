## Prerequisites
1) This pipeline is written in nextflow. Install nextflow by following the steps at [Nextflow Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2) The processes run on Docker containers. Install docker engine by following steps at [Docker Installation Guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/)

## How to use CoPi
1) Download this github repository using `git clone`.
```
git clone https://github.com/lanadelrea/CoPi
```

2) Navigate to the directory where `CoPi` is downloaded. Run the pipeline by indicating the path to the input and output folders. The input directory should contain fasta, fastq, bam, and bam index files. Remove the slash "/" after the input and output directory. You can run the pipeline using docker or conda by indicating the `-profile`. 

```
nextflow run CoPi/main.nf \
--in_dir <input directory> \
--out_dir <output directory> \
-profile <docker or conda>
```

For example:
```
nextflow run CoPi/main.nf --in_dir Samples --out_dir Results -profile docker

or 

nextflow run CoPi/main.nf --in_dir Samples --out_dir Results -profile conda
```

## Results
Output directory will have 8 folders containing results from each step used to test for SARS-CoV-2 variant coinfection. The summary report is in the final folder named `08-Report`. A sample [summary-report.html](https://github.com/lanadelrea/CoPi/blob/main/assets/summary-report.html) can be downloaded in assets. 

## Test the pipeline using sample files
To test the pipeline, use the files from [CoPi Sample Files](https://tinyurl.com/CoPi-Samples). Decompress `samples.tar.xz` by:
```
tar -xvf samples.tar.xz
```
