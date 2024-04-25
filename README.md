![Alt text](https://github.com/lanadelrea/CoPi/blob/main/assets/CoPi_logo.png)

## Prerequisites
1) This pipeline is written in nextflow. Install nextflow by following the steps at [Nextflow Installation Guide](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2) The processes run on Docker containers. Install docker engine by following steps at [Docker Installation Guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/)

## How to use CoPi
1) Download this github repository using `git clone`.
```
cd ~
git clone https://github.com/lanadelrea/CoPi
```

2) Navigate to the directory where `CoPi` is downloaded. Run the pipeline by indicating the path to the input and output folders. The input directory should contain fasta, fastq, bam, and bam index files.

```
cd ~
nextflow run CoPi --in_dir <input directory> --out_dir <output directory>
```

For example:
```
cd ~
nextflow run CoPi --in_dir /mnt/c/Samples --out_dir /mnt/c/Results
```

## Results
Output directory will have 8 folders containing results from each step used to test for SARS-CoV-2 variant coinfection. The summary report is in the final folder named `08-Report`. 

## Test the pipeline using sample files
To test the pipeline, use the files from [CoPi Sample Files](https://tinyurl.com/CoPi-Samples). Decompress `samples.tar.xz` by:
```
tar -xvf samples.tar.xz
```
