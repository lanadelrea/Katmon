## How to use CoPi
1) Download this github repository by 
```git clone https://github.com/lanadelrea/CoPi```

2) The processes are run through docker containers. Here is the list of the docker images to pull:
```
docker pull staphb/pangolin:latest
docker pull nextstrain/nextclade:latest
docker pull ufuomababatunde/bammix:v1.0.0
docker pull staphb/freyja:latest
docker pull pegi3s/samtools_bcftools:latest
docker pull vandhanak/bcftools:1.3.1
```

3) This pipeline is written in nextflow, so make sure to install nextflow first by following the steps at https://www.nextflow.io/docs/latest/getstarted.html#installation

4) If you already have nextflow, run the pipeline by indicating path to the input and output directories. Input directory should contain fasta, fastq, bam, and bam index files. 
```
nextflow run CoPi --in_dir <input directory> --out_dir <output directory>
```

For example:
```
nextflow run CoPi --in_dir /mnt/c/Samples --out_dir /mnt/c/Results
```

## Results
Once the pipeline finished, output directory will have 8 folders containing results from each of the step used to test for SARS-CoV-2 coinfection. Summary report is in the final folder. 

## Test the pipeline using sample files
To test the pipeline, use the `Samples` directory containing necessary files from https://tinyurl.com/CoPi-Samples. Decompress the `samples.tar.xz` by
```
xz -d samples.tar.xz
```