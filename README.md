To run the pipeline:
**nextflow run CoPi --in_dir `<input directory>` --out_dir `<output directory>`**

Input directory should contain the fasta, fastq, bam, and bam index files.

This pipeline is mainly run through docker containers. Here is the list of the docker images to pull:
```
docker pull staphb/pangolin:latest
docker pull nextstrain/nextclade:latest
docker pull ufuomababatunde/bammix:v1.0.0
docker pull staphb/freyja:latest
docker pull pegi3s/samtools_bcftools:latest
docker pull vandhanak/bcftools:1.3.1
```