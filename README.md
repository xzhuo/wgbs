# Wang lab WGBS processing pipeline
Inital pipeline written by Hyung Joo Lee (https://github.com/hyungjoo-lee/wgbs). The first snakemake version written by Xiaoyu Zhuo (https://github.com/xzhuo/wgbs)

## main updates:
* pipeline works on Wang lab server, htcf, or with a docker.
* split Xiaoyu's Snakefile into ss.smk and aggr.smk
* split bismark into multiple steps

## manual:
### ss.smk: snakefile to process single samples
example usage:
```
$ cat ss.sh
#!/bin/bash
sample=$1
cd ~/storage1/wgbs/full/
threads=16
mode="full"
genome="mm10_conv"

echo "processing sample:" $sample

source activate snakemake
rm -rf ss.smk
ln -s ~/software/wgbs/ss.smk ./

snakemake ${sample}_mode_${mode}.txt -p --snakefile ss.smk --config genome=${genome} --jobs ${threads} 
```
* full mode: the entire pipeline
* shallow mode: For shallow sequenced libraries, only run certain parts of the pipeline for the purpose of qc. 
* <mark>note: genome assembly: use conventional chromosomes only. chrM and other unplaced contigs tend to cause error in bismark</mark>
* memory usage is hard coded. Need at least 100G to run smoothly
### aggr.smk: aggregate qc stats from multiple samples
example usage:
```
$ cat aggr.sh
#!/bin/bash
cd ~/storage1/wgbs/full
threads=4
genome="mm10_conv"

source activate snakemake
rm -rf aggr.smk
ln -s ~/software/wgbs/aggr.smk ./

snakemake -p --snakefile aggr.smk --config genome=${genome} --jobs ${threads} 
```
## Docker
example usage for ss.smk:
```
#!/bin/bash
for sample in Dpos Hpos DHpos DHneg
do
	LSF_DOCKER_PRESERVE_ENVIRONMENT=false LSF_DOCKER_VOLUMES='/storage1/fs1/hprc:/storage1/fs1/hprc' bsub \
	-q general -n 16 -G compute-hprc \
	-R 'span[hosts=1] select[mem>150G] rusage[mem=150G]' -M 150G \
	-a 'docker(funchansu/wgbs:v2.0.3)' -oo ~/storage1/wgbs/full/job_${sample}_v2.log \
	/bin/bash ~/storage1/wgbs/full/ss.sh ${sample}
done
```
example usage for aggr.smk:
```
LSF_DOCKER_PRESERVE_ENVIRONMENT=false LSF_DOCKER_VOLUMES='/storage1/fs1/hprc:/storage1/fs1/hprc' bsub -q general-interactive -G compute-hprc \
-Is -R 'span[hosts=1] select[mem>50G] rusage[mem=50G]' -M 50G \
-a 'docker(funchansu/wgbs:v2.0.3)' /bin/bash 
```
* aggr.smk is very fast and can be done with an interactive session.
## Known issues
* cpg bias calculation is hard coded and only supports hg38.