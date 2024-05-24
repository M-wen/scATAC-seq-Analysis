# Alinment to genome using cellranger-atac

For the fastq file, those steps can be skipped to step 3.
## 1.Convert the SRA data form to the fastq from
```shell
for x in `ls *.sra`
do
name=$(echo $x|sed 's/\.sra//g')
echo fastq-dump --split-files --gzip ${x} > ${name}_sra2fastq.sh
nohup sh ${name}_sra2fastq.sh >${name}_sra2fq.log 2>&1 &
done 
```
### 1.1 Check the data integrity

#### 1.1.1 Check if the spots in the log file are equal to the metadata of SRR data.
cat ${name}_sra2fq.log
```shell
#Read 64090285 spots for xx.sra
#Writen xx spots for xx.sra
```
#### 1.1.2 Check if the length matches the metadata of SRR data.
```shell
less xx_1.fastq.gz
less xx_2.fastq.gz
```
## 2.rename the fataq files for cellranger-atac
```shell
for x in `ls *.sra`
do
name=$(echo $x|sed 's/\.sra//g')
mv ${name}_1.fastq.gz ${name}_S1_L001_R1_001.fastq.gz
mv ${name}_2.fastq.gz ${name}_S1_L001_R2_001.fastq.gz
done
```
## 3.align to genome using cellranger-atac
```shell
for x in `cat list`
do
echo cellranger-atac count --id ${x} --reference path/to/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --fastqs /path/to/00.data/01.scATAC/ --sample ${x} --localcores 4 --localmem 100 > cellranger_atac_${x}.sh
nohup bash cellranger_atac_${x}.sh &
done
```
