# Detection of homotypic doublets using Amulet

## 1.download resource database and install Amulet
### 1.1 install Amulet
```shell
pip3 install numpy pandas scipy statsmodels
wget -O Amulet_v1.1.tar.gz "https://github.com/UcarLab/AMULET/archive/refs/tags/v1.1.tar.gz"
tar -xzvf Amulet_v1.1.tar.gz
cd AMULET-1.1
chmod +x AMULET.sh
```
### 1.2 Create repetitive elements file required for Amulet
```shell
mkdir -p $resource_path/repetitive_elements
cd $resource_path/repetitive_elements
# download genomicSuperDups from UCSC 
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz
# download simpleRepeats from UCSC
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
# download the exclusion list from ENCODE
wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
zcat simpleRepeat.txt.gz | awk -v OFS="\t" '{print $2,$3,$4}' > hg38_simpleRepeat.bed
zcat genomicSuperDups.txt.gz | awk -v OFS="\t" '{print $1,$3,$4}' > hg38_genomicSuperDups_1.bed
zcat genomicSuperDups.txt.gz | awk -v OFS="\t" '{print $7,$8,$9}' > hg38_genomicSuperDups_2.bed
gunzip ENCFF356LFX.bed.gz 
mv ENCFF356LFX.bed GRCh38_unified_blacklist.bed
cat hg38_simpleRepeat.bed hg38_genomicSuperDups_1.bed hg38_genomicSuperDups_2.bed GRCh38_unified_blacklist.bed > blacklist_repeats_segdups_rmsk_hg38.bed
```
## 2. Detect homotypic doublets with Amulet.

```shell
AMULET.sh /path/to/fragments.tsv.gz /path/to/singlecell.csv /path/to/human_autosomes.txt /path/to/blacklist_repeats_segdups_rmsk_hg38.bed /path/to/output/ /path/to/shellscript/
```
