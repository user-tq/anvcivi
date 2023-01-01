
一个根据annovar结果关联用药的测试流程

#### MANE基因与转录本
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz

#### 来自civic的突变与用药关系
https://civicdb.org/downloads/01-Dec-2022/01-Dec-2022-ClinicalEvidenceSummaries.tsv

### 目前已知有如下问题
1.初步判断，civic是根据hgvs标准命名，但是annovar的氨基酸水平命名并不遵守此规则（使用-hgvs参数，可以使annovar cDNA突变符合hgvs标准）。
目前设想是:先根据annovar -hgvs 获得满足hgvs命名规则的cDNA突变命名，再借由python hgvs包，解析出氨基酸水平突变命名。最终将其与civic用药信息关联。
由于hgvs依赖与pysam，pysam在windows下无法安装，暂时搁置。

2.civic涵盖突变类型很多，要用于实际生产与科研，必须人工逐一核对。本流程仅作测试用途，未逐一核对，目前仅测试点突变。
### 使用方法
```
#安装依赖
pip install yaml pandas
#下载
git clone .../.git
#运行
python anvcivi.py -i anovar_add_af.tsv  -o outdir/
```
anvcivi 需要与anvcivi.config.yml 同一目录

anvcivi.config.yml需要配置civic.tsv与MANE.txt.gz的绝对路径
### 输入文件格式 
输入文件为调整后的annovar的输出，实际为一个tsv文件
至少需要有AAChange.refGene与af作为表头