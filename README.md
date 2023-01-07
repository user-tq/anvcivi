一个根据annovar结果关联用药的测试流程

#### MANE基因与转录本
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz

#### 来自civic的突变与用药关系
https://civicdb.org/downloads/01-Dec-2022/01-Dec-2022-ClinicalEvidenceSummaries.tsv

### 目前已知有如下问题
1.初步判断，civic是根据hgvs标准命名，但是annovar的氨基酸水平命名并不完全遵守此规则，因此会有些位点遗漏（使用-hgvs参数，可以使annovar cDNA突变符合hgvs标准）。
例如：
    终止密码子突变时，hgvs: p.C2546* annovar: p.C2546X  
    同义突变时，hgvs:p.C123=   annovar: p.C123C
    
关于终止密码子[突变符号描述](http://varnomen.hgvs.org/bg-material/standards/),hgvs的原话是In addition HGVS nomenclature uses “Ter” (three-letter amino acid code) and “*” (three- and one-letter amino acid code) to indicate a translation termination (stop) codon. NOTE: in older versions the “X was used instead.

同义突变的描述的[更改提案](http://www.hgvs.org/mutnomen/accepted001.html)在2015年10月6日接受
总之annovar对hgvs格式不友好，似乎也不能选择氨基酸是简写还是全称。

~~目前设想是:先根据annovar -hgvs 获得满足hgvs命名规则的cDNA突变命名，再借由python hgvs包，解析出氨基酸水平突变命名。最终将其与civic用药信息关联。~~
发现由于转录本版本问题，并不是所有cDAN——to——protein都能成功，见https://github.com/biocommons/hgvs/issues/649#issue-1518426891

要实现hgvs标准化，需要保证annovar与hgvs所用的uta内的转录本版本一致.并且测试发现hgvs性能很差（可能是未本地化原因），一个cDNA突变映射到氨基酸水平需要10秒左右。
如此，考虑VEP对hgvs的支持程度，还不如直接以VEP的结果去关联civic的药物靶点信息.(也许后面会整个vepcivic,还能抄抄pcgr的代码)

2.civic涵盖突变类型很多，要用于实际生产与科研，必须人工逐一核对。本流程仅作测试用途，未逐一核对，目前仅测试点突变。
### 使用方法
```
#安装依赖
pip install yaml pandas
#下载
git clone https://github.com/user-tq/anvcivi.git
#运行
python anvcivi.py -i anovar_add_af.tsv  -o out_dir/
```
anvcivi 需要与anvcivi.config.yml 同一目录

anvcivi.config.yml需要配置civic.tsv与MANE.txt.gz的绝对路径
### 输入文件格式 
输入文件为调整后的annovar的输出，实际为一个tsv文件
至少需要有AAChange.refGene与af作为表头
