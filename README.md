一个根据annovar结果关联用药的测试流程

#### MANE基因与转录本
https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.summary.txt.gz

#### 来自civic的突变与用药关系
https://civicdb.org/downloads/01-Dec-2022/01-Dec-2022-ClinicalEvidenceSummaries.tsv

### 使用须知
1.civic中的蛋白质与cDNA命名满足hgvs标准命名规范，但是annovar的氨基酸水平命名并不完全遵守此规则，因此会有些位点遗漏（使用-hgvs参数，可以使annovar cDNA突变符合hgvs标准）。
例如：
    终止密码子突变时，hgvs: p.C2546* annovar: p.C2546X  
    同义突变时，hgvs:p.C123=   annovar: p.C123C
    
关于终止密码子[突变符号描述](http://varnomen.hgvs.org/bg-material/standards/),hgvs的原话是In addition HGVS nomenclature uses “Ter” (three-letter amino acid code) and “*” (three- and one-letter amino acid code) to indicate a translation termination (stop) codon. NOTE: in older versions the “X was used instead.

同义突变的描述的[更改提案](http://www.hgvs.org/mutnomen/accepted001.html)在2015年10月6日接受
总之annovar对hgvs格式不友好，疑似遵循老版HGVS标准

~~目前设想是:先根据annovar -hgvs 获得满足hgvs命名规则的cDNA突变命名，再借由python hgvs包，解析出氨基酸水平突变命名。最终将其与civic用药信息关联。~~
发现由于转录本版本问题，并不是所有cDNA_to_protein都能成功，见https://github.com/biocommons/hgvs/issues/649#issue-1518426891

要实现hgvs标准化，需要保证作为输入的annovar的注释结果与hgvs所用的uta内的转录本版本一致，测试发现hgvs性能很差（可能是未本地化原因），一个cDNA突变映射到氨基酸水平需要10秒左右。

考虑以上种种原因，显然，基于VEP的结果去关联civic的药物靶点信息，是更好的选择(后续有机会仔细学习VEP，暂定名vepcivic)
目前，只能将*替换为X以实现匹配

2.civic涵盖突变类型很多，要用于实际生产与科研，必须人工逐一核对。本流程仅作测试用途，部分变异暂未能完全解析（见civic_clean/unaccept_var.tsv）

注释层级与优先顺序为 ：
    蛋白质层面准确变异位点，如：EGFR p.L858R
    蛋白质层面位置无准确位点，如：EGFR p.G719
    cDNA准确位点，如：JAK2 c.1641+1dup （注意，这取决于civic内是否用cDNA层级表示）
    任何突变，如KRAS MUTATION

3.civic本身存在变体未细分的问题，如EGFR EXON 19 DEL，EGFR EXON 20 INS，后续我打算将这类突变归类起来，创建一个json文件，直接根据annovar结果映射

### 使用方法
```
#安装依赖
pip install yaml pandas
#下载
git clone https://github.com/user-tq/anvcivi.git
#运行
python anvcivi.py -i example.tsv  -o test/  -d all
```
anvcivi.py 需要与anvcivi.config.yml 同一目录
anvcivi.config.yml需要配置accept_var_info.tsv与MANE.GRCh38.v1.0.summary.txt.gz的绝对路径
### 输入文件格式 
输入文件格式：见example.tsv
