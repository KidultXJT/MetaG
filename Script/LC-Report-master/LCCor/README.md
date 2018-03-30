- [LCCor](# LCCor)
 - [Introduction](## Introduction)
 - [Installation](## Installation)
 - [Quick Start!](## Quick Start!)
 - [Help Page](## Help Page)
 - [Contact](## Contact)
 
 
LCCor
===============
LC == **Lancet**, means solve problem in a **SHORTCUT**. LC is a SERIES, such as LCORD, LCNT, LCGDF, LCBDIV, LCDA, LCMVA...ect
 
Introduction
---------------
In LCCor, LC == Lancet is a SERIES, Cor == **Correlation Analysis**. LCCor is a Software, that solve the Correlation requirment in a SHORTCUT. 

LCCor Based on base[R-Package] and vegan[R-Package], AND jinja2[Python]. 

1. **Correlation Analysis [pearson, spearman, kandell]**

2. **Correlation TEST**

The three methods each estimate the association between paired samples and compute a test of the value being zero. They use different measures of association, all in the range [-1, 1] with 0 indicating no association. These are sometimes referred to as tests of no correlation, but that term is often confined to the default method.

If method is "pearson", the test statistic is based on Pearson's product moment correlation coefficient cor(x, y) and follows a t distribution with length(x)-2 degrees of freedom if the samples follow independent normal distributions. If there are at least 4 complete pairs of observation, an asymptotic confidence interval is given based on Fisher's Z transform.

If method is "kendall" or "spearman", Kendall's tau or Spearman's rho statistic is used to estimate a rank-based measure of association. These tests may be used if the data do not necessarily come from a bivariate normal distribution.

For **Kendall's** test, by default (if exact is NULL), an exact p-value is computed if there are less than 50 paired samples containing finite values and there are no ties. Otherwise, the test statistic is the estimate scaled to zero mean and unit variance, and is approximately normally distributed.

For **Spearman's** test, p-values are computed using algorithm AS 89 for n < 1290 and exact = TRUE, otherwise via the asymptotic t approximation. Note that these are ‘exact’ for n < 10, and use an Edgeworth series approximation for larger sample sizes (the cutoff has been changed from the original paper).

3. **Mantel TEST**

Mantel statistic is simply a correlation between entries of **two dissimilarity matrices** (some use cross products, but these are linearly related). However, the significance cannot be directly assessed, because there are N(N-1)/2 entries for just N observations. Mantel developed asymptotic test, but here we use permutations of N rows and columns of dissimilarity matrix. See permutations for additional details on **permutation tests** in Vegan.

**Partial Mantel** statistic uses **partial correlation** conditioned on the **third matrix**. Only the first matrix is permuted so that the correlation structure between second and first matrices is kept constant. Although mantel.partial silently accepts other methods than "pearson", partial correlations will probably be wrong with other methods.

The function uses cor, which should accept alternatives pearson for product moment correlations and spearman or kendall for rank correlations.

Installation
---------------

Just git clone this Repository, AND install the **R-Package** and **jinja2**, the Software require.

Require R Package::
```r
```

Install Jinja2
```python
```


Quick Start!
---------------
Try This Script, AND Have Fun~ !
```bash
# workdir LCCor

cd LCCor
cd Example
cat work.sh

## python ../MultiLCCor.py -a data/X.xls -b data/Y.xls -o MultiLCCor -k no
## python ../LCCor.py      -a data/X.xls -b data/Y.xls -o LCCor      -k no
```

Help Page
---------------
You can Find the Help Page by:
```bash
python ../LCCor.py --help
# or
python ../MultiLCCor.py --help
```


```text
# LCCor.py
Usage: 
copyright: 
	Guangzhou Sagene Biotech Co.,Ltd

contact e-mail:
	junting.xie@sagene.com.cn

description:
	

options:
	-h	--help

	-a	--infile1	(necessary): The Infile1 (str)
	-b	--infile2	(necessary): The Infile2 (str)
	-o	--outDir	(optional):  The Output Direction (str, Default = "./PAnalysis")
	-m	--method	(optional):  Correlation Method(str, can be a list, Default = "spearman")
	-d	--dist		(optional):  Distance Method(str, can be a list, Default = "euclidean")
	-P,	--sigthrd	(optional):  Correlation TEST Pvalue Result Table (float, Defalult = 0.05)
	-R,	--corthrd	(optional):  Correlation Analysis R statistics Result Tabule (float, Default = 0.5)
	
	-L,	--language	(optional):  Report Language (str, Default = EN)

	Project Information:
	-c	--client	(optional):  Client Name(str, Default = "Kidult")
	-i	--institute	(optional):  Institute Name(str, Default = "Sagene")
	-s	--samples	(optional):  Samples Name(str, Default = "")
	-t	--pjectypes	(optional):  Types of This Project(str, Default = "SGA PAnalysis")
	-p	--project	(optional):  Project ID(str, Default = "SGA0000")

usage e.g.:
	python LCCor.py -h 
	python LCCor.py -a infile1 -b infile2
	python LCCor.py --infile1 infile1 --infile2 infile2

Use -h/--help option to see the usage.


Options:
  -h, --help            show this help message and exit
  -a INFILE1, --infile1=INFILE1
                        Input file1
  -b INFILE2, --infile2=INFILE2
                        Input file2
  -o OUTDIR, --outDir=OUTDIR
                        Working Direction
  -m CORMETHOD, --method=CORMETHOD
                        Correlation Method
  -d DISTMETHOD, --dist=DISTMETHOD
                        distance Method
  -P PTHRESHOLD, --sigthrd=PTHRESHOLD
                        p threshold
  -R RTHRESHOLD, --corthrd=RTHRESHOLD
                        r threshold
  -L LANGUAGE, --language=LANGUAGE
                        Report Language
  -c CLIENT, --client=CLIENT
                        Client Name
  -i INSTITUTE, --institu=INSTITUTE
                        Institute Name
  -s SAMPLES, --samples=SAMPLES
                        Samples Name
  -t TYPES, --types=TYPES
                        Project Types
  -p PROJECT, --project=PROJECT
                        Project ID
  -C CONTRACT, --contract=CONTRACT
                        Contract ID
  -k KEEPRMD, --keepRMD=KEEPRMD
                        KeepRMD or Not
```

```text
Usage: 
copyright: 
	Guangzhou Sagene Biotech Co.,Ltd

contact e-mail:
	junting.xie@sagene.com.cn

description:
	

options:
	-h	--help

	-a	--infile1	(necessary): The Infile1 (str)
	-b	--infile2	(necessary): The Infile2 (str)
	-o	--outDir	(optional):  The Output Direction (str, Default = "./PAnalysis")
	-m	--method	(optional):  Correlation Method(str, can be a list, Default = "pearson,spearman" Can be a List)
	-d	--dist		(optional):  Distance Method(str, can be a list, Default = "euclidean,bray" Can be a List)
	-P,	--sigthrd	(optional):  Correlation TEST Pvalue Result Table (float, Defalult = 0.05)
	-R,	--corthrd	(optional):  Correlation Analysis R statistics Result Tabule (float, Default = 0.5)
	
	-L,	--language	(optional):  Report Language (str, Default = EN)
	
	Project Information:
	-c	--client	(optional):  Client Name(str, Default = "Kidult")
	-i	--institute	(optional):  Institute Name(str, Default = "Sagene")
	-s	--samples	(optional):  Samples Name(str, Default = "")
	-t	--pjectypes	(optional):  Types of This Project(str, Default = "SGA PAnalysis")
	-p	--project	(optional):  Project ID(str, Default = "SGA0000")

usage e.g.:
	python LCCor.py -h 
	python LCCor.py -a infile1 -b infile2
	python LCCor.py --infile1 infile1 --infile2 infile2

Use -h/--help option to see the usage.


Options:
  -h, --help            show this help message and exit
  -a INFILE1, --infile1=INFILE1
                        Input file1
  -b INFILE2, --infile2=INFILE2
                        Input file2
  -o OUTDIR, --outDir=OUTDIR
                        Working Direction
  -m CORMETHOD, --method=CORMETHOD
                        Correlation Method
  -d DISTMETHOD, --dist=DISTMETHOD
                        distance Method
  -P PTHRESHOLD, --sigthrd=PTHRESHOLD
                        p threshold
  -R RTHRESHOLD, --corthrd=RTHRESHOLD
                        r threshold
  -L LANGUAGE, --language=LANGUAGE
                        Report Language
  -c CLIENT, --client=CLIENT
                        Client Name
  -i INSTITUTE, --institu=INSTITUTE
                        Institute Name
  -s SAMPLES, --samples=SAMPLES
                        Samples Name
  -t TYPES, --types=TYPES
                        Project Types
  -p PROJECT, --project=PROJECT
                        Project ID
  -C CONTRACT, --contract=CONTRACT
                        Contract ID
  -k KEEPRMD, --keepRMD=KEEPRMD
                        KeepRMD or Not
```
