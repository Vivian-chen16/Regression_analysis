# Regression_analysis
## Motivation
To explore the unveiling data in whole genome sequence and to see if any significant pattern exists between East asian and Taiwanese, we apply a statistical (linear equation) and visualized (scatter plot) manner in the linear regression model of allele frequencies.

## Features
* With scatter plot and regression line, the correlation of East Asisan(gnomAD) and Taiwanese(TWB1496/official data) can be visulaized.
* Take calling quality as one of metrics, in order to maximize the confidentiality when comparing different dataset.
* Provide comprehensive pharmacogenomics table for TWB1496 cohort.

## Technologies
This project is created with :
* Python version: Python 3.8
  * Pandas 3.3.2
  * numpy 3.3.2
  * Scipy 0.23.2
  * Matplotlib 3.3.2
  * Seaborn 0.11.1

## Installation
Use the package manager [pip](http://pip.pypa.io/en/stable/) to install.
```bash
$ pip install python
```

## How to use
* Subest gnomAD information we interested in annotated VCF 
  * Please check tbb_subset.sh
* Set up conditon for nonsynonymous variants
  * Please check tbb_cond.sh
  * Define nonsynonymous variant is qualify when
    (1) gene function is related to splicing
    (2) exonic variant function is nonsynonymous SNV or related to deletion, insertion, stoploss and stopgain.
* Data clean, define population difference and specific in afterwards regression analysis
  * Please check 

* Regression analysis
  * Please check Regression.sh
* PharGKB tbale
