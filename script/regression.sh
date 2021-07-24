#!/bin/bash
#PBS -q ngs192G
#PBS -P MST109178
#PBS -W group_list=MST109178
#PBS -N Regression
#PBS -o /work1/viviane1695/tbb1496/updates/run.log/regression_out
#PBS -e /work1/viviane1695/tbb1496/updates/run.log/regression_err
#PBS -M vivianchen1695@gmail.com
#PBS -m a

export PATH=/pkg/biology/Python/Python3_default/bin:$PATH
cd /work1/viviane1695/tbb1496/updates/script

python3 TWB1496_EAS_regression.py
python3 TWB1496_official_regression.py
python3 TWB1496_EAS_official_subplot.py
python3 4subplot_by_filter.py


