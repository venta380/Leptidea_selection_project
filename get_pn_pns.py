import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
from pandas.tools.plotting import scatter_matrix
import Bio.Data.CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
import argparse

table=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/"+str(sys.argv[1]))

table['Pin']=table['Pi_nonsyn']/table['nonsyn_sites']
table['Pis']=table['Pi_syn']/table['syn_sites']

table=table.drop_duplicates()

table.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/'+str(sys.argv[1])+'.pnps',index=False, header=True, sep=' ')


