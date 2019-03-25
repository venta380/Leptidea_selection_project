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


goods_irish_juvernica=personal_popgen.output_good_sites_fequency("/proj/b2014034/nobackup/POPULATION_RESEQ/dnds/temp_files/irish_juvernica.frq",20)
goods_kazak_juvernica=personal_popgen.output_good_sites_fequency("/proj/b2014034/nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_juvernica.frq",20)
goods_kazak_sinapis=personal_popgen.output_good_sites_fequency("/proj/b2014034/nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_sinapis.frq",20)
goods_spanish_reali=personal_popgen.output_good_sites_fequency("/proj/b2014034/nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_reali.frq",20)
goods_spanish_sinapis=personal_popgen.output_good_sites_fequency("/proj/b2014034/nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_sinapis.frq",20)
goods_Swedish_sinapis=personal_popgen.output_good_sites_fequency("/proj/b2014034/nobackup/POPULATION_RESEQ/dnds/temp_files/Swedish_sinapis.frq",20)

goods_irish_juvernica=goods_irish_juvernica[['CHROM','POS','minor_freq']]
goods_kazak_juvernica=goods_kazak_juvernica[['CHROM','POS','minor_freq']]
goods_kazak_sinapis=goods_kazak_sinapis[['CHROM','POS','minor_freq']]
goods_spanish_reali=goods_spanish_reali[['CHROM','POS','minor_freq']]
goods_spanish_sinapis=goods_spanish_sinapis[['CHROM','POS','minor_freq']]
goods_Swedish_sinapis=goods_Swedish_sinapis[['CHROM','POS','minor_freq']]


goods_irish_juvernica=goods_irish_juvernica.rename(columns={'minor_freq': 'minor_freq_irish_juvernica'})
goods_kazak_juvernica=goods_kazak_juvernica.rename(columns={'minor_freq': 'minor_freq_kazak_juvernica'})
goods_kazak_sinapis=goods_kazak_sinapis.rename(columns={'minor_freq': 'minor_freq_kazak_sinapis'})
goods_spanish_reali=goods_spanish_reali.rename(columns={'minor_freq': 'minor_freq_spanish_reali'})
goods_spanish_sinapis=goods_spanish_sinapis.rename(columns={'minor_freq': 'minor_freq_spanish_sinapis'})
goods_Swedish_sinapis=goods_Swedish_sinapis.rename(columns={'minor_freq': 'minor_freq_Swedish_sinapis'})



def join_raw_data_base(lists):
	for i in range(1,len(lists)):
		j=i-1
		if i == 1:
			new_2=pandas.merge(lists[j], lists[i], on=['CHROM','POS'], how='outer', suffixes=('_1', '_2'))
		else:
			new_2=pandas.merge(nextone, lists[i], on=['CHROM','POS'], how='outer')

		nextone=new_2
	return nextone

df4=join_raw_data_base([goods_irish_juvernica,goods_kazak_juvernica,goods_kazak_sinapis,goods_spanish_reali,goods_spanish_sinapis,goods_Swedish_sinapis])



#df4[["minor_freq_irish_juvernica","minor_freq_kazak_juvernica","minor_freq_kazak_sinapis","minor_freq_spanish_reali","minor_freq_spanish_sinapis","minor_freq_Swedish_sinapis"]].plot.hist(bins=10, alpha=0.5)


result = df4[["minor_freq_irish_juvernica","minor_freq_kazak_juvernica","minor_freq_kazak_sinapis","minor_freq_spanish_reali","minor_freq_spanish_sinapis","minor_freq_Swedish_sinapis"]].apply(pandas.value_counts)
result.plot(legend=True, color=['#C8C800','#006400','#FF8C00','#0000FF','#FF0000','#C04000'])


