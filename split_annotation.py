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




head_stuff=['scaffold','source','feature','start','end','extra','-','indexing','infor']
#CDSs=pandas.read_table("/home/venkat/bin/snpgenie/test/",skiprows=1, names=head_stuff)

CDS_fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/fasta/cds.fa')



annotation=pandas.read_table("/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/gff/gene-builds/leptidea_sinapis_rc1.gff",skiprows=1, names=head_stuff)
annotation=annotation[['scaffold','source','feature','start','end','-','indexing','infor']]

annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='+'), 'start']-1
annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='+'), 'start']-1

annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='-'), 'start']-1
annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='-'), 'start']-1

annotation.loc[(annotation['indexing']=='0')  & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='0') & (annotation['-']=='-'), 'start']-1
annotation.loc[(annotation['indexing']=='0')  & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='0') & (annotation['-']=='+'), 'start']-1

#annotation_gene=annotation[(annotation.feature=='gene') & (annotation.infor.str[31:37]==';Name=')]
#annotation_gene=annotation[(annotation.infor.str.split(';').str[0].str[3:].isin(CDS_fasta.keys()))]
#
#annotation_gene=annotation[(annotation.infor.str.split(';').str[0].str[3:].isin(CDS_fasta.keys()))].infor.str.split(';').str[1].str[7:]
#
#annotation_gene=annotation[(annotation.infor.str.split(';').str[0].str[3:].isin(list(annotation_gene)))]

annotation_gene=annotation[(annotation.feature=='gene')]


def index_marks(nrows, chunk_size):
    return range(1 * chunk_size, (nrows // chunk_size + 1) * chunk_size, chunk_size)

def split(dfm, chunk_size):
    indices = index_marks(dfm.shape[0], chunk_size)
    return np.split(dfm, indices)



chunks = split(annotation_gene, 2600)

chunks[0].to_csv('part_1.csv')
chunks[1].to_csv('part_2.csv')
chunks[2].to_csv('part_3.csv')
chunks[3].to_csv('part_4.csv')
chunks[4].to_csv('part_5.csv')
chunks[5].to_csv('part_6.csv')


cat codon_df_3_part_1.csv >>codon_df_3.csv
cat codon_df_4d_part_1.csv >>codon_df_4d.csv
cat codon_df_1_part_1.csv >>codon_df_1.csv
cat codon_df_2_part_1.csv >>codon_df_2.csv
cat introns_part_1.csv >>introns.csv

cat codon_df_3_part_2.csv >>codon_df_3.csv
cat codon_df_4d_part_2.csv >>codon_df_4d.csv
cat codon_df_1_part_2.csv >>codon_df_1.csv
cat codon_df_2_part_2.csv >>codon_df_2.csv
cat introns_part_2.csv >>introns.csv

cat codon_df_3_part_3.csv >>codon_df_3.csv
cat codon_df_4d_part_3.csv >>codon_df_4d.csv
cat codon_df_1_part_3.csv >>codon_df_1.csv
cat codon_df_2_part_3.csv >>codon_df_2.csv
cat introns_part_3.csv >>introns.csv

cat codon_df_3_part_4.csv >>codon_df_3.csv
cat codon_df_4d_part_4.csv >>codon_df_4d.csv
cat codon_df_1_part_4.csv >>codon_df_1.csv
cat codon_df_2_part_4.csv >>codon_df_2.csv
cat introns_part_4.csv >>introns.csv

cat codon_df_3_part_5.csv >>codon_df_3.csv
cat codon_df_4d_part_5.csv >>codon_df_4d.csv
cat codon_df_1_part_5.csv >>codon_df_1.csv
cat codon_df_2_part_5.csv >>codon_df_2.csv
cat introns_part_5.csv >>introns.csv

cat codon_df_3_part_6.csv >>codon_df_3.csv
cat codon_df_4d_part_6.csv >>codon_df_4d.csv
cat codon_df_1_part_6.csv >>codon_df_1.csv
cat codon_df_2_part_6.csv >>codon_df_2.csv
cat introns_part_6.csv >>introns.csv


-rw-rw-r-- 1 venkat uppstore2017185 3.1G Mar 20 09:35 introns.csv
-rw-rw-r-- 1 venkat uppstore2017185 186K Mar 20 09:35 start_stuff.csv
-rw-rw-r-- 1 venkat uppstore2017185  84M Mar 20 09:35 codon_df_1.csv
-rw-rw-r-- 1 venkat uppstore2017185  84M Mar 20 09:35 codon_df_2.csv
-rw-rw-r-- 1 venkat uppstore2017185  84M Mar 20 09:35 codon_df_3.csv
-rw-rw-r-- 1 venkat uppstore2017185  35M Mar 20 09:35 codon_df_4d.csv


-rw-rw-r-- 1 venkat uppstore2017185 121M May 22 19:10 codon_df_3.csv
-rw-rw-r-- 1 venkat uppstore2017185 50M May 22 19:10 codon_df_4d.csv
-rw-rw-r-- 1 venkat uppstore2017185 121M May 22 19:10 codon_df_1.csv
-rw-rw-r-- 1 venkat uppstore2017185 121M May 22 19:10 codon_df_2.csv
-rw-rw-r-- 1 venkat uppstore2017185 4.6G May 22 19:10 introns.csv
