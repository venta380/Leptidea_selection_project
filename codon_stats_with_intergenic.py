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
import argparse
from Bio.SeqUtils import GC


part=int(sys.argv[1])
#part=int(4)


dnds_table_windows_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/dnds_table_windows.csv')

dfList = np.array_split(dnds_table_windows_2, 6)

dnds_table_windows_2=dfList[part]


fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/GENOME_ASSEMBLY/assembly_updates/v1.4/N.Backstrom_leptidea.scf.1.4.fasta')

introns=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/introns.csv',names=['POS','CHROM'], sep=' ',header=None)
codon_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_1.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_2.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_3.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_4D=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_4d.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)




dnds_table_windows_2['GC_reference']=0.0
dnds_table_windows_2['GC_intron']=0.0
dnds_table_windows_2['GC_codon_1']=0.0
dnds_table_windows_2['GC_codon_2']=0.0
dnds_table_windows_2['GC_codon_3']=0.0
dnds_table_windows_2['GC_codon_4D']=0.0
dnds_table_windows_2['GC_intergenic']=0.0


for window_index, window in dnds_table_windows_2.iterrows():
    gc_fraction=GC(str(fasta[window.CHROM][int(window.BIN_START-1):int(window.BIN_END-1)]))
    sequnce=np.array(fasta[window.CHROM])
    introns_seq_1=introns[(introns['CHROM']==window.CHROM) & (introns['POS']>=window.BIN_START) & (introns['POS']<=window.BIN_END) ]
    introns_seq=list(introns_seq_1['POS']-1)
    introns_seq=sequnce[introns_seq]
    introns_seq=float(GC(''.join(introns_seq)))
    CD_1_seq_1=codon_df_1[(codon_df_1['CHROM']==window.CHROM) & (codon_df_1['POS']>=window.BIN_START) & (codon_df_1['POS']<=window.BIN_END) ]
    CD_1_seq=list(CD_1_seq_1['POS']-1)
    if len(CD_1_seq) != 0 and max(CD_1_seq) <= len(sequnce):
        CD_1_seq=sequnce[CD_1_seq]
        CD_1_seq=float(GC(''.join(CD_1_seq)))
    else:
        CD_1_seq=np.nan
    CD_2_seq_1=codon_df_2[(codon_df_2['CHROM']==window.CHROM) & (codon_df_2['POS']>=window.BIN_START) & (codon_df_2['POS']<=window.BIN_END) ]
    CD_2_seq=list(CD_2_seq_1['POS']-1)
    if len(CD_2_seq) != 0 and max(CD_2_seq) <= len(sequnce):
        CD_2_seq=sequnce[CD_2_seq]
        CD_2_seq=float(GC(''.join(CD_2_seq)))
    else:
        CD_2_seq=np.nan
    CD_3_seq_1=codon_df_3[(codon_df_3['CHROM']==window.CHROM) & (codon_df_3['POS']>=window.BIN_START) & (codon_df_3['POS']<=window.BIN_END) ]
    CD_3_seq=list(CD_3_seq_1['POS']-1)
    if len(CD_3_seq) != 0 and max(CD_3_seq) < len(sequnce):
        CD_3_seq=sequnce[CD_3_seq]
        CD_3_seq=float(GC(''.join(CD_3_seq)))
    else:
        CD_3_seq=np.nan
    CD_4D_seq_1=codon_df_4D[(codon_df_4D['CHROM']==window.CHROM) & (codon_df_4D['POS']>=window.BIN_START) & (codon_df_4D['POS']<=window.BIN_END) ]
    CD_4D_seq=list(CD_4D_seq_1['POS']-1)
    if len(CD_4D_seq) != 0 and max(CD_4D_seq) < len(sequnce):
        CD_4D_seq=sequnce[CD_4D_seq]
        CD_4D_seq=float(GC(''.join(CD_4D_seq)))
    else:
        CD_4D_seq=np.nan
    list_pos=range(int(window.BIN_START-1), int(window.BIN_END-1))
    merge=pandas.concat([introns_seq_1,CD_1_seq_1,CD_2_seq_1,CD_3_seq_1,CD_4D_seq_1],ignore_index=True)
    intergenic=main_list = np.setdiff1d(list_pos,list(merge.POS))
    intergenic_seq=list(intergenic-1)
    if len(intergenic_seq) != 0 and max(intergenic_seq) < len(sequnce):
        intergenic_seq=sequnce[intergenic_seq]
        intergenic_seq=float(GC(''.join(intergenic_seq)))
    else:
        intergenic_seq=np.nan
    #print str(gc_fraction)+ '   '+ window.CHROM + '     '+ str(window.BIN_START) + '        '+ str(window.BIN_END)
    dnds_table_windows_2['GC_reference'].loc[window_index]=gc_fraction
    dnds_table_windows_2['GC_intron'].loc[window_index]=introns_seq
    dnds_table_windows_2['GC_codon_1'].loc[window_index]=CD_1_seq
    dnds_table_windows_2['GC_codon_2'].loc[window_index]=CD_2_seq
    dnds_table_windows_2['GC_codon_3'].loc[window_index]=CD_3_seq
    dnds_table_windows_2['GC_codon_4D'].loc[window_index]=CD_4D_seq
    dnds_table_windows_2['GC_intergenic'].loc[window_index]=intergenic_seq
    print window_index
    sys.stdout.flush()





dnds_table_windows_2.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_'+str(part)+'.csv')






#! /bin/bash -l
#SBATCH -A snic2017-1-615
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00:00
#SBATCH -J call_GC
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL

cd /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test

python callculate_GC_per_window.py 0 &
python callculate_GC_per_window.py 1 &
python callculate_GC_per_window.py 2 &
python callculate_GC_per_window.py 3 &
python callculate_GC_per_window.py 4 &
python callculate_GC_per_window.py 5 &
wait

