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
import argparse

def join_raw_data_base(lists, on_what, how_is_it):
        for i in range(1,len(lists)):
                j=i-1
                if i == 1:
                        new_2=pandas.merge(lists[j], lists[i], on=on_what, how=how_is_it)
                else:
                        new_2=pandas.merge(nextone, lists[i], on=on_what, how=how_is_it)
                nextone=new_2
        return nextone



#

dnds_table=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/final_dnds_result_species")

dnds_table['Ds_L_sin']=dnds_table['ds'].str.split(':').str[1].str.split(',').str[0]
dnds_table['Ds_L_rea']=dnds_table['ds'].str.split(':').str[2].str.split(')').str[0]
dnds_table['Ds_L_juv']=dnds_table['ds'].str.split(':').str[4].str.split(')').str[0]
dnds_table['Dn_L_sin']=dnds_table['dn'].str.split(':').str[1].str.split(',').str[0]
dnds_table['Dn_L_rea']=dnds_table['dn'].str.split(':').str[2].str.split(')').str[0]
dnds_table['Dn_L_juv']=dnds_table['dn'].str.split(':').str[4].str.split(')').str[0]

dnds_table=dnds_table[['gene_name','Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']]

dnds_table[['Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']]=dnds_table[['Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']].astype('float').round(4)


dnds_table['W_L_sin']=dnds_table['Dn_L_sin']/dnds_table['Ds_L_sin']
dnds_table['W_L_rea']=dnds_table['Dn_L_rea']/dnds_table['Ds_L_rea']
dnds_table['W_L_juv']=dnds_table['Dn_L_juv']/dnds_table['Ds_L_juv']

dnds_table['W_L_sin']=dnds_table['Dn_L_sin']/dnds_table['Ds_L_sin']
dnds_table['W_L_rea']=dnds_table['Dn_L_rea']/dnds_table['Ds_L_rea']
dnds_table['W_L_juv']=dnds_table['Dn_L_juv']/dnds_table['Ds_L_juv']


dnds_table['W_L_sin']=dnds_table['W_L_sin'].replace(np.inf, dnds_table.Dn_L_sin)
dnds_table['W_L_rea']=dnds_table['W_L_rea'].replace(np.inf, dnds_table.Dn_L_rea)
dnds_table['W_L_juv']=dnds_table['W_L_juv'].replace(np.inf, dnds_table.Dn_L_juv)


dnds_table[(dnds_table['W_L_sin'] >1.0) & (dnds_table['W_L_sin'] < 10.0)] 
dnds_table[(dnds_table['W_L_rea'] >1.0) & (dnds_table['W_L_rea'] < 10.0)] 
dnds_table[(dnds_table['W_L_juv'] >1.0) & (dnds_table['W_L_juv'] < 10.0)] 

print dnds_table['W_L_sin'].mean()
print dnds_table['W_L_rea'].mean()
print dnds_table['W_L_juv'].mean()

### load PNPS

fianl_df_1=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_irish_juvernica.csv")
fianl_df_2=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_juvernica.csv")
fianl_df_3=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_kazak_juvernica.csv")
fianl_df_4=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_kazak_sinapis.csv")
fianl_df_5=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_sinapis.csv")
fianl_df_6=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_spanish_reali.csv")
fianl_df_7=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_spanish_sinapis.csv")
fianl_df_8=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Pin_Pis_Swedish_sinapis.csv")



#Dn_DS_for inuvedual comparisions



dn_ds=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/dn_ds.csv")

dn_ds=dn_ds.rename(columns={'gene': 'geneID'})

#######psoitive selection ###############
dn_ds[(dn_ds['Dn_ds_sin_reali'] >1.0) & (dn_ds['Dn_ds_sin_reali'] < 10.0)] 
dn_ds[(dn_ds['Dn_ds_sin_juv'] >1.0) & (dn_ds['Dn_ds_sin_juv'] < 10.0)] 
dn_ds[(dn_ds['Dn_ds_juv_reali'] >1.0) & (dn_ds['Dn_ds_juv_reali'] < 10.0)] 



join_raw_data_base([fianl_df_1,fianl_df_2,fianl_df_3,fianl_df_4,fianl_df_5,fianl_df_6,fianl_df_7,fianl_df_8, dn_ds], 'geneID', 'outer')

