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




sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})

pwd='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/'
os.chdir(pwd)



def load_files(pop):
    Pop_1=pop
    fianl_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf1/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf2/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf3/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_4=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf4/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_5=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf5/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_6=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf6/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_7=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf7/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_8=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf8/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_9=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf9/'+Pop_1+'_1_final_df.csv',index_col=False)
    fianl_df_10=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf10/'+Pop_1+'_1_final_df.csv',index_col=False)
    final_df_final=pandas.concat([fianl_df_1,fianl_df_2,fianl_df_3,fianl_df_4,fianl_df_5,fianl_df_6,fianl_df_7,fianl_df_8,fianl_df_9,fianl_df_10],ignore_index=True)
    final_df_final=final_df_final.rename(columns={'site_pi_codon_1': Pop_1+'_site_pi_codon_1', 'site_pi_codon_2': Pop_1+'_site_pi_codon_2', 'site_pi_codon_3': Pop_1+'_site_pi_codon_3', 'site_pi_codon_4d': Pop_1+'_site_pi_codon_4d', 'site_pi_introns': Pop_1+'_site_pi_introns', 'site_pi_global': Pop_1+'_site_pi_global','sites_codon_1': Pop_1+'_sites_codon_1','sites_codon_2': Pop_1+'_sites_codon_2','sites_codon_3': Pop_1+'_sites_codon_3','sites_codon_4d': Pop_1+'_sites_codon_4d','sites_introns': Pop_1+'_sites_introns','sites_global': Pop_1+'_sites_global'})
    return final_df_final



pop_df_1=load_files('irish_juvernica')
pop_df_2=load_files('kazak_juvernica')
pop_df_3=load_files('spanish_reali')
pop_df_4=load_files('swe_sin_allele')
pop_df_5=load_files('kaz_sin')
pop_df_6=load_files('spanish_sinapis')

pop_df_7=load_files('sinapis')
pop_df_8=load_files('juvernica')
pop_df_9=load_files('reali')


final_stats=personal_popgen.join_raw_data_base([pop_df_1,pop_df_2,pop_df_3,pop_df_4,pop_df_5,pop_df_6,pop_df_7,pop_df_8,pop_df_9])


sns.boxplot(pop_df_1[["irish_juvernica_site_pi_codon_1","irish_juvernica_site_pi_codon_2","irish_juvernica_site_pi_codon_3","irish_juvernica_site_pi_introns","irish_juvernica_site_pi_codon_4d","irish_juvernica_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_2[["kazak_juvernica_site_pi_codon_1","kazak_juvernica_site_pi_codon_2","kazak_juvernica_site_pi_codon_3","kazak_juvernica_site_pi_introns","kazak_juvernica_site_pi_codon_4d","kazak_juvernica_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_3[["spanish_reali_site_pi_codon_1","spanish_reali_site_pi_codon_2","spanish_reali_site_pi_codon_3","spanish_reali_site_pi_introns","spanish_reali_site_pi_codon_4d","spanish_reali_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_4[["swe_sin_allele_site_pi_codon_1","swe_sin_allele_site_pi_codon_2","swe_sin_allele_site_pi_codon_3","swe_sin_allele_site_pi_introns","swe_sin_allele_site_pi_codon_4d","swe_sin_allele_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_5[["kaz_sin_site_pi_codon_1","kaz_sin_site_pi_codon_2","kaz_sin_site_pi_codon_3","kaz_sin_site_pi_introns","kaz_sin_site_pi_codon_4d","kaz_sin_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_6[["spanish_sinapis_site_pi_codon_1","spanish_sinapis_site_pi_codon_2","spanish_sinapis_site_pi_codon_3","spanish_sinapis_site_pi_introns","spanish_sinapis_site_pi_codon_4d","spanish_sinapis_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))

sns.boxplot(pop_df_7[["sinapis_site_pi_codon_1","sinapis_site_pi_codon_2","sinapis_site_pi_codon_3","sinapis_site_pi_introns","sinapis_site_pi_codon_4d","sinapis_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_8[["juvernica_site_pi_codon_1","juvernica_site_pi_codon_2","juvernica_site_pi_codon_3","juvernica_site_pi_introns","juvernica_site_pi_codon_4d","juvernica_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))
sns.boxplot(pop_df_9[["reali_site_pi_codon_1","reali_site_pi_codon_2","reali_site_pi_codon_3","reali_site_pi_introns","reali_site_pi_codon_4d","reali_site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))


pi_colours=sns.color_palette(['#C8C800','#006400','#0000FF','#C04000','#FF8C00','#FF0000'])

sns.boxplot(final_stats[["irish_juvernica_site_pi_codon_1","kazak_juvernica_site_pi_codon_1","spanish_reali_site_pi_codon_1","swe_sin_allele_site_pi_codon_1","kaz_sin_site_pi_codon_1","spanish_sinapis_site_pi_codon_1","irish_juvernica_site_pi_codon_2","kazak_juvernica_site_pi_codon_2","spanish_reali_site_pi_codon_2","swe_sin_allele_site_pi_codon_2","kaz_sin_site_pi_codon_2","spanish_sinapis_site_pi_codon_2","irish_juvernica_site_pi_codon_3","kazak_juvernica_site_pi_codon_3","spanish_reali_site_pi_codon_3","swe_sin_allele_site_pi_codon_3","kaz_sin_site_pi_codon_3","spanish_sinapis_site_pi_codon_3","irish_juvernica_site_pi_introns","kazak_juvernica_site_pi_introns","spanish_reali_site_pi_introns","swe_sin_allele_site_pi_introns","kaz_sin_site_pi_introns","spanish_sinapis_site_pi_introns"]],showfliers=False,palette=pi_colours).set(ylim=(0,0.0105))




pop_df_4.mean()

print final_pi_db_codon1.site_pi.mean()
print final_pi_db_codon2.site_pi.mean()
print final_pi_db_codon3.site_pi.mean()
print final_pi_db_introns.site_pi.mean()


#fianl_df=pandas.merge(final_pi_db_codon1,final_pi_db_codon2, on=['CHROM','BIN_START','BIN_END'], how='inner',suffixes=('_x', '_y'))
#fianl_df=fianl_df.rename(columns={'site_pi_x': 'site_pi_codon_1', 'sites_x':'sites_codon_1','site_pi_y': 'site_pi_codon_2', 'sites_y':'sites_codon_2','site_pi_x': 'site_pi_codon_1', 'sites_x':'sites_codon_1'})
#fianl_df=pandas.merge(fianl_df,final_pi_db_codon3, on=['CHROM','BIN_START','BIN_END'], how='inner',suffixes=('_x', '_y'))
#fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon_3', 'sites':'sites_codon_3'})
#fianl_df=pandas.merge(fianl_df,final_pi_db_introns, on=['CHROM','BIN_START','BIN_END'], how='inner',suffixes=('_x', '_y'))
#fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_introns', 'sites':'sites_introns'})
#
#
#sns.boxplot(fianl_df[["site_pi_codon_1","site_pi_codon_2","site_pi_codon_3","site_pi_introns"]],showfliers=False).set(ylim=(0,0.0105))





/home/venkat/bin/snpgenie/snpgenie.pl --sepfiles --minfreq=0.01 --vcfformat 4 --snpreport='test_vcf2.vcf' --fastafile='test.fasta' --gtffile='test_gff_converted.gtf'


print pop_df_1.mean()
print pop_df_2.mean()
print pop_df_3.mean()
print pop_df_4.mean()
print pop_df_5.mean()
print pop_df_6.mean()

