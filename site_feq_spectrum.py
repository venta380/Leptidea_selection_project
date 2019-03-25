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


def join_data_bases(lists):
    for i in range(1,len(lists)):
        j=i-1
        if i == 1:
            new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        else:
            new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        nextone=new_2
    return nextone


def get_site_freq_spectrum(goods,codon_df_1):
    codon_df_1_temp=pandas.merge(goods,codon_df_1, on=['CHROM','POS'], how='inner')
    codon_df_1_temp.minor_freq=codon_df_1_temp.minor_freq.round(2)
    spectrum=codon_df_1_temp.groupby('minor_freq')['minor_freq'].count()
    return spectrum


def output_good_sites_fequency_new(freq_file1,chr_filter):
    header=['CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'al_1_', 'al_2_','al_3_','al_4_']
    a=pandas.read_table(freq_file1, names=header, engine='python',skiprows=1)
    merged=a[(a['N_CHR'] >= int(chr_filter))]
    bads=[]
    goods=[]
    new_output=pandas.DataFrame(columns=['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq'])
    temp2=merged[merged['N_ALLELES']<=2]
    temp2['minor_allele']=''
    temp2['minor_freq']=0
    temp2['minor_allele'][temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float)]=temp2.al_2_.str[0]
    temp2['minor_allele'][temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float)]=temp2.al_1_.str[0]
    temp2['minor_freq'][temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float)]=temp2.al_2_.str[2:].astype(float)
    temp2['minor_freq'][temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float)]=temp2.al_1_.str[2:].astype(float)
    #
    temp2['major_allele']=''
    temp2['major_freq']=0
    #temp2.assign(major_allele=(temp2.al_1_.str[0]).where(float(temp2.al_1_.str[2:)) > float(temp2.al_2_.str[2:))), 0))
    temp2['major_allele'][temp2.al_2_.str[0] == temp2.minor_allele ]=temp2.al_1_.str[0]
    temp2['major_freq'][temp2.al_2_.str[0] == temp2.minor_allele ]=temp2.al_1_.str[2:].astype(float)
    temp2['major_allele'][temp2.al_1_.str[0] == temp2.minor_allele ]=temp2.al_2_.str[0]
    temp2['major_freq'][temp2.al_1_.str[0] == temp2.minor_allele ]=temp2.al_2_.str[2:].astype(float)
    temp2=temp2[['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq']]
    temp2=temp2[temp2.minor_freq>=0]
    new_output=new_output.append(temp2, ignore_index=True)
    temp2=merged[merged['N_ALLELES']==3]
    temp2['minor_allele']=''
    temp2['minor_freq']=0
    temp2['major_allele']=''
    temp2['major_freq']=0
    temp2['minor_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[0]
    temp2['major_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['minor_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['major_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[0]
    #
    temp2['minor_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['minor_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    #
    temp2['minor_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    temp2['major_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['minor_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['major_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    #
    temp2['minor_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['minor_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    #
    temp2['minor_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    temp2['major_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[0]
    temp2['minor_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[0]
    temp2['major_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    #
    temp2['minor_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    temp2['major_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    temp2['minor_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    temp2['major_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    #
    temp2=temp2[(temp2['minor_allele']!='') & (temp2['major_allele']!='')]
    temp2=temp2[['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq']]
    temp2=temp2[temp2.minor_freq>=0]
    new_output=new_output.append(temp2, ignore_index=True)
    #temp2=merged[merged['N_ALLELES'] ==4 ]
    #if len(temp2) > 0:
    #   for row in temp2.itertuples():
    #       allele={}
    #       freq_list=[float(row.al_1_[2:]), float(row.al_2_[2:]),float(row.al_3_[2:]), float(row.al_4_[2:])]
    #       if 0.0 in freq_list:
    #           if float(row.al_1_[2:]) not in [0.0]: allele[row.al_1_[0]]=float(row.al_1_[2:])
    #           if float(row.al_2_[2:]) not in [0.0]: allele[row.al_2_[0]]=float(row.al_2_[2:])
    #           if float(row.al_3_[2:]) not in [0.0]: allele[row.al_3_[0]]=float(row.al_3_[2:])
    #           if float(row.al_4_[2:]) not in [0.0]: allele[row.al_4_[0]]=float(row.al_4_[2:])
    #       else:
    #           bads.append(row)
    #       if  ((len(allele.values()) == 2) and min(allele.values()) != 0.0):
    #           temp={}
    #           temp['CHROM']=[row.CHROM]
    #           temp['POS']=[row.POS]
    #           temp['minor_allele']=[min(allele, key=lambda k: allele[k])]
    #           temp['minor_freq']=[allele[min(allele, key=lambda k: allele[k])]]
    #           temp['major_allele']=[max(allele, key=lambda k: allele[k])]
    #           temp['major_freq']=[allele[max(allele, key=lambda k: allele[k])]]
    #           temp_df = pandas.DataFrame(temp)
    #           new_output=new_output.append(temp_df, ignore_index=True)
    #new_output['BIN_START']=(np.floor(new_output['POS']/100000)*100000)+1
    #new_output['BIN_END']=new_output['BIN_START']+(100000-1)
    new_output['site_pi']=0.0
    total_comb=(int(chr_filter)*(int(chr_filter)-1.0))/2.0
    new_output['site_pi']=((new_output.minor_freq*int(chr_filter))*(new_output.major_freq*int(chr_filter)))/total_comb
    new_output['BIN_START']=0
    new_output['BIN_END']=0
    new_output['BIN_START']=(np.floor(new_output['POS']/100000)*100000)+1
    new_output['BIN_END']=new_output['BIN_START']+(100000-1)
    return new_output




introns=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/introns.csv',names=['POS','CHROM'], sep=' ',header=None)
codon_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_1.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_2.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_3.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_4d=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_4d.csv',names=['CHROM','POS'], sep=' ',header=None)

#goods1=output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica_freq.frq",40)
#goods2=output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/reali_freq.frq",20)
#goods3=output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis_freq.frq",60)



irish_juvernica_freq = output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_freq.frq", 10)
kazak_juvernica_freq = output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_freq.frq", 10)
spanish_reali_freq = output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_freq.frq", 10)
swe_sin_allele_freq = output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/swe_sin_allele_freq.frq", 10)
kaz_sin_freq = output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq", 10)
spanish_sinapis_freq = output_good_sites_fequency_new("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_freq.frq", 10)

hm=3.54774

######################ireJuv

irish_juvernica_1_sf=get_site_freq_spectrum(irish_juvernica_freq, codon_df_1)
irish_juvernica_2_sf=get_site_freq_spectrum(irish_juvernica_freq, codon_df_2)
irish_juvernica_3_sf=get_site_freq_spectrum(irish_juvernica_freq, codon_df_3)
irish_juvernica_4d_sf=get_site_freq_spectrum(irish_juvernica_freq, codon_df_4d)


irish_juvernica_theta=[(i/hm) for i in list(irish_juvernica_4d_sf)[1:]]
irish_juvernica_expected_numbers=[(list(irish_juvernica_4d_sf)[1:][i-1])/float(i)  for i in range(1,11) ]
irish_juvernica_expected_total=sum(irish_juvernica_expected_numbers)
irish_juvernica_expected_proportions=[(i/irish_juvernica_expected_total) for i in list(irish_juvernica_expected_numbers)]
irish_juvernica_expected_proportions.insert(0,0)

irish_juvernica=pandas.DataFrame()
irish_juvernica['LJire_1']= [i/float(sum(list(irish_juvernica_1_sf)[1:])) for i in list(irish_juvernica_1_sf)]
irish_juvernica['LJire_2']= [i/float(sum(list(irish_juvernica_2_sf)[1:])) for i in list(irish_juvernica_2_sf)]
irish_juvernica['LJire_3']= [i/float(sum(list(irish_juvernica_3_sf)[1:])) for i in list(irish_juvernica_3_sf)]
irish_juvernica['LJire_4D']=[i/float(sum(list(irish_juvernica_4d_sf)[1:])) for i in list(irish_juvernica_4d_sf)]
irish_juvernica['minor_freq']=irish_juvernica_1_sf.index
irish_juvernica['LJire_4D_expcted']=irish_juvernica_expected_proportions


######################KazJuv


kazak_juvernica_1_sf=get_site_freq_spectrum(kazak_juvernica_freq, codon_df_1)
kazak_juvernica_2_sf=get_site_freq_spectrum(kazak_juvernica_freq, codon_df_2)
kazak_juvernica_3_sf=get_site_freq_spectrum(kazak_juvernica_freq, codon_df_3)
kazak_juvernica_4d_sf=get_site_freq_spectrum(kazak_juvernica_freq, codon_df_4d)


kazak_juvernica_theta=[(i/hm) for i in list(kazak_juvernica_4d_sf)[1:]]
kazak_juvernica_expected_numbers=[(list(kazak_juvernica_4d_sf)[1:][i-1])/float(i)  for i in range(1,11) ]
kazak_juvernica_expected_total=sum(kazak_juvernica_expected_numbers)
kazak_juvernica_expected_proportions=[(i/kazak_juvernica_expected_total) for i in list(kazak_juvernica_expected_numbers)]
kazak_juvernica_expected_proportions.insert(0,0)




kazak_juvernica=pandas.DataFrame()
kazak_juvernica['LJkaz_1']= [i/float(sum(list(kazak_juvernica_1_sf)[1:])) for i in list(kazak_juvernica_1_sf)]
kazak_juvernica['LJkaz_2']= [i/float(sum(list(kazak_juvernica_2_sf)[1:])) for i in list(kazak_juvernica_2_sf)]
kazak_juvernica['LJkaz_3']= [i/float(sum(list(kazak_juvernica_3_sf)[1:])) for i in list(kazak_juvernica_3_sf)]
kazak_juvernica['LJkaz_4D']=[i/float(sum(list(kazak_juvernica_4d_sf)[1:])) for i in list(kazak_juvernica_4d_sf)]
kazak_juvernica['minor_freq']=kazak_juvernica_1_sf.index
kazak_juvernica['LJkaz_4D_expcted']=kazak_juvernica_expected_proportions


######################Sparea

spanish_reali_1_sf=get_site_freq_spectrum(spanish_reali_freq, codon_df_1)
spanish_reali_2_sf=get_site_freq_spectrum(spanish_reali_freq, codon_df_2)
spanish_reali_3_sf=get_site_freq_spectrum(spanish_reali_freq, codon_df_3)
spanish_reali_4d_sf=get_site_freq_spectrum(spanish_reali_freq, codon_df_4d)

spanish_reali_theta=[(i/hm) for i in list(spanish_reali_4d_sf)[1:]]
spanish_reali_expected_numbers=[(list(spanish_reali_4d_sf)[1:][i-1])/float(i)  for i in range(1,11) ]
spanish_reali_expected_total=sum(spanish_reali_expected_numbers)
spanish_reali_expected_proportions=[(i/spanish_reali_expected_total) for i in list(spanish_reali_expected_numbers)]
spanish_reali_expected_proportions.insert(0,0)


spanish_reali=pandas.DataFrame()
spanish_reali['Lrspa_1']= [i/float(sum(list(spanish_reali_1_sf)[1:])) for i in list(spanish_reali_1_sf)]
spanish_reali['Lrspa_2']= [i/float(sum(list(spanish_reali_2_sf)[1:])) for i in list(spanish_reali_2_sf)]
spanish_reali['Lrspa_3']= [i/float(sum(list(spanish_reali_3_sf)[1:])) for i in list(spanish_reali_3_sf)]
spanish_reali['Lrspa_4D']=[i/float(sum(list(spanish_reali_4d_sf)[1:])) for i in list(spanish_reali_4d_sf)]
spanish_reali['minor_freq']=spanish_reali_1_sf.index
spanish_reali['Lrspa_4D_expcted']=spanish_reali_expected_proportions

######################Swesin


#group sites into codon positions
swe_sin_allele_1_sf=get_site_freq_spectrum(swe_sin_allele_freq, codon_df_1)
swe_sin_allele_2_sf=get_site_freq_spectrum(swe_sin_allele_freq, codon_df_2)
swe_sin_allele_3_sf=get_site_freq_spectrum(swe_sin_allele_freq, codon_df_3)
swe_sin_allele_4d_sf=get_site_freq_spectrum(swe_sin_allele_freq, codon_df_4d)


#callculate theta
swe_sin_allele_theta=[(i/hm) for i in list(swe_sin_allele_4d_sf)[1:]]

#use theta to callculate the expected prortion of 4D sites under nutrality in each catogory of minor allele freqncies (based on Fu, 1995)
swe_sin_allele_expected_numbers=[(list(swe_sin_allele_4d_sf)[1:][i-1])/float(i)  for i in range(1,11) ]
swe_sin_allele_expected_total=sum(swe_sin_allele_expected_numbers)
swe_sin_allele_expected_proportions=[(i/swe_sin_allele_expected_total) for i in list(swe_sin_allele_expected_numbers)]
swe_sin_allele_expected_proportions.insert(0,0)

swe_sin_allele=pandas.DataFrame()
swe_sin_allele['Lsswe_1']= [i/float(sum(list(swe_sin_allele_1_sf)[1:])) for i in list(swe_sin_allele_1_sf)]
swe_sin_allele['Lsswe_2']= [i/float(sum(list(swe_sin_allele_2_sf)[1:])) for i in list(swe_sin_allele_2_sf)]
swe_sin_allele['Lsswe_3']= [i/float(sum(list(swe_sin_allele_3_sf)[1:])) for i in list(swe_sin_allele_3_sf)]
swe_sin_allele['Lsswe_4D']=[i/float(sum(list(swe_sin_allele_4d_sf)[1:])) for i in list(swe_sin_allele_4d_sf)]
swe_sin_allele['minor_freq']=swe_sin_allele_1_sf.index
swe_sin_allele['Lsswe_4D_expcted']=swe_sin_allele_expected_proportions


######################Kazsin


kaz_sin_1_sf=get_site_freq_spectrum(kaz_sin_freq, codon_df_1)
kaz_sin_2_sf=get_site_freq_spectrum(kaz_sin_freq, codon_df_2)
kaz_sin_3_sf=get_site_freq_spectrum(kaz_sin_freq, codon_df_3)
kaz_sin_4d_sf=get_site_freq_spectrum(kaz_sin_freq, codon_df_4d)

kaz_sin_theta=[(i/hm) for i in list(kaz_sin_4d_sf)[1:]]
kaz_sin_expected_numbers=[(list(kaz_sin_4d_sf)[1:][i-1])/float(i)  for i in range(1,11) ]
kaz_sin_expected_total=sum(kaz_sin_expected_numbers)
kaz_sin_expected_proportions=[(i/kaz_sin_expected_total) for i in list(kaz_sin_expected_numbers)]
kaz_sin_expected_proportions.insert(0,0)



kaz_sin=pandas.DataFrame()
kaz_sin['Lskaz_1']= [i/float(sum(list(kaz_sin_1_sf)[1:])) for i in list(kaz_sin_1_sf)]
kaz_sin['Lskaz_2']= [i/float(sum(list(kaz_sin_2_sf)[1:])) for i in list(kaz_sin_2_sf)]
kaz_sin['Lskaz_3']= [i/float(sum(list(kaz_sin_3_sf)[1:])) for i in list(kaz_sin_3_sf)]
kaz_sin['Lskaz_4D']=[i/float(sum(list(kaz_sin_4d_sf)[1:])) for i in list(kaz_sin_4d_sf)]
kaz_sin['minor_freq']=kaz_sin_1_sf.index
kaz_sin['Lskaz_4D_expcted']=kaz_sin_expected_proportions


######################Spasin


spanish_sinapis_1_sf=get_site_freq_spectrum(spanish_sinapis_freq, codon_df_1)
spanish_sinapis_2_sf=get_site_freq_spectrum(spanish_sinapis_freq, codon_df_2)
spanish_sinapis_3_sf=get_site_freq_spectrum(spanish_sinapis_freq, codon_df_3)
spanish_sinapis_4d_sf=get_site_freq_spectrum(spanish_sinapis_freq, codon_df_4d)


spanish_sinapis_theta=[(i/hm) for i in list(spanish_sinapis_4d_sf)[1:]]
spanish_sinapis_expected_numbers=[(list(spanish_sinapis_4d_sf)[1:][i-1])/float(i)  for i in range(1,11) ]
spanish_sinapis_expected_total=sum(spanish_sinapis_expected_numbers)
spanish_sinapis_expected_proportions=[(i/spanish_sinapis_expected_total) for i in list(spanish_sinapis_expected_numbers)]
spanish_sinapis_expected_proportions.insert(0,0)


spanish_sinapis=pandas.DataFrame()
spanish_sinapis['Lsspa_1']= [i/float(sum(list(spanish_sinapis_1_sf)[1:])) for i in list(spanish_sinapis_1_sf)]
spanish_sinapis['Lsspa_2']= [i/float(sum(list(spanish_sinapis_2_sf)[1:])) for i in list(spanish_sinapis_2_sf)]
spanish_sinapis['Lsspa_3']= [i/float(sum(list(spanish_sinapis_3_sf)[1:])) for i in list(spanish_sinapis_3_sf)]
spanish_sinapis['Lsspa_4D']=[i/float(sum(list(spanish_sinapis_4d_sf)[1:])) for i in list(spanish_sinapis_4d_sf)]
spanish_sinapis['minor_freq']=spanish_sinapis_1_sf.index
spanish_sinapis['Lsspa_4D_expcted']=spanish_sinapis_expected_proportions




one_=pandas.merge(irish_juvernica, kazak_juvernica, on='minor_freq', how='inner')
one_=pandas.merge(one_, spanish_reali, on='minor_freq', how='inner')
one_=pandas.merge(one_, swe_sin_allele, on='minor_freq', how='inner')
one_=pandas.merge(one_, kaz_sin, on='minor_freq', how='inner')
one_=pandas.merge(one_, spanish_sinapis, on='minor_freq', how='inner')
final=one_.iloc[1:]

final_4D=final[['LJire_4D_expcted','LJire_4D','LJkaz_4D_expcted','LJkaz_4D','Lrspa_4D_expcted','Lrspa_4D','Lsswe_4D_expcted','Lsswe_4D','Lskaz_4D_expcted','Lskaz_4D','Lsspa_4D_expcted','Lsspa_4D']]

final_4D[['LJire_4D','LJkaz_4D','Lrspa_4D','Lsswe_4D','Lskaz_4D','Lsspa_4D']].plot.bar()



#final_4D[['LJire_4D_expcted','LJkaz_4D_expcted','Lrspa_4D_expcted','Lsswe_4D_expcted','Lskaz_4D_expcted','Lsspa_4D_expcted']].plot.bar()

y_pos = np.arange(len(final_4D[['LJire_4D_expcted','LJkaz_4D_expcted','Lrspa_4D_expcted','Lsswe_4D_expcted','Lskaz_4D_expcted','Lsspa_4D_expcted']]))

plt.plot(final_4D[['LJire_4D_expcted','LJkaz_4D_expcted','Lrspa_4D_expcted','Lsswe_4D_expcted','Lskaz_4D_expcted','Lsspa_4D_expcted'], marker="D", linestyle="", alpha=0.8, color="r")

#plt.xticks(y_pos, final_4D[['LJire_4D_expcted','LJkaz_4D_expcted','Lrspa_4D_expcted','Lsswe_4D_expcted','Lskaz_4D_expcted','Lsspa_4D_expcted'])




#######################################################################



juvernica_1_sf=get_site_freq_spectrum(goods1, codon_df_1)
juvernica_2_sf=get_site_freq_spectrum(goods1, codon_df_2)
juvernica_3_sf=get_site_freq_spectrum(goods1, codon_df_3)
juvernica_4d_sf=get_site_freq_spectrum(goods1, codon_df_4d)

juvernica=pandas.DataFrame()
juvernica['J_1']= [i/float(len(goods1)) for i in list(juvernica_1_sf)]
juvernica['J_2']= [i/float(len(goods1)) for i in list(juvernica_2_sf)]
juvernica['J_3']= [i/float(len(goods1)) for i in list(juvernica_3_sf)]
juvernica['J_4D']=[i/float(len(goods1)) for i in list(juvernica_4d_sf)]
juvernica['minor_freq']=juvernica_1_sf.index



reali_1_sf=get_site_freq_spectrum(goods2, codon_df_1)
reali_2_sf=get_site_freq_spectrum(goods2, codon_df_2)
reali_3_sf=get_site_freq_spectrum(goods2, codon_df_3)
reali_4d_sf=get_site_freq_spectrum(goods2, codon_df_4d)

reali=pandas.DataFrame()
reali['R_1']=[i/float(len(goods2)) for i in list(reali_1_sf)]
reali['R_2']=[i/float(len(goods2)) for i in list(reali_2_sf)]
reali['R_3']=[i/float(len(goods2)) for i in list(reali_3_sf)]
reali['R_4D']=[i/float(len(goods2)) for i in list(reali_4d_sf)]
reali['minor_freq']=reali_1_sf.index



sinapis_1_sf=get_site_freq_spectrum(goods3, codon_df_1)
sinapis_2_sf=get_site_freq_spectrum(goods3, codon_df_2)
sinapis_3_sf=get_site_freq_spectrum(goods3, codon_df_3)
sinapis_4d_sf=get_site_freq_spectrum(goods3, codon_df_4d)

sinapis=pandas.DataFrame()
sinapis['S_1']=[i/float(len(goods3)) for i in list(sinapis_1_sf)]
sinapis['S_2']=[i/float(len(goods3)) for i in list(sinapis_2_sf)]
sinapis['S_3']=[i/float(len(goods3)) for i in list(sinapis_3_sf)]
sinapis['S_4D']=[i/float(len(goods3)) for i in list(sinapis_4d_sf)]
sinapis['minor_freq']=sinapis_1_sf.index


one_=pandas.merge(juvernica, reali, on='minor_freq', how='inner')
final=pandas.merge(one_, sinapis, on='minor_freq', how='inner')


