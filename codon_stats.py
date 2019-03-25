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



def  assess_coverage_filter_7X_induv(list_1, cov, ind):
        headera=['CHROM', 'POS', 'COV0']
        keys=['COV0']
        list_all=[]
        for i in range(0,len(list_1)):
                if i == 0:
                        n=pandas.read_table(list_1[i], header=None, engine='c')
                        n.columns = headera
                        list_all=[n]
                else:
                        n=pandas.read_table(list_1[i], header=None, engine='c')
                        headera=['CHROM', 'POS', 'COV'+str(i)]
                        n.columns = headera
                        list_all.append(n)
                        keys.append('COV'+str(i))
                        nextone=personal_popgen.join_bam_coverage(list_all)
                        del list_all
                        list_all=[nextone]
        del list_all
        temp=nextone[keys]
        nextone=nextone[temp[temp >= cov].count(axis=1) >= ind]
        del temp
        return nextone[['CHROM','POS']]




def get_args():
    parser = argparse.ArgumentParser(description='''outputs nucleotide diversity from 1st,2nd,3rd codon and intronic sites for 100 kb windows

        -f freq file
        -p which part of the scafolds scaf1 scaf2 scaf3
        -pop pop file
        -c coverage filter for 10 induveduals
        -ind number of chromosomes (induvuduals in the population *2)
         ''')
    parser.add_argument('-f', '--freq', type=str, required=True)
    parser.add_argument('-p', '--part', type=str, required=True)
    parser.add_argument('-pop', '--population', type=str, required=True)
    parser.add_argument('-c', '--coverage', type=int, required=True)
    parser.add_argument('-ind', '--chr', type=float, required=True)
    return parser.parse_args()

def join_data_bases(lists):
    for i in range(1,len(lists)):
        j=i-1
        if i == 1:
            new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        else:
            new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        nextone=new_2
    return nextone


args = get_args()


print args.freq
print args.part
print args.population
print args.coverage
print args.chr


introns=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/introns.csv',names=['POS','CHROM'], sep=' ',header=None)
codon_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_1.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_2.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_3.csv',names=['CHROM','POS','major_allele'], sep=' ',header=None)
codon_df_4d=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_4d.csv',names=['CHROM','POS'], sep=' ',header=None)

introns=introns.drop_duplicates(subset=['POS','CHROM'])
codon_df_1=codon_df_1.drop_duplicates(subset=['POS','CHROM'])
codon_df_2=codon_df_2.drop_duplicates(subset=['POS','CHROM'])
codon_df_3=codon_df_3.drop_duplicates(subset=['POS','CHROM'])
codon_df_4d=codon_df_4d.drop_duplicates(subset=['POS','CHROM'])

#codon_df_1['BIN_START']=0
#codon_df_1['BIN_END']=0
#codon_df_2['BIN_START']=0
#codon_df_2['BIN_END']=0
#codon_df_3['BIN_START']=0
#codon_df_3['BIN_END']=0
#codon_df_4d['BIN_START']=0
#codon_df_4d['BIN_END']=0
#introns['BIN_START']=0
#introns['BIN_END']=0
#codon_df_1['BIN_START']=(np.floor(codon_df_1['POS']/100000)*100000)+1
#codon_df_1['BIN_END']=codon_df_1['BIN_START']+(100000-1)
#codon_df_2['BIN_START']=(np.floor(codon_df_2['POS']/100000)*100000)+1
#codon_df_2['BIN_END']=codon_df_2['BIN_START']+(100000-1)
#codon_df_3['BIN_START']=(np.floor(codon_df_3['POS']/100000)*100000)+1
#codon_df_3['BIN_END']=codon_df_3['BIN_START']+(100000-1)
#codon_df_4d['BIN_START']=(np.floor(codon_df_4d['POS']/100000)*100000)+1
#codon_df_4d['BIN_END']=codon_df_4d['BIN_START']+(100000-1)
#introns['BIN_START']=(np.floor(introns['POS']/100000)*100000)+1
#introns['BIN_END']=introns['BIN_START']+(100000-1)

codon_df_1=codon_df_1[['CHROM','POS']]
codon_df_2=codon_df_2[['CHROM','POS']]
codon_df_3=codon_df_3[['CHROM','POS']]
codon_df_4d=codon_df_4d[['CHROM','POS']]
introns=introns[['CHROM','POS']]


os.chdir('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/')
#coverage_list_pop1='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kaz_sin'
coverage_list_pop1=args.population
coverage_list_pop1=[stuff.strip() for stuff in open(coverage_list_pop1)]

stats_dir="/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/"

list_1=[]
for i in coverage_list_pop1:
        n=[stats_dir+"/"+args.part+"/"+stuff.strip()+".cov" for stuff in open(i)]
        list_1.append(assess_coverage_filter_7X_induv(n, args.coverage, args.chr/2))

#list_1=[]
#for i in coverage_list_pop1:
#        n=[stats_dir+"/"+'scaf7'+"/"+stuff.strip()+".cov" for stuff in open(i)]
#        list_1.append(assess_coverage_filter_7X_induv(n[1:3], 5, 2))



df_pass=reduce(lambda x, y: pandas.merge(x, y, on=['CHROM','POS'], how='inner'), list_1)
print "merged all"

del list_1


#df_merge = pandas.merge(df_pass, chromosome, on=['CHROM','BIN_START','BIN_END'],how='inner')
#df_merge = df_merge.query('BIN_START <= POS and BIN_END >= POS')
#
#print "computed windows"
#del df_pass
#
#out=pandas.DataFrame({output+'_N_sites' :df_merge.groupby(['CHROM','BIN_START','BIN_END']).size()}).reset_index()
#print "computed output"

introns=pandas.merge(df_pass, introns, on=['CHROM','POS'],how='inner')
codon_df_1=pandas.merge(df_pass, codon_df_1, on=['CHROM','POS'],how='inner')
codon_df_2=pandas.merge(df_pass, codon_df_2, on=['CHROM','POS'],how='inner')
codon_df_3=pandas.merge(df_pass, codon_df_3, on=['CHROM','POS'],how='inner')
codon_df_4d=pandas.merge(df_pass, codon_df_4d, on=['CHROM','POS'],how='inner')




GENE=pandas.merge(introns,codon_df_1, on=['CHROM','POS'], how='outer',)
GENE=pandas.merge(GENE,codon_df_2, on=['CHROM','POS'], how='outer')
GENE=pandas.merge(GENE,codon_df_3, on=['CHROM','POS'], how='outer')
GENE=pandas.merge(GENE,codon_df_4d, on=['CHROM','POS'], how='outer')


temp=pandas.merge(df_pass, GENE, right_index=True ,on=['CHROM','POS'], how='inner', )


intergenic_sites=df_pass.drop(np.intersect1d(df_pass.index, temp.index))

del GENE
del temp
#intergenic_sites['BIN_START']=(np.floor(intergenic_sites['POS']/100000)*100000)+1
#intergenic_sites['BIN_END']=intergenic_sites['BIN_START']+(100000-1)
#
#df_pass['BIN_START']=(np.floor(df_pass['POS']/100000)*100000)+1
#df_pass['BIN_END']=df_pass['BIN_START']+(100000-1)




######callculate site PI#######
goods=personal_popgen.output_good_sites_fequency(args.freq,args.chr)

#goods=personal_popgen.output_good_sites_fequency("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq",20)

#goods['site_pi']=0.0
#total_comb=(20*(20-1.0))/2.0
#goods['site_pi']=((goods.minor_freq*20)*(goods.major_freq*20))/total_comb


goods['site_pi']=0.0
total_comb=(args.chr*(args.chr-1.0))/2.0
goods['site_pi']=((goods.minor_freq*args.chr)*(goods.major_freq*args.chr))/total_comb


goods['BIN_START']=0
goods['BIN_END']=0
goods['BIN_START']=(np.floor(goods['POS']/100000)*100000)+1
goods['BIN_END']=goods['BIN_START']+(100000-1)


goods_s_s=goods[(goods.minor_allele.isin(['G','C'])) & (goods.major_allele.isin(['G','C']))]
goods_w_w=goods[(goods.minor_allele.isin(['A','T'])) & (goods.major_allele.isin(['A','T']))]
goods_w_w_s_s=goods_s_s.append(goods_w_w)


goods_w_s_1=goods[(goods.minor_allele.isin(['A','T'])) & (goods.major_allele.isin(['G','C']))]
goods_w_s_2=goods[(goods.minor_allele.isin(['G','C'])) & (goods.major_allele.isin(['A','T']))]
goods_w_s=goods_w_s_1.append(goods_w_s_2)


#start_stuff=start_stuff.rename(columns={'scaffold': 'CHROM', 'start':'POS'})
#start_stuff=pandas.merge(goods, start_stuff, on=['CHROM','POS'], how='inner',right_index=True)
#goods.loc[start_stuff.index, 'POS']=goods.loc[start_stuff.index, 'POS']-1


############test################to check the refrence allele

#goods['POS']=goods['POS'].astype(int)
#
#
#fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/GENOME_ASSEMBLY/assembly_updates/v1.4/N.Backstrom_leptidea.scf.1.4.fasta')
#
#goods['site_from_ref']=''
#
#for POS_index, pos in goods.iterrows():
#   site_from_ref=str(fasta[pos.CHROM][pos.POS-1:pos.POS])
#   if site_from_ref in [pos.major_allele,pos.minor_allele]:
#       goods.loc[POS_index,'site_from_ref']=site_from_ref
#   elif str(fasta[pos.CHROM][pos.POS-1:pos.POS]) in [pos.major_allele,pos.minor_allele]:
#       goods.loc[POS_index,'site_from_ref']=str(fasta[pos.CHROM][pos.POS-1:pos.POS])
#

######merge DB


def get_windowed_pi(goods,codon_df_1,filters):
    codon_df_1['BIN_START']=(np.floor(codon_df_1['POS']/100000)*100000)+1
    codon_df_1['BIN_END']=codon_df_1['BIN_START']+(100000-1)
    codon_df_1_counts= codon_df_1.groupby(['CHROM','BIN_START','BIN_END'],as_index=False).agg({'POS': 'count'})
    codon_df_1_counts=codon_df_1_counts[codon_df_1_counts['POS']>=filters]
    codon_df_1_counts=codon_df_1_counts.rename(columns={'POS': 'sites'})
    codon_df_1_temp=pandas.merge(codon_df_1_counts,codon_df_1, on=['CHROM','BIN_START','BIN_END'], how='inner')
    new_df=pandas.merge(goods,codon_df_1_temp, on=['CHROM','BIN_START','BIN_END','POS'], how='right')
    new_df=new_df.groupby(['CHROM','BIN_START','BIN_END'],as_index=False).agg({'sites': np.mean, 'site_pi': np.sum})
    new_df=new_df[new_df.site_pi> 0.0]
    new_df['site_pi']=new_df['site_pi']/new_df['sites']
    return new_df


def get_site_freq_spectrum(goods,codon_df_1):
    codon_df_1_temp=pandas.merge(goods,codon_df_1, on=['CHROM','POS'], how='inner')
    codon_df_1_temp.minor_freq=codon_df_1_temp.minor_freq.round(2)
    spectrum=codon_df_1_temp.groupby('minor_freq')['minor_freq'].count()
    return spectrum






final_global_pi=get_windowed_pi(goods,df_pass,10000)
final_pi_db_codon1=get_windowed_pi(goods,codon_df_1,500)
final_pi_db_codon2=get_windowed_pi(goods,codon_df_2,500)
final_pi_db_codon3=get_windowed_pi(goods,codon_df_3,500)
final_pi_db_codon4d=get_windowed_pi(goods,codon_df_4d,250)
final_pi_db_introns=get_windowed_pi(goods,introns,  5000)
final_pi_db_intergenic=get_windowed_pi(goods,intergenic_sites,  6000)



final_global_pi_corrected=get_windowed_pi(goods,df_pass,10000)
final_pi_db_codon1_corrected=get_windowed_pi(goods_w_w_s_s,codon_df_1,500)
final_pi_db_codon2_corrected=get_windowed_pi(goods_w_w_s_s,codon_df_2,500)
final_pi_db_codon3_corrected=get_windowed_pi(goods_w_w_s_s,codon_df_3,500)
final_pi_db_codon4d_corrected=get_windowed_pi(goods_w_w_s_s,codon_df_4d,250)
final_pi_db_introns_corrected=get_windowed_pi(goods_w_w_s_s,introns,  5000)
final_pi_db_intergenic_corrected=get_windowed_pi(goods_w_w_s_s,intergenic_sites,  6000)


final_global_pi_w_s=get_windowed_pi(goods,df_pass,10000)
final_pi_db_codon1_w_s=get_windowed_pi(goods_w_s,codon_df_1,500)
final_pi_db_codon2_w_s=get_windowed_pi(goods_w_s,codon_df_2,500)
final_pi_db_codon3_w_s=get_windowed_pi(goods_w_s,codon_df_3,500)
final_pi_db_codon4d_w_s=get_windowed_pi(goods_w_s,codon_df_4d,250)
final_pi_db_introns_w_s=get_windowed_pi(goods_w_s,introns,  5000)
final_pi_db_intergenic_w_s=get_windowed_pi(goods_w_s,intergenic_sites,  6000)

#



out_put='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/'

pop_name=str(args.population).split('/')[-1]


#final_pi_db_codon1.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/' +args.part+'/'+pop_name+'_'+str(args.coverage)+'_final_pi_db_codon1.csv')
#final_pi_db_codon2.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/' +args.part+'/'+pop_name+'_'+str(args.coverage)+'_final_pi_db_codon2.csv')
#final_pi_db_codon3.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/' +args.part+'/'+pop_name+'_'+str(args.coverage)+'_final_pi_db_codon3.csv')
#final_pi_db_introns.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/'+args.part+'/'+pop_name+'_'+str(args.coverage)+'_final_pi_db_introns.csv')

print final_pi_db_codon1.site_pi.mean()
print final_pi_db_codon2.site_pi.mean()
print final_pi_db_codon3.site_pi.mean()
print final_pi_db_codon4d.site_pi.mean()
print final_pi_db_introns.site_pi.mean()
#print final_global_pi.site_pi.mean()
print final_pi_db_intergenic.site_pi.mean()






print final_pi_db_codon1_w_s.site_pi.mean()
print final_pi_db_codon2_w_s.site_pi.mean()
print final_pi_db_codon3_w_s.site_pi.mean()
print final_pi_db_codon4d_w_s.site_pi.mean()
print final_pi_db_introns_w_s.site_pi.mean()
#print final_global_pi.site_pi.mean()
print final_pi_db_intergenic_w_s.site_pi.mean()













fianl_df=pandas.merge(final_pi_db_codon1,final_pi_db_codon2, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi_x': 'site_pi_codon_1', 'sites_x':'sites_codon_1','site_pi_y': 'site_pi_codon_2', 'sites_y':'sites_codon_2','site_pi_x': 'site_pi_codon_1', 'sites_x':'sites_codon_1'})
fianl_df=pandas.merge(fianl_df,final_pi_db_codon3, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon_3', 'sites':'sites_codon_3'})
fianl_df=pandas.merge(fianl_df,final_pi_db_introns, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_introns', 'sites':'sites_introns'})
fianl_df=pandas.merge(fianl_df,final_global_pi, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_global', 'sites':'sites_global'})
fianl_df=pandas.merge(fianl_df,final_pi_db_codon4d, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon_4d', 'sites':'sites_codon_4d'})
fianl_df=pandas.merge(fianl_df,final_pi_db_intergenic, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_intergenic', 'sites':'sites_intergenic'})
fianl_df=pandas.merge(final_pi_db_codon1,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi_x': 'site_pi_codon_1', 'sites_x':'sites_codon_1','site_pi_y': 'site_pi_codon_2', 'sites_y':'sites_codon_2','site_pi_x': 'site_pi_codon_1', 'sites_x':'sites_codon_1'})

fianl_df=pandas.merge(final_pi_db_codon1_corrected,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon1_corrected', 'sites':'sites_codon1_corrected'})
fianl_df=pandas.merge(final_pi_db_codon2_corrected,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon2_corrected', 'sites':'sites_codon2_corrected'})
fianl_df=pandas.merge(final_pi_db_codon3_corrected,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon3_corrected', 'sites':'sites_codon3_corrected'})
fianl_df=pandas.merge(final_pi_db_codon4d_corrected,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_codon4d_corrected', 'sites':'sites_codon4d_corrected'})
fianl_df=pandas.merge(final_pi_db_introns_corrected,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_introns_corrected', 'sites':'sites_introns_corrected'})
fianl_df=pandas.merge(final_pi_db_intergenic_corrected,fianl_df, on=['CHROM','BIN_START','BIN_END'], how='outer',suffixes=('_x', '_y'))
fianl_df=fianl_df.rename(columns={'site_pi': 'site_pi_intergenic_corrected', 'sites':'sites_intergenic_corrected'})



#fianl_df=fianl_df[fianl_df.sites_global>20000]
#sns.boxplot(fianl_df[["site_pi_codon_1","site_pi_codon_2","site_pi_codon_3","site_pi_introns","site_pi_global"]],showfliers=False).set(ylim=(0,0.0105))


fianl_df.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/'+args.part+'/'+pop_name+'_'+str(args.coverage)+'_final_df_with_inter_geneic.csv')




###############################################################
