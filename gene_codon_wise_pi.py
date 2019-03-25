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

pwd='/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/'
os.chdir(pwd)



def get_args():
    parser = argparse.ArgumentParser(description='''outputs sites covered 7X in 7Induveduals in population pairs

        -f freq file
        -pop pop file
        -c coverage file
        ''')
    parser.add_argument('-f', '--freq', type=str, required=True)
    parser.add_argument('-pop', '--population', type=str, required=True)
    parser.add_argument('-c', '--cov', type=str, required=True)
    return parser.parse_args()


args = get_args()


print args.freq
print args.population

introns=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/introns.csv',names=['CHROM','POS'], sep=' ',header=None)
codon_df_1=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_1.csv',names=['CHROM','POS'], sep=' ',header=None)
codon_df_2=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_2.csv',names=['CHROM','POS'], sep=' ',header=None)
codon_df_3=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_3.csv',names=['CHROM','POS'], sep=' ',header=None)
start_stuff=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/start_stuff.csv',names=['CHROM','POS'], sep=' ',header=None)
#coverage=pandas.read_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kaz_sin_1_X.csv',names=['.','CHROM','POS'],index_col=False)
coverage=pandas.read_csv(args.cov,names=['.','CHROM','POS'],index_col=False)
coverage=coverage[['CHROM','POS']]




introns.POS=introns.POS+1
codon_df_1.POS=codon_df_1.POS+1
codon_df_2.POS=codon_df_2.POS+1
codon_df_3.POS=codon_df_3.POS+1



introns=pandas.merge(introns, coverage, on=['CHROM','POS'])
codon_df_1=pandas.merge(codon_df_1, coverage, on=['CHROM','POS'])
codon_df_2=pandas.merge(codon_df_2, coverage, on=['CHROM','POS'])
codon_df_3=pandas.merge(codon_df_3, coverage, on=['CHROM','POS'])
coverage=0



head_stuff=['CHROM','source','feature','BIN_START','BIN_END','.','-','..','infor']
annotation=pandas.read_table("/proj/b2014034/NBIS_annotation_leptidea/gff/gene-builds/leptidea_sinapis_rc1.gff",skiprows=1, names=head_stuff)
annotation=annotation[['CHROM','source','feature','BIN_START','BIN_END','-','infor']]
#annotation_gene=annotation[(annotation.feature=='gene') & (annotation.infor.str[31:37]==';Name=')]
annotation_gene=annotation[(annotation.feature=='gene')]
annotation_gene=annotation_gene[['CHROM','BIN_START','BIN_END','-','infor']]

#goods=personal_popgen.output_good_sites_fequency('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq',20)
goods=personal_popgen.output_good_sites_fequency(args.freq,20)
goods.POS=goods.POS-1
goods.POS=goods.POS.astype('int64')

goods['site_pi']=0.0
total_comb=(20.0*(20.0-1.0))/2.0
goods['site_pi']=((goods.minor_freq*20)*(goods.major_freq*20))/total_comb




def gene_wise_codon_pi(goods, introns, codon_df_1, codon_df_2, codon_df_3, annotation_gene,coverage):
    gene_pi=pandas.DataFrame(columns=['GENE','introns_pi','codon_1_pi','codon_2_pi','codon_3_pi','sites_int','sites_cd1','sites_cd2','sites_cd3','+/-'])
    for scaffold in set(list(annotation_gene['CHROM'])):
        print scaffold
        scaffold_gene=annotation_gene[(annotation_gene.CHROM == scaffold)]
        for gene_index, gene in scaffold_gene.iterrows():
            if gene['-']=='+':
                temp_introns=introns[(gene['CHROM']==introns.CHROM) & (gene['BIN_START']<=introns.POS) & (gene['BIN_END']>=introns.POS)]
                temp_codon_1=codon_df_1[(gene['CHROM']==codon_df_1.CHROM) & (gene['BIN_START']<=codon_df_1.POS) & (gene['BIN_END']>=codon_df_1.POS)]
                temp_codon_2=codon_df_2[(gene['CHROM']==codon_df_2.CHROM) & (gene['BIN_START']<=codon_df_2.POS) & (gene['BIN_END']>=codon_df_2.POS)]
                temp_codon_3=codon_df_3[(gene['CHROM']==codon_df_3.CHROM) & (gene['BIN_START']<=codon_df_3.POS) & (gene['BIN_END']>=codon_df_3.POS)]
                if len(temp_introns) > 0:
                    sites_int=len(temp_introns)
                    sites_cd1=len(temp_codon_1)
                    sites_cd2=len(temp_codon_2)
                    sites_cd3=len(temp_codon_3)
                    temp_introns.POS=temp_introns.POS-1
                    temp_codon_1.POS=temp_codon_1.POS-1
                    temp_codon_2.POS=temp_codon_2.POS-1
                    temp_codon_3.POS=temp_codon_3.POS-1
                    temp_introns=pandas.merge(temp_introns, goods, on=['CHROM','POS'])
                    temp_codon_1=pandas.merge(temp_codon_1, goods, on=['CHROM','POS'])
                    temp_codon_2=pandas.merge(temp_codon_2, goods, on=['CHROM','POS'])
                    temp_codon_3=pandas.merge(temp_codon_3, goods, on=['CHROM','POS'])
                    pi_introns=temp_introns.site_pi.sum()/float(sites_int)
                    pi_codon_1=temp_codon_1.site_pi.sum()/float(sites_cd1)
                    pi_codon_2=temp_codon_2.site_pi.sum()/float(sites_cd2)
                    pi_codon_3=temp_codon_3.site_pi.sum()/float(sites_cd3)
                    print pi_introns
                    print pi_codon_1
                    print pi_codon_2
                    print pi_codon_3
                    temp={}
                    temp['GENE']=pandas.Series(gene.infor)
                    temp['introns_pi']=pandas.Series(pi_introns)
                    temp['codon_1_pi']=pandas.Series(pi_codon_1)
                    temp['codon_2_pi']=pandas.Series(pi_codon_2)
                    temp['codon_3_pi']=pandas.Series(pi_codon_3)
                    temp['sites_int']=pandas.Series(sites_int)
                    temp['sites_cd1']=pandas.Series(sites_cd1)
                    temp['sites_cd2']=pandas.Series(sites_cd2)
                    temp['sites_cd3']=pandas.Series(sites_cd3)
                    temp['+/-']=pandas.Series('+')
                    temp_df = pandas.DataFrame(temp)
                    gene_pi=gene_pi.append(temp_df[['GENE','introns_pi','codon_1_pi','codon_2_pi','codon_3_pi','sites_int', 'sites_cd1', 'sites_cd2', 'sites_cd3','+/-']],ignore_index=True)
            elif gene['-']=='-':
                temp_introns=introns[(gene['CHROM']==introns.CHROM) & (gene['BIN_START']>=introns.POS) & (gene['BIN_END']<=introns.POS)]
                temp_codon_1=codon_df_1[(gene['CHROM']==codon_df_1.CHROM) & (gene['BIN_START']>=codon_df_1.POS) & (gene['BIN_END']<=codon_df_1.POS)]
                temp_codon_2=codon_df_2[(gene['CHROM']==codon_df_2.CHROM) & (gene['BIN_START']>=codon_df_2.POS) & (gene['BIN_END']<=codon_df_2.POS)]
                temp_codon_3=codon_df_3[(gene['CHROM']==codon_df_3.CHROM) & (gene['BIN_START']>=codon_df_3.POS) & (gene['BIN_END']<=codon_df_3.POS)]              
                if len(temp_introns) > 0:
                    sites_int=len(temp_introns)
                    sites_cd1=len(temp_codon_1)
                    sites_cd2=len(temp_codon_2)
                    sites_cd3=len(temp_codon_3)
                    temp_introns.POS=temp_introns.POS-1
                    temp_codon_1.POS=temp_codon_1.POS-1
                    temp_codon_2.POS=temp_codon_2.POS-1
                    temp_codon_3.POS=temp_codon_3.POS-1
                    temp_introns=pandas.merge(temp_introns, goods, on=['CHROM','POS'])
                    temp_codon_1=pandas.merge(temp_codon_1, goods, on=['CHROM','POS'])
                    temp_codon_2=pandas.merge(temp_codon_2, goods, on=['CHROM','POS'])
                    temp_codon_3=pandas.merge(temp_codon_3, goods, on=['CHROM','POS'])
                    pi_introns=temp_introns.site_pi.sum()/float(sites_int)
                    pi_codon_1=temp_codon_1.site_pi.sum()/float(sites_cd1)
                    pi_codon_2=temp_codon_2.site_pi.sum()/float(sites_cd2)
                    pi_codon_3=temp_codon_3.site_pi.sum()/float(sites_cd3)
                    print pi_introns
                    print pi_codon_1
                    print pi_codon_2
                    print pi_codon_3
                    temp={}
                    temp['GENE']=pandas.Series(gene.infor)
                    temp['introns_pi']=pandas.Series(pi_introns)
                    temp['codon_1_pi']=pandas.Series(pi_codon_1)
                    temp['codon_2_pi']=pandas.Series(pi_codon_2)
                    temp['codon_3_pi']=pandas.Series(pi_codon_3)
                    temp['sites_int']=pandas.Series(sites_int)
                    temp['sites_cd1']=pandas.Series(sites_cd1)
                    temp['sites_cd2']=pandas.Series(sites_cd2)
                    temp['sites_cd3']=pandas.Series(sites_cd3)
                    temp['+/-']=pandas.Series('-')
                    temp_df = pandas.DataFrame(temp)
                    gene_pi=gene_pi.append(temp_df[['GENE','introns_pi','codon_1_pi','codon_2_pi','codon_3_pi','sites_int', 'sites_cd1', 'sites_cd2', 'sites_cd3','+/-']],ignore_index=True)
    return gene_pi



output=gene_wise_codon_pi(goods, introns, codon_df_1, codon_df_2, codon_df_3, annotation_gene, coverage)
output.to_csv('/proj/b2014034/nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/'+pop_name+'_gene_codon_stats_df.csv')


introns 0.00131874852048    0.00273215131676    0.00349089030617    0.00305261140973    0.00329677643596    0.00385841689386    
codon_1 0.00497555363819    0.00101927353595    0.000858369098712   0.0060595567867     0.000650920587688   0.0                 
codon_2 0.00468794938165    0.00313194959229    0.000813191777728   0.00429362880886    0.00515157150828    0.000975609756098   
codon_3 0.00135174000575    0.00257598220904    0.00291393720352    0.00512465373961    0.000176991150442   0.00310654685494    




