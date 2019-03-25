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

def join_raw_data_base(lists, on_what, how_is_it):
        for i in range(1,len(lists)):
                j=i-1
                if i == 1:
                        new_2=pandas.merge(lists[j], lists[i], on=on_what, how=how_is_it)
                else:
                        new_2=pandas.merge(nextone, lists[i], on=on_what, how=how_is_it)
                nextone=new_2
        return nextone


def load_files(pop):
    Pop_1=pop
    fianl_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf1/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf2/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf3/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_4=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf4/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_5=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf5/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_6=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf6/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_7=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf7/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_8=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf8/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_9=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf9/'+Pop_1+'_2_final_df.csv',index_col=False)
    fianl_df_10=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf10/'+Pop_1+'_2_final_df.csv',index_col=False)
    final_df_final=pandas.concat([fianl_df_1,fianl_df_2,fianl_df_3,fianl_df_4,fianl_df_5,fianl_df_6,fianl_df_7,fianl_df_8,fianl_df_9,fianl_df_10],ignore_index=True)
    final_df_final=final_df_final.rename(columns={'site_pi_codon_1': Pop_1+'_site_pi_codon_1', 'site_pi_codon_2': Pop_1+'_site_pi_codon_2', 'site_pi_codon_3': Pop_1+'_site_pi_codon_3', 'site_pi_codon_4d': Pop_1+'_site_pi_codon_4d', 'site_pi_introns': Pop_1+'_site_pi_introns', 'site_pi_global': Pop_1+'_site_pi_global','sites_codon_1': Pop_1+'_sites_codon_1','sites_codon_2': Pop_1+'_sites_codon_2','sites_codon_3': Pop_1+'_sites_codon_3','sites_codon_4d': Pop_1+'_sites_codon_4d','sites_introns': Pop_1+'_sites_introns','sites_global': Pop_1+'_sites_global'})
    #final_df_final=final_df_final.rename(columns={'site_pi_codon_1_corrected': Pop_1+'_site_pi_codon_1_corrected', 'site_pi_codon_2_corrected': Pop_1+'_site_pi_codon_2_corrected', 'site_pi_codon_3_corrected': Pop_1+'_site_pi_codon_3_corrected', 'site_pi_codon_4d_corrected': Pop_1+'_site_pi_codon_4d_corrected', 'site_pi_introns_corrected': Pop_1+'_site_pi_introns_corrected', 'sites_codon_1_corrected': Pop_1+'_sites_codon_1_corrected','sites_codon_2_corrected': Pop_1+'_sites_codon_2_corrected','sites_codon_3_corrected': Pop_1+'_sites_codon_3_corrected','sites_codon_4d_corrected': Pop_1+'_sites_codon_4d_corrected','sites_introns_corrected': Pop_1+'_sites_introns_corrected'})
    return final_df_final



pop_df_1=load_files('irish_juvernica')
pop_df_2=load_files('juvernica')
pop_df_3=load_files('kazak_juvernica')
pop_df_4=load_files('kaz_sin')
pop_df_5=load_files('sinapis')
pop_df_6=load_files('spanish_reali')
pop_df_7=load_files('spanish_sinapis')
pop_df_8=load_files('swe_sin_allele')

#pop_df_1=pop_df_1.dropna()
#pop_df_2=pop_df_2.dropna()
#pop_df_3=pop_df_3.dropna()
#pop_df_4=pop_df_4.dropna()
#pop_df_5=pop_df_5.dropna()
#pop_df_6=pop_df_6.dropna()
#pop_df_7=pop_df_7.dropna()
#pop_df_8=pop_df_8.dropna()


#############filter_for_ns##########
lengths=[]
linst_wrong=[]
CDS_fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/fasta/cds.fa')

head_stuff=['CHROM','source','feature','start','end','extra','-','indexing','infor']
annotation=pandas.read_table("/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/gff/gene-builds/leptidea_sinapis_rc1.gff",skiprows=1, names=head_stuff)
annotation=annotation[['CHROM','source','feature','start','end','-','indexing','infor']]
annotation=annotation[(annotation.feature=='mRNA')]
annotation['geneID']=annotation.infor.str[3:31].iloc[0:]

for i in CDS_fasta.keys():
	if str(CDS_fasta[i]).count('N')/float(len(CDS_fasta[i])) > 0.05:
		linst_wrong.append(i)

chromosome_file=[stuff.split() for stuff in open('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/chromosomes.txt')]
chromosomes_dict={}
for j in chromosome_file:
	if j[1] not in chromosomes_dict.keys():
		chromosomes_dict[j[1]]=[]
		chromosomes_dict[j[1]].append(j[0])
	else:
		chromosomes_dict[j[1]].append(j[0])

annotation['CHROM_ID']='.'
annotation['auto_allo']='auto'
for n in chromosomes_dict.keys():
	for nj in chromosomes_dict[n]:
		annotation.ix[(annotation['CHROM'] == nj), 'CHROM_ID']=str(n)

annotation.ix[(annotation['CHROM_ID'] == '21'), 'auto_allo']='Z'

z=annotation[(annotation['auto_allo'] == 'Z')]
a=annotation[(annotation['auto_allo'] != 'Z')]

annotation=annotation[['CHROM','CHROM_ID','auto_allo', 'geneID', 'start', 'end']].drop_duplicates(keep='first')
annotation['MID']=(annotation.start+annotation.end)/2.0
annotation['BIN_START']=0
annotation['BIN_END']=0
annotation['BIN_START']=(np.floor(annotation['MID']/100000)*100000)+1
annotation['BIN_END']=annotation['BIN_START']+(100000-1)


###############test_gene density ###############

gene_density=annotation.groupby(['CHROM','BIN_START','BIN_END'],as_index=False).agg({'geneID': 'count'})
gene_density=gene_density.rename(columns={'geneID': 'gene_counts'})


###################################

pop_df_1=join_raw_data_base([pop_df_1,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_2=join_raw_data_base([pop_df_2,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_3=join_raw_data_base([pop_df_3,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_4=join_raw_data_base([pop_df_4,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_5=join_raw_data_base([pop_df_5,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_6=join_raw_data_base([pop_df_6,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_7=join_raw_data_base([pop_df_7,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')
pop_df_8=join_raw_data_base([pop_df_8,annotation[['CHROM','CHROM_ID','auto_allo']].drop_duplicates(keep='first')],['CHROM'],'inner')


pop_df_1['Relative_irish_juvernica_site_pi']=(pop_df_1.irish_juvernica_site_pi_global - pop_df_1.irish_juvernica_site_pi_global.mean())/pop_df_1.irish_juvernica_site_pi_global.std()
pop_df_2['Relative_juvernica_site_pi']=(pop_df_2.juvernica_site_pi_global - pop_df_2.juvernica_site_pi_global.mean())/pop_df_2.juvernica_site_pi_global.std()
pop_df_3['Relative_kazak_juvernica_site_pi']=(pop_df_3.kazak_juvernica_site_pi_global - pop_df_3.kazak_juvernica_site_pi_global.mean())/pop_df_3.kazak_juvernica_site_pi_global.std()
pop_df_4['Relative_kaz_sin_site_pi']=(pop_df_4.kaz_sin_site_pi_global - pop_df_4.kaz_sin_site_pi_global.mean())/pop_df_4.kaz_sin_site_pi_global.std()
pop_df_5['Relative_sinapis_site_pi']=(pop_df_5.sinapis_site_pi_global - pop_df_5.sinapis_site_pi_global.mean())/pop_df_5.sinapis_site_pi_global.std()
pop_df_6['Relative_spanish_reali_site_pi']=(pop_df_6.spanish_reali_site_pi_global - pop_df_6.spanish_reali_site_pi_global.mean())/pop_df_6.spanish_reali_site_pi_global.std()
pop_df_7['Relative_spanish_sinapis_site_pi']=(pop_df_7.spanish_sinapis_site_pi_global - pop_df_7.spanish_sinapis_site_pi_global.mean())/pop_df_7.spanish_sinapis_site_pi_global.std()
pop_df_8['Relative_swe_sin_allele_site_pi']=(pop_df_8.swe_sin_allele_site_pi_global - pop_df_8.swe_sin_allele_site_pi_global.mean())/pop_df_8.swe_sin_allele_site_pi_global.std()



pop_df_1=join_raw_data_base([pop_df_1,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_2=join_raw_data_base([pop_df_2,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_3=join_raw_data_base([pop_df_3,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_4=join_raw_data_base([pop_df_4,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_5=join_raw_data_base([pop_df_5,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_6=join_raw_data_base([pop_df_6,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_7=join_raw_data_base([pop_df_7,gene_density],['BIN_START','BIN_END','CHROM'],'inner')
pop_df_8=join_raw_data_base([pop_df_8,gene_density],['BIN_START','BIN_END','CHROM'],'inner')

sns.jointplot("irish_juvernica_site_pi_global", "gene_counts", data=pop_df_1, kind="reg")
sns.jointplot("juvernica_site_pi_global", "gene_counts", data=pop_df_2, kind="reg")
sns.jointplot("kazak_juvernica_site_pi_global", "gene_counts", data=pop_df_3, kind="reg")
sns.jointplot("kaz_sin_site_pi_global", "gene_counts", data=pop_df_4, kind="reg")
sns.jointplot("sinapis_site_pi_global", "gene_counts", data=pop_df_5, kind="reg")
sns.jointplot("spanish_reali_site_pi_global", "gene_counts", data=pop_df_6, kind="reg")
sns.jointplot("spanish_sinapis_site_pi_global", "gene_counts", data=pop_df_7, kind="reg")
sns.jointplot("swe_sin_allele_site_pi_global", "gene_counts", data=pop_df_8, kind="reg")











sns.boxplot(pop_df_1[['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_introns','irish_juvernica_site_pi_global','irish_juvernica_site_pi_codon_4d']],showfliers=False)



sns.boxplot(pop_df_8[['swe_sin_allele_site_pi_codon_1','swe_sin_allele_site_pi_codon_2','swe_sin_allele_site_pi_codon_3','swe_sin_allele_site_pi_introns','swe_sin_allele_site_pi_global','swe_sin_allele_site_pi_codon_4d']],showfliers=False)



fianl_df_1=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/irish_juvernica.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_2=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/juvernica.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_3=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_juvernica.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_4=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_sinapis.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_5=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/sinapis.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_6=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_reali.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_7=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_sinapis.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_8=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Swedish_sinapis.table.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]


fianl_df_1=join_raw_data_base([fianl_df_1,annotation],['geneID'],'inner')
fianl_df_2=join_raw_data_base([fianl_df_2,annotation],['geneID'],'inner')
fianl_df_3=join_raw_data_base([fianl_df_3,annotation],['geneID'],'inner')
fianl_df_4=join_raw_data_base([fianl_df_4,annotation],['geneID'],'inner')
fianl_df_5=join_raw_data_base([fianl_df_5,annotation],['geneID'],'inner')
fianl_df_6=join_raw_data_base([fianl_df_6,annotation],['geneID'],'inner')
fianl_df_7=join_raw_data_base([fianl_df_7,annotation],['geneID'],'inner')
fianl_df_8=join_raw_data_base([fianl_df_8,annotation],['geneID'],'inner')


fianl_df_1=fianl_df_1.rename(columns={'Pin': 'Pin_irish_juvernica', 'Pis': 'Pis_irish_juvernica', "Pi_nonsyn": "Pi_nonsyn_irish_juvernica", "Pi_syn": "Pi_syn_irish_juvernica","nonsyn_sites":"nonsyn_sites_irish_juvernica", "syn_sites": "syn_sites_irish_juvernica"})
fianl_df_2=fianl_df_2.rename(columns={'Pin': 'Pin_juvernica', 'Pis': 'Pis_juvernica', "Pi_nonsyn": "Pi_nonsyn_juvernica", "Pi_syn": "Pi_syn_juvernica","nonsyn_sites":"nonsyn_sites_juvernica", "syn_sites": "syn_sites_juvernica"})
fianl_df_3=fianl_df_3.rename(columns={'Pin': 'Pin_kazak_juvernica', 'Pis': 'Pis_kazak_juvernica', "Pi_nonsyn": "Pi_nonsyn_kazak_juvernica", "Pi_syn": "Pi_syn_kazak_juvernica","nonsyn_sites":"nonsyn_sites_kazak_juvernica", "syn_sites": "syn_sites_kazak_juvernica"})
fianl_df_4=fianl_df_4.rename(columns={'Pin': 'Pin_kazak_sinapis', 'Pis': 'Pis_kazak_sinapis', "Pi_nonsyn": "Pi_nonsyn_kazak_sinapis", "Pi_syn": "Pi_syn_kazak_sinapis","nonsyn_sites":"nonsyn_sites_kazak_sinapis", "syn_sites": "syn_sites_kazak_sinapis"})
fianl_df_5=fianl_df_5.rename(columns={'Pin': 'Pin_sinapis', 'Pis': 'Pis_sinapis', "Pi_nonsyn": "Pi_nonsyn_sinapis", "Pi_syn": "Pi_syn_sinapis","nonsyn_sites":"nonsyn_sites_sinapis", "syn_sites": "syn_sites_sinapis"})
fianl_df_6=fianl_df_6.rename(columns={'Pin': 'Pin_spanish_reali', 'Pis': 'Pis_spanish_reali', "Pi_nonsyn": "Pi_nonsyn_spanish_reali", "Pi_syn": "Pi_syn_spanish_reali","nonsyn_sites":"nonsyn_sites_spanish_reali", "syn_sites": "syn_sites_spanish_reali"})
fianl_df_7=fianl_df_7.rename(columns={'Pin': 'Pin_spanish_sinapis', 'Pis': 'Pis_spanish_sinapis', "Pi_nonsyn": "Pi_nonsyn_spanish_sinapis", "Pi_syn": "Pi_syn_spanish_sinapis","nonsyn_sites":"nonsyn_sites_spanish_sinapis", "syn_sites": "syn_sites_spanish_sinapis"})
fianl_df_8=fianl_df_8.rename(columns={'Pin': 'Pin_Swedish_sinapis', 'Pis': 'Pis_Swedish_sinapis', "Pi_nonsyn": "Pi_nonsyn_Swedish_sinapis", "Pi_syn": "Pi_syn_Swedish_sinapis","nonsyn_sites":"nonsyn_sites_Swedish_sinapis", "syn_sites": "syn_sites_Swedish_sinapis"})






fianl_df_1=fianl_df_1[~fianl_df_1['geneID'].isin(linst_wrong)]
fianl_df_2=fianl_df_2[~fianl_df_2['geneID'].isin(linst_wrong)]
fianl_df_3=fianl_df_3[~fianl_df_3['geneID'].isin(linst_wrong)]
fianl_df_4=fianl_df_4[~fianl_df_4['geneID'].isin(linst_wrong)]
fianl_df_5=fianl_df_5[~fianl_df_5['geneID'].isin(linst_wrong)]
fianl_df_6=fianl_df_6[~fianl_df_6['geneID'].isin(linst_wrong)]
fianl_df_7=fianl_df_7[~fianl_df_7['geneID'].isin(linst_wrong)]
fianl_df_8=fianl_df_8[~fianl_df_8['geneID'].isin(linst_wrong)]



fianl_df_1['Pin_Pis_irish_juvernica']=fianl_df_1["Pin_irish_juvernica"]/fianl_df_1["Pis_irish_juvernica"]
fianl_df_2['Pin_Pis_juvernica']=fianl_df_2["Pin_juvernica"]/fianl_df_2["Pis_juvernica"]
fianl_df_3['Pin_Pis_kazak_juvernica']=fianl_df_3["Pin_kazak_juvernica"]/fianl_df_3["Pis_kazak_juvernica"]
fianl_df_4['Pin_Pis_kazak_sinapis']=fianl_df_4["Pin_kazak_sinapis"]/fianl_df_4["Pis_kazak_sinapis"]
fianl_df_5['Pin_Pis_sinapis']=fianl_df_5["Pin_sinapis"]/fianl_df_5["Pis_sinapis"]
fianl_df_6['Pin_Pis_spanish_reali']=fianl_df_6["Pin_spanish_reali"]/fianl_df_6["Pis_spanish_reali"]
fianl_df_7['Pin_Pis_spanish_sinapis']=fianl_df_7["Pin_spanish_sinapis"]/fianl_df_7["Pis_spanish_sinapis"]
fianl_df_8['Pin_Pis_Swedish_sinapis']=fianl_df_8["Pin_Swedish_sinapis"]/fianl_df_8["Pis_Swedish_sinapis"]

fianl_df_1['Pin_Pis_irish_juvernica']=fianl_df_1['Pin_Pis_irish_juvernica'].replace(np.inf, np.nan)
fianl_df_2['Pin_Pis_juvernica']=fianl_df_2['Pin_Pis_juvernica'].replace(np.inf, np.nan)
fianl_df_3['Pin_Pis_kazak_juvernica']=fianl_df_3['Pin_Pis_kazak_juvernica'].replace(np.inf, np.nan)
fianl_df_4['Pin_Pis_kazak_sinapis']=fianl_df_4['Pin_Pis_kazak_sinapis'].replace(np.inf, np.nan)
fianl_df_5['Pin_Pis_sinapis']=fianl_df_5['Pin_Pis_sinapis'].replace(np.inf, np.nan)
fianl_df_6['Pin_Pis_spanish_reali']=fianl_df_6['Pin_Pis_spanish_reali'].replace(np.inf, np.nan)
fianl_df_7['Pin_Pis_spanish_sinapis']=fianl_df_7['Pin_Pis_spanish_sinapis'].replace(np.inf, np.nan)
fianl_df_8['Pin_Pis_Swedish_sinapis']=fianl_df_8['Pin_Pis_Swedish_sinapis'].replace(np.inf, np.nan)

print fianl_df_1["Pin_Pis_irish_juvernica"].mean()
print fianl_df_2["Pin_Pis_juvernica"].mean()
print fianl_df_3["Pin_Pis_kazak_juvernica"].mean()
print fianl_df_4["Pin_Pis_kazak_sinapis"].mean()
print fianl_df_5["Pin_Pis_sinapis"].mean()
print fianl_df_6["Pin_Pis_spanish_reali"].mean()
print fianl_df_7["Pin_Pis_spanish_sinapis"].mean()
print fianl_df_8["Pin_Pis_Swedish_sinapis"].mean()


fianl_df_1=join_raw_data_base([fianl_df_1,pop_df_1[['BIN_START','BIN_END','CHROM','irish_juvernica_site_pi_codon_4d','irish_juvernica_site_pi_global', 'irish_juvernica_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_2=join_raw_data_base([fianl_df_2,pop_df_2[['BIN_START','BIN_END','CHROM','juvernica_site_pi_codon_4d','juvernica_site_pi_global', 'juvernica_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_3=join_raw_data_base([fianl_df_3,pop_df_3[['BIN_START','BIN_END','CHROM','kazak_juvernica_site_pi_codon_4d','kazak_juvernica_site_pi_global', 'kazak_juvernica_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_4=join_raw_data_base([fianl_df_4,pop_df_4[['BIN_START','BIN_END','CHROM','kaz_sin_site_pi_codon_4d','kaz_sin_site_pi_global', 'kaz_sin_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_5=join_raw_data_base([fianl_df_5,pop_df_5[['BIN_START','BIN_END','CHROM','sinapis_site_pi_codon_4d','sinapis_site_pi_global', 'sinapis_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_6=join_raw_data_base([fianl_df_6,pop_df_6[['BIN_START','BIN_END','CHROM','spanish_reali_site_pi_codon_4d','spanish_reali_site_pi_global', 'spanish_reali_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_7=join_raw_data_base([fianl_df_7,pop_df_7[['BIN_START','BIN_END','CHROM','spanish_sinapis_site_pi_codon_4d','spanish_sinapis_site_pi_global', 'spanish_sinapis_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')
fianl_df_8=join_raw_data_base([fianl_df_8,pop_df_8[['BIN_START','BIN_END','CHROM','swe_sin_allele_site_pi_codon_4d','swe_sin_allele_site_pi_global', 'swe_sin_allele_sites_global']]],['BIN_START','BIN_END','CHROM'],'inner')



########################################relative pi plot agains pin/pis

fianl_df_1['Relative_irish_juvernica_site_pi']=(fianl_df_1.irish_juvernica_site_pi_global - fianl_df_1.irish_juvernica_site_pi_global.mean())/fianl_df_1.irish_juvernica_site_pi_global.std()
fianl_df_2['Relative_juvernica_site_pi']=(fianl_df_2.juvernica_site_pi_global - fianl_df_2.juvernica_site_pi_global.mean())/fianl_df_2.juvernica_site_pi_global.std()
fianl_df_3['Relative_kazak_juvernica_site_pi']=(fianl_df_3.kazak_juvernica_site_pi_global - fianl_df_3.kazak_juvernica_site_pi_global.mean())/fianl_df_3.kazak_juvernica_site_pi_global.std()
fianl_df_4['Relative_kaz_sin_site_pi']=(fianl_df_4.kaz_sin_site_pi_global - fianl_df_4.kaz_sin_site_pi_global.mean())/fianl_df_4.kaz_sin_site_pi_global.std()
fianl_df_5['Relative_sinapis_site_pi']=(fianl_df_5.sinapis_site_pi_global - fianl_df_5.sinapis_site_pi_global.mean())/fianl_df_5.sinapis_site_pi_global.std()
fianl_df_6['Relative_spanish_reali_site_pi']=(fianl_df_6.spanish_reali_site_pi_global - fianl_df_6.spanish_reali_site_pi_global.mean())/fianl_df_6.spanish_reali_site_pi_global.std()
fianl_df_7['Relative_spanish_sinapis_site_pi']=(fianl_df_7.spanish_sinapis_site_pi_global - fianl_df_7.spanish_sinapis_site_pi_global.mean())/fianl_df_7.spanish_sinapis_site_pi_global.std()
fianl_df_8['Relative_swe_sin_allele_site_pi']=(fianl_df_8.swe_sin_allele_site_pi_global - fianl_df_8.swe_sin_allele_site_pi_global.mean())/fianl_df_8.swe_sin_allele_site_pi_global.std()



fianl_df_1_new=fianl_df_1[['CHROM','BIN_START','BIN_END','Pin_Pis_irish_juvernica',"irish_juvernica_site_pi_global",'irish_juvernica_sites_global', 'Relative_irish_juvernica_site_pi']].dropna(subset=['Pin_Pis_irish_juvernica']).groupby(['CHROM','BIN_START','BIN_END',"irish_juvernica_site_pi_global",'Relative_irish_juvernica_site_pi','irish_juvernica_sites_global'],as_index=False).agg({'Pin_Pis_irish_juvernica': 'mean'})
fianl_df_2_new=fianl_df_2[['CHROM','BIN_START','BIN_END','Pin_Pis_juvernica',"juvernica_site_pi_global",'juvernica_sites_global', 'Relative_juvernica_site_pi']].dropna(subset=['Pin_Pis_juvernica']).groupby(['CHROM','BIN_START','BIN_END',"juvernica_site_pi_global",'Relative_juvernica_site_pi','juvernica_sites_global'],as_index=False).agg({'Pin_Pis_juvernica': 'mean'})
fianl_df_3_new=fianl_df_3[['CHROM','BIN_START','BIN_END','Pin_Pis_kazak_juvernica',"kazak_juvernica_site_pi_global",'kazak_juvernica_sites_global', 'Relative_kazak_juvernica_site_pi']].dropna(subset=['Pin_Pis_kazak_juvernica']).groupby(['CHROM','BIN_START','BIN_END',"kazak_juvernica_site_pi_global",'Relative_kazak_juvernica_site_pi','kazak_juvernica_sites_global'],as_index=False).agg({'Pin_Pis_kazak_juvernica': 'mean'})
fianl_df_4_new=fianl_df_4[['CHROM','BIN_START','BIN_END','Pin_Pis_kazak_sinapis',"kaz_sin_site_pi_global",'kaz_sin_sites_global', 'Relative_kaz_sin_site_pi']].dropna(subset=['Pin_Pis_kazak_sinapis']).groupby(['CHROM','BIN_START','BIN_END',"kaz_sin_site_pi_global",'Relative_kaz_sin_site_pi','kaz_sin_sites_global'],as_index=False).agg({'Pin_Pis_kazak_sinapis': 'mean'})
fianl_df_5_new=fianl_df_5[['CHROM','BIN_START','BIN_END','Pin_Pis_sinapis',"sinapis_site_pi_global",'sinapis_sites_global', 'Relative_sinapis_site_pi']].dropna(subset=['Pin_Pis_sinapis']).groupby(['CHROM','BIN_START','BIN_END',"sinapis_site_pi_global",'Relative_sinapis_site_pi','sinapis_sites_global'],as_index=False).agg({'Pin_Pis_sinapis': 'mean'})
fianl_df_6_new=fianl_df_6[['CHROM','BIN_START','BIN_END','Pin_Pis_spanish_reali',"spanish_reali_site_pi_global",'spanish_reali_sites_global', 'Relative_spanish_reali_site_pi']].dropna(subset=['Pin_Pis_spanish_reali']).groupby(['CHROM','BIN_START','BIN_END',"spanish_reali_site_pi_global",'Relative_spanish_reali_site_pi','spanish_reali_sites_global'],as_index=False).agg({'Pin_Pis_spanish_reali': 'mean'})
fianl_df_7_new=fianl_df_7[['CHROM','BIN_START','BIN_END','Pin_Pis_spanish_sinapis',"spanish_sinapis_site_pi_global",'spanish_sinapis_sites_global', 'Relative_spanish_sinapis_site_pi']].dropna(subset=['Pin_Pis_spanish_sinapis']).groupby(['CHROM','BIN_START','BIN_END',"spanish_sinapis_site_pi_global",'Relative_spanish_sinapis_site_pi','spanish_sinapis_sites_global'],as_index=False).agg({'Pin_Pis_spanish_sinapis': 'mean'})
fianl_df_8_new=fianl_df_8[['CHROM','BIN_START','BIN_END','Pin_Pis_Swedish_sinapis',"swe_sin_allele_site_pi_global",'swe_sin_allele_sites_global', 'Relative_swe_sin_allele_site_pi']].dropna(subset=['Pin_Pis_Swedish_sinapis']).groupby(['CHROM','BIN_START','BIN_END',"swe_sin_allele_site_pi_global",'Relative_swe_sin_allele_site_pi','swe_sin_allele_sites_global'],as_index=False).agg({'Pin_Pis_Swedish_sinapis': 'mean'})

fianl_df_1_new=fianl_df_1_new[fianl_df_1_new['irish_juvernica_sites_global']>50000]
fianl_df_2_new=fianl_df_2_new[fianl_df_2_new['juvernica_sites_global']>50000]
fianl_df_3_new=fianl_df_3_new[fianl_df_3_new['kazak_juvernica_sites_global']>50000]
fianl_df_4_new=fianl_df_4_new[fianl_df_4_new['kaz_sin_sites_global']>50000]
fianl_df_5_new=fianl_df_5_new[fianl_df_5_new['sinapis_sites_global']>50000]
fianl_df_6_new=fianl_df_6_new[fianl_df_6_new['spanish_reali_sites_global']>50000]
fianl_df_7_new=fianl_df_7_new[fianl_df_7_new['spanish_sinapis_sites_global']>50000]
fianl_df_8_new=fianl_df_8_new[fianl_df_8_new['swe_sin_allele_sites_global']>50000]

print fianl_df_1[fianl_df_1.Relative_irish_juvernica_site_pi > 0].Pin_Pis_irish_juvernica.mean()
print fianl_df_2[fianl_df_2.Relative_juvernica_site_pi > 0].Pin_Pis_juvernica.mean()
print fianl_df_3[fianl_df_3.Relative_kazak_juvernica_site_pi > 0].Pin_Pis_kazak_juvernica.mean()
print fianl_df_4[fianl_df_4.Relative_kaz_sin_site_pi > 0].Pin_Pis_kazak_sinapis.mean()
print fianl_df_5[fianl_df_5.Relative_sinapis_site_pi > 0].Pin_Pis_sinapis.mean()
print fianl_df_6[fianl_df_6.Relative_spanish_reali_site_pi > 0].Pin_Pis_spanish_reali.mean()
print fianl_df_7[fianl_df_7.Relative_spanish_sinapis_site_pi > 0].Pin_Pis_spanish_sinapis.mean()
print fianl_df_8[fianl_df_8.Relative_swe_sin_allele_site_pi > 0].Pin_Pis_Swedish_sinapis.mean()


print fianl_df_2_new[fianl_df_2.Relative_juvernica_site_pi < 0].Pin_Pis_juvernica.mean()
print fianl_df_5_new[fianl_df_5.Relative_sinapis_site_pi < 0].Pin_Pis_sinapis.mean()
print fianl_df_6_new[fianl_df_6.Relative_spanish_reali_site_pi < 0].Pin_Pis_spanish_reali.mean()

sns.jointplot("Relative_irish_juvernica_site_pi", "Pin_Pis_irish_juvernica", data=fianl_df_1_new, kind="reg")
sns.jointplot("Relative_juvernica_site_pi", "Pin_Pis_juvernica", data=fianl_df_2_new, kind="reg")
sns.jointplot("Relative_kazak_juvernica_site_pi", "Pin_Pis_kazak_juvernica", data=fianl_df_3_new, kind="reg")
sns.jointplot("Relative_kaz_sin_site_pi", "Pin_Pis_kazak_sinapis", data=fianl_df_4_new, kind="reg")
sns.jointplot("Relative_sinapis_site_pi", "Pin_Pis_sinapis", data=fianl_df_5_new, kind="reg")
sns.jointplot("Relative_spanish_reali_site_pi", "Pin_Pis_spanish_reali", data=fianl_df_6_new, kind="reg")
sns.jointplot("Relative_spanish_sinapis_site_pi", "Pin_Pis_spanish_sinapis", data=fianl_df_7_new, kind="reg")
sns.jointplot("Relative_swe_sin_allele_site_pi", "Pin_Pis_Swedish_sinapis", data=fianl_df_8_new, kind="reg")


fianl_df_5_new[fianl_df_5_new.Relative_sinapis_site_pi < 0].Pin_Pis_sinapis.mean()
fianl_df_5_new[fianl_df_5_new.Relative_sinapis_site_pi > 0].Pin_Pis_sinapis.mean()

fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
fianl_df_5_new[fianl_df_5_new.Relative_sinapis_site_pi > 0].plot( ax=axes[0,0], kind='scatter', x='Relative_sinapis_site_pi', y='Pin_Pis_sinapis',  color='#72B7A2', alpha=0.2, marker='^', xlim=(-2.0,+2.0),ylim=(0,0.5))
fianl_df_5_new[fianl_df_5_new.Relative_sinapis_site_pi < 0].plot(ax=axes[0,0],kind='scatter', x='Relative_sinapis_site_pi', y='Pin_Pis_sinapis', xlim=(-2.0,+2.0),ylim=(0,0.5), color='#E79776',alpha=0.2,marker="v")
axes[0,0].plot([0, 0], [0, 3], color='#777574', linestyle='--')
axes[0,0].axhline(y=fianl_df_5_new[fianl_df_5_new.Relative_sinapis_site_pi < 0].Pin_Pis_sinapis.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,0].axhline(y=fianl_df_5_new[fianl_df_5_new.Relative_sinapis_site_pi > 0].Pin_Pis_sinapis.mean(), xmin=0.5, xmax=1.0, color='b')

fianl_df_6_new[fianl_df_6_new.Relative_spanish_reali_site_pi > 0].plot( ax=axes[0,1], kind='scatter', x='Relative_spanish_reali_site_pi', y='Pin_Pis_spanish_reali',  color='#72B7A2', alpha=0.2, marker='^', xlim=(-2.0,+2.0),ylim=(0,0.5))
fianl_df_6_new[fianl_df_6_new.Relative_spanish_reali_site_pi < 0].plot(ax=axes[0,1],kind='scatter', x='Relative_spanish_reali_site_pi', y='Pin_Pis_spanish_reali', xlim=(-2.0,+2.0),ylim=(0,0.5), color='#E79776',alpha=0.2,marker="v")
axes[0,1].plot([0, 0], [0, 3], color='#777574', linestyle='--')
axes[0,1].axhline(y=fianl_df_6_new[fianl_df_6_new.Relative_spanish_reali_site_pi < 0].Pin_Pis_spanish_reali.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,1].axhline(y=fianl_df_6_new[fianl_df_6_new.Relative_spanish_reali_site_pi > 0].Pin_Pis_spanish_reali.mean(), xmin=0.5, xmax=1.0, color='b')

fianl_df_2_new[fianl_df_2_new.Relative_juvernica_site_pi > 0].plot( ax=axes[1,0], kind='scatter', x='Relative_juvernica_site_pi', y='Pin_Pis_juvernica',  color='#72B7A2', alpha=0.2, marker='^', xlim=(-2.0,+2.0),ylim=(0,0.5))
fianl_df_2_new[fianl_df_2_new.Relative_juvernica_site_pi < 0].plot(ax=axes[1,0],kind='scatter', x='Relative_juvernica_site_pi', y='Pin_Pis_juvernica', xlim=(-2.0,+2.0),ylim=(0,0.5), color='#E79776',alpha=0.2,marker="v")
axes[1,0].plot([0, 0], [0, 3], color='#777574', linestyle='--')
axes[1,0].axhline(y=fianl_df_2_new[fianl_df_2_new.Relative_juvernica_site_pi < 0].Pin_Pis_juvernica.mean(), xmax=0.0, xmin=0.5, color='r')
axes[1,0].axhline(y=fianl_df_2_new[fianl_df_2_new.Relative_juvernica_site_pi > 0].Pin_Pis_juvernica.mean(), xmin=0.5, xmax=1.0, color='b')



#fianl_df_1['Pin_Pis_irish_juvernica']=fianl_df_1['Pin_Pis_irish_juvernica'].fillna(0.0)

###############################################dn/ds

dnds_table=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/final_dnds_result")
dnds_table_pop=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/final_dnds_populations_result")

dnds_table=dnds_table.drop_duplicates(subset='gene_name')
dnds_table_pop=dnds_table_pop.drop_duplicates(subset='gene_name')

dnds_table=dnds_table[dnds_table.gene_name!='leptidea_sinapisT00000000201.fasta']
dnds_table_pop=dnds_table_pop[dnds_table_pop.gene_name!='leptidea_sinapisT00000000201.fasta']


dnds_table['Ds_L_sin']=dnds_table['ds'].str.split(':').str[1].str.split(',').str[0]
dnds_table['Ds_L_rea']=dnds_table['ds'].str.split(':').str[2].str.split(')').str[0]
dnds_table['Ds_L_juv']=dnds_table['ds'].str.split(':').str[4].str.split(')').str[0]
dnds_table['Dn_L_sin']=dnds_table['dn'].str.split(':').str[1].str.split(',').str[0]
dnds_table['Dn_L_rea']=dnds_table['dn'].str.split(':').str[2].str.split(')').str[0]
dnds_table['Dn_L_juv']=dnds_table['dn'].str.split(':').str[4].str.split(')').str[0]

dnds_table=dnds_table[['gene_name','Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv', 'N_sites', 'S_sites']]
dnds_table[['Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']]=dnds_table[['Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']].astype('float').round(4)

dnds_table['s_subs_L_sin']=dnds_table['Ds_L_sin']*dnds_table['S_sites']
dnds_table['s_subs_L_rea']=dnds_table['Ds_L_rea']*dnds_table['S_sites']
dnds_table['s_subs_L_juv']=dnds_table['Ds_L_juv']*dnds_table['S_sites']
dnds_table['n_subs_L_sin']=dnds_table['Dn_L_sin']*dnds_table['N_sites']
dnds_table['n_subs_L_rea']=dnds_table['Dn_L_rea']*dnds_table['N_sites']
dnds_table['n_subs_L_juv']=dnds_table['Dn_L_juv']*dnds_table['N_sites']
dnds_table['s_subs_L_sin']=dnds_table['s_subs_L_sin'].round(0)
dnds_table['s_subs_L_rea']=dnds_table['s_subs_L_rea'].round(0)
dnds_table['s_subs_L_juv']=dnds_table['s_subs_L_juv'].round(0)
dnds_table['n_subs_L_sin']=dnds_table['n_subs_L_sin'].round(0)
dnds_table['n_subs_L_rea']=dnds_table['n_subs_L_rea'].round(0)
dnds_table['n_subs_L_juv']=dnds_table['n_subs_L_juv'].round(0)

dnds_table['W_L_sin']=dnds_table['Dn_L_sin']/dnds_table['Ds_L_sin']
dnds_table['W_L_rea']=dnds_table['Dn_L_rea']/dnds_table['Ds_L_rea']
dnds_table['W_L_juv']=dnds_table['Dn_L_juv']/dnds_table['Ds_L_juv']


dnds_table['W_L_sin']=dnds_table['W_L_sin'].replace(np.inf, np.nan)
dnds_table['W_L_rea']=dnds_table['W_L_rea'].replace(np.inf, np.nan)
dnds_table['W_L_juv']=dnds_table['W_L_juv'].replace(np.inf, np.nan)


dnds_table=dnds_table.rename(columns={'gene_name': 'geneID'})
dnds_table['geneID']=dnds_table['geneID'].str.split('.').str[0]

dnds_table_gene=join_raw_data_base([annotation,dnds_table,],['geneID'],'inner')

dnds_table_gene=join_raw_data_base([dnds_table_gene,fianl_df_2[['geneID','Pin_juvernica','Pis_juvernica','Pin_Pis_juvernica', "Pi_syn_juvernica","Pi_nonsyn_juvernica", "syn_sites_juvernica", "nonsyn_sites_juvernica"]],fianl_df_5[['geneID','Pin_sinapis','Pis_sinapis','Pin_Pis_sinapis', 'Pi_syn_sinapis','Pi_nonsyn_sinapis',"syn_sites_sinapis", "nonsyn_sites_sinapis"]],fianl_df_6[['geneID','Pin_spanish_reali','Pis_spanish_reali','Pin_Pis_spanish_reali', 'Pi_syn_spanish_reali', 'Pi_nonsyn_spanish_reali', "syn_sites_spanish_reali", "nonsyn_sites_spanish_reali"]]],['geneID'],'inner')

dnds_table_gene['alpha_Sin']=1-(dnds_table_gene['Pin_Pis_sinapis']/dnds_table_gene['W_L_sin'])
dnds_table_gene['alpha_Rea']=1-(dnds_table_gene['Pin_Pis_spanish_reali']/dnds_table_gene['W_L_rea'])
dnds_table_gene['alpha_Juv']=1-(dnds_table_gene['Pin_Pis_juvernica']/dnds_table_gene['W_L_juv'])

dnds_table_gene['alpha_Sin']=dnds_table_gene['alpha_Sin'].replace(-(np.inf), np.nan)
dnds_table_gene['alpha_Rea']=dnds_table_gene['alpha_Rea'].replace(-(np.inf), np.nan)
dnds_table_gene['alpha_Juv']=dnds_table_gene['alpha_Juv'].replace(-(np.inf), np.nan)
dnds_table_gene['alpha_Sin']=dnds_table_gene['alpha_Sin'].replace((np.inf), np.nan)
dnds_table_gene['alpha_Rea']=dnds_table_gene['alpha_Rea'].replace((np.inf), np.nan)
dnds_table_gene['alpha_Juv']=dnds_table_gene['alpha_Juv'].replace((np.inf), np.nan)



dnds_table_gene['W_alpha_Sin']=dnds_table_gene['alpha_Sin']*dnds_table_gene['W_L_sin']
dnds_table_gene['W_alpha_Rea']=dnds_table_gene['alpha_Rea']*dnds_table_gene['W_L_rea']
dnds_table_gene['W_alpha_Juv']=dnds_table_gene['alpha_Juv']*dnds_table_gene['W_L_juv']


dnds_table_gene['DOS_Sin']=(dnds_table_gene['Dn_L_sin']/(dnds_table_gene['Dn_L_sin']+dnds_table_gene['Ds_L_sin'])) - (dnds_table_gene['Pin_sinapis']/(dnds_table_gene['Pin_sinapis']+dnds_table_gene['Pis_sinapis']))
dnds_table_gene['DOS_Rea']=(dnds_table_gene['Dn_L_rea']/(dnds_table_gene['Dn_L_rea']+dnds_table_gene['Ds_L_rea'])) - (dnds_table_gene['Pin_spanish_reali']/(dnds_table_gene['Pin_spanish_reali']+dnds_table_gene['Pis_spanish_reali']))
dnds_table_gene['DOS_Juv']=(dnds_table_gene['Dn_L_juv']/(dnds_table_gene['Dn_L_juv']+dnds_table_gene['Ds_L_juv'])) - (dnds_table_gene['Pin_juvernica']/(dnds_table_gene['Pin_juvernica']+dnds_table_gene['Pis_juvernica']))


dnds_table_gene.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/dnds_table_gene.csv')
#######################dnds_table_windows is for genes #############################

dnds_table_windows=dnds_table_gene.groupby(['CHROM','BIN_START','BIN_END', 'auto_allo'],as_index=False).agg({'N_sites': 'sum','S_sites': 'sum', 's_subs_L_sin': 'sum','s_subs_L_rea': 'sum','s_subs_L_juv': 'sum','n_subs_L_sin': 'sum','n_subs_L_rea': 'sum','n_subs_L_juv': 'sum',"syn_sites_juvernica": 'sum',"nonsyn_sites_juvernica": 'sum',"syn_sites_sinapis": 'sum',"nonsyn_sites_sinapis": 'sum',"syn_sites_spanish_reali" : 'sum',"nonsyn_sites_spanish_reali": 'sum','Pi_syn_juvernica': 'sum','Pi_syn_sinapis': 'sum','Pi_syn_spanish_reali': 'sum','Pi_nonsyn_juvernica': 'sum','Pi_nonsyn_sinapis': 'sum','Pi_nonsyn_spanish_reali': 'sum'})

dnds_table_windows['Dn_L_sin']=dnds_table_windows['n_subs_L_sin']/dnds_table_windows['N_sites']
dnds_table_windows['Dn_L_rea']=dnds_table_windows['n_subs_L_rea']/dnds_table_windows['N_sites']
dnds_table_windows['Dn_L_juv']=dnds_table_windows['n_subs_L_juv']/dnds_table_windows['N_sites']


dnds_table_windows['Ds_L_sin']=dnds_table_windows['s_subs_L_sin']/dnds_table_windows['S_sites']
dnds_table_windows['Ds_L_rea']=dnds_table_windows['s_subs_L_rea']/dnds_table_windows['S_sites']
dnds_table_windows['Ds_L_juv']=dnds_table_windows['s_subs_L_juv']/dnds_table_windows['S_sites']



dnds_table_windows['W_L_sin']=dnds_table_windows['Dn_L_sin']/dnds_table_windows['Ds_L_sin']
dnds_table_windows['W_L_rea']=dnds_table_windows['Dn_L_rea']/dnds_table_windows['Ds_L_rea']
dnds_table_windows['W_L_juv']=dnds_table_windows['Dn_L_juv']/dnds_table_windows['Ds_L_juv']

dnds_table_windows['W_L_sin']=dnds_table_windows['W_L_sin'].replace(np.inf, np.nan)
dnds_table_windows['W_L_rea']=dnds_table_windows['W_L_rea'].replace(np.inf, np.nan)
dnds_table_windows['W_L_juv']=dnds_table_windows['W_L_juv'].replace(np.inf, np.nan)

dnds_table_windows.W_L_sin[dnds_table_windows['W_L_sin'] > 3.0]=np.nan
dnds_table_windows.W_L_rea[dnds_table_windows['W_L_rea'] > 3.0]=np.nan
dnds_table_windows.W_L_juv[dnds_table_windows['W_L_juv'] > 3.0]=np.nan


dnds_table_windows.W_L_sin[dnds_table_windows['W_L_sin'] < -3.0]=np.nan
dnds_table_windows.W_L_rea[dnds_table_windows['W_L_rea'] < -3.0]=np.nan
dnds_table_windows.W_L_juv[dnds_table_windows['W_L_juv'] < -3.0]=np.nan

dnds_table_windows[['Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']]=dnds_table_windows[['Ds_L_sin','Ds_L_rea','Ds_L_juv','Dn_L_sin','Dn_L_rea','Dn_L_juv']].astype('float').round(3)

dnds_table_windows['sites_sinapis']=dnds_table_windows['syn_sites_sinapis']+dnds_table_windows['nonsyn_sites_sinapis']
dnds_table_windows['sites_spanish_reali']=dnds_table_windows['syn_sites_spanish_reali']+dnds_table_windows['nonsyn_sites_spanish_reali']
dnds_table_windows['sites_juvernica']=dnds_table_windows['syn_sites_juvernica']+dnds_table_windows['nonsyn_sites_juvernica']


dnds_table_windows=dnds_table_windows[(dnds_table_windows['sites_sinapis']>1000)]
dnds_table_windows=dnds_table_windows[(dnds_table_windows['sites_spanish_reali']>1000)]
dnds_table_windows=dnds_table_windows[(dnds_table_windows['sites_juvernica']>1000)]


########error  in the past
dnds_table_windows['Pin_L_sin']=dnds_table_windows['Pi_nonsyn_sinapis']/dnds_table_windows['nonsyn_sites_sinapis']
dnds_table_windows['Pin_L_rea']=dnds_table_windows['Pi_nonsyn_spanish_reali']/dnds_table_windows['nonsyn_sites_spanish_reali']
dnds_table_windows['Pin_L_juv']=dnds_table_windows['Pi_nonsyn_juvernica']/dnds_table_windows['nonsyn_sites_juvernica']




dnds_table_windows['Pis_L_sin']=dnds_table_windows['Pi_syn_sinapis']/dnds_table_windows['syn_sites_sinapis']
dnds_table_windows['Pis_L_rea']=dnds_table_windows['Pi_syn_spanish_reali']/dnds_table_windows['syn_sites_spanish_reali']
dnds_table_windows['Pis_L_juv']=dnds_table_windows['Pi_syn_juvernica']/dnds_table_windows['syn_sites_juvernica']

dnds_table_windows['Pin_Pis_sinapis']=dnds_table_windows['Pin_L_sin']/dnds_table_windows['Pis_L_sin']
dnds_table_windows['Pin_Pis_spanish_reali']=dnds_table_windows['Pin_L_rea']/dnds_table_windows['Pis_L_rea']
dnds_table_windows['Pin_Pis_juvernica']=dnds_table_windows['Pin_L_juv']/dnds_table_windows['Pis_L_juv']

dnds_table_windows['Pin_Pis_sinapis']=dnds_table_windows['Pin_Pis_sinapis'].replace(-(np.inf), np.nan)
dnds_table_windows['Pin_Pis_spanish_reali']=dnds_table_windows['Pin_Pis_spanish_reali'].replace(-(np.inf), np.nan)
dnds_table_windows['Pin_Pis_juvernica']=dnds_table_windows['Pin_Pis_juvernica'].replace(-(np.inf), np.nan)
dnds_table_windows['Pin_Pis_sinapis']=dnds_table_windows['Pin_Pis_sinapis'].replace((np.inf), np.nan)
dnds_table_windows['Pin_Pis_spanish_reali']=dnds_table_windows['Pin_Pis_spanish_reali'].replace((np.inf), np.nan)
dnds_table_windows['Pin_Pis_juvernica']=dnds_table_windows['Pin_Pis_juvernica'].replace((np.inf), np.nan)

dnds_table_windows.Pin_Pis_sinapis[dnds_table_windows['Pin_Pis_sinapis'] > 3.0]=np.nan
dnds_table_windows.Pin_Pis_spanish_reali[dnds_table_windows['Pin_Pis_spanish_reali'] > 3.0]=np.nan
dnds_table_windows.Pin_Pis_juvernica[dnds_table_windows['Pin_Pis_juvernica'] > 3.0]=np.nan


dnds_table_windows.Pin_Pis_sinapis[dnds_table_windows['Pin_Pis_sinapis'] < -3.0]=np.nan
dnds_table_windows.Pin_Pis_spanish_reali[dnds_table_windows['Pin_Pis_spanish_reali'] < -3.0]=np.nan
dnds_table_windows.Pin_Pis_juvernica[dnds_table_windows['Pin_Pis_juvernica'] < -3.0]=np.nan


dnds_table_windows['alpha_Sin']=1-(dnds_table_windows['Pin_Pis_sinapis']/dnds_table_windows['W_L_sin'])
dnds_table_windows['alpha_Rea']=1-(dnds_table_windows['Pin_Pis_spanish_reali']/dnds_table_windows['W_L_rea'])
dnds_table_windows['alpha_Juv']=1-(dnds_table_windows['Pin_Pis_juvernica']/dnds_table_windows['W_L_juv'])

dnds_table_windows['alpha_Sin']=dnds_table_windows['alpha_Sin'].replace(-(np.inf), np.nan)
dnds_table_windows['alpha_Rea']=dnds_table_windows['alpha_Rea'].replace(-(np.inf), np.nan)
dnds_table_windows['alpha_Juv']=dnds_table_windows['alpha_Juv'].replace(-(np.inf), np.nan)
dnds_table_windows['alpha_Sin']=dnds_table_windows['alpha_Sin'].replace((np.inf), np.nan)
dnds_table_windows['alpha_Rea']=dnds_table_windows['alpha_Rea'].replace((np.inf), np.nan)
dnds_table_windows['alpha_Juv']=dnds_table_windows['alpha_Juv'].replace((np.inf), np.nan)



dnds_table_windows['W_alpha_Sin']=dnds_table_windows['alpha_Sin']*dnds_table_windows['W_L_sin']
dnds_table_windows['W_alpha_Rea']=dnds_table_windows['alpha_Rea']*dnds_table_windows['W_L_rea']
dnds_table_windows['W_alpha_Juv']=dnds_table_windows['alpha_Juv']*dnds_table_windows['W_L_juv']

dnds_table_windows.W_alpha_Sin[dnds_table_windows['W_alpha_Sin'] > 10.0]=np.nan
dnds_table_windows.W_alpha_Rea[dnds_table_windows['W_alpha_Rea'] > 10.0]=np.nan
dnds_table_windows.W_alpha_Juv[dnds_table_windows['W_alpha_Juv'] > 10.0]=np.nan

dnds_table_windows.W_alpha_Sin[dnds_table_windows['W_alpha_Sin'] < -3.0]=np.nan
dnds_table_windows.W_alpha_Rea[dnds_table_windows['W_alpha_Rea'] < -3.0]=np.nan
dnds_table_windows.W_alpha_Juv[dnds_table_windows['W_alpha_Juv'] < -3.0]=np.nan




dnds_table_windows['DOS_Sin']=(dnds_table_windows['Dn_L_sin']/(dnds_table_windows['Dn_L_sin']+dnds_table_windows['Ds_L_sin'])) - (dnds_table_windows['Pin_L_sin']/(dnds_table_windows['Pin_L_sin']+dnds_table_windows['Pis_L_sin']))
dnds_table_windows['DOS_Rea']=(dnds_table_windows['Dn_L_rea']/(dnds_table_windows['Dn_L_rea']+dnds_table_windows['Ds_L_rea'])) - (dnds_table_windows['Pin_L_rea']/(dnds_table_windows['Pin_L_rea']+dnds_table_windows['Pis_L_rea']))
dnds_table_windows['DOS_Juv']=(dnds_table_windows['Dn_L_juv']/(dnds_table_windows['Dn_L_juv']+dnds_table_windows['Ds_L_juv'])) - (dnds_table_windows['Pin_L_juv']/(dnds_table_windows['Pin_L_juv']+dnds_table_windows['Pis_L_juv']))


pop_df_2['gene_density_juvernica']=((pop_df_2['juvernica_sites_codon_3']+pop_df_2['juvernica_sites_codon_2']+pop_df_2['juvernica_sites_codon_1'])/pop_df_2['juvernica_sites_global'])*100
pop_df_5['gene_density_sinapis']=((pop_df_5['sinapis_sites_codon_3']+pop_df_5['sinapis_sites_codon_2']+pop_df_5['sinapis_sites_codon_1'])/pop_df_5['sinapis_sites_global'])*100
pop_df_6['gene_density_reali']=((pop_df_6['spanish_reali_sites_codon_3']+pop_df_6['spanish_reali_sites_codon_2']+pop_df_6['spanish_reali_sites_codon_1'])/pop_df_6['spanish_reali_sites_global'])*100

pop_df_2['coding_sites']=((pop_df_2['juvernica_sites_codon_3']+pop_df_2['juvernica_sites_codon_2']+pop_df_2['juvernica_sites_codon_1']))
pop_df_5['coding_sites']=((pop_df_5['sinapis_sites_codon_3']+pop_df_5['sinapis_sites_codon_2']+pop_df_5['sinapis_sites_codon_1']))
pop_df_6['coding_sites']=((pop_df_6['spanish_reali_sites_codon_3']+pop_df_6['spanish_reali_sites_codon_2']+pop_df_6['spanish_reali_sites_codon_1']))


pop_df_2_1=pop_df_2[['CHROM','BIN_START','BIN_END','juvernica_site_pi_global','Relative_juvernica_site_pi', 'juvernica_site_pi_codon_4d','gene_density_juvernica','coding_sites']]
pop_df_5_1=pop_df_5[['CHROM','BIN_START','BIN_END','sinapis_site_pi_global','Relative_sinapis_site_pi', 'sinapis_site_pi_codon_4d','gene_density_sinapis','coding_sites']]
pop_df_6_1=pop_df_6[['CHROM','BIN_START','BIN_END','spanish_reali_site_pi_global','Relative_spanish_reali_site_pi', 'spanish_reali_site_pi_codon_4d','gene_density_reali','coding_sites']]


pop_df_2_1['Relative_juvernica_site_pi_log_transform']=np.log10(pop_df_2_1.Relative_juvernica_site_pi)
pop_df_5_1['Relative_sinapis_site_pi_log_transform']=np.log10(pop_df_5_1.Relative_sinapis_site_pi)
pop_df_6_1['Relative_spanish_reali_site_pi_log_transform']=np.log10(pop_df_6_1.Relative_spanish_reali_site_pi)

pi_colours=sns.color_palette(['#C8C800','#006400','#0000FF','#FF8C00','#FF0000','#C04000'])
colour_all=['#C8C800','#006400','#FF8C00','#0000FF','#FF0000','#C04000']
pi_colours_species=sns.color_palette(['#E79676', '#72B7A1'])
global_c=[sns.color_palette(sns.color_palette("Set2"))[2]]

sns.jointplot("gene_density", "juvernica_site_pi_codon_4d", data=pop_df_2[pop_df_2['juvernica_sites_global']>40000], kind="reg", color='#24cc24')
sns.jointplot("gene_density", "sinapis_site_pi_codon_4d", data=pop_df_5[pop_df_5['sinapis_sites_global']>40000], kind="reg", color='#00e1ff')
sns.jointplot("gene_density", "spanish_reali_site_pi_codon_4d", data=pop_df_6[pop_df_6['spanish_reali_sites_global']>40000], kind="reg", color='#a80000')




dnds_table_windows_2=join_raw_data_base([dnds_table_windows,pop_df_2_1,pop_df_5_1,pop_df_6_1],['CHROM','BIN_START','BIN_END'],'outer')



sns.jointplot("Ds_L_juv", "juvernica_site_pi_global", data=dnds_table_windows_2, kind="reg", color='#24cc24')
sns.jointplot("Ds_L_sin", "sinapis_site_pi_global", data=dnds_table_windows_2, kind="reg", color='#00e1ff')
sns.jointplot("Ds_L_rea", "spanish_reali_site_pi_global", data=dnds_table_windows_2, kind="reg", color='#a80000')


#dnds_table_windows_2.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/dnds_table_windows.csv')
#dnds_table_windows_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/dnds_table_windows.csv')

dnds_table_windows_2['Relative_juvernica_site_pi_direction']=''
dnds_table_windows_2['Relative_juvernica_site_pi_direction'][dnds_table_windows_2.Relative_juvernica_site_pi > 0]='high'
dnds_table_windows_2['Relative_juvernica_site_pi_direction'][dnds_table_windows_2.Relative_juvernica_site_pi < 0]='low'

dnds_table_windows_2['Relative_sinapis_site_pi_direction']=''
dnds_table_windows_2['Relative_sinapis_site_pi_direction'][dnds_table_windows_2.Relative_sinapis_site_pi > 0]='high'
dnds_table_windows_2['Relative_sinapis_site_pi_direction'][dnds_table_windows_2.Relative_sinapis_site_pi < 0]='low'

dnds_table_windows_2['Relative_spanish_reali_site_pi_direction']=''
dnds_table_windows_2['Relative_spanish_reali_site_pi_direction'][dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]='high'
dnds_table_windows_2['Relative_spanish_reali_site_pi_direction'][dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]='low'




pi_ps=list(dnds_table_windows_2.Pin_Pis_sinapis)+list(dnds_table_windows_2.Pin_Pis_spanish_reali)+list(dnds_table_windows_2.Pin_Pis_juvernica)
W=list(dnds_table_windows_2.W_L_sin)+list(dnds_table_windows_2.W_L_rea)+list(dnds_table_windows_2.W_L_juv)
alpha=list(dnds_table_windows_2.alpha_Sin)+list(dnds_table_windows_2.alpha_Rea)+list(dnds_table_windows_2.alpha_Juv)
W_alpha=list(dnds_table_windows_2.W_alpha_Sin)+list(dnds_table_windows_2.W_alpha_Rea)+list(dnds_table_windows_2.W_alpha_Juv)
N_sites=list(dnds_table_windows_2.N_sites)+list(dnds_table_windows_2.N_sites)+list(dnds_table_windows_2.N_sites)
S_sites=list(dnds_table_windows_2.S_sites)+list(dnds_table_windows_2.S_sites)+list(dnds_table_windows_2.S_sites)
gene_density=list(dnds_table_windows_2.gene_density_sinapis)+list(dnds_table_windows_2.gene_density_reali)+list(dnds_table_windows_2.gene_density_juvernica)
DOS=list(dnds_table_windows_2.DOS_Sin)+list(dnds_table_windows_2.DOS_Rea)+list(dnds_table_windows_2.DOS_Juv)

relative=list(dnds_table_windows_2.Relative_sinapis_site_pi_direction)+list(dnds_table_windows_2.Relative_spanish_reali_site_pi_direction)+list(dnds_table_windows_2.Relative_juvernica_site_pi_direction)
species=['L_sin']*len(list(dnds_table_windows_2.Pin_Pis_sinapis))+['L_Rea']*len(list(dnds_table_windows_2.Pin_Pis_spanish_reali))+['L_Juv']*len(list(dnds_table_windows_2.Pin_Pis_juvernica))
violin_df = pandas.DataFrame()
violin_df['pn_ps']=pi_ps
violin_df['alpha']=alpha
violin_df['gene_density']=np.array(gene_density)
violin_df['W']=W
violin_df['W_alpha']=W_alpha
violin_df['DOS']=DOS

violin_df['N_sites']=N_sites
violin_df['S_sites']=S_sites
violin_df['relative']=relative
violin_df['species']=species
violin_df=violin_df[violin_df.relative!='']
violin_df['sites']=violin_df['S_sites']+violin_df['N_sites']


ax=sns.boxplot(x="species", y="pn_ps", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
ax.legend_.remove()
ax.set(ylim=(0,0.7))

sns.boxplot(x="species", y="W", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'] ).set(ylim=(-1,1))

ax=sns.boxplot(x="species", y="W_alpha", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
ax.legend_.remove()
ax.set(ylim=(-1,1))

sns.boxplot(x="species", y="DOS", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high']).set(ylim=(-3,3))

sns.boxplot(x="species", y="N_sites", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
sns.boxplot(x="species", y="S_sites", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])

sns.boxplot(x="species", y="alpha", data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high']).set(ylim=(-1.5,1.2))

fig, (axe1, axe2) = plt.subplots(ncols=2, sharey=True)
a=sns.boxplot(x="species", y="alpha", ax=axe1,hue="relative" ,data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
a.legend_.remove()
a.set(ylim=(-2.0,1.2))
b=sns.boxplot(x="species", y="alpha",ax=axe2,data=violin_df, palette=global_c, showfliers=False, hue_order=['low','high'])
b.legend_.remove()
b.set(ylim=(-2.0,1.2))

fig, (axe1, axe2) = plt.subplots(ncols=2, sharey=True)
a=sns.boxplot(x="species", y="W_alpha", ax=axe1,hue="relative" ,data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
a.legend_.remove()
a.set(ylim=(-2.0,1.2))
b=sns.boxplot(x="species", y="W_alpha",ax=axe2,data=violin_df, palette=global_c, showfliers=False, hue_order=['low','high'])
b.legend_.remove()
b.set(ylim=(-2.0,1.2))




ax=sns.boxplot(x="species", y="gene_density", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
ax.legend_.remove()


import scipy

print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['Pin_Pis_sinapis'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['Pin_Pis_sinapis'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['Pin_Pis_spanish_reali'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['Pin_Pis_spanish_reali'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['Pin_Pis_juvernica'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['Pin_Pis_juvernica'].dropna(how='all')))



print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['W_L_sin'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['W_L_sin'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['W_L_rea'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['W_L_rea'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['W_L_juv'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['W_L_juv'].dropna(how='all')))

print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['Pin_Pis_sinapis'].median()
print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['Pin_Pis_sinapis'].median()
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['Pin_Pis_spanish_reali'].median()
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['Pin_Pis_spanish_reali'].median()
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['Pin_Pis_juvernica'].median()
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['Pin_Pis_juvernica'].median()

print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['alpha_Sin'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['alpha_Sin'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['alpha_Rea'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['alpha_Rea'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['alpha_Juv'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['alpha_Juv'].dropna(how='all')))


print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['W_alpha_Sin'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['W_alpha_Sin'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['W_alpha_Rea'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['W_alpha_Rea'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['W_alpha_Juv'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['W_alpha_Juv'].dropna(how='all')))


print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['N_sites'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['N_sites'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['N_sites'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['N_sites'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['N_sites'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['N_sites'].dropna(how='all')))


print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['S_sites'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['S_sites'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['S_sites'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['S_sites'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['S_sites'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['S_sites'].dropna(how='all')))


dnds_table_windows_2['sites']=dnds_table_windows_2['N_sites']+dnds_table_windows_2['S_sites']

print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['gene_density_sinapis'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['gene_density_sinapis'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['gene_density_reali'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['gene_density_reali'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['gene_density_juvernica'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['gene_density_juvernica'].dropna(how='all')))






print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0.0][['Pin_Pis_sinapis']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0.0][['Pin_Pis_sinapis']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0.0][['Pin_Pis_spanish_reali']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0.0][['Pin_Pis_spanish_reali']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0.0][['Pin_Pis_juvernica']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0.0][['Pin_Pis_juvernica']].dropna().mean().round(3)


print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0.0][['W_L_sin']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0.0][['W_L_sin']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0.0][['W_L_rea']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0.0][['W_L_rea']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0.0][['W_L_juv']].dropna().mean().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0.0][['W_L_juv']].dropna().mean().round(3)







print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0.0][['alpha_Sin']].dropna().median().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0.0][['alpha_Sin']].dropna().median().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0.0][['alpha_Rea']].dropna().median().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0.0][['alpha_Rea']].dropna().median().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0.0][['alpha_Juv']].dropna().median().round(3)
print dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0.0][['alpha_Juv']].dropna().median().round(3)

print dnds_table_windows_2[['alpha_Sin']].dropna().median().round(3)
print dnds_table_windows_2[['alpha_Rea']].dropna().median().round(3)
print dnds_table_windows_2[['alpha_Juv']].dropna().median().round(3)


dnds_table_windows_2[['alpha_Sin','alpha_Rea','alpha_Juv']].plot.box(color=color, showmeans=True)

color = dict(boxes='DarkGreen', whiskers='DarkOrange', medians='DarkBlue', caps='Gray')
dnds_table_windows_2[['W_alpha_Sin','W_alpha_Rea','W_alpha_Juv']].plot.box(color=color, showmeans=True)










sns.jointplot("sinapis_site_pi_global", "gene_counts", data=pop_df_5, kind="reg")
sns.jointplot("sinapis_site_pi_global", "gene_counts", data=pop_df_5, kind="reg")



sns.jointplot("W_alpha_Sin", "W_alpha_Rea", data=dnds_table_gene, kind="reg")
sns.jointplot("W_alpha_Juv", "W_alpha_Rea", data=dnds_table_gene, kind="reg")
sns.jointplot("W_alpha_Juv", "W_alpha_Sin", data=dnds_table_gene, kind="reg")


#######################window stats pnps#############################

fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].plot( ax=axes[0,0], kind='scatter', x='Relative_sinapis_site_pi', y='Pin_Pis_sinapis',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-5.0,+5.0),ylim=(0,10.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].plot(ax=axes[0,0],kind='scatter', x='Relative_sinapis_site_pi', y='Pin_Pis_sinapis', xlim=(-5.0,+5.0),ylim=(0,10.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,0].plot([0, 0], [0, 3], color='#777574', linestyle='--')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].Pin_Pis_sinapis.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].Pin_Pis_sinapis.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].plot( ax=axes[0,1], kind='scatter', x='Relative_spanish_reali_site_pi', y='Pin_Pis_spanish_reali',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-5.0,+5.0),ylim=(0,10.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].plot(ax=axes[0,1],kind='scatter', x='Relative_spanish_reali_site_pi', y='Pin_Pis_spanish_reali', xlim=(-5.0,+5.0),ylim=(0,10.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,1].plot([0, 0], [0, 3], color='#777574', linestyle='--')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].Pin_Pis_spanish_reali.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].Pin_Pis_spanish_reali.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].plot( ax=axes[1,0], kind='scatter', x='Relative_juvernica_site_pi', y='Pin_Pis_juvernica',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-5.0,+5.0),ylim=(0,10.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].plot(ax=axes[1,0],kind='scatter', x='Relative_juvernica_site_pi', y='Pin_Pis_juvernica', xlim=(-5.0,+5.0),ylim=(0,10.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[1,0].plot([0, 0], [0, 3], color='#777574', linestyle='--')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].Pin_Pis_juvernica.mean(), xmax=0.0, xmin=0.5, color='r')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].Pin_Pis_juvernica.mean(), xmin=0.5, xmax=1.0, color='b')

##########################window stats dnds plot######################################

fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].plot( ax=axes[0,0], kind='scatter', x='Relative_sinapis_site_pi', y='W_L_sin',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].plot(ax=axes[0,0],kind='scatter', x='Relative_sinapis_site_pi', y='W_L_sin', xlim=(-2.0,+2.0),ylim=(0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,0].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].W_L_sin.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].W_L_sin.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].plot( ax=axes[0,1], kind='scatter', x='Relative_spanish_reali_site_pi', y='W_L_rea',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].plot(ax=axes[0,1],kind='scatter', x='Relative_spanish_reali_site_pi', y='W_L_rea', xlim=(-2.0,+2.0),ylim=(0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,1].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].W_L_rea.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].W_L_rea.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].plot( ax=axes[1,0], kind='scatter', x='Relative_juvernica_site_pi', y='W_L_juv',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].plot(ax=axes[1,0],kind='scatter', x='Relative_juvernica_site_pi', y='W_L_juv', xlim=(-2.0,+2.0),ylim=(0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[1,0].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].W_L_juv.mean(), xmax=0.0, xmin=0.5, color='r')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].W_L_juv.mean(), xmin=0.5, xmax=1.0, color='b')


#############################W_alpha_plot
fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].plot( ax=axes[0,0], kind='scatter', x='Relative_sinapis_site_pi', y='W_alpha_Sin',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(-3.0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].plot(ax=axes[0,0],kind='scatter', x='Relative_sinapis_site_pi', y='W_alpha_Sin', xlim=(-2.0,+2.0),ylim=(-3.0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,0].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].W_alpha_Sin.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].W_alpha_Sin.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].plot( ax=axes[0,1], kind='scatter', x='Relative_spanish_reali_site_pi', y='W_alpha_Rea',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(-3.0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].plot(ax=axes[0,1],kind='scatter', x='Relative_spanish_reali_site_pi', y='W_alpha_Rea', xlim=(-2.0,+2.0),ylim=(-3.0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,1].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].W_alpha_Rea.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].W_alpha_Rea.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].plot( ax=axes[1,0], kind='scatter', x='Relative_juvernica_site_pi', y='W_alpha_Juv',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(-3.0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].plot(ax=axes[1,0],kind='scatter', x='Relative_juvernica_site_pi', y='W_alpha_Juv', xlim=(-2.0,+2.0),ylim=(-3.0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[1,0].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].W_alpha_Juv.mean(), xmax=0.0, xmin=0.5, color='r')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].W_alpha_Juv.mean(), xmin=0.5, xmax=1.0, color='b')

#############################DOS_alpha_plot
fig, axes = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True)
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].plot( ax=axes[0,0], kind='scatter', x='Relative_sinapis_site_pi', y='DOS_Sin',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(-3.0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].plot(ax=axes[0,0],kind='scatter', x='Relative_sinapis_site_pi', y='DOS_Sin', xlim=(-2.0,+2.0),ylim=(-3.0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,0].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0].DOS_Sin.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0].DOS_Sin.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].plot( ax=axes[0,1], kind='scatter', x='Relative_spanish_reali_site_pi', y='DOS_Rea',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(-3.0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].plot(ax=axes[0,1],kind='scatter', x='Relative_spanish_reali_site_pi', y='DOS_Rea', xlim=(-2.0,+2.0),ylim=(-3.0,3.0), color='#E79776',alpha=0.2,s=60, marker="v")
axes[0,1].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0].DOS_Rea.mean(), xmax=0.0, xmin=0.5, color='r')
axes[0,1].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0].DOS_Rea.mean(), xmin=0.5, xmax=1.0, color='b')

dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].plot( ax=axes[1,0], kind='scatter', x='Relative_juvernica_site_pi', y='DOS_Juv',  color='#72B7A2', alpha=0.2, s=60, marker='^', xlim=(-2.0,+2.0),ylim=(-3.0,3.0))
dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].plot(ax=axes[1,0],kind='scatter', x='Relative_juvernica_site_pi', y='DOS_Juv', xlim=(-2.0,+2.0),ylim=(-3.0,3.0), color='#E79776',alpha=0.2, s=60, marker="v")
axes[1,0].plot([0, 0], [-3, 3], color='#777574', linestyle='--')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0].DOS_Juv.mean(), xmax=0.0, xmin=0.5, color='r')
axes[1,0].axhline(y=dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0].DOS_Juv.mean(), xmin=0.5, xmax=1.0, color='b')

########################################################################

FOR_DOS=join_raw_data_base([dnds_table_gene,pop_df_2_1,pop_df_5_1,pop_df_6_1],['CHROM','BIN_START','BIN_END'],'inner') #15217
DM_genes=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/gene_ortholog_flybase/overlap_new.csv", names=['DM_gene','LS_gene'], sep=' ')
DM_genes=DM_genes.rename(columns={'LS_gene': 'geneID'})

sns.jointplot("Relative_sinapis_site_pi", "DOS_Sin", data=FOR_DOS[FOR_DOS.s_subs_L_sin>=3], kind="reg")
sns.jointplot("Relative_spanish_reali_site_pi", "DOS_Rea", data=FOR_DOS[FOR_DOS.s_subs_L_rea>=3], kind="reg")
sns.jointplot("Relative_juvernica_site_pi", "DOS_Juv", data=FOR_DOS[FOR_DOS.s_subs_L_juv>=3], kind="reg")



#############################McDonald–Kreitman test#####################

FOR_MK=dnds_table_gene[['CHROM','geneID','s_subs_L_sin','s_subs_L_rea','s_subs_L_juv','n_subs_L_sin','n_subs_L_rea','n_subs_L_juv', 'Pi_syn_sinapis','Pi_nonsyn_sinapis','Pi_syn_spanish_reali','Pi_nonsyn_spanish_reali', 'Pi_syn_juvernica', 'Pi_nonsyn_juvernica', 'N_sites','S_sites', 'DOS_Sin', 'DOS_Rea', 'DOS_Juv']]

FOR_MK_from_R=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/dn_ds_species.csv")
FOR_MK_from_R=FOR_MK_from_R.rename(columns={'gene': 'geneID'})
FOR_MK_from_R['geneID']=FOR_MK_from_R['geneID'].str.split('.').str[0]

FOR_MK=join_raw_data_base([FOR_MK, FOR_MK_from_R],['geneID'],'inner')

FOR_MK["Dn_sin_reali"]=FOR_MK["Dn_sin_reali"]*FOR_MK.N_sites
FOR_MK["Dn_sin_juv"]=FOR_MK["Dn_sin_juv"]*FOR_MK.N_sites
FOR_MK["Dn_juv_reali"]=FOR_MK["Dn_juv_reali"]*FOR_MK.N_sites
FOR_MK["Ds_sin_reali"]=FOR_MK["Ds_sin_reali"]*FOR_MK.S_sites
FOR_MK["Ds_sin_juv"]=FOR_MK["Ds_sin_juv"]*FOR_MK.S_sites
FOR_MK["Ds_juv_reali"]=FOR_MK["Ds_juv_reali"]*FOR_MK.S_sites

FOR_MK["Dn_sin_reali"]=FOR_MK["Dn_sin_reali"].round(0)
FOR_MK["Dn_sin_juv"]=FOR_MK["Dn_sin_juv"].round(0)
FOR_MK["Dn_juv_reali"]=FOR_MK["Dn_juv_reali"].round(0)
FOR_MK["Ds_sin_reali"]=FOR_MK["Ds_sin_reali"].round(0)
FOR_MK["Ds_sin_juv"]=FOR_MK["Ds_sin_juv"].round(0)
FOR_MK["Ds_juv_reali"]=FOR_MK["Ds_juv_reali"].round(0)

FOR_MK=join_raw_data_base([FOR_MK,DM_genes],['geneID'],'outer')

FOR_MK=FOR_MK.drop_duplicates(subset=['geneID'])

adaptive_evolu_L_sin = dnds_table_gene[dnds_table_gene.W_L_sin > dnds_table_gene.Pin_Pis_sinapis]
adaptive_evolu_L_rea = dnds_table_gene[dnds_table_gene.W_L_rea > dnds_table_gene.Pin_Pis_spanish_reali]
adaptive_evolu_L_juv = dnds_table_gene[dnds_table_gene.W_L_juv > dnds_table_gene.Pin_Pis_juvernica]

def McDonald_Kreitman_fisher_exact(FOR_MK, list_1, DOS_list,relative_ID):
	import scipy.stats as stats
	#np.array([[ps, ds], [pn, dn]])

	significat={}
	non_significat={}
	for i_index, i in FOR_MK.iterrows():
		if i.Pi_syn_sinapis > 0 and i.Ds_sin_juv > 0 and i.Pi_nonsyn_sinapis >0 and i.Dn_sin_juv >0 and i[DOS_list] > 0:
			relative=float(FOR_DOS[FOR_DOS.geneID==i.geneID][relative_ID])
			obs=np.array([[i[list_1[0]],i[list_1[1]]], [i[list_1[2]], i[list_1[3]]]])
			fisher, p_val=stats.fisher_exact(obs)

			if p_val < 0.05:
				significat[i.geneID]=[p_val, i.DM_gene, obs, i[DOS_list],relative]
			else:
				non_significat[i.geneID]=[p_val, i.DM_gene ,obs, i[DOS_list],relative]
	return significat

sinapis_juv=McDonald_Kreitman_fisher_exact(FOR_MK, ['Pi_syn_sinapis','Ds_sin_juv','Pi_nonsyn_sinapis','Dn_sin_juv'], 'DOS_Sin', 'Relative_sinapis_site_pi' )
reali_juv=McDonald_Kreitman_fisher_exact(FOR_MK, ['Pi_syn_spanish_reali','Ds_juv_reali','Pi_nonsyn_spanish_reali','Dn_juv_reali'], 'DOS_Rea', 'Relative_spanish_reali_site_pi')

sinapis_reali=McDonald_Kreitman_fisher_exact(FOR_MK, ['Pi_syn_sinapis','Ds_sin_reali','Pi_nonsyn_sinapis','Dn_sin_reali'], 'DOS_Sin', 'Relative_sinapis_site_pi')
reali_sinapis=McDonald_Kreitman_fisher_exact(FOR_MK, ['Pi_syn_spanish_reali','Ds_sin_reali','Pi_nonsyn_spanish_reali','Dn_sin_reali'], 'DOS_Rea', 'Relative_spanish_reali_site_pi')

juvernica_sinapis=McDonald_Kreitman_fisher_exact(FOR_MK, ['Pi_syn_juvernica','Ds_sin_juv','Pi_nonsyn_juvernica','Dn_sin_juv'], 'DOS_Juv', 'Relative_juvernica_site_pi')
juvernica_reali=McDonald_Kreitman_fisher_exact(FOR_MK, ['Pi_syn_juvernica','Ds_juv_reali','Pi_nonsyn_juvernica','Dn_juv_reali'], 'DOS_Juv', 'Relative_juvernica_site_pi')



FOR_DOS_L_sin=FOR_DOS[(FOR_DOS.s_subs_L_sin >= 3) & (FOR_DOS.DOS_Sin > 0.0)]
FOR_DOS_L_rea=FOR_DOS[(FOR_DOS.s_subs_L_rea >= 3) & (FOR_DOS.DOS_Rea > 0.0)]
FOR_DOS_L_juv=FOR_DOS[(FOR_DOS.s_subs_L_juv >= 3) & (FOR_DOS.DOS_Juv > 0.0)]

FOR_DOS_L_sin_DM_genes=join_raw_data_base([FOR_DOS_L_sin,DM_genes],['geneID'],'inner')
FOR_DOS_L_rea_DM_genes=join_raw_data_base([FOR_DOS_L_rea,DM_genes],['geneID'],'inner')
FOR_DOS_L_juv_DM_genes=join_raw_data_base([FOR_DOS_L_juv,DM_genes],['geneID'],'inner')


###########################################################


fianl_df_1_a=fianl_df_1[fianl_df_1['auto_allo']!='Z']
fianl_df_2_a=fianl_df_2[fianl_df_2['auto_allo']!='Z']
fianl_df_3_a=fianl_df_3[fianl_df_3['auto_allo']!='Z']
fianl_df_4_a=fianl_df_4[fianl_df_4['auto_allo']!='Z']
fianl_df_5_a=fianl_df_5[fianl_df_5['auto_allo']!='Z']
fianl_df_6_a=fianl_df_6[fianl_df_6['auto_allo']!='Z']
fianl_df_7_a=fianl_df_7[fianl_df_7['auto_allo']!='Z']
fianl_df_8_a=fianl_df_8[fianl_df_8['auto_allo']!='Z']

fianl_df_1_Z=fianl_df_1[fianl_df_1['auto_allo']=='Z']
fianl_df_2_Z=fianl_df_2[fianl_df_2['auto_allo']=='Z']
fianl_df_3_Z=fianl_df_3[fianl_df_3['auto_allo']=='Z']
fianl_df_4_Z=fianl_df_4[fianl_df_4['auto_allo']=='Z']
fianl_df_5_Z=fianl_df_5[fianl_df_5['auto_allo']=='Z']
fianl_df_6_Z=fianl_df_6[fianl_df_6['auto_allo']=='Z']
fianl_df_7_Z=fianl_df_7[fianl_df_7['auto_allo']=='Z']
fianl_df_8_Z=fianl_df_8[fianl_df_8['auto_allo']=='Z']

final_df=join_raw_data_base([fianl_df_1,fianl_df_2,fianl_df_3,fianl_df_4,fianl_df_5,fianl_df_6,fianl_df_7,fianl_df_8],['geneID','length_orf'],'outer')



sns.boxplot(final_df[["Pin_juvernica","Pis_juvernica","Pin_spanish_reali","Pis_spanish_reali","Pin_sinapis","Pis_sinapis"]],palette="PRGn",showfliers=False)


sns.boxplot(final_df[["Pin_juvernica","Pis_juvernica","Pin_spanish_reali","Pis_spanish_reali","Pin_sinapis","Pis_sinapis"]],palette="PRGn",showfliers=False)




sns.boxplot(final_df[["Pin_Pis_juvernica","Pin_Pis_spanish_reali","Pin_Pis_sinapis"]],palette="PRGn",showfliers=False).set(ylim=(0,3.5))

final_df[["Pin_Pis_juvernica","Pin_Pis_spanish_reali","Pin_Pis_sinapis"]].replace([np.inf, -np.inf], np.nan).mean()


fig, ax = plt.subplots(ncols=2, sharey=True)
box = pop_df_1[pop_df_1.auto_allo_x=='auto'][['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_introns','irish_juvernica_site_pi_global','irish_juvernica_site_pi_codon_4d']].boxplot(ax=ax[0], sym='')
box2 = fianl_df_1_a[['Pin_irish_juvernica','Pis_irish_juvernica']].boxplot(ax=ax[1], sym='')
ax[0].grid(False)
ax[1].grid(False)



fig, ax = plt.subplots(ncols=2, sharey=True)
box = pop_df_1[pop_df_1.auto_allo_x=='Z'][['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_introns','irish_juvernica_site_pi_global','irish_juvernica_site_pi_codon_4d']].boxplot(ax=ax[0], sym='')
box2 = fianl_df_1_Z[['Pin_irish_juvernica','Pis_irish_juvernica']].boxplot(ax=ax[1], sym='')

fig, ax = plt.subplots(ncols=2, sharey=True)
box = pop_df_1[['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_introns','irish_juvernica_site_pi_global','irish_juvernica_site_pi_codon_4d']].boxplot(ax=ax[0], sym='')
box2 = fianl_df_1[['Pin_irish_juvernica','Pis_irish_juvernica']].boxplot(ax=ax[1], sym='')



#############################callculate DN DS for all pairs#######################


print FOR_DOS.s_subs_L_rea.sum()
print FOR_DOS.s_subs_L_sin.sum()
print FOR_DOS.s_subs_L_juv.sum()

print FOR_DOS.n_subs_L_rea.sum()
print FOR_DOS.n_subs_L_sin.sum()
print FOR_DOS.n_subs_L_juv.sum()

###############################################GC contenets of genes

from Bio.SeqUtils import GC

fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034/NBIS_annotation_leptidea/fasta/cds.fa')



dnds_table_gene['GC_gene']=0


for window_index, window in dnds_table_gene.iterrows():
    gc_fraction=GC(str(fasta[window.geneID]))
    #print str(gc_fraction)+ '   '+ window.CHROM + '     '+ str(window.BIN_START) + '        '+ str(window.BIN_END)
    dnds_table_gene['GC_gene'].iloc[window_index]=gc_fraction
    print window_index

dnds_table_gene['GC_gene']=dnds_table_gene['GC_gene'].astype('float')


sns.jointplot("Pin_Pis_juvernica", "GC_gene", data=dnds_table_gene, kind="reg")
sns.jointplot("Pin_Pis_spanish_reali", "GC_gene", data=dnds_table_gene, kind="reg")
sns.jointplot("Pin_Pis_sinapis", "GC_gene", data=dnds_table_gene, kind="reg")


sns.jointplot("Pin_Pis_juvernica", "GC_gene", data=dnds_table_gene, kind="reg")
sns.jointplot("Pin_Pis_spanish_reali", "GC_gene", data=dnds_table_gene, kind="reg")
sns.jointplot("Pin_Pis_sinapis", "GC_gene", data=dnds_table_gene, kind="reg")



dnds_table_gene[['GC_gene']].plot.hist(alpha=0.5, bins=20)


################################################GC contenets of windows


from Bio.SeqUtils import GC

fasta=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/GENOME_ASSEMBLY/assembly_updates/v1.4/N.Backstrom_leptidea.scf.1.4.fasta')




dnds_table_windows_2['GC_reference']=0.0

for window_index, window in dnds_table_windows_2.iterrows():
    gc_fraction=GC(str(fasta[window.CHROM][int(window.BIN_START-1):int(window.BIN_END-1)]))
    #print str(gc_fraction)+ '   '+ window.CHROM + '     '+ str(window.BIN_START) + '        '+ str(window.BIN_END)
    dnds_table_windows_2['GC_reference'].iloc[window_index]=gc_fraction
    print window_index


dnds_table_windows_2['GC_reference']=dnds_table_windows_2['GC_reference'].astype('float')


a=list(dnds_table_windows_2['GC_reference'])
b=list(dnds_table_gene['GC_gene'])

fig, axes  = plt.subplots(nrows=3, sharex=True, gridspec_kw={"height_ratios": (.10,.10 ,.80)})

sns.boxplot(a, ax=axes[0], color='#A1BAD9', showfliers=False)
sns.boxplot(b, ax=axes[1], color='#AAD3B3', showfliers=False)
sns.distplot(a, ax=axes[2], color='#A1BAD9')
sns.distplot(b, ax=axes[2], color='#AAD3B3')



sns.jointplot("GC_reference", "gene_density", data=dnds_table_windows_2, kind="reg", color='#24cc24')

sns.jointplot("GC_reference", "juvernica_site_pi_global", data=dnds_table_windows_2, kind="reg", color='#24cc24')
sns.jointplot("GC_reference", "sinapis_site_pi_global", data=dnds_table_windows_2, kind="reg", color='#24cc24')
sns.jointplot("GC_reference", "spanish_reali_site_pi_global", data=dnds_table_windows_2, kind="reg", color='#24cc24')




pi_ps=list(dnds_table_windows_2.Pin_Pis_sinapis)+list(dnds_table_windows_2.Pin_Pis_spanish_reali)+list(dnds_table_windows_2.Pin_Pis_juvernica)
W=list(dnds_table_windows_2.W_L_sin)+list(dnds_table_windows_2.W_L_rea)+list(dnds_table_windows_2.W_L_juv)
W_alpha=list(dnds_table_windows_2.W_alpha_Sin)+list(dnds_table_windows_2.W_alpha_Rea)+list(dnds_table_windows_2.W_alpha_Juv)
N_sites=list(dnds_table_windows_2.N_sites)+list(dnds_table_windows_2.N_sites)+list(dnds_table_windows_2.N_sites)
S_sites=list(dnds_table_windows_2.S_sites)+list(dnds_table_windows_2.S_sites)+list(dnds_table_windows_2.S_sites)
coding_sites=list(dnds_table_windows_2.S_sites)+list(dnds_table_windows_2.S_sites)+list(dnds_table_windows_2.S_sites)
DOS=list(dnds_table_windows_2.DOS_Sin)+list(dnds_table_windows_2.DOS_Rea)+list(dnds_table_windows_2.DOS_Juv)
GC_reference=list(dnds_table_windows_2.GC_reference)+list(dnds_table_windows_2.GC_reference)+list(dnds_table_windows_2.GC_reference)


relative=list(dnds_table_windows_2.Relative_sinapis_site_pi_direction)+list(dnds_table_windows_2.Relative_spanish_reali_site_pi_direction)+list(dnds_table_windows_2.Relative_juvernica_site_pi_direction)
species=['L_sin']*len(list(dnds_table_windows_2.Pin_Pis_sinapis))+['L_Rea']*len(list(dnds_table_windows_2.Pin_Pis_spanish_reali))+['L_Juv']*len(list(dnds_table_windows_2.Pin_Pis_juvernica))
violin_df = pandas.DataFrame()
violin_df['pn_ps']=pi_ps
violin_df['W']=W
violin_df['W_alpha']=W_alpha
violin_df['DOS']=DOS
violin_df['GC_reference']=GC_reference

violin_df['N_sites']=N_sites
violin_df['S_sites']=S_sites
violin_df['relative']=relative
violin_df['species']=species
violin_df=violin_df[violin_df.relative!='']
violin_df['sites']=violin_df['S_sites']+violin_df['N_sites']


sns.boxplot(x="species", y="pn_ps", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high']).set(ylim=(0,10))
sns.boxplot(x="species", y="W", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'] ).set(ylim=(0,1))
sns.boxplot(x="species", y="W_alpha", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high']).set(ylim=(-3,3))
sns.boxplot(x="species", y="DOS", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high']).set(ylim=(-3,3))

sns.boxplot(x="species", y="N_sites", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
sns.boxplot(x="species", y="S_sites", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])

new=sns.boxplot(x="species", y="GC_reference", hue="relative",data=violin_df, palette=pi_colours_species, showfliers=False, hue_order=['low','high'])
new.legend_.remove()


print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['GC_reference'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['GC_reference'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['GC_reference'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['GC_reference'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['GC_reference'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['GC_reference'].dropna(how='all')))


print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['Pin_Pis_sinapis'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['Pin_Pis_sinapis'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['Pin_Pis_spanish_reali'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['Pin_Pis_spanish_reali'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['Pin_Pis_juvernica'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['Pin_Pis_juvernica'].dropna(how='all')))

print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi > 0]['W_L_sin'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_sinapis_site_pi < 0]['W_L_sin'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi > 0]['W_L_rea'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_spanish_reali_site_pi < 0]['W_L_rea'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi > 0]['W_L_juv'].dropna(how='all')), list(dnds_table_windows_2[dnds_table_windows_2.Relative_juvernica_site_pi < 0]['W_L_juv'].dropna(how='all')))


################################codon GC content#########################


fianl_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_0.csv',index_col=False)
fianl_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_1.csv',index_col=False)
fianl_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_2.csv',index_col=False)
fianl_4=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_3.csv',index_col=False)
fianl_5=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_4.csv',index_col=False)
fianl_6=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/test/dnds_table_windows_5.csv',index_col=False)

dnds_table_windows_2_2=pandas.concat([fianl_1,fianl_2,fianl_3,fianl_4, fianl_5, fianl_6],ignore_index=True)

dnds_table_windows_2_2[["GC_codon_1","GC_codon_2","GC_codon_3","GC_codon_4D","GC_intron","GC_intergenic"]]

dnds_table_windows_2_2['CDS_sites']=dnds_table_windows_2_2['N_sites']+dnds_table_windows_2_2['S_sites']

fig, axes  = plt.subplots(nrows=7, sharex=True, gridspec_kw={"height_ratios": (.04,.04 ,.04,.04, .04, .04 ,.76)})

sns.boxplot(dnds_table_windows_2_2["GC_codon_1"], ax=axes[0], color='#5975A4', showfliers=False)
sns.boxplot(dnds_table_windows_2_2["GC_codon_2"], ax=axes[1], color='#5F9E6F', showfliers=False)
sns.boxplot(dnds_table_windows_2_2["GC_codon_3"], ax=axes[2], color='#B65D61', showfliers=False)
sns.boxplot(dnds_table_windows_2_2["GC_codon_4D"], ax=axes[3], color='#72AEC0', showfliers=False)
sns.boxplot(dnds_table_windows_2_2["GC_intron"], ax=axes[4], color='#867AAA', showfliers=False)
sns.boxplot(dnds_table_windows_2_2["GC_intergenic"], ax=axes[5], color='#C0B37E', showfliers=False)

sns.distplot(dnds_table_windows_2_2["GC_codon_1"].dropna(), ax=axes[6], color='#5975A4')
sns.distplot(dnds_table_windows_2_2["GC_codon_2"].dropna(), ax=axes[6], color='#5F9E6F')
sns.distplot(dnds_table_windows_2_2["GC_codon_3"].dropna(), ax=axes[6], color='#B65D61')
sns.distplot(dnds_table_windows_2_2["GC_codon_4D"].dropna(), ax=axes[6], color='#72AEC0')
sns.distplot(dnds_table_windows_2_2["GC_intron"].dropna(), ax=axes[6], color='#867AAA')
sns.distplot(dnds_table_windows_2_2["GC_intergenic"].dropna(), ax=axes[6], color='#C0B37E')

axes[0].set_ylabel('')    
axes[0].set_xlabel('')
axes[1].set_ylabel('')    
axes[1].set_xlabel('')
axes[2].set_xlabel('')
axes[2].set_ylabel('')
axes[3].set_xlabel('')
axes[3].set_ylabel('')
axes[4].set_xlabel('')
axes[4].set_ylabel('')
axes[5].set_xlabel('')
axes[5].set_ylabel('')
axes[6].set_xlabel('')
axes[6].set_ylabel('')








print dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_sinapis_site_pi < 0]['GC_reference'].std()
print dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_sinapis_site_pi > 0]['GC_reference'].std()
print dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_spanish_reali_site_pi < 0]['GC_reference'].std()
print dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_spanish_reali_site_pi > 0]['GC_reference'].std()
print dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_juvernica_site_pi < 0]['GC_reference'].std()
print dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_juvernica_site_pi > 0]['GC_reference'].std()



print scipy.stats.mannwhitneyu(list(dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_sinapis_site_pi > 0]['GC_reference'].dropna(how='all')), list(dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_sinapis_site_pi < 0]['GC_reference'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_spanish_reali_site_pi > 0]['GC_reference'].dropna(how='all')), list(dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_spanish_reali_site_pi < 0]['GC_reference'].dropna(how='all')))
print scipy.stats.mannwhitneyu(list(dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_juvernica_site_pi > 0]['GC_reference'].dropna(how='all')), list(dnds_table_windows_2_2[dnds_table_windows_2_2.Relative_juvernica_site_pi < 0]['GC_reference'].dropna(how='all')))


     
      

sns.regplot(x="GC_reference", y="CDS_sites", data=dnds_table_windows_2_2, scatter_kws={"s": 80}, order=2, ci=None, truncate=True)

sns.lmplot(x="GC_reference", y="CDS_sites", data=dnds_table_windows_2_2 )

sns.jointplot("GC_reference", "CDS_sites", data=dnds_table_windows_2_2, kind="reg", color='#24cc24')


##################################################rho###recombination read_table


irish_juvernica_rho=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/irish_juvernica.rho")
kazak_juvernica_rho=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/kazak_juvernica.rho")
kazak_sinapis_rho=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/kazak_sinapis.rho")
spanish_reali_rho=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/spanish_reali.rho")
spanish_sinapis_rho=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/spanish_sinapis.rho")
swedish_sinapis_rho=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private/swedish_sinapis.rho")

irish_juvernica_rho.BIN_START=irish_juvernica_rho.BIN_START+1
kazak_juvernica_rho.BIN_START=kazak_juvernica_rho.BIN_START+1
kazak_sinapis_rho.BIN_START=kazak_sinapis_rho.BIN_START+1
spanish_reali_rho.BIN_START=spanish_reali_rho.BIN_START+1
spanish_sinapis_rho.BIN_START=spanish_sinapis_rho.BIN_START+1
swedish_sinapis_rho.BIN_START=swedish_sinapis_rho.BIN_START+1



#irish_juvernica_rho.irish_juvernica_rho=irish_juvernica_rho.irish_juvernica_rho/10000
#kazak_juvernica_rho.kazak_juvernica_rho=kazak_juvernica_rho.kazak_juvernica_rho/10000
#kazak_sinapis_rho.kazak_sinapis_rho=kazak_sinapis_rho.kazak_sinapis_rho/10000
#spanish_reali_rho.spanish_reali_rho=spanish_reali_rho.spanish_reali_rho/10000
#spanish_sinapis_rho.spanish_sinapis_rho=spanish_sinapis_rho.spanish_sinapis_rho/10000
#swedish_sinapis_rho.swedish_sinapis_rho=swedish_sinapis_rho.swedish_sinapis_rho/10000


rho=join_raw_data_base([irish_juvernica_rho,kazak_juvernica_rho, kazak_sinapis_rho, spanish_reali_rho, spanish_sinapis_rho, swedish_sinapis_rho],['CHROM', 'BIN_START', 'BIN_END'],'inner')
rho['POS']=rho.BIN_START
rho['BIN_START']=0
rho['BIN_END']=0
rho['BIN_START']=(np.floor(rho['POS']/100000)*100000)+1
rho['BIN_END']=rho['BIN_START']+(100000-1)

rho=rho.groupby(['CHROM','BIN_START','BIN_END'],as_index=False).agg({'irish_juvernica_rho': 'mean','kazak_juvernica_rho': 'mean','kazak_sinapis_rho': 'mean','spanish_reali_rho': 'mean','spanish_sinapis_rho': 'mean','swedish_sinapis_rho': 'mean'})

rho=join_raw_data_base([dnds_table_windows_2_2,rho],['CHROM', 'BIN_START', 'BIN_END'],'inner')

temp_df_1=pop_df_1[['CHROM','BIN_START','BIN_END', 'irish_juvernica_site_pi_global', 'irish_juvernica_site_pi_codon_4d']]
temp_df_3=pop_df_3[['CHROM','BIN_START','BIN_END', 'kazak_juvernica_site_pi_global','kazak_juvernica_site_pi_codon_4d']]
temp_df_4=pop_df_4[['CHROM','BIN_START','BIN_END', 'kaz_sin_site_pi_global','kaz_sin_site_pi_codon_4d']]
temp_df_6=pop_df_6[['CHROM','BIN_START','BIN_END', 'spanish_reali_site_pi_global','spanish_reali_site_pi_codon_4d']]
temp_df_7=pop_df_7[['CHROM','BIN_START','BIN_END', 'spanish_sinapis_site_pi_global','spanish_sinapis_site_pi_codon_4d']]
temp_df_8=pop_df_8[['CHROM','BIN_START','BIN_END', 'swe_sin_allele_site_pi_global','swe_sin_allele_site_pi_codon_4d']]


temp=join_raw_data_base([temp_df_1,temp_df_3,temp_df_4,temp_df_6,temp_df_7,temp_df_8], ['CHROM', 'BIN_START', 'BIN_END'],'inner')


rho=join_raw_data_base([rho,temp],['CHROM', 'BIN_START', 'BIN_END'],'inner')
rho_2=join_raw_data_base([rho,dnds_table_windows_2],['CHROM', 'BIN_START', 'BIN_END'],'inner')




rho_list=list(rho.irish_juvernica_rho)+list(rho.kazak_juvernica_rho)+list(rho.kazak_sinapis_rho)+list(rho.spanish_reali_rho)+list(rho.spanish_sinapis_rho)+list(rho.swedish_sinapis_rho)
pi=list(rho.irish_juvernica_site_pi_global)+list(rho.kazak_juvernica_site_pi_global)+list(rho.kaz_sin_site_pi_global)+list(rho.spanish_reali_site_pi_global)+list(rho.spanish_sinapis_site_pi_global)+list(rho.swe_sin_allele_site_pi_global)
populations=['LjIre']*len(list(rho.irish_juvernica_site_pi_global))+['LjKaz']*len(list(rho.kazak_juvernica_site_pi_global))+['LsKaz']*len(list(rho.kaz_sin_site_pi_global))+['LrSpa']*len(list(rho.spanish_reali_site_pi_global))+['LsSpa']*len(list(rho.spanish_sinapis_site_pi_global))+['LsSwe']*len(list(rho.swe_sin_allele_site_pi_global))

violin= pandas.DataFrame()
violin['rho']=rho_list
violin['pi']=pi
violin['populations']=populations



sns.boxplot(x="populations", y="rho",data=violin, palette=pi_colours, showfliers=False, order=['LjIre','LjKaz', 'LrSpa', 'LsKaz', 'LsSpa', 'LsSwe'] )
plt.savefig('/home/venkat/box_rho.png')
plt.close()


sns.jointplot("irish_juvernica_rho", "irish_juvernica_site_pi_global", data=rho, kind="reg", color='#C8C800')
plt.savefig('/home/venkat/LjIre_rho_pi.png')
plt.close()
sns.jointplot("kazak_juvernica_rho", "kazak_juvernica_site_pi_global", data=rho, kind="reg", color='#006400')
plt.savefig('/home/venkat/LjKaz_rho_pi.png')
plt.close()
sns.jointplot("kazak_sinapis_rho", "kaz_sin_site_pi_global", data=rho, kind="reg", color='#FF8C00')
plt.savefig('/home/venkat/LsKaz_rho_pi.png')
plt.close()
sns.jointplot("spanish_reali_rho", "spanish_reali_site_pi_global", data=rho, kind="reg", color='#0000FF')
plt.savefig('/home/venkat/LrSpa_rho_pi.png')
plt.close()
sns.jointplot("spanish_sinapis_rho", "spanish_sinapis_site_pi_global", data=rho, kind="reg", color='#FF0000')
plt.savefig('/home/venkat/LsSpa_rho_pi.png')
plt.close()
sns.jointplot("swedish_sinapis_rho", "swe_sin_allele_site_pi_global", data=rho, kind="reg", color='#C04000')
plt.savefig('/home/venkat/LsSwe_rho_pi.png')
plt.close()



sns.jointplot("GC_reference", "irish_juvernica_site_pi_codon_4d", data=rho, kind="reg", color='#C8C800')
plt.savefig('/home/venkat/LjIre_rho_pi.png')
plt.close()
sns.jointplot("GC_reference", "kazak_juvernica_site_pi_codon_4d", data=rho, kind="reg", color='#006400')
plt.savefig('/home/venkat/LjKaz_rho_pi.png')
plt.close()
sns.jointplot("GC_reference", "kaz_sin_site_pi_codon_4d", data=rho, kind="reg", color='#FF8C00')
plt.savefig('/home/venkat/LsKaz_rho_pi.png')
plt.close()
sns.jointplot("GC_reference", "spanish_reali_site_pi_codon_4d", data=rho, kind="reg", color='#0000FF')
plt.savefig('/home/venkat/LrSpa_rho_pi.png')
plt.close()
sns.jointplot("GC_reference", "spanish_sinapis_site_pi_codon_4d", data=rho, kind="reg", color='#FF0000')
plt.savefig('/home/venkat/LsSpa_rho_pi.png')
plt.close()
sns.jointplot("GC_reference", "swe_sin_allele_site_pi_codon_4d", data=rho, kind="reg", color='#C04000')
plt.savefig('/home/venkat/LsSwe_rho_pi.png')
plt.close()

 


def corr_plot(df):
        axes= scatter_matrix(df,diagonal='none')
        corr = df.corr(method='pearson').as_matrix()
        mask = np.zeros_like(corr, dtype=np.bool)
        mask[np.tril_indices_from(mask)] = True

        for i, j in zip(*plt.np.triu_indices_from(axes, k=1)):
            axes[i, j].set_visible(False)
        
        for i, j in zip(*plt.np.diag_indices_from(axes)):
            axes[i, j].set_visible(False)

        for i, j in zip(*plt.np.tril_indices_from(axes, k=1)):
            axes[i, j].set()
            axes[i, j].set()
        plt.show()


def corr_pvalue(df):
        from scipy.stats import pearsonr
        import numpy as np
        import pandas as pd
        numeric_df = df.dropna()._get_numeric_data()
        cols = numeric_df.columns
        mat = numeric_df.values
        arr = np.zeros((len(cols),len(cols)), dtype=object)
        for xi, x in enumerate(mat.T):
                for yi, y in enumerate(mat.T[xi:]):
                    arr[xi, yi+xi] = map(lambda _: round(_,3), pearsonr(x,y))
                    arr[yi+xi, xi] = arr[xi, yi+xi]
        
        return pandas.DataFrame(arr, index=cols, columns=cols)



rho[['spanish_sinapis_rho','swedish_sinapis_rho','kazak_juvernica_rho','kazak_sinapis_rho','irish_juvernica_rho','spanish_reali_rho']]

valuess=rho[['swedish_sinapis_rho','kazak_sinapis_rho','spanish_sinapis_rho','spanish_reali_rho','kazak_juvernica_rho','irish_juvernica_rho']]





order=['LjIre','LjKaz', 'LrSpa', 'LsSpa', 'LsKaz', 'LsSwe']

pi_colours=sns.color_palette(['#C8C800','#006400','#0000FF','#FF0000','#FF8C00','#C04000'])




fig, axes  = plt.subplots(nrows=2, sharex=True, gridspec_kw={"height_ratios": (.30,.70)})

sns.boxplot(y="populations", x="rho", orient='h', ax=axes[0], data=violin,palette=pi_colours[::-1], showfliers=False,order=['LjIre','LjKaz', 'LrSpa', 'LsSpa', 'LsKaz', 'LsSwe'][::-1] )

sns.distplot(violin[violin.populations == 'LjIre'].rho, ax=axes[1], hist=False,color='#C8C800')
sns.distplot(violin[violin.populations == 'LjKaz'].rho, ax=axes[1], hist=False,color='#006400')
sns.distplot(violin[violin.populations == 'LrSpa'].rho, ax=axes[1], hist=False,color='#0000FF')
sns.distplot(violin[violin.populations == 'LsKaz'].rho, ax=axes[1], hist=False,color='#FF8C00')
sns.distplot(violin[violin.populations == 'LsSpa'].rho, ax=axes[1], hist=False,color='#FF0000')
sns.distplot(violin[violin.populations == 'LsSwe'].rho, ax=axes[1], hist=False,color='#C04000')



axes[0].set_ylabel('')    
axes[0].set_xlabel('')
axes[1].set_ylabel('')    
axes[1].set_xlabel('')

corr_plot(valuess[::-1])

corr_pvalue(valuess[::-1])



