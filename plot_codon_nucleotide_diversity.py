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
    fianl_df_1=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf1/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_2=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf2/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_3=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf3/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_4=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf4/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_5=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf5/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_6=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf6/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_7=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf7/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_8=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf8/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_9=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf9/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    fianl_df_10=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/scaf10/'+Pop_1+'_2_final_df_with_inter_geneic.csv',index_col=False)
    final_df_final=pandas.concat([fianl_df_1,fianl_df_2,fianl_df_3,fianl_df_4,fianl_df_5, fianl_df_6, fianl_df_7, fianl_df_8, fianl_df_9, fianl_df_10],ignore_index=True)
    final_df_final=final_df_final.rename(columns={'site_pi_codon_1': Pop_1+'_site_pi_codon_1', 'site_pi_codon_2': Pop_1+'_site_pi_codon_2', 'site_pi_codon_3': Pop_1+'_site_pi_codon_3', 'site_pi_codon_4d': Pop_1+'_site_pi_codon_4d', 'site_pi_introns': Pop_1+'_site_pi_introns', 'site_pi_global': Pop_1+'_site_pi_global','sites_codon_1': Pop_1+'_sites_codon_1','sites_codon_2': Pop_1+'_sites_codon_2','sites_codon_3': Pop_1+'_sites_codon_3','sites_codon_4d': Pop_1+'_sites_codon_4d','sites_introns': Pop_1+'_sites_introns','sites_global': Pop_1+'_sites_global'})
    #final_df_final=final_df_final.rename(columns={'site_pi_codon_1_corrected': Pop_1+'_site_pi_codon_1_corrected', 'site_pi_codon_2_corrected': Pop_1+'_site_pi_codon_2_corrected', 'site_pi_codon_3_corrected': Pop_1+'_site_pi_codon_3_corrected', 'site_pi_codon_4d_corrected': Pop_1+'_site_pi_codon_4d_corrected', 'site_pi_introns_corrected': Pop_1+'_site_pi_introns_corrected', 'sites_codon_1_corrected': Pop_1+'_sites_codon_1_corrected','sites_codon_2_corrected': Pop_1+'_sites_codon_2_corrected','sites_codon_3_corrected': Pop_1+'_sites_codon_3_corrected','sites_codon_4d_corrected': Pop_1+'_sites_codon_4d_corrected','sites_introns_corrected': Pop_1+'_sites_introns_corrected'})
    final_df_final=final_df_final[final_df_final[Pop_1+'_sites_global']>20000]
    return final_df_final



pop_df_1=load_files('irish_juvernica')
pop_df_3=load_files('kazak_juvernica')
pop_df_4=load_files('kaz_sin')
pop_df_6=load_files('spanish_reali')
pop_df_7=load_files('spanish_sinapis')
pop_df_8=load_files('swe_sin_allele')



#Nucleotide_diversity_LjIre=list(pop_df_1.site_pi_global)+list(pop_df_1.site_pi_intergenic)+list(pop_df_1.site_pi_introns)+list(pop_df_1.site_pi_codon_1)+list(pop_df_1.site_pi_codon_2)+list(pop_df_1.site_pi_codon_3)+list(pop_df_1.site_pi_codon_4d)
#Nucleotide_diversity_LjKaz=list(pop_df_3.site_pi_global)+list(pop_df_3.site_pi_intergenic)+list(pop_df_3.site_pi_introns)+list(pop_df_3.site_pi_codon_1)+list(pop_df_3.site_pi_codon_2)+list(pop_df_3.site_pi_codon_3)+list(pop_df_3.site_pi_codon_4d)
#Nucleotide_diversity_LsKaz=list(pop_df_4.site_pi_global)+list(pop_df_4.site_pi_intergenic)+list(pop_df_4.site_pi_introns)+list(pop_df_4.site_pi_codon_1)+list(pop_df_4.site_pi_codon_2)+list(pop_df_4.site_pi_codon_3)+list(pop_df_4.site_pi_codon_4d)
#Nucleotide_diversity_LrSpa=list(pop_df_6.site_pi_global)+list(pop_df_6.site_pi_intergenic)+list(pop_df_6.site_pi_introns)+list(pop_df_6.site_pi_codon_1)+list(pop_df_6.site_pi_codon_2)+list(pop_df_6.site_pi_codon_3)+list(pop_df_6.site_pi_codon_4d)
#Nucleotide_diversity_LsSpa=list(pop_df_7.site_pi_global)+list(pop_df_7.site_pi_intergenic)+list(pop_df_7.site_pi_introns)+list(pop_df_7.site_pi_codon_1)+list(pop_df_7.site_pi_codon_2)+list(pop_df_7.site_pi_codon_3)+list(pop_df_7.site_pi_codon_4d)
#Nucleotide_diversity_Lsswe=list(pop_df_8.site_pi_global)+list(pop_df_8.site_pi_intergenic)+list(pop_df_8.site_pi_introns)+list(pop_df_8.site_pi_codon_1)+list(pop_df_8.site_pi_codon_2)+list(pop_df_8.site_pi_codon_3)+list(pop_df_8.site_pi_codon_4d)
#
#Population_LjIre=['LjIre']*len(list(pop_df_1.site_pi_global))+['LjIre']*len(list(pop_df_1.site_pi_intergenic))+['LjIre']*len(list(pop_df_1.site_pi_introns))+['LjIre']*len(list(pop_df_1.site_pi_codon_1))+['LjIre']*len(list(pop_df_1.site_pi_codon_2))+['LjIre']*len(list(pop_df_1.site_pi_codon_3))+['LjIre']*len(list(pop_df_1.site_pi_codon_4d))
#Population_LjKaz=["LjKaz"]*len(list(pop_df_3.site_pi_global))+["LjKaz"]*len(list(pop_df_3.site_pi_intergenic))+["LjKaz"]*len(list(pop_df_3.site_pi_introns))+["LjKaz"]*len(list(pop_df_3.site_pi_codon_1))+["LjKaz"]*len(list(pop_df_3.site_pi_codon_2))+["LjKaz"]*len(list(pop_df_3.site_pi_codon_3))+["LjKaz"]*len(list(pop_df_3.site_pi_codon_4d))
#Population_LsKaz=["LsKaz"]*len(list(pop_df_4.site_pi_global))+["LsKaz"]*len(list(pop_df_4.site_pi_intergenic))+["LsKaz"]*len(list(pop_df_4.site_pi_introns))+["LsKaz"]*len(list(pop_df_4.site_pi_codon_1))+["LsKaz"]*len(list(pop_df_4.site_pi_codon_2))+["LsKaz"]*len(list(pop_df_4.site_pi_codon_3))+["LsKaz"]*len(list(pop_df_4.site_pi_codon_4d))
#Population_LrSpa=["LrSpa"]*len(list(pop_df_6.site_pi_global))+["LrSpa"]*len(list(pop_df_6.site_pi_intergenic))+["LrSpa"]*len(list(pop_df_6.site_pi_introns))+["LrSpa"]*len(list(pop_df_6.site_pi_codon_1))+["LrSpa"]*len(list(pop_df_6.site_pi_codon_2))+["LrSpa"]*len(list(pop_df_6.site_pi_codon_3))+["LrSpa"]*len(list(pop_df_6.site_pi_codon_4d))
#Population_LsSpa=["LsSpa"]*len(list(pop_df_7.site_pi_global))+["LsSpa"]*len(list(pop_df_7.site_pi_intergenic))+["LsSpa"]*len(list(pop_df_7.site_pi_introns))+["LsSpa"]*len(list(pop_df_7.site_pi_codon_1))+["LsSpa"]*len(list(pop_df_7.site_pi_codon_2))+["LsSpa"]*len(list(pop_df_7.site_pi_codon_3))+["LsSpa"]*len(list(pop_df_7.site_pi_codon_4d))
#Population_Lsswe=["Lsswe"]*len(list(pop_df_8.site_pi_global))+["Lsswe"]*len(list(pop_df_8.site_pi_intergenic))+["Lsswe"]*len(list(pop_df_8.site_pi_introns))+["Lsswe"]*len(list(pop_df_8.site_pi_codon_1))+["Lsswe"]*len(list(pop_df_8.site_pi_codon_2))+["Lsswe"]*len(list(pop_df_8.site_pi_codon_3))+["Lsswe"]*len(list(pop_df_8.site_pi_codon_4d))
#
#
#category_LjIre=['all']*len(list(pop_df_1.site_pi_global))+['intergenic']*len(list(pop_df_1.site_pi_intergenic))+['introns']*len(list(pop_df_1.site_pi_introns))+['codon_1']*len(list(pop_df_1.site_pi_codon_1))+['codon_2']*len(list(pop_df_1.site_pi_codon_2))+['codon_3']*len(list(pop_df_1.site_pi_codon_3))+['codon_4d']*len(list(pop_df_1.site_pi_codon_4d))
#category_LjKaz=["all"]*len(list(pop_df_3.site_pi_global))+["intergenic"]*len(list(pop_df_3.site_pi_intergenic))+["introns"]*len(list(pop_df_3.site_pi_introns))+["codon_1"]*len(list(pop_df_3.site_pi_codon_1))+["codon_2"]*len(list(pop_df_3.site_pi_codon_2))+["codon_3"]*len(list(pop_df_3.site_pi_codon_3))+["codon_4d"]*len(list(pop_df_3.site_pi_codon_4d))
#category_LsKaz=["all"]*len(list(pop_df_4.site_pi_global))+["intergenic"]*len(list(pop_df_4.site_pi_intergenic))+["introns"]*len(list(pop_df_4.site_pi_introns))+["codon_1"]*len(list(pop_df_4.site_pi_codon_1))+["codon_2"]*len(list(pop_df_4.site_pi_codon_2))+["codon_3"]*len(list(pop_df_4.site_pi_codon_3))+["codon_4d"]*len(list(pop_df_4.site_pi_codon_4d))
#category_LrSpa=["all"]*len(list(pop_df_6.site_pi_global))+["intergenic"]*len(list(pop_df_6.site_pi_intergenic))+["introns"]*len(list(pop_df_6.site_pi_introns))+["codon_1"]*len(list(pop_df_6.site_pi_codon_1))+["codon_2"]*len(list(pop_df_6.site_pi_codon_2))+["codon_3"]*len(list(pop_df_6.site_pi_codon_3))+["codon_4d"]*len(list(pop_df_6.site_pi_codon_4d))
#category_LsSpa=["all"]*len(list(pop_df_7.site_pi_global))+["intergenic"]*len(list(pop_df_7.site_pi_intergenic))+["introns"]*len(list(pop_df_7.site_pi_introns))+["codon_1"]*len(list(pop_df_7.site_pi_codon_1))+["codon_2"]*len(list(pop_df_7.site_pi_codon_2))+["codon_3"]*len(list(pop_df_7.site_pi_codon_3))+["codon_4d"]*len(list(pop_df_7.site_pi_codon_4d))
#category_Lsswe=["all"]*len(list(pop_df_8.site_pi_global))+["intergenic"]*len(list(pop_df_8.site_pi_intergenic))+["introns"]*len(list(pop_df_8.site_pi_introns))+["codon_1"]*len(list(pop_df_8.site_pi_codon_1))+["codon_2"]*len(list(pop_df_8.site_pi_codon_2))+["codon_3"]*len(list(pop_df_8.site_pi_codon_3))+["codon_4d"]*len(list(pop_df_8.site_pi_codon_4d))
#
#
#Nucleotide_diversity=Nucleotide_diversity_LjIre+Nucleotide_diversity_LjKaz+Nucleotide_diversity_LsKaz+Nucleotide_diversity_LrSpa+Nucleotide_diversity_LsSpa+Nucleotide_diversity_Lsswe
#Population=Population_LjIre+Population_LjKaz+Population_LsKaz+Population_LrSpa+Population_LsSpa+Population_Lsswe
#Category=category_LjIre+category_LjKaz+category_LsKaz+category_LrSpa+category_LsSpa+category_Lsswe
#
#
#violin_df = pandas.DataFrame()
#violin_df['Nucleotide_diversity']=Nucleotide_diversity
#violin_df['Population']=Population
#violin_df['Category']=Category
#pi_colours=sns.color_palette(['#C8C800','#006400','#0000FF','#FF8C00','#FF0000','#C04000'])
#
#
#sns.boxplot(x="Category", y="Nucleotide_diversity", hue="Population",data=violin_df, palette=pi_colours, hue_order=['low','high'], showfliers=False,)



sns.boxplot(pop_df_1[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LjIre_corrected.png')
plt.close()
sns.boxplot(pop_df_3[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LjKaz_corrected.png')
plt.close()
sns.boxplot(pop_df_4[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LsKaz_corrected.png')
plt.close()
sns.boxplot(pop_df_6[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LrSpa_corrected.png')
plt.close()
sns.boxplot(pop_df_7[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LsSpa_corrected.png')
plt.close()
sns.boxplot(pop_df_8[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LsSwe_corrected.png')
plt.close()


sns.boxplot(pop_df_1[['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_codon_4d','irish_juvernica_site_pi_introns','site_pi_intergenic']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LjIre_all.png')
plt.close()
sns.boxplot(pop_df_3[['kazak_juvernica_site_pi_codon_1','kazak_juvernica_site_pi_codon_2','kazak_juvernica_site_pi_codon_3','kazak_juvernica_site_pi_codon_4d','kazak_juvernica_site_pi_introns','site_pi_intergenic']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LjKaz_all.png')
plt.close()
sns.boxplot(pop_df_4[['kaz_sin_site_pi_codon_1','kaz_sin_site_pi_codon_2','kaz_sin_site_pi_codon_3','kaz_sin_site_pi_codon_4d','kaz_sin_site_pi_introns','site_pi_intergenic']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LsKaz_all.png')
plt.close()
sns.boxplot(pop_df_6[['spanish_reali_site_pi_codon_1','spanish_reali_site_pi_codon_2','spanish_reali_site_pi_codon_3','spanish_reali_site_pi_codon_4d','spanish_reali_site_pi_introns','site_pi_intergenic']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LrSpa_all.png')
plt.close()
sns.boxplot(pop_df_7[['spanish_sinapis_site_pi_codon_1','spanish_sinapis_site_pi_codon_2','spanish_sinapis_site_pi_codon_3','spanish_sinapis_site_pi_codon_4d','spanish_sinapis_site_pi_introns','site_pi_intergenic']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LsSpa_all.png')
plt.close()
sns.boxplot(pop_df_8[['swe_sin_allele_site_pi_codon_1','swe_sin_allele_site_pi_codon_2','swe_sin_allele_site_pi_codon_3','swe_sin_allele_site_pi_codon_4d','swe_sin_allele_site_pi_introns','site_pi_intergenic']],showfliers=False)
plt.ylim([0.0,0.015])
plt.savefig('/home/venkat/glob/LsSwe_all.png')
plt.close()





fig, axes  = plt.subplots(nrows=6, ncols=2, sharex=True)
figsize=(10,60)


sns.boxplot(pop_df_1[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']], ax=axes[5,1],showfliers=False)
sns.boxplot(pop_df_3[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']], ax=axes[4,1],showfliers=False)
sns.boxplot(pop_df_4[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']], ax=axes[1,1],showfliers=False)
sns.boxplot(pop_df_6[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']], ax=axes[3,1],showfliers=False)
sns.boxplot(pop_df_7[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']], ax=axes[2,1],showfliers=False)
sns.boxplot(pop_df_8[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_codon4d_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected']], ax=axes[0,1],showfliers=False)


sns.boxplot(pop_df_1[['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_codon_4d','irish_juvernica_site_pi_introns','site_pi_intergenic']], ax=axes[5,0],showfliers=False)
sns.boxplot(pop_df_3[['kazak_juvernica_site_pi_codon_1','kazak_juvernica_site_pi_codon_2','kazak_juvernica_site_pi_codon_3','kazak_juvernica_site_pi_codon_4d','kazak_juvernica_site_pi_introns','site_pi_intergenic']], ax=axes[4,0],showfliers=False)
sns.boxplot(pop_df_4[['kaz_sin_site_pi_codon_1','kaz_sin_site_pi_codon_2','kaz_sin_site_pi_codon_3','kaz_sin_site_pi_codon_4d','kaz_sin_site_pi_introns','site_pi_intergenic']], ax=axes[1,0],showfliers=False)
sns.boxplot(pop_df_6[['spanish_reali_site_pi_codon_1','spanish_reali_site_pi_codon_2','spanish_reali_site_pi_codon_3','spanish_reali_site_pi_codon_4d','spanish_reali_site_pi_introns','site_pi_intergenic']], ax=axes[3,0],showfliers=False)
sns.boxplot(pop_df_7[['spanish_sinapis_site_pi_codon_1','spanish_sinapis_site_pi_codon_2','spanish_sinapis_site_pi_codon_3','spanish_sinapis_site_pi_codon_4d','spanish_sinapis_site_pi_introns','site_pi_intergenic']], ax=axes[2,0],showfliers=False)
sns.boxplot(pop_df_8[['swe_sin_allele_site_pi_codon_1','swe_sin_allele_site_pi_codon_2','swe_sin_allele_site_pi_codon_3','swe_sin_allele_site_pi_codon_4d','swe_sin_allele_site_pi_introns','site_pi_intergenic']], ax=axes[0,0],showfliers=False)


axes[0,1].set_ylim([0.0,0.0033])
axes[1,1].set_ylim([0.0,0.0033])
axes[2,1].set_ylim([0.0,0.0033])
axes[3,1].set_ylim([0.0,0.0033])
axes[4,1].set_ylim([0.0,0.0033])
axes[5,1].set_ylim([0.0,0.0033])
axes[0,0].set_ylim([0.0,0.015])
axes[1,0].set_ylim([0.0,0.015])
axes[2,0].set_ylim([0.0,0.015])
axes[3,0].set_ylim([0.0,0.015])
axes[4,0].set_ylim([0.0,0.015])
axes[5,0].set_ylim([0.0,0.015])

plt.savefig('/home/venkat/glob/codon_stats_all.png')
plt.close()










print pop_df_1[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected','site_pi_codon4d_corrected']].std()
print pop_df_3[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected','site_pi_codon4d_corrected']].std()
print pop_df_6[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected','site_pi_codon4d_corrected']].std()
print pop_df_4[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected','site_pi_codon4d_corrected']].std()
print pop_df_7[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected','site_pi_codon4d_corrected']].std()
print pop_df_8[['site_pi_x','site_pi_codon2_corrected','site_pi_codon3_corrected','site_pi_introns_corrected','site_pi_intergenic_corrected','site_pi_codon4d_corrected']].std()


print pop_df_1[['irish_juvernica_site_pi_codon_1','irish_juvernica_site_pi_codon_2','irish_juvernica_site_pi_codon_3','irish_juvernica_site_pi_introns','site_pi_intergenic','irish_juvernica_site_pi_codon_4d']].std()
print pop_df_3[['kazak_juvernica_site_pi_codon_1','kazak_juvernica_site_pi_codon_2','kazak_juvernica_site_pi_codon_3','kazak_juvernica_site_pi_introns','site_pi_intergenic','kazak_juvernica_site_pi_codon_4d']].std()
print pop_df_4[['kaz_sin_site_pi_codon_1','kaz_sin_site_pi_codon_2','kaz_sin_site_pi_codon_3','kaz_sin_site_pi_introns','site_pi_intergenic','kaz_sin_site_pi_codon_4d']].std()
print pop_df_6[['spanish_reali_site_pi_codon_1','spanish_reali_site_pi_codon_2','spanish_reali_site_pi_codon_3','spanish_reali_site_pi_introns','site_pi_intergenic','spanish_reali_site_pi_codon_4d']].std()
print pop_df_7[['spanish_sinapis_site_pi_codon_1','spanish_sinapis_site_pi_codon_2','spanish_sinapis_site_pi_codon_3','spanish_sinapis_site_pi_introns','site_pi_intergenic','spanish_sinapis_site_pi_codon_4d']].std()
print pop_df_8[['swe_sin_allele_site_pi_codon_1','swe_sin_allele_site_pi_codon_2','swe_sin_allele_site_pi_codon_3','swe_sin_allele_site_pi_introns','site_pi_intergenic','swe_sin_allele_site_pi_codon_4d']].std()


