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
import Bio.Data.CodonTable
from Bio.Seq import Seq
import argparse
import threading 
import multiprocessing
from multiprocessing import  Pool
from Bio import SeqIO

#conda activate /proj/uppstore2017185/b2014034_nobackup/Venkat/conda_envs

sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})

pwd='/proj/uppstore2017185/b2014034_nobackup/Venkat/Lepidoptera_dnds_project/asymptotic_alpha/'
os.chdir(pwd)

fasta=personal_popgen.fasta_dict('/proj/uppoff2020002/private/result_files/Leptidea/Genome_assembly/HiC_assemblies/LsinapisSweM.fasta')

CDS_fasta_1=personal_popgen.fasta_dict('/lutra/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/fasta/LsinapisSweM.all.maker.proteins.fasta')

CDS_fasta={}

for i in CDS_fasta_1: 
	if CDS_fasta_1[i][0] == 'M':
		if '*' not in CDS_fasta_1[i][1:]:
			CDS_fasta[i]=CDS_fasta_1[i]




def get_intron(CDS):
	df1=CDS
	df1.index = range(1, 2*len(df1)+1, 2)
	df2 = df1.reindex(index=range(2*len(df1)))
	df2 = df2.iloc[1:]
	df2.scaffold = df2.scaffold.ffill()
	df2.source = df2.source.ffill()
	df2.feature = df2.feature.fillna('intron')
	df2.infor = df2.infor.ffill()
	int_starts=list(df2.end.dropna())[:-1]
	int_end=list(df2.start.dropna())[1:]
	df2.loc[df2.start.isnull(), 'start'] = int_starts
	df2.loc[df2.end.isnull(), 'end'] = int_end
	return df2

def join_data_bases(lists):
    for i in range(1,len(lists)):
        j=i-1
        if i == 1:
            new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        else:
            new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        nextone=new_2
    return nextone



def join_raw_data_base(lists):
	for i in range(1,len(lists)):
		j=i-1
		if i == 1:
			new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END'], how='inner', suffixes=('_1', '_2'))
		else:
			new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END'], how='inner')

		nextone=new_2
	return nextone

def join_CDS_seq(lists):
	seq=''
	if lists[0] == '-':
		list_to_use = lists[1:]
		for i in list_to_use:
			my_seq = str(i)
			seq += Seq(my_seq).translate()	
	elif lists[0] == '+':
		list_to_use = lists[1:]
		for i in list_to_use:
			seq += Seq(str(i)).translate()
	return seq


four_d=["TC","CT","CC","CG","GT","GC","GG","AC"]

head_stuff=['scaffold','source','feature','start','end','extra','-','indexing','infor']
#CDSs=pandas.read_table("/home/venkat/bin/snpgenie/test/",skiprows=1, names=head_stuff)

annotation=pandas.read_table("/lutra/proj/uppoff2020002/private/result_files/Leptidea/Annotation/Lsin_swe/gff/LsinapisSweM.all.maker.genes.gff",skiprows=1, names=head_stuff)
annotation=annotation[['scaffold','source','feature','start','end','-','indexing','infor']]

annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='+'), 'start']-1
annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='+'), 'start']-1

annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='-'), 'start']-1
annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='-'), 'start']-1

annotation.loc[(annotation['indexing']=='0')  & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='0') & (annotation['-']=='-'), 'start']-1
annotation.loc[(annotation['indexing']=='0')  & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='0') & (annotation['-']=='+'), 'start']-1

annotation_gene=annotation[(annotation.feature=='mRNA')]
annotation_gene['marker_name']=annotation_gene.infor.str.split(';').str.get(0).str.split('=').str.get(-1)
annotation_gene=annotation_gene.sort_values(by=['marker_name'])
annotation_gene=annotation_gene[annotation_gene.marker_name.isin(CDS_fasta.keys())]

annotation=annotation[annotation.feature=='CDS']
annotation['marker_name']=annotation.infor.str.split(';').str.get(0).str.split('=').str.get(-1).str.split(':').str.get(0)
annotation=annotation[annotation.marker_name.isin(CDS_fasta.keys())]


def main(part):
	match_list_Goods=[]
	match_list_Bad=[]
	list_genes=np.array_split(annotation_gene, 20)[part]
	codon_df_1=pandas.DataFrame(columns=['scaffold','codon1','ref1'])
	codon_df_2=pandas.DataFrame(columns=['scaffold','codon2','ref2'])
	codon_df_3=pandas.DataFrame(columns=['scaffold','codon3','ref3'])
	codon_df_4d=pandas.DataFrame(columns=['CHROM','POS'])
	introns=pandas.DataFrame(columns=['scaffold','POS',])
	for gene_index, gene in list_genes.iterrows():
		#print scaffold
		scaffold=gene.scaffold
		sequence=str(fasta[scaffold]).upper()
		scaffold_gene=annotation_gene[(annotation_gene.scaffold == scaffold)]
		marker_name=gene.marker_name
		gene_all_temp=annotation[(annotation.start >= int(gene['start']-1)) & (annotation.end <= int(gene['end'])) & (gene['scaffold'] == annotation.scaffold)]
		gene_all=gene_all_temp[gene_all_temp.feature=='CDS']
		if len(gene_all) > 0:
			gene_all=gene_all.drop_duplicates(subset=['scaffold', 'source', 'feature', 'start', 'end', '-'], keep='first') 
			CDS_name=gene_all[gene_all.feature=='CDS'].infor.str.split(';').str.get(0).str.split('=').str.get(-1).str.split(':').str.get(0).iloc[0]
			CDS=gene_all[(gene_all.feature=='CDS')]
			CDS=CDS[(CDS.infor.str[3:20]==CDS_name) ]
			exon_stuff=get_intron(CDS)
			exon_stuff.loc[exon_stuff.feature=='intron','start']=exon_stuff.loc[exon_stuff.feature=='intron','start']+1
			exon_stuff.loc[exon_stuff.feature=='intron','end']=exon_stuff.loc[exon_stuff.feature=='intron','end']-1
			#CDS=CDS.append(exon_stuff[exon_stuff.feature=='intron'], ignore_index=True)
			gene_sequnce_cds=CDS_fasta[CDS_name]
			gene_sequnce_extract_list=[gene['-']]
			if CDS.dropna().empty != 'True':
					#CDS.start=CDS.start-1
					#CDS.loc[1,'start'] == '-':
					#CDS.loc[1,'start']=CDS.iloc[0]['start']+1
					#CDS.loc[(CDS.index==len(CDS))&(CDS['-']=='-'),'end']=CDS.loc[(CDS.index==len(CDS))&(CDS['-']=='-'),'end']+1
					for codon_index ,codon in CDS.iterrows():
							if codon.feature=='CDS' and codon['-']=='+':
								start=int(codon.start)
								end=int(codon.end)
								codon1=[i-int(codon.indexing) for i in range(start,end,3)]
								codon2=[i+1 for i in codon1]
								codon3=[i+1 for i in codon2]
								ref1 = [sequence[value:value+1] for value in codon1]
								ref2 = [sequence[value:value+1] for value in codon2]
								ref3 = [sequence[value:value+1] for value in codon3]
								scaffold_1=[codon.scaffold]*len(codon1)
								temp={}
								temp['scaffold']=pandas.Series(scaffold_1)
								temp['codon1']=pandas.Series(codon1)
								temp['codon2']=pandas.Series(codon2)
								temp['codon3']=pandas.Series(codon3)
								temp['ref1']=pandas.Series(ref1)
								temp['ref2']=pandas.Series(ref2)
								temp['ref3']=pandas.Series(ref3)
								temp_df = pandas.DataFrame(temp)
								temp_df[['codon1','codon2','codon3']] = temp_df[(temp_df[['codon1','codon2','codon3']] >= start) & (temp_df[['codon1','codon2','codon3']] <= end)][['codon1','codon2','codon3']]
								#print codon_df
								codon_df_1=codon_df_1.append(temp_df[['scaffold','codon1','ref1']].dropna(),ignore_index=True)
								codon_df_2=codon_df_2.append(temp_df[['scaffold','codon2','ref2']].dropna(),ignore_index=True)
								codon_df_3=codon_df_3.append(temp_df[['scaffold','codon3','ref3']].dropna(),ignore_index=True)
								temp_seq=sequence[start+int(codon.indexing):end]
								gene_sequnce_extract_list.append(temp_seq)
								ref4d = [str(ref1[value]+ref2[value]) for value in range(0,len(ref1))]
								ref4d = [value[0] for value in enumerate(ref4d) if value[1] in four_d ]
								ref4d = [codon3[i] for i in ref4d]
								temp_4d={}
								scaffold_4d=[codon.scaffold]*len(ref4d)
								temp_4d['CHROM']=scaffold_4d
								temp_4d['POS']=ref4d
								temp_df_4d = pandas.DataFrame(temp_4d)
								codon_df_4d=codon_df_4d.append(temp_df_4d,ignore_index=True)
								sys.stdout.flush()
							#elif codon.feature=='CDS' and codon['-']=='-':
							#	#print CDS
							#	#print gene_all[gene_all.feature=='mRNA']['infor']
							#	start=int(codon.start)
							#	end=int(codon.end)
							#	codon3=[i-int(codon.indexing) for i in range(start,end,3)]
							#	codon2=[i+1 for i in codon3]
							#	codon1=[i+1 for i in codon2]
							#	ref1 = [sequence[value:value+1] for value in codon1]
							#	ref2 = [sequence[value:value+1] for value in codon2]
							#	ref3 = [sequence[value:value+1] for value in codon3]
							#	ref4d = [str(ref1[value]+ref2[value]) for value in range(0,len(ref1))]
							#	ref4d = [value[0] for value in enumerate(ref4d) if value[1] in four_d ]
							#	ref4d = [codon3[i] for i in ref4d]
							#	scaffold_1=[codon.scaffold]*len(codon1)
							#	temp={}
							#	temp['scaffold']=pandas.Series(codon.scaffold)
							#	temp['codon1']=pandas.Series(codon1)
							#	temp['codon2']=pandas.Series(codon2)
							#	temp['codon3']=pandas.Series(codon3)
							#	temp['ref1']=pandas.Series(ref1)
							#	temp['ref2']=pandas.Series(ref2)
							#	temp['ref3']=pandas.Series(ref3)
							#	temp_df = pandas.DataFrame(temp)
							#	temp_df[['codon1','codon2','codon3']] = temp_df[(temp_df[['codon1','codon2','codon3']] >= start) & (temp_df[['codon1','codon2','codon3']] <= end)][['codon1','codon2','codon3']]
							#	#print codon_df
							#	codon_df_1=codon_df_1.append(temp_df[['scaffold','codon1','ref1']],ignore_index=True)
							#	codon_df_2=codon_df_2.append(temp_df[['scaffold','codon2','ref2']],ignore_index=True)
							#	codon_df_3=codon_df_3.append(temp_df[['scaffold','codon3','ref3']],ignore_index=True)
							#	temp_seq=str(Seq(sequence[start:end-int(codon.indexing)]).reverse_complement())
							#	gene_sequnce_extract_list.append(temp_seq)
							#	ref4d = [str(ref1[value]+ref2[value]) for value in range(0,len(ref1))]
							#	ref4d = [value[0] for value in enumerate(ref4d) if value[1] in four_d ]
							#	ref4d = [codon3[i] for i in ref4d]
							#	temp_4d={}
							#	scaffold_4d=[codon.scaffold]*len(ref4d)
							#	temp_4d['CHROM']=scaffold_4d
							#	temp_4d['POS']=ref4d
							#	temp_df_4d = pandas.DataFrame(temp_4d)
							#	codon_df_4d=codon_df_4d.append(temp_df_4d,ignore_index=True)
							#	sys.stdout.flush()
							#elif codon.feature=='intron':
							#	temp_int={}
							#	positions=range(int(codon.start), int(codon.end))
							#	final=[]
							#	scaffold_int=[codon.scaffold]*len(positions)
							#	temp_int['scaffold']=scaffold_int
							#	temp_int['POS']=positions
							#	temp_df_int = pandas.DataFrame(temp_int)
							#	introns=introns.append(temp_df_int,ignore_index=True)
							#	sys.stdout.flush()
				#print gene_sequnce_cds
				#print "\n"
				#print join_CDS_seq(gene_sequnce_extract_list)
				#print "\n"
				#print CDS
				#print (CDS_name+'    '+str(gene_sequnce_cds) == join_CDS_seq(gene_sequnce_extract_list))
			final_seq=str(join_CDS_seq(gene_sequnce_extract_list))[0:-1]
			#if final_seq == str(gene_sequnce_cds):
				#print (gene.infor)
				#print ('Match  '+marker_name+' '+CDS_name+' '+str(scaffold)+ ' '+str(codon['-']))
				#match_list_Goods.append(marker_name)
			#else:
				#print ('No Match  '+marker_name+' '+CDS_name+' '+str(scaffold)+ ' '+str(codon['-']))
				#match_list_Bad.append(marker_name)
			print (CDS)
			print ('\n')
			print (final_seq)
			print (gene_sequnce_cds)
			print ('End of gene'+'\n')
	return [codon_df_1, codon_df_2, codon_df_3, codon_df_4d, introns, match_list_Goods, match_list_Bad]



p = Pool(processes=20)
data = p.map(main, [i for i in range(20)])
p.close()



codon_df_1 = [i[0] for i in data] 
codon_df_2 = [i[1] for i in data] 
codon_df_3 = [i[2] for i in data] 
codon_df_4d = [i[3] for i in data] 
introns = [i[4] for i in data] 
goods = [j  for i in data for j in i[-2]] 
bads = [j for i in data for j in i[-1]] 


codon_df_1=pandas.concat(codon_df_1, ignore_index=True)
codon_df_2=pandas.concat(codon_df_2, ignore_index=True)
codon_df_3=pandas.concat(codon_df_3, ignore_index=True)
codon_df_4d=pandas.concat(codon_df_4d, ignore_index=True)
introns=pandas.concat(introns, ignore_index=True)

codon_df_1=codon_df_1.dropna()
codon_df_2=codon_df_2.dropna()
codon_df_3=codon_df_3.dropna()
codon_df_4d=codon_df_4d.dropna()

codon_df_1=codon_df_1.rename(columns={'codon1': 'POS'})
codon_df_2=codon_df_2.rename(columns={'codon2': 'POS'})
codon_df_3=codon_df_3.rename(columns={'codon3': 'POS'})
introns=introns.rename(columns={'codon3': 'POS'})

codon_df_4d=codon_df_4d.rename(columns={'CHROM': 'scaffold'})


codon_df_1.POS=codon_df_1.POS.astype(int)
codon_df_2.POS=codon_df_2.POS.astype(int)
codon_df_3.POS=codon_df_3.POS.astype(int)
introns.POS=introns.POS.astype(int)
codon_df_4d.POS=codon_df_4d.POS.astype(int)



codon_df_1.POS=codon_df_1.POS+1
codon_df_2.POS=codon_df_2.POS+1
codon_df_3.POS=codon_df_3.POS+1
codon_df_4d.POS=codon_df_4d.POS+1



codon_df_1=codon_df_1.rename(columns={'scaffold': 'CHROM'})
codon_df_2=codon_df_2.rename(columns={'scaffold': 'CHROM'})
codon_df_3=codon_df_3.rename(columns={'scaffold': 'CHROM'})
codon_df_4d=codon_df_4d.rename(columns={'scaffold': 'CHROM'})

introns=introns.rename(columns={'scaffold': 'CHROM'})




introns.POS=introns.POS+1


#start_stuff.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/start_stuff_'+part+'.csv', sep=' ', header=False, index= False)
introns.to_csv('introns_new.csv', sep='\t', header=False, index= False)
codon_df_1[['CHROM', 'POS']].drop_duplicates(subset=['CHROM','POS']).to_csv('codon_df_1_New.csv', sep='\t', header=False, index= False)
codon_df_2[['CHROM', 'POS']].drop_duplicates(subset=['CHROM','POS']).to_csv('codon_df_2_New.csv', sep='\t', header=False, index= False)
codon_df_3[['CHROM', 'POS']].drop_duplicates(subset=['CHROM','POS']).to_csv('codon_df_3_New.csv', sep='\t', header=False, index= False)
codon_df_4d[['CHROM', 'POS']].drop_duplicates(subset=['CHROM','POS']).to_csv('codon_df_4d_New.csv', sep='\t', header=False, index= False)







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
#start_stuff['BIN_START']=0
#start_stuff['BIN_END']=0
#
#
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
#start_stuff['BIN_START']=(np.floor(start_stuff['POS']/100000)*100000)+1
#start_stuff['BIN_END']=start_stuff['BIN_START']+(100000-1)
#
#
#start_stuff.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/start_stuff.csv', sep=' ', header=False, index= False)
#
#introns.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/introns.csv', sep=' ', header=False, index= False)
#codon_df_1.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_1.csv', sep=' ', header=False, index= False)
#codon_df_2.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_2.csv', sep=' ', header=False, index= False)
#codon_df_3.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_3.csv', sep=' ', header=False, index= False)
#codon_df_4d.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/codon_df_4d.csv', sep=' ', header=False, index= False)

