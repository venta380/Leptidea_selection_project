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
from Bio.Seq import MutableSeq
import argparse
from Bio import codonalign
from Bio.codonalign.codonseq import cal_dn_ds
from Bio.codonalign.codonseq import default_codon_table
from Bio.codonalign import CodonSeq
from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna


pwd='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/'
os.chdir(pwd)



fasta_1=personal_popgen.fasta_dict('sinapis.CDS.fasta')
fasta_2=personal_popgen.fasta_dict('spanish_reali.CDS.fasta')
fasta_3=personal_popgen.fasta_dict('juvernica.CDS.fasta')




#DM_genes=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/gene_ortholog_flybase/overlap_new.csv", names=['DM_gene','LS_gene'], sep=' ')


for i in fasta_1.keys():
		F = open("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments_species/"+str(i)+'.fasta',"w") 
		#F.write(">"+str(gene_name)+'_L_Sin'+'\n')
		F.write(">"+'L_Sin'+'\n')
		F.write(str(fasta_1[i])+'\n')
		#F.write(">"+str(gene_name)+'_L_Rea'+'\n')
		F.write(">"+'L_Rea'+'\n')
		F.write(str(fasta_2[i])+'\n')
		#F.write(">"+str(gene_name)+'_L_Juv'+'\n')
		F.write(">"+'L_Juv'+'\n')
		F.write(str(fasta_3[i])+'\n')
		F.close()
		


fasta_4=personal_popgen.fasta_dict("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/irish_juvernica.CDS.fasta")
fasta_5=personal_popgen.fasta_dict("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_juvernica.CDS.fasta")
fasta_6=personal_popgen.fasta_dict("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_reali.CDS.fasta")
fasta_7=personal_popgen.fasta_dict("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_sinapis.CDS.fasta")
fasta_8=personal_popgen.fasta_dict("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_sinapis.CDS.fasta")
fasta_9=personal_popgen.fasta_dict("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/Swedish_sinapis.CDS.fasta")


for i in fasta_4.keys():
		F = open("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments/"+str(i)+'.fasta',"w") 
		F.write(">"+'irish_juvernica'+'\n')
		F.write(str(fasta_4[i])+'\n')
		F.write(">"+'kazak_juvernica'+'\n')
		F.write(str(fasta_5[i])+'\n')
		F.write(">"+'spanish_reali'+'\n')
		F.write(str(fasta_6[i])+'\n')
		F.write(">"+'kazak_sinapis'+'\n')
		F.write(str(fasta_7[i])+'\n')
		F.write(">"+'spanish_sinapis'+'\n')
		F.write(str(fasta_8[i])+'\n')
		F.write(">"+'Swedish_sinapis'+'\n')
		F.write(str(fasta_9[i])+'\n')
		F.close()
	else:
		F = open("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments/"+str(i)+'.fasta',"w") 
		F.write(">"+'irish_juvernica'+'\n')
		F.write(str(fasta_4[i])+'\n')
		F.write(">"+'kazak_juvernica'+'\n')
		F.write(str(fasta_5[i])+'\n')
		F.write(">"+'spanish_reali'+'\n')
		F.write(str(fasta_6[i])+'\n')
		F.write(">"+'kazak_sinapis'+'\n')
		F.write(str(fasta_7[i])+'\n')
		F.write(">"+'spanish_sinapis'+'\n')
		F.write(str(fasta_8[i])+'\n')
		F.write(">"+'Swedish_sinapis'+'\n')
		F.write(str(fasta_9[i])+'\n')
		F.close()




