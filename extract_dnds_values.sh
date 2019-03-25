#!/bin/bash -l

#SBATCH -A snic2017-1-615
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J paml
#SBATCH -o paml.out
#SBATCH -e paml.error
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL




module  load bioinfo-tools paml
gene_seqs='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments/'


echo 'gene_name' > gene_name_temp
echo 'ds' > ds_temp
echo 'dn' > dn_temp
echo 'w' > w_temp

paste   gene_name_temp   ds_temp   dn_temp   w_temp > final_dnds_result_species

for i in $(ls $gene_seqs);
do
perl catfasta2phyml.pl -c -s -v '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments/'$i > /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_analysis/temp.nuc

codeml codeml.ctl

echo $i > gene_name_temp
less /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_analysis/temp.out  | tail -10 | head -2 | tail -1 > ds_temp
less /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_analysis/temp.out  | tail -10 | head -4 | tail -1 > dn_temp
less /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_analysis/temp.out  | tail -10 | head -7 | tail -1 > w_temp

paste   gene_name_temp   ds_temp   dn_temp   w_temp >> final_dnds_result_species

less /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_analysis/temp.out | tail -20 | head -6 > '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_out_puts/'$i'.out'

done



