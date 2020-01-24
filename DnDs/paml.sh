#!/bin/bash -l

#SBATCH -A snic2017-1-615
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH -J paml_1
#SBATCH -o paml_1.out
#SBATCH -e paml_1.error
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL




module  load bioinfo-tools paml
gene_seqs='/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments_species/'

cd /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/_1/

echo 'gene_name' > gene_name_temp
echo 'dn' > dn_temp
echo 'ds' > ds_temp
echo 'N_sites' > N_temp
echo 'S_sites' > S_temp


paste   gene_name_temp   dn_temp   ds_temp   N_temp   S_temp > final_dnds_result

for i in $(ls $gene_seqs | head -1560);
do
perl /proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/catfasta2phyml.pl -c -s -v '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/gene_allignments_species/'$i > temp_1.nuc

codeml codeml.ctl

echo $i > gene_name_temp
less temp_1.out  | tail -10 | head -2 | tail -1 > ds_temp
less temp_1.out  | tail -10 | head -4 | tail -1 > dn_temp
less temp_1.out  | tail -20 | less temp_1.out | tail -18 | head -1  |  tr -s [:blank:] | cut -f 4 -d " " > N_temp
less temp_1.out  | tail -20 | less temp_1.out | tail -18 | head -1  |  tr -s [:blank:] | cut -f 5 -d " " > S_temp
paste   gene_name_temp   dn_temp   ds_temp   N_temp   S_temp >> final_dnds_result

#less temp_1.out | tail -20 | head -6 > '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/paml_out_puts/'$i'.out'

done


