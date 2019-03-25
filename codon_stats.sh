#! /bin/bash -l
#SBATCH -A snic2017-1-615
#SBATCH -p core  -C mem1TB
#SBATCH -n 16
#SBATCH -t 100:00:00
#SBATCH -J get_codon_stats
#SBATCH --mail-user venkat.talla@ebc.uu.se
#SBATCH --mail-type=ALL



for i in scaf1 scaf2 scaf3 scaf4 scaf5 scaf6 scaf7 scaf8 scaf9 scaf10
do
echo '/'$i'/'
python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/irish_juvernica"  -p $i -c 1 -ind 20 &
python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kazak_juvernica"  -p $i -c 1 -ind 20 &
wait
python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_freq.frq"   -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_reali"    -p $i -c 1 -ind 20 &
python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/swe_sin_allele_freq.frq"  -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/swe_sin_allele"   -p $i -c 1 -ind 20 &
wait
python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq"         -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kaz_sin"          -p $i -c 1 -ind 20 &
python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_sinapis"  -p $i -c 1 -ind 20 &
wait
python codon_stats_with_intergenic.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica_freq.frq'       -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/juvernica'        -p $i -c 1 -ind 40 &
python codon_stats_with_intergenic.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/reali_freq.frq'           -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/reali'            -p $i -c 1 -ind 20 &
python codon_stats_with_intergenic.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis_freq.frq'         -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/sinapis'          -p $i -c 1 -ind 60 &
wait
done



#for i in scaf1 scaf2 scaf3 scaf4 scaf5 scaf6 scaf7 scaf8 scaf9 scaf10
#do
#echo '/'$i'/'
#python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/irish_juvernica"  -p $i -c 2 -ind 20 &
#python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kazak_juvernica"  -p $i -c 2 -ind 20 &
#python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_freq.frq"   -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_reali"    -p $i -c 2 -ind 20 &
#python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/swe_sin_allele_freq.frq"  -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/swe_sin_allele"   -p $i -c 2 -ind 20 &
#python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq"         -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kaz_sin"          -p $i -c 2 -ind 20 &
#wait
#python codon_stats_with_intergenic.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_sinapis"  -p $i -c 2 -ind 20 &
#python codon_stats_with_intergenic.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica_freq.frq'       -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/juvernica'        -p $i -c 2 -ind 40 &
#python codon_stats_with_intergenic.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/reali_freq.frq'           -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/reali'            -p $i -c 2 -ind 20 &
#python codon_stats_with_intergenic.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis_freq.frq'         -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/sinapis'          -p $i -c 2 -ind 60 &
#wait
#done




#for i in scaf1 scaf2 scaf3 scaf4 scaf5 scaf6 scaf7 scaf8 scaf9 scaf10
#do
#echo '/'$i'/'
#wait
#python codon_stats.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/irish_juvernica_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/irish_juvernica"  -p $i -c 5 -ind 20 &
#python codon_stats.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kazak_juvernica_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kazak_juvernica"  -p $i -c 5 -ind 20 &
#python codon_stats.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_reali_freq.frq"   -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_reali"    -p $i -c 5 -ind 20 &
#python codon_stats.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/swe_sin_allele_freq.frq"  -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/swe_sin_allele"   -p $i -c 5 -ind 20 &
#python codon_stats.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/kaz_sin_freq.frq"         -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/kaz_sin"          -p $i -c 5 -ind 20 &
#wait
#python codon_stats.py -f "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/spanish_sinapis_freq.frq" -pop "/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/spanish_sinapis"  -p $i -c 5 -ind 20 &
#python codon_stats.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/juvernica_freq.frq'       -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/juvernica'        -p $i -c 5 -ind 40 &
#python codon_stats.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/reali_freq.frq'           -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/reali'            -p $i -c 5 -ind 20 &
#python codon_stats.py -f '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/pop_stats_shared_private_2/sinapis_freq.frq'         -pop '/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/final_VCF/bam_covrage_files/sinapis'          -p $i -c 5 -ind 60 &
#wait
#done

