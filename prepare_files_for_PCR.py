dnds_table_pop=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/final_dnds_populations_result")
dnds_table_pop=dnds_table_pop.drop_duplicates(subset='gene_name')

dnds_table_pop=dnds_table_pop[dnds_table_pop.gene_name!='leptidea_sinapisT00000000201.fasta']


dnds_table_pop['Ds_LsSwe']=dnds_table_pop['ds'].str.split(':').str[6].str.split(',').str[0]
dnds_table_pop['Ds_LsKaz']=dnds_table_pop['ds'].str.split(':').str[5].str.split(',').str[0]
dnds_table_pop['Ds_LsSpa']=dnds_table_pop['ds'].str.split(':').str[4].str.split(',').str[0]
dnds_table_pop['Ds_LrSpa']=dnds_table_pop['ds'].str.split(':').str[3].str.split(',').str[0]
dnds_table_pop['Ds_LjKaz']=dnds_table_pop['ds'].str.split(':').str[2].str.split(',').str[0].str[:9]
dnds_table_pop['Ds_LjIre']=dnds_table_pop['ds'].str.split(':').str[1].str.split(',').str[0]


dnds_table_pop['Dn_LsSwe']=dnds_table_pop['dn'].str.split(':').str[6].str.split(',').str[0]
dnds_table_pop['Dn_LsKaz']=dnds_table_pop['dn'].str.split(':').str[5].str.split(',').str[0]
dnds_table_pop['Dn_LsSpa']=dnds_table_pop['dn'].str.split(':').str[4].str.split(',').str[0]
dnds_table_pop['Dn_LrSpa']=dnds_table_pop['dn'].str.split(':').str[3].str.split(',').str[0]
dnds_table_pop['Dn_LjKaz']=dnds_table_pop['dn'].str.split(':').str[2].str.split(',').str[0].str[:9]
dnds_table_pop['Dn_LjIre']=dnds_table_pop['dn'].str.split(':').str[1].str.split(',').str[0]

dnds_table_pop['geneID']=dnds_table_pop['gene_name'].str.split('.').str[0]


dnds_table_pop=join_raw_data_base([dnds_table[['geneID','N_sites', 'S_sites']], dnds_table_pop,],['geneID'],'inner')

dnds_table_pop=dnds_table_pop[['geneID',"Ds_LsSwe","Ds_LsKaz","Ds_LsSpa","Ds_LrSpa","Ds_LjKaz","Ds_LjIre","Dn_LsSwe","Dn_LsKaz","Dn_LsSpa","Dn_LrSpa","Dn_LjKaz","Dn_LjIre", 'N_sites', 'S_sites']]
dnds_table_pop=dnds_table_pop.drop_duplicates(subset='geneID')

dnds_table_pop[["Ds_LsSwe","Ds_LsKaz","Ds_LsSpa","Ds_LrSpa","Ds_LjKaz","Ds_LjIre","Dn_LsSwe","Dn_LsKaz","Dn_LsSpa","Dn_LrSpa","Dn_LjKaz","Dn_LjIre"]]=dnds_table_pop[["Ds_LsSwe","Ds_LsKaz","Ds_LsSpa","Ds_LrSpa","Ds_LjKaz","Ds_LjIre","Dn_LsSwe","Dn_LsKaz","Dn_LsSpa","Dn_LrSpa","Dn_LjKaz","Dn_LjIre"]].astype('float')


dnds_table_pop['s_subs_LsSwe']=(dnds_table_pop['Ds_LsSwe']*dnds_table_pop['S_sites']).round(0)
dnds_table_pop['s_subs_LsKaz']=(dnds_table_pop['Ds_LsKaz']*dnds_table_pop['S_sites']).round(0)
dnds_table_pop['s_subs_LsSpa']=(dnds_table_pop['Ds_LsSpa']*dnds_table_pop['S_sites']).round(0)
dnds_table_pop['s_subs_LrSpa']=(dnds_table_pop['Ds_LrSpa']*dnds_table_pop['S_sites']).round(0)
dnds_table_pop['s_subs_LjKaz']=(dnds_table_pop['Ds_LjKaz']*dnds_table_pop['S_sites']).round(0)
dnds_table_pop['s_subs_LjIre']=(dnds_table_pop['Ds_LjIre']*dnds_table_pop['S_sites']).round(0)


dnds_table_pop['n_subs_LsSwe']=(dnds_table_pop['Dn_LsSwe']*dnds_table_pop['N_sites']).round(0)
dnds_table_pop['n_subs_LsKaz']=(dnds_table_pop['Dn_LsKaz']*dnds_table_pop['N_sites']).round(0)
dnds_table_pop['n_subs_LsSpa']=(dnds_table_pop['Dn_LsSpa']*dnds_table_pop['N_sites']).round(0)
dnds_table_pop['n_subs_LrSpa']=(dnds_table_pop['Dn_LrSpa']*dnds_table_pop['N_sites']).round(0)
dnds_table_pop['n_subs_LjKaz']=(dnds_table_pop['Dn_LjKaz']*dnds_table_pop['N_sites']).round(0)
dnds_table_pop['n_subs_LjIre']=(dnds_table_pop['Dn_LjIre']*dnds_table_pop['N_sites']).round(0)

dnds_table_gene_pop=join_raw_data_base([annotation,dnds_table_pop,],['geneID'],'inner')

dnds_table_windows_pop=dnds_table_gene_pop.groupby(['CHROM','BIN_START','BIN_END', 'auto_allo'],as_index=False).agg({'N_sites': 'sum','S_sites': 'sum', 's_subs_LsSwe': 'sum','s_subs_LsKaz': 'sum','s_subs_LsSpa': 'sum','s_subs_LrSpa': 'sum','s_subs_LjKaz': 'sum','s_subs_LjIre': 'sum','n_subs_LsSwe': 'sum','n_subs_LsKaz': 'sum','n_subs_LsSpa': 'sum','n_subs_LrSpa': 'sum','n_subs_LjKaz': 'sum','n_subs_LjIre': 'sum'})

dnds_table_windows_pop['Dn_LsSwe']=dnds_table_windows_pop['n_subs_LsSwe']/dnds_table_windows_pop['N_sites']
dnds_table_windows_pop['Dn_LsKaz']=dnds_table_windows_pop['n_subs_LsKaz']/dnds_table_windows_pop['N_sites']
dnds_table_windows_pop['Dn_LsSpa']=dnds_table_windows_pop['n_subs_LsSpa']/dnds_table_windows_pop['N_sites']
dnds_table_windows_pop['Dn_LrSpa']=dnds_table_windows_pop['n_subs_LrSpa']/dnds_table_windows_pop['N_sites']
dnds_table_windows_pop['Dn_LjKaz']=dnds_table_windows_pop['n_subs_LjKaz']/dnds_table_windows_pop['N_sites']
dnds_table_windows_pop['Dn_LjIre']=dnds_table_windows_pop['n_subs_LjIre']/dnds_table_windows_pop['N_sites']

dnds_table_windows_pop['Ds_LsSwe']=dnds_table_windows_pop['s_subs_LsSwe']/dnds_table_windows_pop['S_sites']
dnds_table_windows_pop['Ds_LsKaz']=dnds_table_windows_pop['s_subs_LsKaz']/dnds_table_windows_pop['S_sites']
dnds_table_windows_pop['Ds_LsSpa']=dnds_table_windows_pop['s_subs_LsSpa']/dnds_table_windows_pop['S_sites']
dnds_table_windows_pop['Ds_LrSpa']=dnds_table_windows_pop['s_subs_LrSpa']/dnds_table_windows_pop['S_sites']
dnds_table_windows_pop['Ds_LjKaz']=dnds_table_windows_pop['s_subs_LjKaz']/dnds_table_windows_pop['S_sites']
dnds_table_windows_pop['Ds_LjIre']=dnds_table_windows_pop['s_subs_LjIre']/dnds_table_windows_pop['S_sites']







pop_df_1['gene_density']=((pop_df_1['irish_juvernica_sites_codon_3']+pop_df_1['irish_juvernica_sites_codon_2']+pop_df_1['irish_juvernica_sites_codon_1'])/pop_df_1['irish_juvernica_sites_global'])*100
pop_df_3['gene_density']=((pop_df_3['kazak_juvernica_sites_codon_3']+pop_df_3['kazak_juvernica_sites_codon_2']+pop_df_3['kazak_juvernica_sites_codon_1'])/pop_df_3['kazak_juvernica_sites_global'])*100 
pop_df_4['gene_density']=((pop_df_4['kaz_sin_sites_codon_3']+pop_df_4['kaz_sin_sites_codon_2']+pop_df_4['kaz_sin_sites_codon_1'])/pop_df_4['kaz_sin_sites_global'])*100 
pop_df_6['gene_density']=((pop_df_6['spanish_reali_sites_codon_3']+pop_df_6['spanish_reali_sites_codon_2']+pop_df_6['spanish_reali_sites_codon_1'])/pop_df_6['spanish_reali_sites_global'])*100 
pop_df_7['gene_density']=((pop_df_7['spanish_sinapis_sites_codon_3']+pop_df_7['spanish_sinapis_sites_codon_2']+pop_df_7['spanish_sinapis_sites_codon_1'])/pop_df_7['spanish_sinapis_sites_global'])*100 
pop_df_8['gene_density']=((pop_df_8['swe_sin_allele_sites_codon_3']+pop_df_8['swe_sin_allele_sites_codon_2']+pop_df_8['swe_sin_allele_sites_codon_1'])/pop_df_8['swe_sin_allele_sites_global'])*100 


pop_df_1=pop_df_1[(pop_df_1['gene_density']>0)]
pop_df_3=pop_df_3[(pop_df_3['gene_density']>0)]
pop_df_4=pop_df_4[(pop_df_4['gene_density']>0)]
pop_df_6=pop_df_6[(pop_df_6['gene_density']>0)]
pop_df_7=pop_df_7[(pop_df_7['gene_density']>0)]
pop_df_8=pop_df_8[(pop_df_8['gene_density']>0)]
#
#
#pop_df_1['gene_density'][~(pop_df_1['gene_density']>0)]=1.0
#pop_df_3['gene_density'][~(pop_df_3['gene_density']>0)]=1.0
#pop_df_4['gene_density'][~(pop_df_4['gene_density']>0)]=1.0
#pop_df_6['gene_density'][~(pop_df_6['gene_density']>0)]=1.0
#pop_df_7['gene_density'][~(pop_df_7['gene_density']>0)]=1.0
#pop_df_8['gene_density'][~(pop_df_8['gene_density']>0)]=1.0

print pop_df_1.irish_juvernica_sites_codon_4d.sum()
print pop_df_3.kazak_juvernica_sites_codon_4d.sum()
print pop_df_4.kaz_sin_sites_codon_4d.sum()
print pop_df_6.spanish_reali_sites_codon_4d.sum()
print pop_df_7.spanish_sinapis_sites_codon_4d.sum()
print pop_df_8.swe_sin_allele_sites_codon_4d.sum()


rho_2=rho[['CHROM','BIN_START','BIN_END','irish_juvernica_rho','kazak_juvernica_rho','kazak_sinapis_rho','spanish_reali_rho','spanish_sinapis_rho','swedish_sinapis_rho']]

GC_ref=dnds_table_windows_2[['CHROM','BIN_START','BIN_END','GC_reference']]

new_1=join_raw_data_base([fianl_df_1_new, pop_df_1, GC_ref, dnds_table_windows_pop , rho_2[['CHROM','BIN_START','BIN_END','irish_juvernica_rho']]],['CHROM','BIN_START','BIN_END'],'inner')
new_3=join_raw_data_base([fianl_df_3_new, pop_df_3, GC_ref, dnds_table_windows_pop , rho_2[['CHROM','BIN_START','BIN_END','kazak_juvernica_rho']]],['CHROM','BIN_START','BIN_END'],'inner')
new_4=join_raw_data_base([fianl_df_4_new, pop_df_4, GC_ref, dnds_table_windows_pop , rho_2[['CHROM','BIN_START','BIN_END','kazak_sinapis_rho']]],['CHROM','BIN_START','BIN_END'],'inner')
new_6=join_raw_data_base([fianl_df_6_new, pop_df_6, GC_ref, dnds_table_windows_pop , rho_2[['CHROM','BIN_START','BIN_END','spanish_reali_rho']]],['CHROM','BIN_START','BIN_END'],'inner')
new_7=join_raw_data_base([fianl_df_7_new, pop_df_7, GC_ref, dnds_table_windows_pop , rho_2[['CHROM','BIN_START','BIN_END','spanish_sinapis_rho']]],['CHROM','BIN_START','BIN_END'],'inner')
new_8=join_raw_data_base([fianl_df_8_new, pop_df_8, GC_ref, dnds_table_windows_pop , rho_2[['CHROM','BIN_START','BIN_END','swedish_sinapis_rho']]],['CHROM','BIN_START','BIN_END'],'inner')







new_1=new_1[['irish_juvernica_site_pi_codon_4d', 'gene_density', 'GC_reference', 'Pin_Pis_irish_juvernica', 'irish_juvernica_rho', "irish_juvernica_sites_global_x", 'n_subs_LjIre','Ds_LjIre']]
new_3=new_3[['kazak_juvernica_site_pi_codon_4d', 'gene_density', 'GC_reference', 'Pin_Pis_kazak_juvernica', 'kazak_juvernica_rho', "kazak_juvernica_sites_global_x", 'n_subs_LjKaz','Ds_LjKaz']]
new_4=new_4[['kaz_sin_site_pi_codon_4d', 'gene_density', 'GC_reference', 'Pin_Pis_kazak_sinapis', 'kazak_sinapis_rho', "kaz_sin_sites_global_x", 'n_subs_LsKaz','Ds_LsKaz']]
new_6=new_6[['spanish_reali_site_pi_codon_4d', 'gene_density', 'GC_reference', 'Pin_Pis_spanish_reali', 'spanish_reali_rho', "spanish_reali_sites_global_x",'n_subs_LrSpa','Ds_LrSpa']]
new_7=new_7[['spanish_sinapis_site_pi_codon_4d', 'gene_density', 'GC_reference', 'Pin_Pis_spanish_sinapis', 'spanish_sinapis_rho', "spanish_sinapis_sites_global_x",'n_subs_LsSpa','Ds_LsSpa']]
new_8=new_8[['swe_sin_allele_site_pi_codon_4d', 'gene_density', 'GC_reference', 'Pin_Pis_Swedish_sinapis', 'swedish_sinapis_rho', "swe_sin_allele_sites_global_x", 'n_subs_LsSwe','Ds_LsSwe']]

X=list(new_1.gene_density)
Y=list(new_1.GC_reference)
Z=list(new_1.irish_juvernica_rho)


new_1=new_1[new_1.irish_juvernica_site_pi_codon_4d > 0]
new_3=new_3[new_3.kazak_juvernica_site_pi_codon_4d > 0]
new_4=new_4[new_4.kaz_sin_site_pi_codon_4d > 0]
new_6=new_6[new_6.spanish_reali_site_pi_codon_4d > 0]
new_7=new_7[new_7.spanish_sinapis_site_pi_codon_4d > 0]
new_8=new_8[new_8.swe_sin_allele_site_pi_codon_4d > 0]

new_1=new_1[['irish_juvernica_site_pi_codon_4d', 'gene_density', 'GC_reference', 'irish_juvernica_rho', "irish_juvernica_sites_global_x",]]
new_3=new_3[['kazak_juvernica_site_pi_codon_4d', 'gene_density', 'GC_reference', 'kazak_juvernica_rho', "kazak_juvernica_sites_global_x",]]
new_4=new_4[['kaz_sin_site_pi_codon_4d', 'gene_density', 'GC_reference', 'kazak_sinapis_rho', "kaz_sin_sites_global_x",]]
new_6=new_6[['spanish_reali_site_pi_codon_4d', 'gene_density', 'GC_reference', 'spanish_reali_rho', "spanish_reali_sites_global_x"]]
new_7=new_7[['spanish_sinapis_site_pi_codon_4d', 'gene_density', 'GC_reference', 'spanish_sinapis_rho', "spanish_sinapis_sites_global_x"]]
new_8=new_8[['swe_sin_allele_site_pi_codon_4d', 'gene_density', 'GC_reference', 'swedish_sinapis_rho', "swe_sin_allele_sites_global_x",]]




new_1=new_1[['irish_juvernica_site_pi_codon_4d', 'gene_density', 'GC_reference', 'irish_juvernica_rho', "irish_juvernica_sites_global_x", 'n_subs_LjIre','Ds_LjIre']]
new_3=new_3[['kazak_juvernica_site_pi_codon_4d', 'gene_density', 'GC_reference', 'kazak_juvernica_rho', "kazak_juvernica_sites_global_x", 'n_subs_LjKaz','Ds_LjKaz']]
new_4=new_4[['kaz_sin_site_pi_codon_4d', 'gene_density', 'GC_reference', 'kazak_sinapis_rho', "kaz_sin_sites_global_x", 'n_subs_LsKaz','Ds_LsKaz']]
new_6=new_6[['spanish_reali_site_pi_codon_4d', 'gene_density', 'GC_reference', 'spanish_reali_rho', "spanish_reali_sites_global_x",'n_subs_LrSpa','Ds_LrSpa']]
new_7=new_7[['spanish_sinapis_site_pi_codon_4d', 'gene_density', 'GC_reference', 'spanish_sinapis_rho', "spanish_sinapis_sites_global_x",'n_subs_LsSpa','Ds_LsSpa']]
new_8=new_8[['swe_sin_allele_site_pi_codon_4d', 'gene_density', 'GC_reference', 'swedish_sinapis_rho', "swe_sin_allele_sites_global_x", 'n_subs_LsSwe','Ds_LsSwe']]


#new_1['gene_density']=(new_1['gene_density']- new_1['gene_density'].mean())/new_1['gene_density'].std()
#new_3['gene_density']=(new_3['gene_density']- new_3['gene_density'].mean())/new_3['gene_density'].std()
#new_4['gene_density']=(new_4['gene_density']- new_4['gene_density'].mean())/new_4['gene_density'].std()
#new_6['gene_density']=(new_6['gene_density']- new_6['gene_density'].mean())/new_6['gene_density'].std()
#new_7['gene_density']=(new_7['gene_density']- new_7['gene_density'].mean())/new_7['gene_density'].std()
#new_8['gene_density']=(new_8['gene_density']- new_8['gene_density'].mean())/new_8['gene_density'].std()
#
#
#new_1['GC_reference']=(new_1['GC_reference']- new_1['GC_reference'].mean())/new_1['GC_reference'].std()
#new_3['GC_reference']=(new_3['GC_reference']- new_3['GC_reference'].mean())/new_3['GC_reference'].std()
#new_4['GC_reference']=(new_4['GC_reference']- new_4['GC_reference'].mean())/new_4['GC_reference'].std()
#new_6['GC_reference']=(new_6['GC_reference']- new_6['GC_reference'].mean())/new_6['GC_reference'].std()
#new_7['GC_reference']=(new_7['GC_reference']- new_7['GC_reference'].mean())/new_7['GC_reference'].std()
#new_8['GC_reference']=(new_8['GC_reference']- new_8['GC_reference'].mean())/new_8['GC_reference'].std()
#
#
#new_1['Pin_Pis_irish_juvernica']=(new_1['Pin_Pis_irish_juvernica']- new_1['Pin_Pis_irish_juvernica'].mean())/new_1['Pin_Pis_irish_juvernica'].std()
#new_3['Pin_Pis_kazak_juvernica']=(new_3['Pin_Pis_kazak_juvernica']- new_3['Pin_Pis_kazak_juvernica'].mean())/new_3['Pin_Pis_kazak_juvernica'].std()
#new_4['Pin_Pis_kazak_sinapis']=(new_4['Pin_Pis_kazak_sinapis']- new_4['Pin_Pis_kazak_sinapis'].mean())/new_4['Pin_Pis_kazak_sinapis'].std()
#new_6['Pin_Pis_spanish_reali']=(new_6['Pin_Pis_spanish_reali']- new_6['Pin_Pis_spanish_reali'].mean())/new_6['Pin_Pis_spanish_reali'].std()
#new_7['Pin_Pis_spanish_sinapis']=(new_7['Pin_Pis_spanish_sinapis']- new_7['Pin_Pis_spanish_sinapis'].mean())/new_7['Pin_Pis_spanish_sinapis'].std()
#new_8['Pin_Pis_Swedish_sinapis']=(new_8['Pin_Pis_Swedish_sinapis']- new_8['Pin_Pis_Swedish_sinapis'].mean())/new_8['Pin_Pis_Swedish_sinapis'].std()
#
#new_1['irish_juvernica_rho']=(new_1['irish_juvernica_rho']- new_1['irish_juvernica_rho'].mean())/new_1['irish_juvernica_rho'].std()
#new_3['kazak_juvernica_rho']=(new_3['kazak_juvernica_rho']- new_3['kazak_juvernica_rho'].mean())/new_3['kazak_juvernica_rho'].std()
#new_4['kazak_sinapis_rho']=(new_4['kazak_sinapis_rho']- new_4['kazak_sinapis_rho'].mean())/new_4['kazak_sinapis_rho'].std()
#new_6['spanish_reali_rho']=(new_6['spanish_reali_rho']- new_6['spanish_reali_rho'].mean())/new_6['spanish_reali_rho'].std()
#new_7['spanish_sinapis_rho']=(new_7['spanish_sinapis_rho']- new_7['spanish_sinapis_rho'].mean())/new_7['spanish_sinapis_rho'].std()
#new_8['swedish_sinapis_rho']=(new_8['swedish_sinapis_rho']- new_8['swedish_sinapis_rho'].mean())/new_8['swedish_sinapis_rho'].std()
#
#
#new_1['irish_juvernica_site_pi_codon_4d']=(new_1['irish_juvernica_site_pi_codon_4d']- new_1['irish_juvernica_site_pi_codon_4d'].mean())/new_1['irish_juvernica_site_pi_codon_4d'].std()
#new_3['kazak_juvernica_site_pi_codon_4d']=(new_3['kazak_juvernica_site_pi_codon_4d']- new_3['kazak_juvernica_site_pi_codon_4d'].mean())/new_3['kazak_juvernica_site_pi_codon_4d'].std()
#new_4['kaz_sin_site_pi_codon_4d']=(new_4['kaz_sin_site_pi_codon_4d']- new_4['kaz_sin_site_pi_codon_4d'].mean())/new_4['kaz_sin_site_pi_codon_4d'].std()
#new_6['spanish_reali_site_pi_codon_4d']=(new_6['spanish_reali_site_pi_codon_4d']- new_6['spanish_reali_site_pi_codon_4d'].mean())/new_6['spanish_reali_site_pi_codon_4d'].std()
#new_7['spanish_sinapis_site_pi_codon_4d']=(new_7['spanish_sinapis_site_pi_codon_4d']- new_7['spanish_sinapis_site_pi_codon_4d'].mean())/new_7['spanish_sinapis_site_pi_codon_4d'].std()
#new_8['swe_sin_allele_site_pi_codon_4d']=(new_8['swe_sin_allele_site_pi_codon_4d']- new_8['swe_sin_allele_site_pi_codon_4d'].mean())/new_8['swe_sin_allele_site_pi_codon_4d'].std()



new_1.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/irish_juvernica_input_for_PCR.csv', index=False)
new_3.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_juvernica_input_for_PCR.csv', index=False)
new_4.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_sinapis_input_for_PCR.csv', index=False)
new_6.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_reali_input_for_PCR.csv', index=False)
new_7.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_sinapis_input_for_PCR.csv', index=False)
new_8.to_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/swedish_sinapis_input_for_PCR.csv', index=False)



corr_pvalue(new_1)
corr_pvalue(new_3)
corr_pvalue(new_4)
corr_pvalue(new_6)
corr_pvalue(new_7)
corr_pvalue(new_8)


sns.residplot(new_1['irish_juvernica_site_pi_codon_4d'], new_1['irish_juvernica_rho'])





pcr=pandas.read_csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/PCR_plot_vaues.csv',index_col=False)


pcr_LjIre=pcr[pcr.Pop=='LjIre'][["Comp","gene_density","GC_reference","rho"]]
pcr_LjKaz=pcr[pcr.Pop=='LjKaz'][["Comp","gene_density","GC_reference","rho"]]
pcr_LrSpa=pcr[pcr.Pop=='LrSpa'][["Comp","gene_density","GC_reference","rho"]]
pcr_LsSpa=pcr[pcr.Pop=='LsSpa'][["Comp","gene_density","GC_reference","rho"]]
pcr_LsKaz=pcr[pcr.Pop=='LsKaz'][["Comp","gene_density","GC_reference","rho"]]
pcr_LsSwe=pcr[pcr.Pop=='LsSwe'][["Comp","gene_density","GC_reference","rho"]]

pcr_LjIre[["gene_density","GC_reference","rho"]]=pcr_LjIre[["gene_density","GC_reference","rho"]].abs()
pcr_LjKaz[["gene_density","GC_reference","rho"]]=pcr_LjKaz[["gene_density","GC_reference","rho"]].abs()
pcr_LrSpa[["gene_density","GC_reference","rho"]]=pcr_LrSpa[["gene_density","GC_reference","rho"]].abs()
pcr_LsSpa[["gene_density","GC_reference","rho"]]=pcr_LsSpa[["gene_density","GC_reference","rho"]].abs()
pcr_LsKaz[["gene_density","GC_reference","rho"]]=pcr_LsKaz[["gene_density","GC_reference","rho"]].abs()
pcr_LsSwe[["gene_density","GC_reference","rho"]]=pcr_LsSwe[["gene_density","GC_reference","rho"]].abs()



pcr_LjIre.plot.bar(x='Comp', stacked=True)
plt.savefig('/home/venkat/pcr_LjIre.png')
plt.close()
pcr_LjKaz.plot.bar(x='Comp', stacked=True)
plt.savefig('/home/venkat/pcr_LjKaz.png')
plt.close()
pcr_LrSpa.plot.bar(x='Comp', stacked=True)
plt.savefig('/home/venkat/pcr_LrSpa.png')
plt.close()
pcr_LsSpa.plot.bar(x='Comp', stacked=True)
plt.savefig('/home/venkat/pcr_LsSpa.png')
plt.close()
pcr_LsKaz.plot.bar(x='Comp', stacked=True)
plt.savefig('/home/venkat/pcr_LsKaz.png')
plt.close()
pcr_LsSwe.plot.bar(x='Comp', stacked=True)
plt.savefig('/home/venkat/pcr_LsSwe.png')
plt.close()

