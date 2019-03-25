library(seqinr)

files=read.table('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/scripts/allignment_list_species.txt')

col.names = c("gene", "Dn_sin_reali", "Dn_sin_juv", "Dn_juv_reali", "Ds_sin_reali", "Ds_sin_juv", "Ds_juv_reali", "Dn_ds_sin_reali", "Dn_ds_sin_juv", "Dn_ds_juv_reali")
df <- read.table(text = "",col.names = col.names)



for (i in files$V1){
x=read.alignment(i,format="fasta")
if(getLength(x$seq[[1]][1]) > 24){
 
	x_new=kaks(x)
	(x_new$a0+(x_new$l0*x_new$b0)
	final=data.frame(c(strsplit(i, "/")[[1]][8], x_new$ka[1], x_new$ka[2], x_new$ka[3], x_new$ks[1], x_new$ks[2], x_new$ks[3], x_new$ka[1]/x_new$ks[1], x_new$ka[2]/x_new$ks[2], x_new$ka[3]/x_new$ks[3]))
	names(final)=1
	final=t(final)
	names(final)=col.names
	df <- rbind(df, final)

}
}


names(df)=col.names

write.csv(df, file = "dn_ds.csv", sep=" ")


((x_new$l2 * x_new$a2) + (x_new$l4 * x_new$a4))/




ka=(x_new$l0 * x_new$b0) + (x_new$l2 * x_new$b2)
ks=


