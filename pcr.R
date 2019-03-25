require(pls)
set.seed (1000)


#install.packages("Hmisc")
#install.packages("corrgram")
#install.packages("ellipse")
#install.packages("lmSupport")
#install.packages("RColorBrewer")


library(car)
library(Hmisc)
library(corrgram)
library(ellipse)
library(lmSupport)
library(RColorBrewer)



#irish_juvernica <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/irish_juvernica_input_for_PCR.csv', header=TRUE)
#kazak_juvernica <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_juvernica_input_for_PCR.csv', header=TRUE)
#kazak_sinapis <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/kazak_sinapis_input_for_PCR.csv', header=TRUE)
#spanish_reali <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_reali_input_for_PCR.csv', header=TRUE)
#spanish_sinapis <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/spanish_sinapis_input_for_PCR.csv', header=TRUE)
#swedish_sinapis <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/swedish_sinapis_input_for_PCR.csv', header=TRUE)

irish_juvernica <- read.csv('/proj/uppstore2017185/b2014034_nobackup/POPULATION_RESEQ/dnds/temp_files/swedish_sinapis_input_for_PCR.csv', header=TRUE)


colnames(irish_juvernica) <- c('site_pi_codon_4d',"gene_density","GC_reference","rho","sites")


irish_juvernica.trans <- within(irish_juvernica, {
  site_pi_codon_4d <- sqrt(site_pi_codon_4d)
  gene_density <- logit(gene_density)
  GC_reference <- logit(GC_reference)
  rho <- sqrt(rho)
  sites <- sites
})





irish_juvernica.ztrans <- irish_juvernica.trans
irish_juvernica.ztrans[c("gene_density","GC_reference","rho")] <- lapply(irish_juvernica.ztrans[c("gene_density","GC_reference","rho")],scale)

cor(irish_juvernica.ztrans[c("gene_density","GC_reference","rho")])

#write.csv(irish_juvernica.ztrans, file = "temp.csv",row.names=TRUE)

#######################################Multiple Regression


#fit_all <- glm(site_pi_codon_4d ~ (GC_reference * rho)+gene_density+GC_reference+rho, data=irish_juvernica.ztrans)
fit_all <- glm(site_pi_codon_4d ~ GC_reference * rho * gene_density - GC_reference:rho:gene_density, data=irish_juvernica.ztrans)


summary(fit_all)
#drop1(fit_all, test="F")






fit_all <- glm(site_pi_codon_4d ~ gene_density+GC_reference+rho+(gene_density * GC_reference) , data=irish_juvernica.ztrans)


#summary(fit_all)
drop1(fit_all, test="F")



fit_all_sep <- glm(site_pi_codon_4d ~ gene_density+GC_reference+rho, data=irish_juvernica.ztrans)

summary(fit_all_sep)
drop1(fit_all_sep, test="F")

##################PCR#######################

### Principal component regression

library(pls)

# dataset preparation
data.ztrans_all <- irish_juvernica.ztrans[,-7]
data.ztrans_all$site_pi_codon_4d <- scale(irish_juvernica.ztrans$site_pi_codon_4d)
head(data.ztrans_all)

#PCR-Principal Component Regression
model.pcr <- pcr(site_pi_codon_4d~as.matrix(data.ztrans_all[,c("gene_density","GC_reference","rho")]), data=data.ztrans_all,  validation="CV", y=T)
summary(model.pcr)
str(summary(model.pcr))

summary(lm(model.pcr$y ~ model.pcr$scores))

# Explained variance of pi by component and variable
cumYexplained <- 100 * drop(R2(model.pcr, estimate = "train", intercept = FALSE)$val) 
Yexplained  <- cumYexplained - c(0,cumYexplained[1:4]); Yexplained
rel.contr <-  model.pcr$projection^2;rownames(rel.contr) <- colnames(data.ztrans_all[,c("gene_density","GC_reference","rho")]); rel.contr
rel.contrtoR2 <- t(t(rel.contr)*Yexplained);rel.contrtoR2
str(rel.contrtoR2)
rel.contrtoR2_ordered <- rel.contrtoR2[c("gene_density","GC_reference","rho"),]

t(rel.contrtoR2_ordered)






#########################################















# Check distribution of z-transformed variables 
pdf("z_transformed.pdf",paper="a4")
par(mfrow=c(3,2))
for (i in 1:ncol(irish_juvernica.ztrans)){
  Sys.sleep(0.1)
  plot(density(irish_juvernica.ztrans[,i]), main=colnames(irish_juvernica.ztrans)[i])
}
dev.off()


## Check correlations in the z-transformed dataset

pdf("pairwise_correlations.pdf",paper="a4")
irish_juvernica.ztrans.corr <- irish_juvernica.ztrans[,1:7]
colnames(irish_juvernica.ztrans.corr) <- c('site_pi_codon_4d',"gene_density","GC_reference","rho","sites","N.dN","dS")
pairs(irish_juvernica.ztrans.corr)
plotcorr(cor(irish_juvernica.ztrans.corr))
corrgram(irish_juvernica.ztrans.corr,lower.panel=panel.shade, upper.panel=panel.pie)
corrgram(irish_juvernica.ztrans.corr,lower.panel=panel.ellipse, upper.panel=panel.pts)
dev.off()
cor(irish_juvernica.ztrans.corr)
rcorr(as.matrix(irish_juvernica.ztrans.corr),type="pearson")






#Check linear relationship between predictors and dependent variables

pdf("residuals_for_each_predictor.pdf",paper="a4", height=9)
par(mfrow=c(5,2), mar=c(4,4,1,1))
model<- lm(site_pi_codon_4d ~ gene_density, weights=sites, data=irish_juvernica.ztrans)
plot(model$fitted.values, model$residuals*sqrt(irish_juvernica.ztrans$sites),xlab="Fitted values", ylab="Weighted Residuals") # weighted residuals vs fitted
abline(h=0)
lo <- loess(model$residuals*sqrt(irish_juvernica.ztrans$sites)~model$fitted.values)
j <- order(model$fitted.values)
lines(model$fitted.values[j],lo$fitted[j],col="red",lwd=1)
qqnorm(model$residuals*sqrt(irish_juvernica.ztrans$sites));qqline(model$residuals*sqrt(irish_juvernica.ztrans$sites))

model<- lm(site_pi_codon_4d ~ GC_reference, weights=sites, data=irish_juvernica.ztrans)
plot(model$fitted.values, model$residuals*sqrt(irish_juvernica.ztrans$sites),xlab="Fitted values", ylab="Weighted Residuals") # weighted residuals vs fitted
abline(h=0)
lo <- loess(model$residuals*sqrt(irish_juvernica.ztrans$sites)~model$fitted.values)
j <- order(model$fitted.values)
lines(model$fitted.values[j],lo$fitted[j],col="red",lwd=1)
qqnorm(model$residuals*sqrt(irish_juvernica.ztrans$sites));qqline(model$residuals*sqrt(irish_juvernica.ztrans$sites))


model<- lm(site_pi_codon_4d ~ rho, weights=sites, data=irish_juvernica.ztrans)
plot(model$fitted.values, model$residuals*sqrt(irish_juvernica.ztrans$sites),xlab="Fitted values", ylab="Weighted Residuals") # weighted residuals vs fitted
abline(h=0)
lo <- loess(model$residuals*sqrt(irish_juvernica.ztrans$sites)~model$fitted.values)
j <- order(model$fitted.values)
lines(model$fitted.values[j],lo$fitted[j],col="red",lwd=1)
qqnorm(model$residuals*sqrt(irish_juvernica.ztrans$sites));qqline(model$residuals*sqrt(irish_juvernica.ztrans$sites))


model<- lm(site_pi_codon_4d ~ N.dN, weights=sites, data=irish_juvernica.ztrans)
plot(model$fitted.values, model$residuals*sqrt(irish_juvernica.ztrans$sites),xlab="Fitted values", ylab="Weighted Residuals") # weighted residuals vs fitted
abline(h=0)
lo <- loess(model$residuals*sqrt(irish_juvernica.ztrans$sites)~model$fitted.values)
j <- order(model$fitted.values)
lines(model$fitted.values[j],lo$fitted[j],col="red",lwd=1)
qqnorm(model$residuals*sqrt(irish_juvernica.ztrans$sites));qqline(model$residuals*sqrt(irish_juvernica.ztrans$sites))

model<- lm(site_pi_codon_4d ~ dS, weights=sites, data=irish_juvernica.ztrans)
plot(model$fitted.values, model$residuals*sqrt(irish_juvernica.ztrans$sites),xlab="Fitted values", ylab="Weighted Residuals") # weighted residuals vs fitted
abline(h=0)
lo <- loess(model$residuals*sqrt(irish_juvernica.ztrans$sites)~model$fitted.values)
j <- order(model$fitted.values)
lines(model$fitted.values[j],lo$fitted[j],col="red",lwd=1)
qqnorm(model$residuals*sqrt(irish_juvernica.ztrans$sites));qqline(model$residuals*sqrt(irish_juvernica.ztrans$sites))

dev.off()




## Linear model

#Transformed reponse and transformed and scaled predictors, residuals weighted by "sites"
model<- lm(site_pi_codon_4d ~ gene_density+GC_reference+rho+N.dN+dS, weights=sites, data=irish_juvernica.ztrans)
summary(model)
par(mfrow=c(2,2))
plot(model)

pdf("residual_plots.pdf",paper="a4",width=4, height=9)
par(mfrow=c(2,1))
plot(model$fitted.values, model$residuals*sqrt(irish_juvernica.ztrans$sites),xlab="Fitted values", ylab="Weighted Residuals", main="Residuals vs Fitted") # weighted residuals vs fitted
abline(h=0)
lo <- loess(model$residuals*sqrt(irish_juvernica.ztrans$sites)~model$fitted.values)
j <- order(model$fitted.values)
lines(model$fitted.values[j],lo$fitted[j],col="red",lwd=1)
qqnorm(model$residuals*sqrt(irish_juvernica.ztrans$sites));qqline(model$residuals*sqrt(irish_juvernica.ztrans$sites))
dev.off()

# global F, F stats, single hypothesis test and R squared
summary(model)
drop1(model, test="F")
step(model, test="F")



### Principal component regression

library(pls)

# dataset preparation
data.ztrans_all <- irish_juvernica.ztrans[,-7]
data.ztrans_all$site_pi_codon_4d <- scale(irish_juvernica.ztrans$site_pi_codon_4d)
head(data.ztrans_all)

#PCR-Principal Component Regression
model.pcr <- pcr(site_pi_codon_4d~as.matrix(data.ztrans_all[,c("gene_density","GC_reference","Pin_Pis","rho")]), data=data.ztrans_all,  validation="CV", y=T)
summary(model.pcr)
str(summary(model.pcr))

summary(lm(model.pcr$y ~ model.pcr$scores))

# Explained variance of pi by component and variable
cumYexplained <- 100 * drop(R2(model.pcr, estimate = "train", intercept = FALSE)$val) 
Yexplained  <- cumYexplained - c(0,cumYexplained[1:4]); Yexplained
rel.contr <-  model.pcr$projection^2;rownames(rel.contr) <- colnames(data.ztrans_all[,c("gene_density","GC_reference","Pin_Pis","rho")]); rel.contr
rel.contrtoR2 <- t(t(rel.contr)*Yexplained);rel.contrtoR2
str(rel.contrtoR2)
rel.contrtoR2_ordered <- rel.contrtoR2[c("gene_density","GC_reference","Pin_Pis","rho"),]

t(rel.contrtoR2_ordered)


pdf("PCR.pdf",paper="a4")
par(mfrow=(c(1,1)),mar=c(5,5,4,3))
barplot(rel.contrtoR2_ordered, col=brewer.pal(4,"Set1"),legend=c("gene_density","GC_reference","Pin_Pis","rho"),las=1, names.arg=c("1","2","3","4"), xlab="Components", ylab=expression(paste("% of variance explained (R"^" 2",")")),ylim=c(-1,1), args.legend=list(bty="n"),axes=F); axis(2,seq(-1,1,5),las=2)
dev.off()

