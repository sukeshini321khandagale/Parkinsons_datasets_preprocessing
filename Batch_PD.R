##Combining illumina and affymetrix datasets using combat##
#merging 4 datasets in form of dataframe#
library(dplyr)

df_datasets2<-merge(df1, df2, by="SYMBOL") #GSE6613 & GSE72267
df_datasets3<-merge(df_datasets2, df3, by="SYMBOL") #df_datasets2 & GSE99039
df_datasets4<-merge(df_datasets3, df4, by="SYMBOL") #df-datasets3 & df_57475
#final merged datasets will be df_datasets4
write.csv(df_datasets4, file = "PD_metadata_sn712.csv")

#after merge check for similar dimentions#
dim(df_datasets4) #sample 712, features=3886

#Specifying the type of batch samples belong to#
pdata<- read.csv("D:/Machine learning/Nikita_machine_learning/PD_Meta-phenodata.csv")
dim(pdata)
#very imp when you have to convert matrix to numeric#no gene symbol #
data_PD<-column_to_rownames(df_datasets4, var = "SYMBOL")
is.numeric(data_PD) #gives false
mutdata<-data_PD %>% mutate_if(is.character,as.numeric)#this also converts matrix to numeric but also gives gene symbol#
#T1D <- apply(as.matrix(T1Ddata), 2,  as.numeric) this converts matrix to numeric but deletes gene symbols#


#Batch correction using GSEID
Batch<- pdata$GSE_ID
modcombat = model.matrix(~1, data=pdata)

#Batch correction with matrix with telling groups depending on GSEID#
batch_correct<- ComBat (mutdata,
                        batch = Batch,
                        mod = modcombat,
                        par.prior = TRUE,
                        prior.plots = FALSE)

library(ggplot2)
color=pdata$GSE_ID
color
autoplot(prcomp(t(mutdata)), scale=T, data = pdata, colour="GSE_ID")+ ggtitle("PCA plot for raw data")
autoplot(prcomp(t(batch_correct)), scale=T, data = pdata, colour="GSE_ID")+ ggtitle("PCA plot for Batch Corrected data")

write.csv(batch_correct, file = "PD_BatchCorrected_metadata.csv")

#using limma for batch correction##commandas below
##PD_Metadata<-removeBatchEffect(mutdata, batch = Batch, covariates = NULL) #function by Limma for batch correction#
##autoplot(prcomp(t(PD_Metadata)), scale=T, data = pdata, colour="GSE_ID")+ ggtitle("PCA plot for raw data")
