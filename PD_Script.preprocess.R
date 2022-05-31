#WGCNA FOR multiple datasets
getwd()
setwd("D:/Machine learning/PD/Parkinsons_disease")
getwd()
install.packages("BiocManager")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("hgu133a.db")
BiocManager::install("hgu133b.db")
BiocManager::install("hgu133plus2.db")
BiocManager::install("a4Base")
BiocManager::install("oligo")
BiocManager::install("EnrichmentBrowser")
BiocManager::install("limma")
BiocManager::install("illuminaio")
BiocManager::install("sva")

getwd()
##calling library for function##
library(BiocManager)
library(affy)
library(oligo) #if both affy and oligo packaages are used specify the function with package name #
library(Biobase)
library(preprocessCore)
library(limma)
library(sva)
library("annotate")
library("a4Base")
library("EnrichmentBrowser")
library(ggfortify)
library(tibble)
library(plyr)
library(ggplot2)
library(limma)
library(illuminaio)
library(lumi)

##############################################################################################################################
#DATASET 1
#preprocessing dataset GSE6613#annotation#Affymetrix Human Genome U133A Array#
celpath = "D:/Machine learning/PD/Parkinsons_disease/PD_Datasets/GSE6613_RAW"
list = list.files(celpath,full.names=TRUE)
GSE6613 = read.celfiles(list)
Norm1<-rma(GSE6613)
par(mar=c(5,5,5,5))
boxplot(GSE6613, las=2)
boxplot(Norm1, las=2)
PROBE1<-EnrichmentBrowser::probe2gene(Norm1, chip = "hgu133plus2", from = "PROBEID",to = "SYMBOL",multi.to = "first",multi.from = "mean")
Expr1<-assay(PROBE1, "exprs")           
STR1<-strsplit(colnames(Expr1), split = ".CEL.gz")
colnames(Expr1)<-STR1
ColS1<-rownames_to_column(as.data.frame(Expr1), var = "SYMBOL")
REP1<-avereps(ColS1, ColS1$SYMBOL)
df1<-as.data.frame(REP1)
write.csv(df1, file="df1_GSE6613.csv")

#DATASET 2
#preprocessing dataset GSE6613#annotation#Affymetrix Human Genome U133A Array#
celpath2 = "D:/Machine learning/PD/Parkinsons_disease/PD_Datasets/GSE72267_RAW"
list2 = list.files(celpath2,full.names=TRUE)
GSE72267 = read.celfiles(list2)
Norm2<-rma(GSE72267)
par(mar=c(5,5,5,5))
boxplot(GSE72267, las=2)
boxplot(Norm2, las=2)
PROBE2<-EnrichmentBrowser::probe2gene(Norm2, chip = "hgu133plus2", from = "PROBEID",to = "SYMBOL",multi.to = "first",multi.from = "mean")
Expr2<-assay(PROBE2, "exprs")           
STR2<-strsplit(colnames(Expr2), split = ".CEL.gz")
colnames(Expr2)<-STR2
ColS2<-rownames_to_column(as.data.frame(Expr2), var = "SYMBOL")
REP2<-avereps(ColS2, ColS2$SYMBOL)
df2<-as.data.frame(REP2)
write.csv(df2, file = "df2_GSE72267.csv")

#DATASET 3
#preprocessing dataset GSE99039#annotation#Affymetrix Human Genome U133Aplus2 Array#
celpath3 = "D:/Machine learning/PD/Parkinsons_disease/PD_Datasets/GSE99039_RAW"
list3 = list.files(celpath3,full.names=TRUE)
GSE99039 = read.celfiles(list3)
Norm3<-rma(GSE99039)
par(mar=c(10,5,5,5))
boxplot(GSE99039, las=2)
boxplot(Norm3, las=2)
PROBE3<-EnrichmentBrowser::probe2gene(Norm3, chip = "hgu133plus2", from = "PROBEID",to = "SYMBOL",multi.to = "first",multi.from = "mean")
Expr3<-assay(PROBE3, "exprs") 
STR3<-strsplit(colnames(Expr3), split = "_HG-U133_Plus_2_.CEL.gz")
#T.Expr3<-t(Expr3)
#rownames(T.Expr3)<- gsub("_.*", "", rownames(T.Expr3))
colnames(Expr3)<-STR3
COLS3<-rownames_to_column(as.data.frame(Expr3), var ="SYMBOL")
REP3<-avereps(COLS3, COLS3$SYMBOL)
df3<-as.data.frame(REP3)
write.csv(df3, file = "df3_GSE99039.csv")

#illumina dataset preproccessing
#DATASET 4
#Seperately processed for GSE57475 as its illumina and file saved as GSE57475_final, script name GSE57475_illumn_script.R
getwd()
mywd57475=setwd("D:/Machine learning/PD/Parkinsons_disease/PD_Datasets/GSE57475_RAW")
getwd()
GSE57475<-lumiR("GSE57475_non-normalized.txt", sep = NULL, na.rm = TRUE, convertNuID = TRUE, lib.mapping = 'lumiHumanIDMapping', dec = '.', parseColumnName = TRUE, checkDupId = TRUE, QC = TRUE,
                columnNameGrepPattern = list(exprs='SAMPLE', se.exprs='NA', detection="Detection", beadNum='NA'))
# filter probes based on detection p value, such that the applied condition is satisfied in atleast 
# >= 3
fil_probes_57475=rowSums(GSE57475@assayData$detection<0.05)>=3 #changing probe to gene
View(fil_probes_57475)#The number of 'TRUE' character in the object represents the actual number
# of probes satisfied the condition. 
GSE57475=GSE57475[fil_probes_57475,]# Add these filtered probes to lumiBatchobject
dim(GSE57475) #16236      142
# Pheno data_View
# the series matrix from the GEO had more sample information compared to the GSEphenodata from the non normalised file
# The information from the series matrix file was merged with the phenodata from GSE38481
GSE57475@phenoData
pData(GSE57475)
# Read series matrix file
Pheno_57475=read.delim("GSE57475_series_matrix_phenodata.txt", colClasses="factor", header = T)
Pheno_57475[]=lapply(Pheno_57475, as.character)
phenoData(GSE57475) <- AnnotatedDataFrame(Pheno_57475) #Add phenotypic info (source and Gender) to the pDATA 
phenoData(GSE57475)
pData(GSE57475)
colnames(GSE57475)=GSE57475@phenoData@data$GSM_ID # Change the sample ID in GSE38481 S4 class object to GSM IDs from phenodata
View(exprs(GSE57475))
dim(GSE57475)

# Quantile normalization
colnames(phenoData(GSE57475))
colnames(GSE57475)
norm57475=lumiN(GSE57475, method = "quantile") # watch for warning message for NANs
enorm57475<-exprs(norm57475)

# Check for negative values and assign 0 
has.neg57475=exprs(norm57475)
has.neg57475 <- apply(has.neg57475, 1, function(row) any(row < 0))
length(which(has.neg57475))#1407 rows with negative values
df57475_exp=exprs(norm57475)
df57475_exp[df57475_exp<0]=0 # Assign 0 to negative values
df57475_exp[df57475_exp<0] # Confirm if there are any negative values

# Gene symbols/EntrezIDs (nuID2EntrezID=for Entrez IDs and nuID2targetID for Gene Symbols)
library(magrittr) #works well w/o this package
library(tibble)

df57475_exp=as.data.frame(df57475_exp)
df57475_exp=rownames_to_column(df57475_exp)
entrezid_df57475=nuID2targetID(df57475_exp$rowname, lib.mapping='lumiHumanIDMapping')
entrezid_df57475=rownames_to_column(data.frame(entrezid_df57475), var = "rowname")
df57475_join=plyr::join(entrezid_df57475,df57475_exp , by= "rowname", type="inner")
View(df57475_join)
df57475_GID=df57475_join[,c(-1)]
View((df57475_GID))
dim(df57475_GID) #16023 143
n_occur57475 <- data.frame(table(df57475_GID$entrezid_df57475))## Check for duplicates (optional)

# gives you a data frame with a list of ids and the number of times they occurred.
n_occur57475[n_occur57475$Freq > 1,] #tells you which ids occurred more than once.
df57475_GID[df57475_GID$entrezid_df57475 %in% n_occur57475$Var1[n_occur57475$Freq > 1],]# returns the records with more than one occurrence.

## Averaging out the intensities for duplicates using avereps function
df57475_final=limma::avereps(df57475_GID, ID= df57475_GID$entrezid_df57475) %>% data.frame
dim(df57475_GID)  #16236 143
dim(df57475_final) # reduced to 14442
# check and remove unannotated genes
table(is.na(df57475_final))# Assuming that unannotated genes would have NA is.na was used
# No NAs were detected. 

# GeneIDs/entrezids as Rownames
df57475_final=column_to_rownames(df57475_final, "entrezid_df57475")
dim(df57475_final) #14442 142
df57475_final[df57475_final<0]

df4<-rownames_to_column(df57475_final)
colnames(df4)[1] <- "SYMBOL"

# df57475_final is a df with log2 transformed values and can be used for meta-analysis
write.csv(df4, "df4_GSE57475_ill.csv")

##############################################################################################################################



