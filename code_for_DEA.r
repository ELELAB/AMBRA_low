#This script downloads both  the primary tumor barcodes with an AMBRA1 gene expression lower than the
#minimum value of the normal samples and the normal samples and perform the Differential Expression Analysis.
#We performed the Differential Expression Analysis on the eight dataset showed in the table:

#Cancer Type,min
#BLCA,705
#COAD,1120 
#KIRC,998
#KIRP,858
#LUAD,533
#LUSC,962
#PRAD,1295
#UCEC,850

# The usage of the script is illustrated for the BLCA dataset, but it is possible to use it on a different dataset replacing the
# name of the cancer type and the specific value of the minimum expression value

#R-packages needed
library(SummarizedExperiment)  
library(TCGAbiolinks) 

#To serch  specific GDC data for download:
queryBLCAtp <- GDCquery(project = "TCGA-BLCA",       
                        data.category = "Gene expression",                       
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        legacy = TRUE,
                        sample.type = c("Primary solid Tumor")) 
samplesTPDown <- queryBLCAtp$results[[1]]$cases
dataSmTP_BLCAtp <- TCGAquery_SampleTypes(barcode = samplesTPDown,
                                         typesample = "TP")    #To select primary solid tumor samples
GDCdownload(queryBLCAtp)
BLCAtp <- GDCprepare(query = queryBLCAtp, save = FALSE)  #To reads the data downloaded and prepare it into an R object.

#To remove samples outliers  using a Pearson correlation cutoff of 0.6.
dataPrepBLCAtp <- TCGAanalyze_Preprocessing(object =BLCAtp , cor.cut = 0.6)  


gene1_ID <- "AMBRA1|55626" #The gene of interest

dataPrepBLCAtp_t <- t (dataPrepBLCAtp) #Transpose matrix
#To select barcodes with the AMBRA1 expression higher than the minimum value in the normal samples and save it in a table
gene1_BLCAtp_t <- dataPrepBLCAtp_t[,gene1_ID, drop=FALSE]
gene1_BLCAtp_t_noambralow <- subset(gene1_BLCAtp_t, gene1_BLCAtp_t[,1] >705) 
write.table(gene1_BLCAtp_t_noambralow, file="BLCAtp_noAMBRAlow.txt", sep = " ")

#To know the number of this samples
nrow(gene1_BLCAtp_t_noambralow)



#To search  specific GDC data for download
queryBLCAnt <- GDCquery(project = "TCGA-BLCA",                 
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        legacy = TRUE,
                        sample.type = c("Solid Tissue Normal"))    #Normal samples
samplesNTDown <- queryBLCAnt$results[[1]]$cases
dataSmNT_BLCAnt <- TCGAquery_SampleTypes(barcode = samplesNTDown,
                                         typesample = "NT")
GDCdownload(queryBLCAnt)
BLCAnt <- GDCprepare(query = queryBLCAnt, save = FALSE) #To read the data downloaded and prepare it into an R object.
#To remove samples outliers  using a Pearson correlation cutoff of 0.6.
dataPrepBLCAnt <- TCGAanalyze_Preprocessing(object =BLCAnt , cor.cut = 0.6)
dataPrepBLCAnt_t <- t (dataPrepBLCAnt) #transpose matrix
gene1_ID <- "AMBRA1|55626"  #The gene of interest
#to select barcodes  of the normal samples and save it in a table
gene1_BLCAnt_t <- dataPrepBLCAnt_t[,gene1_ID, drop=FALSE]
write.table (gene1_BLCAnt_t, file="BLCAnt.csv", sep = " ")
read.csv("BLCAnt.csv")

#To know the number of this samples
nrow(gene1_BLCAnt_t)

#Extract the barcodes of the AMBRAlow

barcodes <- rownames(gene1_BLCAtp_t_noambralow)
barcodes

#Extract the barcodes of the normal sample

barcodesNT <- rownames(gene1_BLCAnt_t)
barcodesNT

#Differential expression analysis:

queryBLCAtp_noAMBRAlow <- GDCquery(project = "TCGA-BLCA",   #to load the AMBRAlow barcodes
                                   data.category = "Gene expression",
                                   data.type = "Gene expression quantification",
                                   platform = "Illumina HiSeq",
                                   file.type = "results",
                                   experimental.strategy = "RNA-Seq",
                                   barcode = barcodes,
                                   legacy = TRUE,
                                   sample.type = c("Primary solid Tumor"))

queryBLCAnt <- GDCquery(project ="TCGA-BLCA",               #to load the normal samples barcodes
                        data.category = "Gene expression",
                        data.type = "Gene expression quantification",
                        platform = "Illumina HiSeq",
                        file.type = "results",
                        experimental.strategy = "RNA-Seq",
                        barcode = barcodesNT,
                        legacy = TRUE,
                        sample.type= c("Solid Tissue Normal"))

dataBLCA_noAMBRAlow <- TCGAquery_SampleTypes(barcode = barcodes, typesample = c("TP"))
dataBLCAnt <- TCGAquery_SampleTypes(barcode = barcodesNT, typesample =c("NT"))     #to return both barcode "TP" and "NT"
querydown_BLCA_noAMBRAlowvsNT <- GDCquery(project ="TCGA-BLCA",                     #to search data with these arguments
                                          legacy = TRUE,
                                          data.category = "Gene expression",
                                          data.type = "Gene expression quantification",
                                          platform = "Illumina HiSeq",
                                          file.type = "results",
                                          experimental.strategy = "RNA-Seq",
                                          barcode = c(dataBLCA_noAMBRAlow, dataBLCAnt),
                                          sample.type = c("Primary solid Tumor", "Solid Tissue Normal"))
GDCdownload(querydown_BLCA_noAMBRAlowvsNT)
BLCA_noAMBRAlowvsNT <- GDCprepare(query = querydown_BLCA_noAMBRAlowvsNT, save = FALSE)  #reads the data downloaded and prepare it into an R object.
dataPrepBLCA <- TCGAanalyze_Preprocessing(object = BLCA_noAMBRAlowvsNT , cor.cut = 0.6) #remove samples outliers  using a Pearson correlation cutoff of 0.6.
write.table (dataPrepBLCA, file= "dataPrepBLCA_noAMBRAlowvsNT.txt", sep=" ", quote= FALSE)
dataNormBLCA <- TCGAanalyze_Normalization (tabDF = dataPrepBLCA, geneInfo = geneInfo, method = "gcContent") #gcContent normalization
write.table (dataNormBLCA, file ="dataNormBLCA_noAMBRAlowvsNT.txt", sep=" ", quote = FALSE)
dataNormBLCA2 <- TCGAanalyze_Normalization (tabDF = dataNormBLCA, geneInfo = geneInfo, method = "geneLength") #geneLength normalization
write.table (dataNormBLCA, file ="dataNormBLCA2_noAMBRAlowvsNT.txt", sep=" ", quote = FALSE)
dataFiltBLCA <- TCGAanalyze_Filtering (tabDF = dataNormBLCA2, method = "quantile", qnt.cut = 0.25) #filter mRNA transcript selecting a threshold (quantile).
write.table (dataFiltBLCA, file="dataFiltBLCA_noAMBRAlowvsNT.txt", sep=" ", quote = FALSE )
#check that the gene of interest is in the dataFiltBLCA data set and if not then check at which step was removed

# to perform Differentially expression analysis (DEA), using edgeR package
dataDEGsFilt <- TCGAanalyze_DEA(mat1= dataFiltBLCA[,dataBLCAnt],   
                                mat2= dataFiltBLCA[,dataBLCA_noAMBRAlow],
                                Cond1type = "dataFilt - NT",
                                Cond2type = "dataFilt - No AMBRAlow",
                                fdr.cut = 0.01,
                                logFC.cut = 1,
                                method = "glmLRT")
write.table(dataDEGsFilt, file ="DEGs_dataFiltBLCA_noAMBRAlowvsNT.txt", sep = " ", quote=FALSE)
dataDEGsFiltUP <- dataDEGsFilt[dataDEGsFilt[,"logFC"]>=1, ]
dataDEGsFiltDOWN <- dataDEGsFilt[dataDEGsFilt[,"logFC"]<=-1, ]
write.table(dataDEGsFiltUP, file="DEGsUP_dataFiltBLCA_noAMBRAlowvsNT.txt",sep= " ", quote=FALSE)
write.table(dataDEGsFiltDOWN, file="DEGsDOWN_dataFiltBLCA_noAMBRAlowvsNT.txt", sep = " ", quote = FALSE)
dataDEGsNorm2 <- TCGAanalyze_DEA(mat1= dataNormBLCA2[,dataBLCAnt],
                                 mat2 = dataNormBLCA2[,dataBLCA_noAMBRAlow],
                                 Cond1type = "dataNorm2 - NT",
                                 Cond2type = "dataNorm2 - No AMBRA low",
                                 fdr.cut = 0.01,
                                 logFC.cut =1)
write.table(dataDEGsNorm2, file ="DEGs_dataNorm2BLCA_noAMBRAlowvsNT.txt", sep = " ", quote = FALSE)
dataDEGsNorm <- TCGAanalyze_DEA(mat1= dataNormBLCA[,dataBLCAnt],
                                mat2 = dataNormBLCA[,dataBLCA_noAMBRAlow],
                                Cond1type = "dataNorm - NT",
                                Cond2type = "dataNorm - No AMBRAlow",
                                fdr.cut = 0.01,
                                logFC.cut = 1)
write.table(dataDEGsNorm, file = "DEGs_dataNormBLCA_noAMBRAlowvsNT.txt", sep = " ", quote = FALSE)
#check if CHEK1 is in the list of DEGs and in which one.
