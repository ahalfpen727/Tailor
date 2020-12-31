
library(cummeRbund)

setwd("/project/umb_triley/cpct/rna-seq/urine1/cuffdiff_results_GRCh38_default/LUTS-over-CTRL")
#setwd("D://LUTS-over-CTRL")
## Read the Cuffdiff Result and map it to sample per genes, using cummerbnund
cuff<-readCufflinks(dir=getwd(), rebuild=F)
cuff

gene_diff_data <- diffData(genes(cuff))
sig_gene_data <- subset(gene_diff_data, significant=="yes")
head(sig_gene_data)


gene.repFpkm<-repFpkm(genes(cuff))
head(gene.repFpkm)




###################################################
### code chunk number 8: Gene symbol by sample
###################################################
gene.rep.matrix<-repFpkmMatrix(genes(cuff))
head(gene.rep.matrix)


############# mapping the gene Symbol into Entrez ID
library(org.Hs.eg.db)
library(netClass)

sampleByGene=t(gene.rep.matrix)

geneName=colnames(sampleByGene)

EntrezID=geneName

#getting the binding data
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

for(i in 1:length(geneName)){

    tmp=EN2SY[geneName[i]]
    if(length(tmp[[1]])>0)
    EntrezID[i]=tmp[[1]]
    else
      EntrezID[i]=NA  ## gene symbols without any Entrez ID
}

colnames(sampleByGene)=EntrezID




## removing those genes without ENtrez ID

df=sampleByGene
colnames(df)=EntrezID

#temprory Matrix to remove the column without Entrez ID
ColnameMatrix=matrix(,nrow=length(EntrezID)- sum(is.na(EntrezID)),ncol=1)
l=0
p=0
for(i in 1:length(EntrezID))
{
  if(is.na(colnames(df)[i]))
   l=l+1
  else
   {p=p+1
     ColnameMatrix[p,1]=colnames(df)[i]
   }
}
# final Matrix that has fpkm of each gene per sample
AnotherMatrix=matrix(,nrow=18,ncol=length(EntrezID)- sum(is.na(EntrezID)))
k=0
for(i in 1:33117)
{
  if(is.na(colnames(df)[i]))
  {
    
  }
  else
  {
   k=k+1
   AnotherMatrix[,k]=df[,i]
  }
  
}
colnames(AnotherMatrix)=ColnameMatrix[,1]
write.table(t(AnotherMatrix),file="removedEntrez.txt")
rownames(AnotherMatrix)=rownames(sampleByGene)



load("a.ppi.rda")

AdjColName=colnames(A.ppi)# a vector of name of adj matrix
EntrezIDColName=colnames(AnotherMatrix) # a vector of name of genes

## Trying to extend the matrix by finding the diferent genes.

NotinNewAdj=setdiff(AdjColName,EntrezIDColName)
SelectedGenes=setdiff(AdjColName, NotinNewAdj)

#ValuesofGenes=subset(AnotherMatrix, select=SelectedGenes)

## New Adj Matrix

AppendGenes=setdiff(EntrezIDColName,SelectedGenes)


nn=subset(A.ppi, select=SelectedGenes)

kk=subset(t(nn), select=SelectedGenes)

for(i in 1:length(AppendGenes))
{
  kk=cbind(kk,0)
  kk=rbind(kk,0)
}
colnames(kk)=c(SelectedGenes,AppendGenes)
rownames(kk)=c(SelectedGenes,AppendGenes)

## Bootstrapping the Matrix
numberOfRows=nrow(ValuesofGenes)

for(i in 1:numberOfRows)
{
BootstrappedData=sample(ValuesofGenes[i,],replace = TRUE)
ValuesofGenes=rbind(ValuesofGenes,BootstrappedData)
}
#Assign the row names for those samples which recently added
rownames(ValuesofGenes)=c("CTRL_0" ,	"CTRL_1" ,	"CTRL_2"	, "CTRL_3"	, "CTRL_4" ,	"CTRL_5"	, "CTRL_6" ,	"CTRL_7","CTRL_8"	 ,	"LUTS_0" ,	"LUTS_1" ,	"LUTS_2" ,	"LUTS_3"	 , "LUTS_4" ,	"LUTS_5" ,	"LUTS_6" ,	"LUTS_7"	,"LUTS_8","CTRL_9" ,	"CTRL_10" ,	"CTRL_11" ,	"CTRL_12"	, "CTRL_13"	, "CTRL_14" ,	"CTRL_15"	, "CTRL_16" ,	"CTRL_17"	, "LUTS_9" ,	"LUTS_10" ,	"LUTS_11" ,	"LUTS_12" ,	"LUTS_13"	 , "LUTS_14" ,	"LUTS_15" ,	"LUTS_16" ,	"LUTS_17"	, "CTRL_18" ,	"CTRL_19"	, "CTRL_20" ,	"CTRL_21"	, "CTRL_22" ,	"CTRL_23"	,"CTRL_24" ,	"CTRL_25" ,	"CTRL_26"	, "LUTS_18" ,	"LUTS_19" ,	"LUTS_20" ,	"LUTS_21" ,	"LUTS_22" ,	"LUTS_23"	 ,"LUTS_24"	, "LUTS_25"	, "LUTS_26"	, "CTRL_27"	, "CTRL_28" ,	"CTRL_29" ,	"CTRL_30" ,	"CTRL_31"	, "CTRL_32" , "CTRL_33"	, "CTRL_34"	, "CTRL_35" , "LUTS_27" ,	"LUTS_28" ,	"LUTS_29" ,	"LUTS_30" ,	"LUTS_31" ,	"LUTS_32" ,	"LUTS_33" , "LUTS_34" , "LUTS_35")

save(ValuesofGenes, file="geneData.rda")
save(kk,file="Adj.rda")
