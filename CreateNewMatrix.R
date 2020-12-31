
library(cummeRbund)

setwd("/project/umb_triley/cpct/rna-seq/urine1/cuffdiff_results_GRCh38_gtf_only/LUTS-over-CTRL")
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

library(org.Hs.eg.db)
sampleByGene=t(gene.rep.matrix)

geneName=colnames(sampleByGene)

EntrezID=geneName

#getting the binding data
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]

for(i in 1:length(geneName)){

    tmp=xx[geneName[i]]
    if(length(tmp[[1]])>0)
    EntrezID[i]=tmp[[1]]
    else
      EntrezID[i]=NA
}

colnames(sampleByGene)=EntrezID



df=sampleByGene
colnames(df)=EntrezID

ColnameMatrix=matrix(,nrow=27466,ncol=1)
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
AnotherMatrix=matrix(,nrow=18,ncol=27466)
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
library(netClass)

load("a.ppi.rda")
AdjColName=colnames(A.ppi)
EntrezIDColName=colnames(AnotherMatrix)
DifferentGenes=setdiff(EntrezIDColName,AdjColName)
SameGenesbyNetwork=intersect(AdjColName,EntrezIDColName)
AddedGenes=setdiff(AdjColName,SameGenesbyNetwork)
tmp2=union(AdjColName,EntrezIDColName)

NewAdjMatrix=matrix(data=0,ncol=length(tmp2),nrow=length(tmp2))
colnames(NewAdjMatrix)=c(tmp2)
rownames(NewAdjMatrix)=c(tmp2)
for(i in 1:dim(A.ppi)[1]){
  for(j in 1:dim(A.ppi)[1])
  {
      if(A.ppi[i,j]>0)
      {
        NewAdjMatrix[AdjColName[i],AdjColName[j]]=A.ppi[i,j]
      }
  }

}
save(NewAdjMatrix,file="NewAdj.RData")