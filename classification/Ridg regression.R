library(glmnet)

#Files directory
setwd("/project/umb_triley/cpct/rna-seq/urine1/ClassificationScripts")
filenames <- list.files("/project/umb_triley/cpct/rna-seq/urine1/leave_out_sig_features", pattern="*.txt", full.names=T)
# two vectors to ma
indexofFiles=c(10,11,12,13,14,15,16,17,18,19,20,3,4,5,6,7,8,9)
indexofsample=c(13,5,14,6,15,7,16,8,17,9,18,1,10,2,11,3,12,4)


#load the test sample list
AllGenesBuiltModel=read.table( file="genebysample.txt")
sample.total=dim(t(AllGenesBuiltModel))[1]

#sample label
y=c(rep(1,9),rep(0,9))

predictions=rep(0,sample.total)

nonzero.genes = {}
nonzero.coeffs = {}
list.model={}
#ListofModels=c()
for(i in 1:sample.total)
{
  
  #perprocessing the data
  dataTable=read.table(filenames[i])
  dataTable=t(dataTable)
  
  tmp=dataTable
  colnames(tmp)= dataTable[2,]
  dataTable=dataTable[3:dim(dataTable)[1],]
  k = dim(dataTable)[1]
  dataMatrix=matrix(,nrow=k, ncol=dim(dataTable)[2])
  for(j in 1: k){
    dataMatrix[j,]=as.numeric(dataTable[j,])
  }
  #Find best lambda by training the data
  cv.fit<-cv.glmnet(dataMatrix, y[-indexofsample[i]], family="binomial", type.measure = "deviance", nfolds = 10)
  best.lambda= cv.fit$lambda.min
  
  lambdas <- cv.fit$lambda
  best.lam.ind <- which(lambdas == best.lambda)
  
  # Call elastic network for test set
  model=glmnet(dataTable, y[-indexofsample[i]], family="binomial", lambda = cv.fit$lambda , alpha=0.02)
  
  list.model=c(list.model,model)
  predictions[i] = c( predict(model, newx=t(AllGenesBuiltModel[,indexofsample[i]]), s=best.lambda ,type = "response", mode = "lambda"))
  nonzero.genes = c(nonzero.genes, colnames(tmp[,which(model$beta[, best.lam.ind] != 0)]))
  nonzero.coeffs = c(nonzero.coeffs, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])
  
}
# plot the performance
plot(predictions,ylab="Performance",main=NULL, col=(c(rep(c("red","blue"),9))))
legend(14.5,.15, c("LUTS","CTRL"),col=c("red","blue"),lty = c(1,1),lwd=c(2.5,2.5))
abline(a=0.5,b=0)

#confustion matrix
truth=c(rep(0:1,9))
table(truth,ifelse((predictions)>0.5,1,0))

library(hash)
#size of nonZero genes
lengthOfNonZeroGenes=length(unique(nonzero.genes))

vote=hash(keys=unique(nonzero.genes), values=1:lengthOfNonZeroGenes)
frequencyLength=c()
for(i in 1:lengthOfNonZeroGenes)
{
  
  index=grep(unique(nonzero.genes)[i], names(nonzero.coeffs))  
  frequencyLength=c(frequencyLength,length(index))
  vote[unique(nonzero.genes)[i]]=nonzero.coeffs[index]
}
names(frequencyLength)=unique(nonzero.genes)
o=order(frequencyLength,decreasing = T)
selectedName=names(frequencyLength[o])


o=nonzero.coeffs>1
sigGenes= unique( names(nonzero.coeffs[o]))
vote[sigGenes[3]]
mean(nonzero.coeffs)
sd(nonzero.coeffs)


####Script 2
#reading saved model

geneName=read.table("noneZerogeneName.txt")
geneCoeficient=read.table("noneZeroCoef.txt")

geneSymbolList=c()
for(i in 1:length(unique(nonzero.genes)))
{
  geneSymbolList=c(geneSymbolList,dataTable[unique(nonzero.genes)[i],2])  
}

row.names(dataTable)[1]
m=intersect(unlist(sig), unlist(dataTable[geneSymbolList,2]))
