library(glmnet)


# Get file name of the input data


filenamesForSigGenes <- list.files("/project/umb_triley/cpct/rna-seq/urine1/leave_out_sig_features/SigGenes", pattern="*.txt", full.names=T)
filenames <- list.files("/project/umb_triley/cpct/rna-seq/urine1/leave_out_sig_features", pattern="*.txt", full.names=T)

#Mapping vector between sample id and LUTS/CTRL id
indexofFiles=c(10,11,12,13,14,15,16,17,18,19,20,3,4,5,6,7,8,9)
indexofsample=c(13,5,14,6,15,7,16,8,17,9,18,1,10,2,11,3,12,4)

#List of Sig Genes between all 17 models
ListOfSigGenes={}

#we have 17 sample in total
# collect all the gene Names
for(i in 1:17)
{
  
  dataTable=read.table(filenames[i])
  ListOfSigGenes=union(dataTable[[2]],ListOfSigGenes )
}


#load all the data for the 
dgene.rep.matrix=read.table( file="genebysample.txt")
sample.total=dim(t(dgene.rep.matrix))[1]
x=dgene.rep.matrix

y=c(rep(1,9),rep(0,9))

predictions=rep(0,sample.total)

nonzero.genes = {}
nonzero.coeffs = {}
list.model={}
ListofModels=c()
for(i in 1:sample.total)
{
  
  dataTable=read.table(filenames[i])
  dataTable=t(dataTable)
  tempDf=dataTable
  
  colnames(dataTable)=dataTable[2,]
  
  dataTable=dataTable[3:dim(dataTable)[1],]
  
  dataTableSigGenes=read.table(filenamesForSigGenes[i])
  geneNames=dataTableSigGenes[[2]]
  GeneSymbol=dataTableSigGenes[[1]]
  tempDf=subset(dataTable[1:2,], select =geneNames )
  
  dataTable=subset(dataTable, select =GeneSymbol )
  k = dim(dataTable)[1]
  dataMatrix=matrix(,nrow=k, ncol=dim(dataTable)[2])
  for(j in 1: k){
    dataMatrix[j,]=as.numeric(dataTable[j,])
  }
  #Find best lambda 
  cv.fit<-cv.glmnet(dataMatrix, y[-indexofsample[i]], family="binomial", type.measure = "deviance", nfolds = 10)
  best.lambda= cv.fit$lambda.min
  
  lambdas <- cv.fit$lambda
  best.lam.ind <- which(lambdas == best.lambda)
  
  
  model=glmnet(dataTable, y[-indexofsample[i]], family="binomial", lambda = cv.fit$lambda , alpha=0.70)
  list.model=c(list.model,model)
  
  
  
  
  
  SelectedGenes=subset(as.data.frame(t( x[,indexofsample[i]])),select =  geneNames )
  colnames(SelectedGenes)=tempDf[2,]
  colnames(SelectedGenes)=GeneSymbol
  predictions[i] = c( predict(model, newx=as.matrix( SelectedGenes), s=best.lambda ,type = "response", mode = "lambda"))
  
  
  
  nonzero.genes = c(nonzero.genes, colnames(dataTable[,which(model$beta[, best.lam.ind] != 0)]))
  nonzero.coeffs = c(nonzero.coeffs, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])
  
}  

truth=c(rep(0:1,9))

table(truth,ifelse((predictions)>0.5,1,0))

# plot the coeffients 
plot(predictions,ylab="Performance",main="Sig Genes classification performance", col=(c(rep(c("red","blue"),9))))
legend("bottomright", c("LUTS","CTRL"),col=c("red","blue"),lty = c(1,1),lwd=c(2.5,2.5))
abline(a=0.5,b=0)



# Counting dictionary

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

hist(frequencyLength, breaks = 18,col="red", main="Histogram of sig-genes based on seen on each class")
o=nonzero.coeffs>1
sigGenes= unique( names(nonzero.coeffs[o]))

