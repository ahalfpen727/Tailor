library(glmnet)

dataTable=read.table("C:/Users/amir/Desktop/netClass/man/leave_out_sig_features/RepFpkmMatrix_leaveout3.genes.txt")
fullPath="C:/Users/amir/Desktop/netClass/man/leave_out_sig_features"
filenames <- list.files("C:/Users/amir/Desktop/netClass/man/leave_out_sig_features", pattern="*.txt", full.names=T)

indexofFiles=c(10,11,12,13,14,15,16,17,18,19,20,3,4,5,6,7,8,9)
indexofsample=c(13,5,14,6,15,7,16,8,17,9,18,1,10,2,11,3,12,4)

#load all the data for the 
gene.rep.matrix=read.table( file="SamplebyGene.txt")
sample.total=dim(t(gene.rep.matrix))[1]
x=gene.rep.matrix

y=c(rep(1,9),rep(0,9))

predictions=rep(0,sample.total)

nonzero.genes = {}
nonzero.coeffs = {}
list.model={}
#ListofModels=c()
for(i in 1:sample.total)
{
  
  dataTable=read.table(filenames[i])
  dataTable=t(dataTable)
  dataTable=dataTable[3:dim(dataTable)[1],]
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
  
  
  model=glmnet(dataTable, y[-indexofsample[i]], family="binomial", lambda = cv.fit$lambda , alpha=0.02)
  list.model=c(list.model,model)
  predictions[i] = c( predict(model, newx=t(x[,indexofsample[i]]), s=best.lambda ,type = "response", mode = "lambda"))
  
  nonzero.genes = c(nonzero.genes, colnames(dataTable[,which(model$beta[, best.lam.ind] != 0)]))
  nonzero.coeffs = c(nonzero.coeffs, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])
  
}
# plot the coeffients 
plot(predictions,ylab="Performance",main=NULL, col=(c(rep(c("red","blue"),9))))
legend(14.5,.15, c("LUTS","CTRL"),col=c("red","blue"),lty = c(1,1),lwd=c(2.5,2.5))
abline(a=0.5,b=0)

#confustion matrix
truth=c(rep(0:1,9))

table(truth,ifelse((predictions)>0.5,1,0))

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
  names(frequencyLength)=unique(nonzero.genes)[i]
o=nonzero.coeffs>1
sigGenes= unique( names(nonzero.coeffs[o]))
vote[sigGenes[3]]
mean(nonzero.coeffs)
sd(nonzero.coeffs)
