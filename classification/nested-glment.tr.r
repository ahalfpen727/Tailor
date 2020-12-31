


library(glmnet)
library(dplyr)


############################################################
### Nested  N-Fold function
########################################
nestedNfold = function(Xtrain, YLabels,alphaGiven, numberIter , nfold ,Output ,verbose = TRUE){
  
  NumberOfSamples = length(YLabels)
  pred = matrix(NA, nrow = NumberOfSamples, ncol = numberIter)
  LambdaMatrix = matrix(0, nrow = (nfold+1), ncol = 2*numberIter)
  
  sizeOfFolds = as.integer(floor(NumberOfSamples/(nfold+1)))
  IndexOfSplitOuter = c(seq.int(0L,by=sizeOfFolds,len=(nfold+1)),NumberOfSamples)
  OuterIndex = lapply(1:(nfold+1), function(x) list()) # predicted responses 		
  OuterIndex = lapply(1:numberIter, function(x) OuterIndex) ## repeate the list for all itterations
  
  
  labs = matrix(NA, nrow = NumberOfSamples, ncol = numberIter)
  NonZeroGenes = {}
  NonZeroCoefficient = {}
  ## generating CV folds.
  for(iter in 1:numberIter){
    if(verbose){
      cat('\n\n')
      cat(paste('running outer iteration', iter))
      cat('\n')
    }
    
    
    while(TRUE){
      Flag = TRUE
      ind = sample(1:NumberOfSamples)
      for (o in 1:(nfold+1))
      {
        IndexOfOuter = ind[(IndexOfSplitOuter[o]+1):IndexOfSplitOuter[o+1]]
        IndexOfInner = setdiff(ind,IndexOfOuter)
        if(sum(YLabels[IndexOfInner]==0) <= 3 || sum(YLabels[IndexOfInner]==1) <= 3){
          Flag = FALSE
          break
        }
      }
      
      if(Flag)
        break
    }
    for (o in 1:(nfold+1))
    {
      if (verbose)
        cat("\n*** OUTER NFOLD", o, "***")
      IndexOfOuter = ind[(IndexOfSplitOuter[o]+1):IndexOfSplitOuter[o+1]]
      IndexOfInner = setdiff(ind,IndexOfOuter)
      OuterIndex[[iter]][[o]] = list(IndexOfOuter)
      if (verbose)
        cat("\n*** Running CRE ***")		
      innerData = list(x = Xtrain[IndexOfInner,], y = YLabels[IndexOfInner])
      
      ## Cross-validation
      if (verbose)
        cat("\n*** Running Cross-Validation ***")
      ## looping through the inner folds ans selecting the best model (lambda)
      cv.fit = NULL
      while(is.null(cv.fit)){
        try(cv.fit=cv.glmnet(innerData$x, innerData$y, family="binomial", type.measure = "deviance", nfolds = nfold))
      }
      #cv.fit <- cv.glmnet(innerData$x, innerData$y, family="binomial", type.measure = "deviance", nfolds = nfold)
      lambdas = cv.fit$lambda
      lambdaMin = cv.fit$lambda.min
      bestLambdaIndex = which(lambdas == lambdaMin)
      LambdaMatrix[o, (2*iter-1):(2*iter)] = c(lambdaMin, bestLambdaIndex)
      cat(alphaGiven)
      ## fir on all folds 
      inner.fit = glmnet(innerData$x, innerData$y, family="binomial", alpha = alphaGiven,lambda = lambdas)
      
      Cutoff=sum(innerData$y)/length(innerData$y)
      ##x.outer = standardize(Xtrain[IndexOfOuter, ,drop = F])$x
      x.outer = Xtrain[IndexOfOuter, ,drop = F]
      pred[IndexOfOuter, iter] = predict(inner.fit, newx = x.outer, s = lambdaMin, type = "response", mode = "lambda")
      labs[IndexOfOuter, iter] = ifelse(pred[IndexOfOuter, iter] > Cutoff, 1, 0)
      
      
      ##
      NonZeroGenes = c(NonZeroGenes, colnames(Xtrain[,which(inner.fit$beta[, bestLambdaIndex] != 0)]))
      NonZeroCoefficient = c(NonZeroCoefficient, inner.fit$beta[, bestLambdaIndex][which(inner.fit$beta[, bestLambdaIndex] != 0)])
    }
    
    fileName=paste0("NonZeroGenes_Alpha_",alphaGiven,"_iteration_",iter)
    
    path=paste0(Output,alphaGiven,"/")
    dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    
    write.csv(NonZeroGenes, file=paste0(path,fileName,"_geneName.txt"))
    write.csv(NonZeroCoefficient, file=paste0(path,fileName,"_geneCoef.txt"))
 
    NonZeroGenes={}
    NonZeroCoefficient={}
  }
  
  L <- list(pred = pred, labs = labs, LambdaMatrix = LambdaMatrix, OuterIndex = OuterIndex)
  
  class(L) = c("nestedNfold")
  
  return(L)
  
} 
######################################################
############## nfold
######################################################
Nfold = function(Xtrain, YLabels,alphaGiven, numberIter , nfold ,Output ,verbose = TRUE){
  
  sampleSize=length(YLabels)
  
  predictions=matrix(0, ncol=numberIter, nrow=sampleSize)
  labs=matrix(0, ncol=numberIter, nrow=sampleSize)
  for(iterationIndex in 1:numberIter){
    
    #Randomly shuffle the data
    ShuffleIndex=sample(sampleSize)
    DataTable<-data[ShuffleIndex,]
    
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(DataTable)),breaks=nfold,labels=FALSE)
    
    #Perform nfold cross validation
    for(foldi in 1:nfold){
      
      #Segement your data by fold using the which() function 
      testIndexes <- which(folds==foldi,arr.ind=TRUE)
      testData <- DataTable[testIndexes, ]
      trainData <- DataTable[-testIndexes, ]
      
      NonZeroGenes = {}
      NonZeroCoefficient = {}
      
      #Find best lambda 
      cv.fit=cv.glmnet(trainData, YLabels[ShuffleIndex[-testIndexes]], family="binomial", type.measure = "deviance", nfolds = 4)
      LambdaMin= cv.fit$lambda.min
      lambdas = cv.fit$lambda
      best.lam.ind = which(lambdas == LambdaMin)
      
      model=glmnet(trainData,YLabels[ShuffleIndex[-testIndexes]], family="binomial", lambda = lambdas , alpha=alphaGiven)
      
      predictions[ShuffleIndex[testIndexes],iterationIndex] = c( predict(model, newx=as.matrix((testData)), s=LambdaMin ,type = "response", mode = "lambda"))
      
      Cutoff=(sum(YLabels[ShuffleIndex[-testIndexes]])/length(ShuffleIndex[-testIndexes]))
      
      labs[ShuffleIndex[testIndexes],iterationIndex]=c(ifelse( predictions[ShuffleIndex[testIndexes],iterationIndex] >= Cutoff, 1, 0))
      
      NonZeroGenes = c(NonZeroGenes, colnames(testData[,which(model$beta[, best.lam.ind] != 0)]))
      NonZeroCoefficient = c(NonZeroCoefficient, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])
    }  
    fileName=paste0("NonZeroGenes_Alpha_",alphaGiven,"_Iteration_",iterationIndex)
    
    path=paste0(Output,alphaGiven,"/")
    dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    
    write.csv(NonZeroGenes, file=paste0(path,fileName,"_geneName.txt"))
    write.csv(NonZeroCoefficient, file=paste0(path,fileName,"_geneCoef.txt"))
    
  }
  
  
  L=list(pred=predictions,labs=labs)
  class(L) = c("Nfold")
  return(L)
}

########################################
### Leave-one out function
########################################
LeaveOneOut = function(Xtrain, YLabels, alphaGiven,numberIter=1 , Output, verbose = TRUE){
  
  
  sampleSize=length(YLabels)
  
  predictions=rep(0,sampleSize)
  labs={}
  
  for(i in 1:sampleSize)
  
    {
    
    
    NonZeroGenes = {}
    NonZeroCoefficient = {}
    
    #Find best lambda 
    cv.fit=cv.glmnet(Xtrain[-i,], YLabels[-i], family="binomial", type.measure = "deviance", nfolds = 3)
    LambdaMin= cv.fit$lambda.min
    lambdas = cv.fit$lambda
    best.lam.ind = which(lambdas == LambdaMin)

    model=glmnet(Xtrain[-i,], YLabels[-i], family="binomial", lambda = lambdas , alpha=alphaGiven)
   
    predictions[i] = c( predict(model, newx=as.matrix(t(Xtrain[i,])), s=LambdaMin ,type = "response", mode = "lambda"))
    
      Cutoff=sum(YLabels[-i])/length(YLabels[-i])
    
    labs=c(labs,ifelse(predictions[i] > Cutoff, 1, 0))
    
    NonZeroGenes = c(NonZeroGenes, colnames(Xtrain[,which(model$beta[, best.lam.ind] != 0)]))
    NonZeroCoefficient = c(NonZeroCoefficient, model$beta[, best.lam.ind][which(model$beta[, best.lam.ind] != 0)])
    
    fileName=paste0("NonZeroGenes_Alpha_",alphaGiven,"_Sample_",i,"_out")
    
    path=paste0(Output,alphaGiven,"/")
    dir.create(path, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    
    write.csv(NonZeroGenes, file=paste0(path,fileName,"_geneName.txt"))
    write.csv(NonZeroCoefficient, file=paste0(path,fileName,"_geneCoef.txt"))
    
    
    
  }
  L=list(pred=predictions,labs=labs)
  class(L) = c("leave-oneOut")
  return(L)
} 

#####################################
#AUC

RocCurveFunction = function(specificity, Sensitivity, fileName){
  thresh=seq(0.01,1,0.01)
  xval=1-specificity
  yval=Sensitivity
  
  xydf= data.frame(xval, yval)
  xydf= xydf[complete.cases(xydf),]
  id = order(xydf$xval)  #find the id of sorted values
  #AUC inside a 1 by 1 square calculated by trapezoid method
  b = -diff(xydf$xval)  #get the consecutive differences
  my.auc = min(xydf$xval)*min(xydf$yval)/2 +  #area below the minimum y value
    sum(b*xydf$yval[-length(xydf$yval)])+(1-max(xydf$xval))*max(xydf$yval) +  #area of rectangles and triangles
    ((1-max(xydf$xval))*(1-max(xydf$yval))/2) #area after the max y value
  pdf(file =paste0(fileName,".pdf"))
  
  plot(yval~xval, xydf, type="l", xlab="1-Specificity", ylab="Sensitivity", main="ROC Curve", col="red")
  abline(0,1,lty=2)
  dev.off()
}
########################################
####################################################################################################
#main Iterator
####################################################################################################
#####################

# Gets arguments that were passed in via command line
args = commandArgs(TRUE)

for (i in 1:length(args)) {
  eval(parse(text=args[[i]]))
}
##set parameters

fold=Sys.getenv("FOLD") #-1 for leave-one out
method=Sys.getenv("METHOD")
nested=Sys.getenv("NESTED") #NO and YES
NumberOfIterations=Sys.getenv("ITERATIONS")
BootStrap=Sys.getenv("BOOTSTRAP_ITERATIONS")
TransformFunction=Sys.getenv("TRANSFORM")   # (NA, log, log2. log10, exp) 
maxBiomarkers=Sys.getenv("MAXBIOMARKERS")
scale=Sys.getenv("SCALE") # none, row, col
fpkmMatrix=Sys.getenv("FPKM_MATRIX")

##############################################
### read the data table
##############################################
data=read.table("/media/amir/DCF8-C228/ClassificationScripts/SamplebyGene.txt")
data=t(data)

if(TransformFunction=="log")
  data=log(data)
if(TransformFunction=="log2")
  data=log2(data)
if(TransformFunction=="log10")
  data=log10(data)
if(TransformFunction=="exp")
  data=exp(data)

# list of different alphas that the elastic net will be iterated on them.
ALphaList={}
if(method=="elasticNet"){
  ALphaList=seq(0,1,by=0.04)
}
if(method=="ridge")
{
  ALphaList=c(0)
}
if(method=="lasso")
{
  ALphaList=c(1)
  
}
resultMatrix=matrix(0, ncol=4, nrow=length(ALphaList))
Output="/home/amir/Desktop/Urine1/"

NumberOfFolds=fold
PositiveClass=rep(1,9)
NegativeClass=rep(0,9)
scale="row" #for heatmap
truth=c(PositiveClass,NegativeClass)

if(NumberOfFolds==-1)
{
  Output=paste0(Output,"Leave-OneOut/")
  dir.create(Output, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
if(NumberOfFolds>0){
  
  Output=paste0(Output,"nestedNfold/")
  
  dir.create(Output, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
ROCCurveMatrix=matrix(0,ncol=4,nrow=101)
###################################################
# call the function
###################################################
Index=1
for(i in ALphaList)
{
  #call the nested
  if(NumberOfFolds>0)
  {
    if(nested=="YES"){
         
         result=nestedNfold(data, truth, alphaGiven=i ,numberIter=NumberOfIterations, nfold=NumberOfFolds, Cutoff=CutOff, Output=Output)
         m=apply(result$pred,1,sum)
         m=m/NumberOfIterations
         
         #calculate the ROC Curve with different cutoff
         for(j in 0:100)
         {
           
           predictions=table(truth, ifelse((m>=(j/100)),1,0))  
           ROCCurveMatrix[j+1,]=predictions[1:4]
           
         }
         SensitivityVector=c()
         for(k in 1:101)
         {
           SensitivityVector=c(SensitivityVector,(ROCCurveMatrix[k,1])/(ROCCurveMatrix[k,1]+ROCCurveMatrix[k,3]))  
         }
         
         specificityVector=c()
         for(k in 1:101)
         {
           specificityVector=c(specificityVector,(ROCCurveMatrix[k,4])/(ROCCurveMatrix[k,4]+ROCCurveMatrix[k,2]))  
         }
         
         RocCurveFunction(specificity=specificityVector,Sensitivity=SensitivityVector, fileName = paste0(Output,"Alpha_",i,"_RocCurve"))
         
         # final Prediction with correct cutoff
         t=table(truth, ifelse((m>(CutOff)),1,0))
         resultMatrix[Index,]=t[1:4]
    }
    else{
      
      
      result=Nfold(data, truth, alphaGiven=i ,numberIter=NumberOfIterations, nfold=NumberOfFolds,  Output=Output)
      m=apply(result$pred,1,sum)
      m=m/NumberOfIterations
      
      #calculate the ROC Curve with different cutoff
      for(j in 0:100)
      {
        
        predictions=table(truth, ifelse((m>=(j/100)),1,0))  
        ROCCurveMatrix[j+1,]=predictions[1:4]
        
      }
      SensitivityVector=c()
      for(k in 1:101)
      {
        SensitivityVector=c(SensitivityVector,(ROCCurveMatrix[k,1])/(ROCCurveMatrix[k,1]+ROCCurveMatrix[k,3]))  
      }
      
      specificityVector=c()
      for(k in 1:101)
      {
        specificityVector=c(specificityVector,(ROCCurveMatrix[k,4])/(ROCCurveMatrix[k,4]+ROCCurveMatrix[k,2]))  
      }
      
      RocCurveFunction(specificity=specificityVector,Sensitivity=SensitivityVector, fileName = paste0(Output,"Alpha_",i,"_RocCurve"))
      
      # final Prediction with correct cutoff
      t=table(truth, ifelse((m>(CutOff)),1,0))
      resultMatrix[Index,]=t[1:4]
      
    }
  }
  #call leave-one out
  else{
   
     result=LeaveOneOut(data,truth,alphaGiven=i ,numberIter=1,Output=Output)
     NumberOfIterations=1
    
    #calculate the ROC Curve with different cutoff
    for(j in 0:100)
    {
      
      predictions=table(truth, ifelse((result$pred>=(j/100)),1,0))  
      ROCCurveMatrix[j+1,]=predictions[1:4]
      
    }
    
    SensitivityVector=c()
    for(k in 1:101)
    {
      SensitivityVector=c(SensitivityVector,(ROCCurveMatrix[k,1])/(ROCCurveMatrix[k,1]+ROCCurveMatrix[k,3]))  
    }
    
    specificityVector=c()
    for(k in 1:101)
    {
      specificityVector=c(specificityVector,(ROCCurveMatrix[k,4])/(ROCCurveMatrix[k,4]+ROCCurveMatrix[k,2]))  
    }
    
    RocCurveFunction(specificity=specificityVector,Sensitivity=SensitivityVector, fileName = paste0(Output,"Alpha_",i,"_RocCurve"))
    
    t=table(truth, result$labs)
    resultMatrix[Index,]=t[1:4]
  }
  
 
  Index=Index+1
}


######################################################
#Calculate the Accuracy
############################################3

FinalResultOfAllAccMeasures=matrix(NA, ncol=length(ALphaList),nrow=7)
colnames(FinalResultOfAllAccMeasures)=c(as.character( ALphaList))
rownames(FinalResultOfAllAccMeasures)=c("Accuracy","Sensitivity","Specificity","FDR","FPR","F1","MCC")


Accuracy=c()

for(i in 1:length( ALphaList))
{
  Accuracy=c(Accuracy,(resultMatrix[i,1]+resultMatrix[i,4])/sum(resultMatrix[i,]))  
}
FinalResultOfAllAccMeasures[1,]=Accuracy

Sensitivity=c()
for(i in 1:length( ALphaList))
{
  Sensitivity=c(Sensitivity,(resultMatrix[i,1])/(resultMatrix[i,1]+resultMatrix[i,3]))  
}
FinalResultOfAllAccMeasures[2,]=Sensitivity

specificity=c()
for(i in 1:length( ALphaList))
{
  specificity=c(specificity,(resultMatrix[i,4])/(resultMatrix[i,4]+resultMatrix[i,2]))  
}
FinalResultOfAllAccMeasures[3,]=specificity

FDR=c()
for(i in 1:length( ALphaList))
{
  FDR=c(FDR,(resultMatrix[i,2])/(resultMatrix[i,1]+resultMatrix[i,2]))  
}

FinalResultOfAllAccMeasures[4,]=FDR

FPR=c()
for(i in 1:length( ALphaList))
{
  FPR=c(FPR,(resultMatrix[i,2])/(resultMatrix[i,4]+resultMatrix[i,2]))  
}
FinalResultOfAllAccMeasures[5,]=FPR

F1=c()
for(i in 1:length( ALphaList))
{
  F1=c(F1,(2*resultMatrix[i,1])/(2*resultMatrix[i,1]+resultMatrix[i,2]+resultMatrix[i,3]))  
}
FinalResultOfAllAccMeasures[6,]=F1

MCC=c()
for(i in 1:length( ALphaList))
{
  MCC=c(MCC,(resultMatrix[i,4]*resultMatrix[i,1]-resultMatrix[i,2]*resultMatrix[i,3])/sqrt((resultMatrix[i,1]+resultMatrix[i,2])*
                                                                                             (resultMatrix[i,1]+resultMatrix[i,3])*
                                                                                             (resultMatrix[i,2]+resultMatrix[i,4])*
                                                                                             (resultMatrix[i,4]+resultMatrix[i,3])))  
}

FinalResultOfAllAccMeasures[7,]=MCC

write.table(FinalResultOfAllAccMeasures, file=paste0(Output,"AllAccuracyMeasures.txt"))

######################################
#Select biomarker genes and generate heatmap for them
##########################################
  for(alpha in ALphaList){
    
    
    # Read the non-zero genes and count them.
    dir=paste0(Output,alpha)
    setwd(dir)
    listofFile=list.files(path=dir, pattern = "_geneName.txt")
    nonZeroGenes={}
    for(file in listofFile)
    {
      tryCatch({
            nonZeroGenes=rbind(nonZeroGenes,read.csv(file))
        
      },error=function(e){})
    }
    #top picked genes
    GroupedByGeneName = group_by(nonZeroGenes,nonZeroGenes$x)
    (CountPerGenes   <- summarise(GroupedByGeneName, nonZeroGenes = n()))
    colnames(CountPerGenes)=c("GeneName","Count")
    o=order(CountPerGenes$Count,decreasing =T)

    markerSize=ifelse(dim(CountPerGenes)[1]<maxBiomarkers,dim(CountPerGenes)[1],maxBiomarkers)  

    MarkerGenes=CountPerGenes[o[1:markerSize],]$GeneName
    MarkerMatrix=subset(data,select=as.character( MarkerGenes))
    filename=paste0("Alpha_",alpha)
    pdf(file =paste0(filename,".pdf"))
    heatmap(t(MarkerMatrix), Rowv=NA, Colv=NA, col = cm.colors(10), scale=scale, margins=c(5,10),na.rm = T)
    
    dev.off()
    write.table(CountPerGenes[o,], file="GroupedByGenes.txt", row.names = F)
  
  }

