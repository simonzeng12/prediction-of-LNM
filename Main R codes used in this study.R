library(dplyr)
library(affy)
library(sva)
library(xlsx)
library(rms)
library(pROC)
library(glmnet)
library(ggplot2)
library(PresenceAbsence)
library(rmda)
### Gene annotation 1
ann1=data.table::fread("GPL6947-13512.txt")
ann1=ann1[,c(1,14)]
ann1=ann1[nchar(ann1$Symbol)>0,] 
## Gene annotation 2
ann2=data.table::fread("GPL570-55999.txt",header = T)
ann2=ann2[,c(1,11)];colnames(ann2)=c("ID","Symbol")
ann2=ann2 %>% tidyr::separate(Symbol,into=c("Symbol","junk"),sep='///');ann2=ann2[,-3]
ann2=ann2[nchar(ann2$Symbol)>0,] 


# Raw data (training set)
rt=read.table("GSE84437_series_matrix.txt.gz",header = T,comment.char = "!")
rt = rt %>% tibble::column_to_rownames("ID_REF")
rt=log2(rt)
# Remove duplicated genes
duplicateANN1 = function(x){
  x = as.data.frame(x) %>% tibble::rownames_to_column("ID_REF")
  x=dplyr::inner_join(ann1,x,by=c("ID"="ID_REF"))
  x$median=apply(x[,-c(1:2)],1,median) 
  table(duplicated(x$Symbol))
  x=x[order(x$median,decreasing = T),]
  x=x[!duplicated(x$Symbol),]
  x=subset(x,select = -c(ID,median))
  x=as.data.frame(x)
}

rt=duplicateANN1(rt) 

# Clinical information
pheno=read.table("pheno.txt",sep="\t",header = T)
# Raw data (testing set 1)
rawdata_GSE62254 <- ReadAffy(celfile.path  = "download_GSE62254")
eset1.rma <- rma(rawdata_GSE62254)
# Raw data (testing set 2)
rawdata_GSE57303 <- ReadAffy(celfile.path  = "download_GSE57303")
eset2.rma <- rma(rawdata_GSE57303)
# column names modified
e1=as.data.frame(exprs(eset1.rma))
colnames(e1)=substr(colnames(e1),1,10)
e2=as.data.frame(exprs(eset2.rma))
colnames(e2)=substr(colnames(e2),1,10)

# Remove duplicated genes
duplicateANN2 = function(x){
  x = as.data.frame(x) %>% tibble::rownames_to_column("ID_REF")
  x=dplyr::inner_join(ann2,x,by=c("ID"="ID_REF"))
  x$median=apply(x[,-c(1:2)],1,median) 
  table(duplicated(x$Symbol))
  x=x[order(x$median,decreasing = T),]
  x=x[!duplicated(x$Symbol),]
  x=subset(x,select = -c(ID,median))
  x=as.data.frame(x)
}
e1=duplicateANN2(e1)
e2=duplicateANN2(e2) 

# Remove batch effects
edata=dplyr::inner_join(rt,e1,by=c("Symbol")) %>% inner_join(y=e2,by=c("Symbol"))
edata = edata %>% tibble::column_to_rownames("Symbol") %>% as.matrix()
combat_edata=ComBat(dat = edata,batch = c(rep(1,433),rep(2,300),rep(3,70)),par.prior = T)


# Training set
pheno.train=pheno %>% filter(batch =="GSE84437") 
pheno.train = pheno.train %>% mutate(N=ifelse(pheno.train$N=="negative",0,1)) 
trainSet = as.data.frame(t(combat_edata)) %>% tibble::rownames_to_column("sample") %>%
  inner_join(x=pheno.train,by=c("sample"="sample")) %>%
  tibble::column_to_rownames("sample")

# testing set 1
pheno.a = pheno %>% filter(batch =="GSE62254") 
pheno.a = pheno.a %>% mutate(N=ifelse(pheno.a$N=="negative",0,1)) 
test.a = as.data.frame(t(combat_edata)) %>% tibble::rownames_to_column("sample") %>%
  inner_join(x=pheno.a,by=c("sample"="sample")) %>%
  tibble::column_to_rownames("sample")

# testing set 2
pheno.b = pheno %>% filter(batch =="GSE57303")
pheno.b = pheno.b %>% mutate(N=ifelse(pheno.b$N=="negative",0,1)) 
test.b = as.data.frame(t(combat_edata)) %>% tibble::rownames_to_column("sample") %>%
  inner_join(x=pheno.b,by=c("sample"="sample")) %>%
  tibble::column_to_rownames("sample")

### Logistic regression
IMM=read.xlsx("ImmuneGenes.xls",sheetName = "Geneappend3")
IMM=as.character(IMM$Symbol)
logisticGene=function(df=null){
  outTab=data.frame()
  for(i in colnames(df[2:ncol(df)])){
    model <- glm(N ~ df[,i], data = df,family=binomial)
    modelSummary = summary(model)
    logisticP=try(modelSummary$coefficients[2,"Pr(>|z|)"])
    if("try-error"%in% class(logisticP)){
      next
    } else
      
      
      if(logisticP<0.05) {
        
        outTab=rbind(outTab,
                     cbind(id=i,
                           Coef=modelSummary$coefficients[2,"Estimate"],
                           OR=exp(coef(model))[2],
                           OR.95L=exp(confint(model))[2,"2.5 %"],
                           OR.95H=exp(confint(model))[2,"97.5 %"],
                           pvalue=logisticP)
        )
      }
    
  }
  outTab=as.data.frame(outTab)
  outTab$id = as.character(outTab$id)
  rownames(outTab)=outTab$id
  outTab[,2:6]=apply(outTab[,2:6],2,as.numeric)
  return(outTab)
}  
uni_output=logisticGene(df=trainSet[,c("N",intersect(IMM,colnames(trainSet)))])

### Lasso regression
set.seed(100001)
x=as.matrix(trainSet[,rownames(uni_output)]) 
y=trainSet$N
fit.lasso <- glmnet(x, y, family = "binomial", maxit = 1000)
plot(fit.lasso, xvar = "lambda", label = TRUE)
cvfit <- cv.glmnet(x, y, family="binomial", maxit = 1000)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
coef <- coef(fit.lasso, s = cvfit$lambda.min)
coef=as.matrix(coef)
index <- which(coef != 0)
lasso_output = rownames(coef)[index][-1]

# stepwise regression
formula = formula(paste(paste("N","~",collapse = " "),paste(lasso_output,collapse = "+")))
multiLogistic=glm(formula, data=trainSet,family=binomial)
multiLogistic=step(multiLogistic,direction = "both")
multiLogSum=summary(multiLogistic)
multiLogSum


calPlot=function(label=NULL,prediction=NULL,N.bins=NULL,lim=NULL,xlab=NULL){
  par(mar=c(5.1 ,4.1, 4.1 ,2.1),mgp=c(2.5, 1, 0))
  calibration_data = data.frame(name = 0, label, prediction) 
  plot_data=calibration.plot(calibration_data, which.model = 1, na.rm = TRUE, N.bins =N.bins)
  plot(x=plot_data$BinPred,y=plot_data$BinObs,xlim = lim,ylim = lim,type = "b",col="blue",
       pch=16,lwd=2,bty="o",cex.lab=2,cex.axis=1.2,xlab =xlab,ylab="Observed Risk",main = "Calibration curve",cex.main=2) ;abline(0,1);box(lwd = 2)
  
}    
plotroc=function(label=NULL, prediction=NULL){
  roc <- roc(label, prediction)
  plot(roc,type="l",lwd = 2.5,cex.lab=2,cex.axis=1.2, col="red",main="ROC curve",cex.main=2)
  text(0.2,0.1, paste0("AUC = ",sprintf("%0.3f",roc$auc)),adj = 0.5,cex = 3)
}

formula = formula(paste(paste("N","~",collapse = " "),paste(c("IRF3",	"IL12A",	"LRSAM1",	"NOV",	"BMPR1A",	"RXRB"),collapse = "+")))
fit=lrm(formula, data=trainSet,x=TRUE,y=TRUE)
fit

# ROC, calibration plots and decision analysis curve (training set)
trainSet$prob <- predict(fit,type = "fitted")
trainSet$lp <- predict(fit)
plotroc(label=trainSet$N, prediction=trainSet$prob)
calPlot(label=trainSet$N, prediction=trainSet$prob,xlab="Predicted Risk",lim=c(0,1),N.bins = 4)
dca.train=decision_curve(formula=N~lp, data=trainSet,family = "binomial",confidence.intervals = F,bootstraps = 10,fitted.risk = F)
plot_decision_curve(dca.train,curve.names = "Immune signature")

# ROC, calibration plots and decision analysis curve (testing set 1)
test.a$prob <- predict(fit,type = "fitted",newdata=test.a)
test.a$lp <- predict(fit,newdata=test.a)
plotroc(label=test.a$N, prediction=test.a$prob)
calPlot(label=test.a$N, prediction=test.a$prob,xlab="Predicted Risk",lim=c(0,1),N.bins = 4)
dca.a=decision_curve(formula=formula, data=test.a,family = "binomial",confidence.intervals = F,bootstraps = 10,fitted.risk = F)
plot_decision_curve(dca.a,curve.names = "Immune signature")

# ROC, calibration plots and decision analysis curve (testing set 2)
test.b$prob <- predict(fit,type = "fitted",newdata=test.b)
test.b$lp <- predict(fit,newdata=test.b)
plotroc(label=test.b$N, prediction=test.b$prob)
calPlot(label=test.b$N, prediction=test.b$prob,xlab="Predicted Risk",lim=c(0.3,1),N.bins = 6)
dca.b=decision_curve(formula=formula, data=test.b,family = "binomial",confidence.intervals = F,bootstraps = 10,fitted.risk = F)
plot_decision_curve(dca.b,curve.names = "Immune signature")
