#######set directory
setwd('D:\\Peptide prediction\\Tumor homing peptides\\Data set')
#######Load package
library(RWeka)
library(caret)
library(randomForest)
library(pls)
library(Interpol)
library(protr)
library(seqinr)
library(Peptides)
library(reshape)
################## customRF
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

################## Kennard stone 
AAindex_main = c("CHAM830101","ANDN920101","CHAM820102","BUNA790103","BROC820102","CHAM830102","BURA740101","BEGF750102","CHAM810101","ARGP820103","BURA740102","BUNA790102","FASG890101","BROC820101",
"BUNA790101","BULH740102","BULH740101","CHAM820101","BIOV880102","BIOV880101","ARGP820102","BHAR880101","BEGF750103","BEGF750101","ARGP820101","BIGC670101","CHOC760103","CHAM830107","CHAM830104",
"CHOP780204")

AAindex_main90 = c("NAKH900101","AURR980120","ONEK900102","AURR980119","GUYH850101","ONEK900101","NAKH900102","AURR980118","NAGK730102","NAGK730103",
"HUTJ700102","LAWE840101","KRIW790101","HUTJ700101","HOPT810101","LEVM760104","KYTJ820101","HUTJ700103","BURA740101","KRIW790103",
"HOPA770101","CHAM830101","ISOY800101","NAKH900105","KLEP840101","CHAM820101","ISOY800103","BUNA790102","KANM800103","KRIW790102")

AAindex_small = c("ANDN920101","CIDH920103","BUNA790103","CIDH920101","FASG890101","CHOP780206","BIOV880101","BIOV880102","CHOP780214","BHAR880101",
"ARGP820103","COHE430101","BULH740101","CIDH920104","BEGF750102","CHOP780216","CHAM830102","BULH740102","CHOC750101","CHOP780212",
"CIDH920105","CIDH920102","BIGC670101","ARGP820102","CHAM830104")


para1 = 3
para2 = 0.6
#######Read data
A <- read.fasta('small.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("label small.csv", header = TRUE) 
m = length(A)
aac<- t(sapply(A, extractAAC))
dpc <- t(sapply(A, extractDC))

paac <- matrix(nrow = m, ncol = 20 + para1)
for(i in 1:m){ 
paac[i, ] = extractPAAC(A[[i]][1],lambda = para1 , w = para2, props = AAindex_small)
}

#### Dividing dataset into internal and external sets
data = data.frame(aac,paac, Class = D[,ncol(D)])

X <- subset(data, data[,ncol(data)] == 'pos')
Y <- subset(data, data[,ncol(data)] == 'neg')
npos = round(nrow(X)*0.75)
nneg = round(nrow(Y)*0.75)
kenX <- kenStone(X[,-ncol(X)],k=npos,metric='euclid',pc=2)
kenY <- kenStone(Y[,-ncol(Y)],k=nneg,metric='euclid',pc=2)
Index_Pos = rbind(cbind(kenX$model,rep('train',length(kenX$model))),cbind(kenX$test,rep('test',length(kenX$test))))
Index_Neg = rbind(cbind(kenY$model,rep('train',length(kenY$model))),cbind(kenY$test,rep('test',length(kenY$test))))

Index_Pos2 = data.frame(index = as.numeric(Index_Pos[,1]),Index_Pos[,2])
Final_Pos = cbind(X, index = Index_Pos2[order(Index_Pos2[,1]),])
ndel = ncol(Final_Pos)
del = c(ndel-1, ndel )
TR_Pos <- subset(Final_Pos, Final_Pos[,ncol(Final_Pos)] == 'train')[,-del]
TS_Pos <- subset(Final_Pos, Final_Pos[,ncol(Final_Pos)] == 'test')[,-del]
Index_Neg2 = data.frame(index = as.numeric(Index_Neg[,1]),Index_Neg[,2])
Final_Neg = cbind(Y, index = Index_Neg2[order(Index_Neg2[,1]),])
TR_Neg <- subset(Final_Neg, Final_Neg[,ncol(Final_Neg)] == 'train')[,-del]
TS_Neg <- subset(Final_Neg, Final_Neg[,ncol(Final_Neg)] == 'test')[,-del]

internal = rbind(TR_Pos,TR_Neg)
external = rbind(TS_Pos,TS_Neg)


###################################################

data = data.frame(aac, Class = D[,ncol(D)])

X <- subset(data, data[,ncol(data)] == 'pos')
Y <- subset(data, data[,ncol(data)] == 'neg')
npos = round(nrow(X)*0.75)
nneg = round(nrow(Y)*0.75)

Index_Pos2 = data.frame(index = as.numeric(Index_Pos[,1]),Index_Pos[,2])
Final_Pos = cbind(X, index = Index_Pos2[order(Index_Pos2[,1]),])
ndel = ncol(Final_Pos)
del = c(ndel-1, ndel )
TR_Pos <- subset(Final_Pos, Final_Pos[,ncol(Final_Pos)] == 'train')[,-del]
TS_Pos <- subset(Final_Pos, Final_Pos[,ncol(Final_Pos)] == 'test')[,-del]
Index_Neg2 = data.frame(index = as.numeric(Index_Neg[,1]),Index_Neg[,2])
Final_Neg = cbind(Y, index = Index_Neg2[order(Index_Neg2[,1]),])
TR_Neg <- subset(Final_Neg, Final_Neg[,ncol(Final_Neg)] == 'train')[,-del]
TS_Neg <- subset(Final_Neg, Final_Neg[,ncol(Final_Neg)] == 'test')[,-del]

internal = rbind(TR_Pos,TR_Neg)
external = rbind(TS_Pos,TS_Neg)

tunegrid <- expand.grid(.mtry=c(1:5), .ntree=seq(100,500,100))
RFmodel <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=trainControl(method = "cv", number=5))

######Loop for 10-fold CV
k <- 5;
folds <- cvsegments(nrow(internal), k);
true <- data.frame()
label <- data.frame()

for (fold in 1:k){
  currentFold <- folds[fold][[1]];
  RF = randomForest(Class ~ ., internal[-currentFold,], ntree= as.numeric(RFmodel$ bestTune[2]) ,mtry = as.numeric(RFmodel$ bestTune[1]),orm.votes=TRUE,keep.forest=TRUE, importance=TRUE) ## Building RF model
  true = rbind(true, predict(RF, internal[currentFold,],type="prob")[,2])
  label = rbind(label, internal[currentFold,]$Class)
  }
probcv = data.frame(melt(true)[,2])[,1]
write.csv(as.matrix(melt(label)$value), "matrix.csv", row.names=TRUE, na="")
label = read.csv("matrix.csv", header = TRUE)
labelcv = label[,2]

################### External validation
RF = randomForest(Class ~ ., internal, ntree= as.numeric(RFmodel$ bestTune[2]) ,mtry = as.numeric(RFmodel$ bestTune[1]),orm.votes=TRUE,keep.forest=TRUE, importance=TRUE) ## Building RF on internal with the optimized parameter
probext = predict(RF, external,type = "prob")[,2]
labelext = external$Class

RF_auc_cv = cbind(probcv,labelcv)
RF_auc_ext= cbind(probext,labelext)

