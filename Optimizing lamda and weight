#######set directory
setwd('D:\\Peptide prediction\\Tumor homing peptides\\Data set')
#######Load package
library(randomForest)
library(protr)
library(seqinr)
library(C50)
library(RWeka)
library(Interpol)

AAindex_main = c("CHAM830101","ANDN920101","CHAM820102","BUNA790103","BROC820102","CHAM830102","BURA740101","BEGF750102","CHAM810101","ARGP820103","BURA740102","BUNA790102","FASG890101","BROC820101",
"BUNA790101","BULH740102","BULH740101","CHAM820101","BIOV880102","BIOV880101","ARGP820102","BHAR880101","BEGF750103","BEGF750101","ARGP820101","BIGC670101","CHOC760103","CHAM830107","CHAM830104",
"CHOP780204")

AAindex_main90 = c("NAKH900101","AURR980120","ONEK900102","AURR980119","GUYH850101","ONEK900101","NAKH900102","AURR980118","NAGK730102","NAGK730103",
"HUTJ700102","LAWE840101","KRIW790101","HUTJ700101","HOPT810101","LEVM760104","KYTJ820101","HUTJ700103","BURA740101","KRIW790103",
"HOPA770101","CHAM830101","ISOY800101","NAKH900105","KLEP840101","CHAM820101","ISOY800103","BUNA790102","KANM800103","KRIW790102")

AAindex_small = c("ANDN920101","CIDH920103","BUNA790103","CIDH920101","FASG890101","CHOP780206","BIOV880101","BIOV880102","CHOP780214","BHAR880101",
"ARGP820103","COHE430101","BULH740101","CIDH920104","BEGF750102","CHOP780216","CHAM830102","BULH740102","CHOC750101","CHOP780212","CIDH920105",
"CIDH920102","BIGC670101","ARGP820102","CHAM830104","CHOP780210","CHOP780203","BURA740101","CHAM830107","BEGF750103")


x <- read.fasta('main.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("label main.csv", header = TRUE) 

lam = 10
w = 10
m = length(x)
A <- x[(sapply(x, protcheck))]
index = seq(0.1,1, 0.1)
error  <- matrix(nrow = length(index), ncol = lam)

for(pse in 1:lam){
paac <- matrix(nrow = m, ncol = 20 + pse)
for(weight in 1:w){
for(i in 1:m){ 
paac[i, ] = extractPAAC(A[[i]][1],lambda = pse, w = index[weight], props = AAindex[1:25])
}
internal= data.frame (paac, Class = label)
ntree <- randomForest(Class ~ ., internal, ntree= 100)
error[weight,pse] <- sum(ntree $ confusion[,3])
}}

t(error)
