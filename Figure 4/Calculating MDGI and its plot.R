setwd('D:\\Peptide prediction\\Tumor homing peptides\\Data set')
library(randomForest)
library(Interpol)
library(protr)
library(seqinr)
library(reshape2)
library(ggplot2)
#######Read data
A <- read.fasta('main.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("label main.csv", header = TRUE) 
aa = read.csv("AA.csv", header = TRUE) 
dp = read.csv("DP.csv", header = TRUE) 
m = length(A)
aac<- t(sapply(A, extractAAC))
dpc <- t(sapply(A, extractDC))

internal = data.frame(aac,Class = label)

###########Prepare MDGI for AAC
ind= c(2,3,5,7,9,11,13,15,17,20)
n = ncol(internal)-1
gini = matrix(nrow = n, ncol = 10)

for (i in 1:10){
RF<-randomForest(Class~.,data=internal,ntree=100,mtry=ind[i],importance=TRUE)
gini[,i] = as.matrix(RF$ importance [,4])
}

meangini = matrix(nrow = nrow(gini), ncol = 1)
for (i in 1:nrow(gini)){
meangini[i,]<-mean(gini[i,])
}

###########Plot MDGI for AAC
data= data.frame(gini, aa)
data_melt1 <- melt(data[1:20,], id.vars = "AA")
data_melt1$feature <- factor(data_melt1$AA)

par(mfrow = c(2,1 ),oma=c(0, 0, 0, 2))

plot1 <- ggplot(data_melt1, aes(x = reorder(feature, value, FUN = mean), y = value)) +
  geom_boxplot(fill = "green", colour = "black", alpha = 0.5) +
  theme_bw() + xlab("") + ylab("Mean decrease of Gini index") + coord_flip() + theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 13, colour = "black"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 13, face = "bold", colour = "black")
  )

###########Prepare MDGI for DPC
internal = data.frame(dpc,Class = label)

############# Feature extraxtion using RF model
ind= c(2,3,5,7,9,11,13,15,17,20)
n = ncol(internal)-1
gini = matrix(nrow = n, ncol = 10)

for (i in 1:10){
RF<-randomForest(Class~.,data=internal,ntree=100,mtry=ind[i],importance=TRUE)
gini[,i] = as.matrix(RF$ importance [,4])
}

meangini = matrix(nrow = nrow(gini), ncol = 1)
for (i in 1:nrow(gini)){
meangini[i,]<-mean(gini[i,])
}

###########Plot MDGI for DPC
data= data.frame(gini, dp,meangini)
data2 = data[order(data$meangini, decreasing=T),]
data_melt1 <- melt(data2[1:20,-ncol(data2)], id.vars = "DP")
data_melt1$feature <- factor(data_melt1$DP)


plot2 <- ggplot(data_melt1, aes(x = reorder(feature, value, FUN = mean), y = value)) +
  geom_boxplot(fill = "blue", colour = "black", alpha = 0.5) +
  theme_bw() + xlab("") + ylab("Mean decrease of Gini index") + coord_flip() + theme(
    axis.text.y = element_text(size = 13, colour = "black"),
    axis.text.x = element_text(size = 13, colour = "black"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 13, face = "bold", colour = "black")
  )
#############################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#################### multiplot 
multiplot(plot1,plot2, cols=2)
