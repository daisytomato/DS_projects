#############load package#########
library(stats)
library(cluster)
library(FPDclustering)
library(MixGHD)
library(factoextra)
#install.packages("dbscan")
library(dbscan)
#install.packages("Spectrum")
#library(Spectrum)
#install.packages("factoextra")
library(factoextra)
#install.packages("ggfortify")
library(ggfortify)
#library(devtools)
library(ggfortify); library(ggplot2)
#install.packages("DataCombine")
library(DataCombine)
##install.packages("caTools")
library(caTools)
#install.packages("corrplot")
library(corrplot)

########################################
############Cellular Data###############
########################################


# df3<-read.csv("ML-MATT-CompetitionQT1920_train.csv")
# save(df3,file="Cellular.Rdata")
#View(df3)
load("Cellular.Rdata")
summary(df3)
dim(df3)
dim(df3[is.na(df3$meanUE_UL)==FALSE,])
#subset dataset 3
X<-as.matrix(df3[,3:10])
dim(X)
dim(df3[df3$Unusual==1,])
dim(df3[df3$Unusual==0,])

#PLOT THE CORRELAITON
M3<-cor(X)
corrplot(M3, method="color", tl.srt=45,tl.col="black",tl.cex=0.5)

p<-ggplot(df3,aes(factor(Unusual)))
p+geom_bar(aes(fill=factor(Unusual)))+xlab("Type")+ylab("Counts")+theme(legend.position="none")+scale_x_discrete(labels=c("Non-Fraud", "Fraud"))
pairs(X)


##plot PCA
pca_res <- prcomp(X, scale. = TRUE)
summary(pca_res)
autoplot(pca_res,data=X,colour=df3$Unusual+1L,main="True Label",shape=21)

################
#####kmeans#####
################


kmeans.out<-kmeans(X,2,nstart=10)

#calculate the silhouette index 
ec2<-eclust(X, FUNcluster = "kmeans",k=2, hc_metric = "euclidean",nstart=25)
ARI(ec2$cluster,label)
fviz_silhouette(ec2)

#find out the clustering from kmeans and compare with true labels
dim(df3[kmeans.out$cluster==1,])
dim(df3[(kmeans.out$cluster==2),])
table(kmeans.out$cluster,df3[,14])
ARI(kmeans.out$cluster,df3[,14])

#plot on the kmeans with pre-PCA
fviz_cluster(kmeans.out,data=X,geom="points")

autoplot(pca_res,data=X,colour=kmeans.out$cluster,main="K-means",shape=21)

################
######PDQ#######
################

pdq.out<-PDQ(X,K=2)

table(df3[,14],pdq.out$label)

Silh(pdq.out$probability) #plot silhouette probability 
ARI(pdq.out$label,df3[,14])

autoplot(pca_res,data=X,colour=pdq.out$label,main="PDQ",shape=21)

################
######MGHD######
################

ghd.out<-MGHD(X,G=2,modelSel="BIC")

plot(ghd.out)
table(ghd.out@map,df3[,14])
ARI(ghd.out@map,df3[,14])
fviz_cluster(ghd.out,data=X)
autoplot(pca_res,data=X,colour=ghd.out@map,main="MGHD",shape=21)

################
#####DBSCAN#####
################

X_scale<-scale(X)
#The optimal value of "eps" parameter can be determined as follow:
{kNNdistplot(X_scale, k =5)
  abline(h=2,lty=2)
}
dbscan.out<-dbscan(X_scale,eps=2,minPts=5)
dbscan.out
pairs(X,col=dbscan.out$cluster+1L)
table(dbscan.out$cluster,df3[,14])
ARI(dbscan.out$cluster,df3[,14])
plot(X[,3:4],col=dbscan.out$cluster)
fviz_cluster(dbscan.out,data=X_scale,geom="point")
autoplot(pca_res,data=X,colour=dbscan.out$cluster+1L,main="DBSCAN",shape=21)


########################################
############Accounting Data#############
########################################

# df<-read.csv("uscecchini28.csv",header=T)
# #View(df)
# dim(df)
# summary(df)

# save(df,file="accouting.Rdata")
##drop the rows with NAs

## loading accounting Rdata
load('accounting.Rdata')

###################################
##############preprocessing########
###################################
df<-DropNA(df[,-c(1:8)])
summary(df)
non_fraud<-nrow(df[df$misstate==0,])#number of non-fraud
fraud<-nrow(df[df$misstate==1,]) #number of fraud
##ratio of fraud/non-fraud
ratio1=non_fraud/(fraud+non_fraud)
ratio2=fraud/(non_fraud+fraud)
###randomly choose the data samples
set.seed(1)
sample1<-sample(10000,round(10000*ratio1,0))
sample2<-sample(909,round(10000*ratio2,0))

df1<-df[df$misstate==0,]
df2<-df[df$misstate==1,]
dim(df1[sample1,])
summary(df2)
summary(df1)
dim(df2[sample2,])

df_merge<-data.frame(t(df1[sample1,]),t(df2[sample2,]))
df_merge <- data.frame(t(df_merge))


summary(df_merge)
dim(df_merge)
summary(df)

df1<-df_merge[,-1]
dim(df1[label==0,])
dim(df1)


label<-df_merge[,1]

#PLOT THE DATASET ON PCA
p<-ggplot(df_merge,aes(factor(misstate)))
p+geom_bar(aes(fill=factor(misstate)))+xlab("Type")+ylab("Counts")+theme(legend.position="none")+scale_x_discrete(labels=c("Non-Fraud", "Fraud"))

pca_df1 <- prcomp(df1[,28:42], scale. = TRUE)
pca_df1$rotation[,1]
pca1<-as.data.frame(pca_df1$x)

##correlation plot
M<-cor(df1)
dim(df1[,28:42])
corrplot(M, method="color", tl.pos="n")
dev.off()
autoplot(pca_df1,data=df1,colour=label+1,main="TureLabel",shape=21)


################
#####kmeans#####
################
kmeans.out<-kmeans(df1[,28:42],2,nstart=10)#no NAS are allowed in kmeans

#try out eclust function and plot out silhouette 
ec1<-eclust(df1[,28:42], FUNcluster = "kmeans",k=2, hc_metric = "euclidean",nstart=25)
ARI(ec1$cluster,label)
fviz_silhouette(ec1)
##find out the clustering from kmeans and compare with true labels
ARI(kmeans.out$cluster,label)
table(kmeans.out$cluster,label)

#plot out the k-means label on PCA
autoplot(pca_df1,data=df1,colour=kmeans.out$cluster,main="K-means",shape=21)

################
#######PDQ######
################

df1<-as.matrix(df1)
pdq.out<-PDQ(df1[,28:42],K=2)
pdq.out$label
table(label,pdq.out$label)
Silh(pdq.out$probability)
ARI(pdq.out$label,label)
pca_df1 <- prcomp(df1, scale. = TRUE)

#plot PDQ label on PCA
autoplot(pca_df1,data=df1,colour=pdq.out$label,main="PDQ",shape=21)

################
#######GHD######
################

ghd.out<-MGHD(df1[,28:42],G=2,modelSel="BIC")
ghd.out_unlimited<-MGHD(df1,G=2:4,modelSel="BIC")
ghd.out
plot(ghd.out)
table(ghd.out@map,label)
ARI(ghd.out@map,label)
ARI(ghd.out_unlimited@map,label)
#try to use another method for initialization BUT DOES NOT WORK
#ghd.out2<-MGHD(df1,G=2,method="modelBased",modelSel="BIC")
#table(ghd.out2@map,label)
#ARI(ghd.out2@map,label)
autoplot(pca_df1,data=df1,colour=ghd.out@map,main="MGHD",shape=21)
contourpl(ghd.out@map)

################
#####DBSCAN#####
################

df1_scale<-scale(df1[,28:42])
{kNNdistplot(df1_scale, k =4)
  abline(h=5,lty=2)
}
dbscan.out<-dbscan(df1_scale,eps=5,minPts=4)
dbscan.out
pairs(df1,col=dbscan.out$cluster+1L)
fviz_cluster(dbscan.out, df1_scale, geom = "point")

table(dbscan.out$cluster,label)
ARI(dbscan.out$cluster,label)
ARI(hdbscan.out$cluster,label)

#plot the labels on PCA
autoplot(pca_df1,data=df1,colour=dbscan.out$cluster+1L,main="DBSCAN",shape=21)

########################################
##############Cancer data###############
########################################

# cancer <- read.csv("data.csv", header = TRUE)
# labels <- read.csv("labels.csv", header = TRUE)

# ##subset the largest class (BRCA) and the smallest class (COAD)
# BRCA <- cancer[labels$Class=='BRCA', ]
# COAD <- cancer[labels$Class=='COAD', ]
# labels.new <- rep(c('BRCA', 'COAD'), c(300,78))
# cancer.new <- rbind(BRCA, COAD)
# #select the first 150 variables
# cancer.new <- as.matrix(cancer.new[, 2:151]) 
# #save as RData files
# save(cancer.new, file = "cancer1.RData")
# save(labels.new, file = "labels1.RData")

###############################################
########## load the RData files################
###############################################
########## use the first 150 columns ##########
###############################################

load('cancer1.RData')
load('labels1.RData')

## kmeans
cancer.kmeans <- kmeans(cancer.new, centers = 2, nstart = 50)
table(cancer.kmeans$cluster, labels.new)
ARI(cancer.kmeans$cluster, labels.new)

## PDQ
library(FPDclustering)
cancer.pdq <- PDQ(cancer.new, K=2)
table(cancer.pdq$label, labels.new)
ARI(cancer.pdq$label, labels.new)

## GHD
## Error in MGHD(cancer.new, G = 2) : No NAs allowed.
cancer.ghd <- MGHD(cancer.new, G= 2) 
ARI(cancer.ghd@map, labels.new)
table(cancer.ghd@map, labels.new)

## DBSCAN
library(dbscan)
library(fpc)
dbscan::kNNdistplot(cancer.new, k = 4)
abline(h = 17.7, lty = 2)
cancer.fpc1 <- fpc::dbscan(cancer.new, eps = 17.15, MinPts = 4) 
ARI(cancer.fpc1$cluster, labels.new)
table(cancer.fpc1$cluster, labels.new)
fviz_cluster(cancer.fpc1, data = cancer.new, geom = "point")

######################################################
########## remove columns if there are less ##########
######### than 11 rows with nonzero values ###########
rm.col <- c()
for (col in 1:50) {
  if (sum(cancer.new[, col] !=0) <11) {
    rm.col <- c(rm.col, col)
  }
}
temp <- cancer.new[, -1*rm.col]
## there are some columns with correlations, need to find them and remove
correlated_cols <- findCorrelation(cor(temp, 
                                       use = "pairwise.complete.obs"), 
                                   cutoff = 0.3, verbose = TRUE)
cancer.new.rm <- temp[, -1*(correlated_cols)]

## kmeans
cancer.rm.kmeans <- kmeans(cancer.new.rm, centers = 2, nstart = 50)
table(cancer.rm.kmeans$cluster, labels.new)
ARI(cancer.rm.kmeans$cluster, labels.new)
fviz_cluster(cancer.rm.kmeans, data = cancer.new.rm, geom = "point")

## PDQ
library(FPDclustering)
cancer.rm.pdq <- PDQ(cancer.new.rm, K=2)
table(cancer.rm.pdq$label, labels.new)
ARI(cancer.pdq$label, labels.new)

## GHD
cancer.rm.ghd <- MGHD(cancer.new.rm, G= 2)
ARI(cancer.rm.ghd@map, labels.new)
table(cancer.rm.ghd@map, labels.new)


## DBSCAN for the new data set removed the correlated columns
dbscan::kNNdistplot(cancer.new.rm, k = 5)
abline(h = 9.5, lty = 2)
cancer.rm.fpc1 <- fpc::dbscan(cancer.new.rm, eps = 9.5, 
                              MinPts = 5) 
ARI(cancer.rm.fpc1$cluster, labels.new)
table(cancer.rm.fpc1$cluster, labels.new)
fviz_cluster(cancer.rm.fpc1, data = cancer.new.rm, geom = "point")


###########################################################
########### cluster plot in the first 2 dimensions ########
################### for the cancer.new.rm #################
###########################################################

# PCA
pca <- prcomp(cancer.new.rm, scale. = TRUE)

my_cols <- c("#110F0F", "#FF6666")
my_cols_1 <- c("#FF6666","#110F0F")
par(mfrow = c(3,2),mai = c(0.55, 0.6, 0.3, 0.3))
plot(pca$x[, 1], pca$x[, 2], col = my_cols_1[as.factor(labels.new)], 
     xlab = "PC1 (7%)", ylab = "PC2 (6%)")
text(3.4, -6, "True clusters", cex = 1.5)

#empty plot
plot.new()

plot(pca$x[, 1], pca$x[, 2], col = my_cols[as.factor(cancer.rm.kmeans$cluster)], 
     xlab = "PC1 (7%)", ylab = "PC2 (6%)")
text(3.4, -6, "K-means", cex = 1.5)
plot(pca$x[, 1], pca$x[, 2], col = my_cols[as.factor(cancer.rm.pdq$label)],
     xlab = "PC1 (7%)", ylab = "PC2 (6%)")
text(3.4, -6, "PDQ", cex = 1.5)
plot(pca$x[, 1], pca$x[, 2], col = my_cols[as.factor(cancer.rm.ghd@map)],
     xlab = "PC1 (7%)", ylab = "PC2 (6%)")
text(3.4, -6, "MGHD", cex = 1.5)
plot(pca$x[, 1], pca$x[, 2], col = my_cols[as.factor(cancer.rm.fpc1$cluster)],
     xlab = "PC1 (7%)", ylab = "PC2 (6%)")
text(3.4, -6, "DBSCAN", cex = 1.5)
