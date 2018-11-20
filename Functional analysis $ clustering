
install.packages('data.table')
install.packages("cluster")
install.packages("shape")
library(cluster)
library(shape)
library('data.table')
library('fda')
#############################################PROGRAM#############################################
data=fread("gene.txt") #READ DATA
gene=data[1:34]
Data <- subset( gene, select = -gene )
Data=t(as.matrix(sapply(Data, as.numeric)))
daybasis10 = create.fourier.basis(c(0,10),35)
harmLfd = vec2Lfd(c(0,(2*pi/10)^2,0), c(0, 10))
tempfdPar = fdPar(daybasis10,harmLfd,1e4)
tempfd = smooth.basis(1:10,Data,tempfdPar)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='Time',ylab='Gene expression',cex.lab=1,cex.axis=1)
lines(tempfdPar$meanfd,lwd=2.5,col=2)
tempfdPar2 = fdPar(daybasis10,harmLfd,1e-2)
tempfdPar3 = fdPar(daybasis10,harmLfd,1e4)
tempfd2 = smooth.basis(1:10,Data,tempfdPar2)

# do FPCA with a roughness penalty on FPCs. 

ptemppca = pca.fd(tempfd2$fd,nharm=4,harmfdPar=tempfdPar3)

# Get FPCs
pharmfd = ptemppca$harmonics
pharmvals = eval.fd(1:10,pharmfd)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:10,pharmvals,xlab='Time',ylab='PCs',
        lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')

legend(0.7,-0.01,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=2)
title('Gene Expression Component Functions')
####################################################################


daybasis10 = create.fourier.basis(c(0,10),10)

# harmonic acceleration differential operator
harmLfd = vec2Lfd(c(0,(2*pi/10)^2,0), c(0, 10))
tempfdPar = fdPar(daybasis10,harmLfd,1e4)
tempfd = smooth.basis(1:10,Data,tempfdPar)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='day',ylab='temperature',cex.lab=1.5,cex.axis=1.5)

tempvar = var.fd(tempfd$fd)
tvvals = eval.bifd(1:10,1:10,tempvar)
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
contour(1:10,1:10,tvvals,xlab='time',ylab='time',cex.lab=1.5,cex.axis=1.5)

quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(cumsum(ptemppca$values[1:10])/sum(ptemppca$values),xlab='component',
     ylab='cumulative variance explained',col=4,cex.lab=1.5,
     cex.axis=1.5,cex=1)
abline(h=0.99)

ptemppca = pca.fd(tempfd2$fd,nharm=4,harmfdPar=tempfdPar3)

# Get FPCs
pharmfd = ptemppca$harmonics
pharmvals = eval.fd(1:10,pharmfd)

# plot the FPCs
quartz()
par(mfrow=c(2,2),mar = c(8, 8, 4, 2))
for (i in 1:4){
plot(1:10,pharmvals[,i],xlab='day',ylab='PCs',
     lwd=2,lty=1,cex.lab=1,cex.axis=1,type='l', col=i)
}
# plot all 4 FPCs
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
matplot(1:10,pharmvals,xlab='day',ylab='PCs',
        lwd=4,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
legend(1,-0.03,c('PC1','PC2','PC3','PC4'),col=1:4,lty=1,lwd=5)
title('Gene Principle Component Functions')
# plot the first FPC scores vs. the second FPC scores 
quartz()
par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(ptemppca$scores[,1:2],xlab='PC Score 1',ylab='PC Score 2',col=11,
     cex.lab=1.5,cex.axis=1.5,cex=1, lwd=4)
#Clustring
set.seed(100)
kmean<-kmeans(ptemppca$scores[,c(1,4)],4,iter.max = 5000,nstart = 10)
clusplot(ptemppca$scores[,c(1,4)],kmean$cluster,line=0,shade = TRUE,color = TRUE,labels = 4 ,plotchar = TRUE,span = TRUE,main = paste('Clusters of the gene'))

kmean$cluster
cluster=kmean$cluster
cluster1 <- which(cluster == 1)
cluster1_data <- Data[, cluster1]
cluster2 <- which(cluster == 2)
cluster2_data <- Data[, cluster2]
cluster3 <- which(cluster == 3)
cluster3_data <- Data[, cluster3]
cluster4 <- which(cluster == 4)
cluster4_data <- Data[, cluster4]
tempc1 = smooth.basis(1:10,cluster1_data,tempfdPar)
tempc2 = smooth.basis(1:10,cluster2_data,tempfdPar)
tempc3 = smooth.basis(1:10,cluster3_data,tempfdPar)
tempc4 = smooth.basis(1:10,cluster4_data,tempfdPar)
par(mfrow=c(2,2),mar = c(8, 8, 4, 2))
plot(tempc1$fd,xlab='Time',ylab='Gene expression',cex.lab=1,cex.axis=1)
title('Gene Expression Component Functions of Cluster #1')
plot(tempc2$fd,xlab='Time',ylab='Gene expression',cex.lab=1,cex.axis=1)
title('Gene Expression Component Functions of Cluster #2')
plot(tempc3$fd,xlab='Time',ylab='Gene expression',cex.lab=1,cex.axis=1)
title('Gene Expression Component Functions of Cluster #3')
plot(tempc4$fd,xlab='Time',ylab='Gene expression',cex.lab=1,cex.axis=1)
title('Gene Expression Component Functions of Cluster #4')

  
