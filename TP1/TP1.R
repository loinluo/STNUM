library(ggplot2)
library(FactoMineR)
library(cluster)
library(corrplot)
setwd("F:\\R\\TP1\\")

# read the data
painting8 <- read.table("painting8.txt",sep=";",dec=".") 
# repeat R 40 times, repeat VG 44 times and combine them
painting8$peintre <- as.factor(c(rep("R",40),rep("VG",44)))  
painting64 <- read.table("painting64.txt",sep=";",dec=".")
painting64$peintre <- as.factor(c(rep("R",40),rep("VG",44)))

# calculate the variance of p8
v <- var(painting8[,1:7]) 

# draw the Boites a moustaches of p8
boxplot(painting8[,1:7],main="Boites a moustaches")

# draw the Nuages des points of p8
plot(painting8[,1:7],main="Nuages des points")

# calculate the covariance of p8
c <- cor(painting8[,1:7])
corrplot(c, method = "number", col = "black", cl.pos = "n")

# calculate acp8n
out_acp8n <- PCA(painting8,scale.unit=TRUE,ncp=7,quali.sup=8,graph=FALSE)

# draw the Eboulis 
val_prop1 <- out_acp8n$eig[,"eigenvalue"]
val_prop_cum1 <- cumsum(val_prop1)/sum(val_prop1)
cp1 <- 1:length(val_prop1)
vp1 <- data.frame(cp1=cp1,val_prop1=val_prop1)
vp_cum1 <- data.frame(cp1=cp1,val_prop_cum1=val_prop_cum1)

ggplot(data=vp1,aes(x=cp1,y=val_prop1))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Eboulis des valeurs propres")+
  xlab("Nombre de composantes principales")+
  ylab("Valeurs propres")+
  scale_x_continuous(breaks=cp1)

ggplot(data=vp_cum1,aes(x=cp1,y=val_prop_cum1))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Part d'inertie expliqu¨¦e en fonction du nombre de CP")+
  xlab("Nombre de composantes principales")+
  ylab("Part d'inertie expliqu¨¦e")+
  scale_x_continuous(breaks=cp1)

# draw the Cercle des correlations
plot.PCA(out_acp8n,shadow=TRUE,cex=0.8,axes=c(1,2),choix="var",new.plot=TRUE,
         title="Cercle des correlations")

# draw the Projections des individus en fonction des peintre
plot.PCA(out_acp8n,shadow=TRUE,cex=2,axes=c(1,2),choix="ind",habillage=8,label="none",new.plot=TRUE,
         title="Projection des individus : en fonction des peintre")

# calculate acp8nn
out_acp8nn <- PCA(painting8,scale.unit=FALSE,ncp=7,quali.sup=8,graph=FALSE)

# draw the Eboulis 
val_prop2 <- out_acp8nn$eig[,"eigenvalue"]
val_prop_cum2 <- cumsum(val_prop2)/sum(val_prop2)
cp2 <- 1:length(val_prop2)
vp2 <- data.frame(cp2=cp2,val_prop2=val_prop2)
vp_cum2 <- data.frame(cp2=cp2,val_prop_cum2=val_prop_cum2)

ggplot(data=vp2,aes(x=cp2,y=val_prop2))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Eboulis des valeurs propres")+
  xlab("Nombre de composantes principales")+
  ylab("Valeurs propres")+
  scale_x_continuous(breaks=cp2)

ggplot(data=vp_cum2,aes(x=cp2,y=val_prop_cum2))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Part d'inertie expliqu¨¦e en fonction du nombre de CP")+
  xlab("Nombre de composantes principales")+
  ylab("Part d'inertie expliqu¨¦e")+
  scale_x_continuous(breaks=cp2)

# draw the Projections des individus en fonction des peintre
plot.PCA(out_acp8nn,shadow=TRUE,cex=2,axes=c(1,2),choix="ind",habillage=8,label="none",new.plot=TRUE,
         title="Projection des individus : en fonction des peintre")

# calculate acp64nn
out_acp64nn <- PCA(painting64,scale.unit=FALSE,ncp=63,quali.sup=64,graph=FALSE)

# draw the Eboulis 
val_prop3 <- out_acp64nn$eig[,"eigenvalue"]
val_prop_cum3 <- cumsum(val_prop3)/sum(val_prop3)
cp3 <- 1:length(val_prop3)
vp3 <- data.frame(cp3=cp3,val_prop3=val_prop3)
vp_cum3 <- data.frame(cp3=cp3,val_prop_cum3=val_prop_cum3)

ggplot(data=vp3,aes(x=cp3,y=val_prop3))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Eboulis des valeurs propres")+
  xlab("Nombre de composantes principales")+
  ylab("Valeurs propres")+
  scale_x_continuous(breaks=cp3)

ggplot(data=vp_cum3,aes(x=cp3,y=val_prop_cum3))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Part d'inertie expliqu¨¦e en fonction du nombre de CP")+
  xlab("Nombre de composantes principales")+
  ylab("Part d'inertie expliqu¨¦e")+
  scale_x_continuous(breaks=cp3)

# draw the Projections des individus en fonction des peintre
plot.PCA(out_acp64nn,shadow=TRUE,cex=2,axes=c(1,2),choix="ind",habillage=64,label="none",new.plot=TRUE,
         title="Projection des individus : en fonction des peintre")

# draw the ward
cah_ward <- agnes(painting8[,1:7],metric="euclidean",stand=TRUE,method="ward")
plot(as.dendrogram(cah_ward),main="Ward")

ei <- data.frame(k=2:dim(painting8)[1],height=sort(cah_ward$height,decreasing=TRUE))

# draw the Gain d'inertie inter-classes lors du passage de (k-1) ¨¤ k classes
ggplot(data=ei,aes(x=k,y=height))+
  geom_bar(stat="identity",fill="steelblue")+
  theme_minimal()+
  ggtitle("Gain d'inertie inter-classes lors du passage de (k-1) ¨¤ k classes")+
  xlab("k")+
  ylab("Indice d'aggr¨¦gation")+
  scale_x_continuous(breaks=2:dim(painting8)[1])

# choose the number of class and draw the Ward 
plot(as.dendrogram(cah_ward),main="Ward")
rect.hclust(cah_ward,k=3)

# do the classification
cluster_ward_3cl <- data.frame(nom=rownames(painting8),classe=cutree(cah_ward,k=3))
cluster_ward_3cl[order(cluster_ward_3cl$classe),]

# K-mean k=2
painting8_nor <- scale(painting8[,1:7],center=TRUE,scale=TRUE)
cluster_kmeans_3cl <- kmeans(painting8_nor,centers=2,nstart=20)
print(cluster_kmeans_3cl)


      
