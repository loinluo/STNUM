library(ggplot2)
library(FactoMineR)
library(cluster)
library(corrplot)
library(leaps)
library(lmtest)
library(car)
setwd("F:\\R\\Projet\\")

# read the data
newdeath <- read.csv(file="coviddeaths.csv",row.names=1,header=T) 

# create a variable containing the price brackets (delimited by the quartiles):
newdeath$nd_cl <- cut(newdeath$nd,breaks=quantile(newdeath$nd),include.lowest=TRUE)
levels(newdeath$nd)<-c("S1","S2","S3","S4")

# Partie 1

# calculate the variance
r <- 3:10
v <- rep(0,8)
for (i in r) {
  v[i-2] <- var(newdeath[,i])
  
}
print(v)

# draw the Diagrammes en etoile
stars(newdeath[3:10],key.loc=c(10,1.8),main="Diagrammes en etoile",cex=0.5,flip.labels=TRUE)
summary(newdeath)

# draw the Boites a moustaches
boxplot(newdeath[3:10],main="Boites a moustaches")

# draw the Nuages des points
plot(newdeath[3:10],main="Nuages des points")

# calculate the covariance
c <- cor(newdeath[3:10])
corrplot(c, method = "number", col = "black", cl.pos = "n")

# calculate acp
out_acpcase <- PCA(newdeath[3:11],scale.unit=TRUE,ncp=8,quali.sup=9,graph=FALSE)

summary(out_acpcase)

out_acpcase

# draw the Eboulis 
val_prop1 <- out_acpcase$eig[,"eigenvalue"]
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
plot.PCA(out_acpcase,shadow=TRUE,cex=0.8,axes=c(1,2),choix="var",new.plot=TRUE,
         title="Cercle des correlations")

# draw the Projections des individus
plot.PCA(out_acpcase,shadow=TRUE,cex=0.8,axes=c(1,2),choix="ind",label="ind",new.plot=TRUE,
         title="Projection des individus : avec les libell¨¦s des individus")

# draw the Projections des individus en fonction des nouveaux cas pour cent
plot.PCA(out_acpcase,shadow=TRUE,cex=2,axes=c(1,2),choix="ind",habillage=9,label="none",new.plot=TRUE,
         title="Projection des individus : en fonction des nouveaux deces pour cent")

c <- matrix(rep(0,62))
for (i in 1:62) {
  c[i] <- newdeath$nd[i]
  
}

# Partie 2
e <- hist(c,breaks = 14,main="Hisgrame of newdeaths")

# test d¡¯adequation a une loi normale
nortest1<-shapiro.test(c)
nortest1

mu_estim <- mean(c[1:62,])

theta2_estim <- sum((c[1:62,]-(mu_estim))^2)/(62-1)

# Partie 3

# On effectue une regression lineaire multiple
reg_multi <- lm(nd~fv+pd+ma+a65+gdp+cda+hb+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~pd+ma+a65+gdp+cda+hb+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~ma+a65+gdp+cda+hb+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~ma+a65+gdp+cda+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~ma+gdp+cda+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~ma+cda+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~ma+le,data=newdeath)
summary(reg_multi)

reg_multi <- lm(nd~ma,data=newdeath)
summary(reg_multi)

# On initialise l'analyse
alpha <- 0.05
n <- 62
p <- 3
analyses <- data.frame(obs=1:n)

# On calcule les leviers
analyses$levier <- hat(model.matrix(reg_multi))
seuil_levier <- 2*p/n

ggplot(data=analyses,aes(x=obs,y=levier))+
  geom_bar(stat="identity",fill="steelblue")+
  geom_hline(yintercept=seuil_levier,col="red")+
  theme_minimal()+
  xlab("Observation")+
  ylab("Leviers")+
  scale_x_continuous(breaks=seq(0,n,by=5))

idl <- abs(analyses$levier)>seuil_levier
analyses$levier[idl]

# On calcule les r¨¦sidus studentis¨¦s
analyses$rstudent <- rstudent(reg_multi)
seuil_rstudent <- qt(1-alpha/2,n-p-1)

ggplot(data=analyses,aes(x=obs,y=rstudent))+
  geom_bar(stat="identity",fill="steelblue")+
  geom_hline(yintercept=-seuil_rstudent,col="red")+
  geom_hline(yintercept=seuil_rstudent,col="red")+
  theme_minimal()+
  xlab("Observation")+
  ylab("R¨¦sidus studentis¨¦s")+
  scale_x_continuous(breaks=seq(0,n,by=5))

ids <- abs(analyses$rstudent)>seuil_rstudent
analyses$rstudent[ids]

# On calcule les distances de Cook
influence <- influence.measures(reg_multi)
names(influence)
colnames(influence$infmat)
analyses$dcook <- influence$infmat[,"cook.d"]
seuil_dcook <- 4/(n-p)

ggplot(data=analyses,aes(x=obs,y=dcook))+
  geom_bar(stat="identity",fill="steelblue")+
  geom_hline(yintercept=seuil_dcook,col="red")+
  theme_minimal()+
  xlab("Observation")+
  ylab("Distance de Cook")+
  scale_x_continuous(breaks=seq(0,n,by=5))

idc <- analyses$dcook>seuil_dcook
analyses$dcook[idc]

vif(reg_multi)

bptest(reg_multi)

shapiro.test(reg_multi$residuals)

# On procede a une selection automatique de modeles

reg_null <- lm(nd~1,data=newdeath)
reg_tot <- lm(nd~fv+pd+ma+a65+gdp+cda+hb+le,data=newdeath)

reg_backward <- step(reg_tot,direction="backward")

# On separer aleatoirement l¡¯echantillon en 2 parties
select <- sample(1:n,round(n*0.8))
newdeath_apprent <- newdeath[select,]
newdeath_test <- newdeath[-select,]
reg_multi_apprent <- lm(nd~ma,data=newdeath_apprent)
summary(reg_multi_apprent)
newdeath_test$nd_prev <- predict(reg_multi_apprent,newdeath_test)
newdeath_test$err_prev <- newdeath_test$nd-newdeath_test$nd_prev

# On estime les RMSE et MAPE
rmse <- sqrt(mean((newdeath_test$err_prev)^2))
round(rmse,digits=2)
