library(ggplot2) # Pour les graphiques ggplot
library(gridExtra) # Pour l'organisation des graphiques ggplot
library(EnvStats) # Pour le test sur les variances
setwd("F:\\R\\TP2\\")

# Q1 Simuler n=500 réalisations indépendantes d’une v.a. de loi N 
x <- matrix(rnorm(500,mean=5,sd=2),1,500)

# Q2 Draw the density function and decide theta
d <- density(x)
plot(d,main="Density function")

# Q3 la consistance des estimateurs
consist_estim <- function(x,theta)
{
  n <- length(x)
  
  theta1 <- cumsum(x)/(1:n) 
  
  theta2 <- 1:n
  for (i in 1:n)
  {
    theta2[i] <- median(x[1:i])
  }
  
  theta3= 1:n
  for (i in 1:n)
  {
    theta3[i] <- (min(x[1:i])+max(x[1:i]))/2
  }
  
  df <- data.frame(x=1:n,theta1=theta1,theta2=theta2,theta3=theta3)
  
  g1 <- ggplot(data=df,aes(x=x,y=theta1))+
    geom_line()+
    geom_hline(yintercept=theta,col="red")+
    labs(x="n",y="theta1")
  
  g2 <- ggplot()+
    geom_line(data=df,aes(x=x,y=theta2))+
    geom_hline(yintercept=theta,col="red")+
    labs(x="n",y="theta2")
  
  g3 <- ggplot()+
    geom_line(data=df,aes(x=x,y=theta3))+
    geom_hline(yintercept=theta,col="red")+
    labs(x="n",y="theta3")
  
  grid.arrange(g1,g2,g3,ncol=1,nrow=3)
  
  return(matrix(c(theta1,theta2,theta3),n,2))
}

# Q4 Lancer la fonction
c1 <- consist_estim(x,5)

# Q5 biais et la variance des estimateurs
compar_estim_norm <- function(N,n)
{
  x <- matrix(rnorm(N*n,mean=5,sd=2),N,n)
  
  theta1 <- apply(x,1,mean)
  theta2 <- apply(x,1,median)
  theta3 <- (apply(x,1,min)+apply(x,1,max))/2
  est <- rep(c("1","2","3"),each=N)
  
  df <- data.frame(est=est,theta=c(theta1,theta2,theta3))
  
  ggplot(df,aes(x=est,y=theta))+
    geom_boxplot()+
    stat_summary(fun=mean,geom="point",shape=20,size=4,col="red")+
    scale_x_discrete(labels=expression(widehat(theta)[1],
                                       widehat(theta)[2],
                                       widehat(theta)[3]))+
    labs(x="Estimateurs",y="")
}

# Q6 Lancer la fonction
compar_estim_norm(200,500)

# Q7 Simuler n=500 réalisations indépendantes d’une v.a. de loi u
y=matrix(runif(500,min=0,max=8),1,500)

# Q8 la consistance des estimateurs
c2 <- consist_estim(y,4)

# Q9 le biais et la variance des estimateurs
compar_estim_unif <- function(N,n)
{
  x <- matrix(runif(N*n,min=0,max=8),N,n)
  
  theta1 <- apply(x,1,mean)
  theta2 <- apply(x,1,median)
  theta3 <- (apply(x,1,min)+apply(x,1,max))/2
  est <- rep(c("1","2","3"),each=N)
  
  df <- data.frame(est=est,theta=c(theta1,theta2,theta3))
  
  ggplot(df,aes(x=est,y=theta))+
    geom_boxplot()+
    stat_summary(fun=mean,geom="point",shape=20,size=4,col="red")+
    scale_x_discrete(labels=expression(widehat(theta)[1],
                                       widehat(theta)[2],
                                       widehat(theta)[3]))+
    labs(x="Estimateurs",y="")
}

compar_estim_unif(200,500)

# read the data 
poulpe1 <- read.table("poulpe1.txt",header=TRUE,sep=";",dec=".")

n <- dim(poulpe1)[1]

# Q10 Estimation d’une moyenne
mu_estim <- mean(poulpe1[1:n,])

# Q11 Estimation d’une variance non biaise
theta2_estim <- sum((poulpe1[1:n,]-(mu_estim))^2)/(n-1)

# Q12 l’intervalle de confiance pour μ
alpha1 <- 0.1
t.test(poulpe1[1:n,],conf.level=1-alpha1)

alpha2 <- 0.05
t.test(poulpe1[1:n,],conf.level=1-alpha2)

# read the data
poulpe2 <- read.table("poulpe2.txt",header=TRUE,sep=";",dec=".")

# Q13 Test de comparaison
femelle <- poulpe2[poulpe2$sexe=="femelle",]$poids
male <- poulpe2[poulpe2$sexe=="male",]$poids

t.test(femelle,male,var.equal=TRUE)

var.test(femelle,male)

# Q14 Déterminer si les naissances sont equiréparties dans la semaine

n <- c(564772,629408,609596,604812,605280,450840,404456)
p0 <- rep(1/7,7)
chisq.test(n,p=p0)

