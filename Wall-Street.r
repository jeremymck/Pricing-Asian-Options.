rm(list=ls())
library(MASS)
library(tidyverse)
library(stargazer)
library(randtoolbox)
"
une id??e g??n??rale: pour le rendu graphique, on peut plotter l'??volution de la moyenne
?? chaque p??riode (k p??riodes) avec des error bars
"

########################################
############# QUESTION 1 ###############
########################################

#on doit simuler plusieurs fois la grandeur dans l'esp??rance
#le seul ??l??ment stochastique est le mouvement Brownien
# pour les variables antith??tiques: le mouvement brownien est tr??s stable



############# SIMULATION DU MB ###############
#dans un premier temps on fait une construction random walk 
k<-20 #longueur du MB
n<-1000 #nombre de simulations
W <- rnorm(n=n*k,mean = 0,sd=1/k)
W<-matrix(W,nrow=n,ncol=k)
W<- t(apply(W,1,cumsum))
dim(W)

plot(W[1,],type='l',xlim = c(0,20),ylim = c(-3,3))
for (i in 2:n){
  lines(W[i,])
}

plot(colMeans(W),type='l')


############# MC STD & ANTITHETIQUE ###############
# concernant le MB, on sait que -W a m??me loi que W. 

monteCarlo<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,mu,antithet=T, randomized = F){
  start<-Sys.time()
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non sp??cifi??, le drift et le taux d'int??r??t sont ??gaux
  t <- (1:k)/k
  W <- rnorm(n=N*(k-1), mean=0,sd=1/sqrt(k))
  W <- matrix(W,nrow = N,ncol = k-1)
  W<-t(apply(W,1,cumsum))
  W1 <- c(rep(0,N),W)
  W1<-matrix(W1,nrow = N,ncol = k)
  S <- rep(K*exp((mu-sigma**2/2)*t),N) 
  S <-matrix(S,nrow=N,ncol=k)
  S1 <-S*exp(sigma*W1)
  if (antithet){
    W2 <- -W1
    S2 <-S*exp(sigma*W2)
    C <- 0.5*exp(-r*T)*(pmax(rowMeans(S1)-K,0)+pmax(rowMeans(S2)-K,0))
  }
  else{
    C <- exp(-r*T)*pmax(rowMeans(S1)-K,0)
  }
  #beta <- cov(W1[,k],C)/var(W1[,k])
  #C <- C-beta*W1[,k]
  #return (list('mean'=mean(C),'sd'=sd(C)))
  end<-Sys.time()
  return (c(mean(C),sd(C),end-start))
}

monteCarlo(antithet=F)
monteCarlo(antithet = T)

############# VARIABLE DE CONTROLE ###############
# id??e (?? voir): prendre le Black Scholes comme variable de contr??le, ou bien juste le BM 

#on change la fa??on de simuler parce qu'on va avoir besoin 
#d'estimateurs de cov et de variances

control<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,antithet=F,mu,randomized=F){
  start<-Sys.time()
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non sp??cifi??, le drift et le taux d'int??r??t sont ??gaux
  t <- (1:k)/k
  W <- rnorm(n=N*(k-1), mean=0,sd=sqrt(1/k))
  W <- matrix(W,nrow = N,ncol = k-1)
  W<-t(apply(W,1,cumsum))
  W1 <- c(rep(0,N),W)
  W1<-matrix(W1,nrow = N,ncol = k)
  S <- rep(K*exp((mu-sigma**2/2)*t),N) 
  S <-matrix(S,nrow=N,ncol=k)
  S1 <-S*exp(sigma*W1)
  if (antithet){
    W2 <- -W1
    S2 <-S*exp(sigma*W2)
    C <- 0.5*exp(-r*T)*(pmax(rowMeans(S1)-K,0)+pmax(rowMeans(S2)-K,0))
  }
  else{
    C <- exp(-r*T)*pmax(rowMeans(S1)-K,0)
  }
  beta <- cov(W1[,k],C)/var(W1[,k])
  C <- C-beta*W1[,k]
  #return (list('mean'=mean(C),'sd'=sd(C)))
  end<-Sys.time()
  return (c(mean(C),sd(C),end-start))
}

control(antithet=T)



############# QUASI MC ###############
#si on devait coder un g??n??rateur QMC nous m??mes on pourrait faire du Box Muller
# avec des uniformes en treillis
#mais on va utiliser la librairie d??j?? faite

# on va utiliser les fonctions halton, torus et sobol pour g??n??rer des gaussiennes
# pour faire le RQMC on utilise l'option scrambling de la fonction sobol

qmc<-function(N=1000,T=1,r=0.05,k=20,K=5,sigma=0.3,mu,funcStr='halton',antithet=T,randomized=F){
  start<-Sys.time()
  if (missing(mu)){
    mu<-r
  } #on suppose que, si non sp??cifi??, le drift et le taux d'int??r??t sont ??gaux
  t <- (1:k)/k
  if (randomized){
    W <- sqrt(1/k)*sobol(normal=T,n = N,dim = k-1,scrambling = 3)
  }
  else{
    switch (funcStr,
            halton={
              W <- sqrt(1/k)*halton(normal=T,n = N,dim = k-1)
            },
            torus={
              W <- sqrt(1/k)*torus(normal=T,n = N,dim = k-1)
            },
            sobol={
              W <- sqrt(1/k)*sobol(normal=T,n = N,dim = k-1)
            },
            print('qmc function not specified')
    )
  }
  
  W<-t(apply(W,1,cumsum))
  W1 <- c(rep(0,N),W)
  W1<-matrix(W1,nrow = N,ncol = k)
  S <- rep(K*exp((mu-sigma**2/2)*t),N) 
  S <-matrix(S,nrow=N,ncol=k)
  S1 <-S*exp(sigma*W1)
  if (antithet){
    W2 <- -W1
    S2 <-S*exp(sigma*W2)
    C <- 0.5*exp(-r*T)*(pmax(rowMeans(S1)-K,0)+pmax(rowMeans(S2)-K,0))
  }
  else{
    C <- exp(-r*T)*pmax(rowMeans(S1)-K,0)
  }
  beta <- cov(W1[,k],C)/var(W1[,k])
  C <- C-beta*W1[,k]
  #return (list('mean'=mean(C),'sd'=sd(C)))
  end<-Sys.time()
  return (c(mean(C),sd(C),end-start))
}
qmc(funcStr = 'sobol')




# par contre les m??thodes antith??tiques donnent de meilleures variances
# on pourrait rajouter l'option antith??tique ou non ?? la m??thode contr??le 

########################################
############# QUESTION 2 ###############
########################################

multiLevel<-function(N,epsilon,T=1,r=0.05,sigma=0.3,K=5,mu,k=20){
  start<-Sys.time()
  if (missing(mu)){mu<-r}
  
  dates<-T/k*1:k
  h0<-T/k #on fait k steps d??s le d??but
  L<- as.integer(-log(epsilon/T)/log(2)) #nombre de niveaux 
  M0 <- as.integer(-log(epsilon)/epsilon**2) #nombre de simulations au niveau 0
  
  Sg <- matrix(K,nrow=N,ncol=M0)
  MSg<- 1/k*Sg
  for (i in 1:k){
    g <- rnorm(n = N*M0)
    g <- matrix(g,nrow=N,ncol=M0)
    Sg <- Sg*(1+mu*h0)+sigma*sqrt(h0)*g
    MSg <- MSg + 1/k*Sg
  }
  estim <- rowMeans(pmax(MSg-K,0)) 
  
  for (l in 1:L){
    
    Ml <- as.integer(M0/2**l) #nombre de simulations au niveau l 
    hl<-h0/(2**l) #pas de discr??tisation du niveau l 
    sig <- sigma*sqrt(hl) #ecart type correspondant
    Sf <- matrix(K,nrow=N,ncol=Ml) #schema d'euler de pas fin
    Sg <- matrix(K,nrow=N,ncol=Ml) #schema d'euler de pas grossier
    
    MSf<-(1/k)*Sf #la moyenne arithm??tique du processus, pas fin 
    MSg<-(1/k)*Sg #la moyenne arithm??tique du processus, pas grossier
    
    for (i in 1:2**(l-1)){
      
      g1 <- rnorm(n=N*Ml)
      g2 <- rnorm(n=N*Ml)
      g1<-matrix(g1,nrow = N,ncol = Ml)
      g2<-matrix(g2,nrow = N,ncol = Ml)
      
      #evolution de deux schemas: pas fin et pas grossier
      Sf <- Sf*exp((mu-sigma**2/2)*hl+sig*g1)
      Sf<-Sf*(1+mu*hl)+sig*g1
      if (i*T*2**l/k %in% dates){MSf<-MSf+(1/k)*Sf}
      Sf <- Sf*exp((mu-sigma**2/2)*hl+sig*g2)
      if ((i+1)*T*2**l/k %in% dates){MSf<-MSf+(1/k)*Sf}
      Sg <- Sg*exp(2*(mu-sigma**2/2)*hl+sig*(g1+g2))
      if (i*T*2**(l-1)/k %in% dates){MSg<-MSg+(1/k)*Sg}
    }
    Cf<-pmax(MSf-K,0)
    Cg<-pmax(MSg-K,0)
    estim<-estim+rowMeans(Cf-Cg)
  }
  e<-exp(-r*T)*estim
  end<-Sys.time()
  return (c(mean(e),sd(e),end-start))
}

multiLevel(N=1000,epsilon = 0.05)




########################################
############# QUESTION 3 ###############
########################################
library(gtools)

multiLevelQMC<-function(N,epsilon,T=1,r=0.05,sigma=0.3,K=5,mu,k=20,func){
  start<-Sys.time()
  if (missing(mu)){mu<-r}
  
  dates<-T/k*1:k
  h0<-T/k #on fait k steps d??s le d??but
  L<- as.integer(-log(epsilon/T)/log(2)) #nombre de niveaux 
  M0 <- as.integer(-log(epsilon)/epsilon**2) #nombre de simulations au niveau 0
  
  Sg <- matrix(K,nrow=N,ncol=M0)
  MSg<- 1/k*Sg
  for (i in 1:k){
    g<-func(normal = T,n=N*M0)
    g<-permute(g)
    g <- matrix(g,nrow=N,ncol=M0)
    Sg <- Sg*(1+mu*h0)+sigma*sqrt(h0)*g
    MSg <- MSg + 1/k*Sg
  }
  estim <- rowMeans(pmax(MSg-K,0)) 
  mean(estim)
  length(estim)
  sd(estim)
  
  for (l in 1:L){
    
    Ml <- as.integer(M0/2**l) #nombre de simulations au niveau l 
    hl<-h0/(2**l) #pas de discr??tisation du niveau l 
    sig <- sigma*sqrt(hl) #ecart type correspondant
    Sf <- matrix(K,nrow=N,ncol=Ml) #schema d'euler de pas fin
    Sg <- matrix(K,nrow=N,ncol=Ml) #schema d'euler de pas grossier
    
    MSf<-(1/k)*Sf #la moyenne arithm??tique du processus, pas fin 
    MSg<-(1/k)*Sg #la moyenne arithm??tique du processus, pas grossier
    
    for (i in 1:2**(l-1)){
      
      g1 <- func(normal=T,n=N*Ml)
      g2 <- func(normal=T,n=N*Ml)
      g1 <- matrix(permute(g1),nrow=N,ncol=Ml)
      g2 <- matrix(permute(g2),nrow=N,ncol=Ml)
      
      
      #evolution de deux schemas: pas fin et pas grossier
      Sf <- Sf*exp((mu-sigma**2/2)*hl+sig*g1)
      Sf<-Sf*(1+mu*hl)+sig*g1
      if (i*T*2**l/k %in% dates){MSf<-MSf+(1/k)*Sf}
      Sf <- Sf*exp((mu-sigma**2/2)*hl+sig*g2)
      if ((i+1)*T*2**l/k %in% dates){MSf<-MSf+(1/k)*Sf}
      Sg <- Sg*exp(2*(mu-sigma**2/2)*hl+sig*(g1+g2))
      if (i*T*2**(l-1)/k %in% dates){MSg<-MSg+(1/k)*Sg}
    }
    Cf<-pmax(MSf-K,0)
    Cg<-pmax(MSg-K,0)
    estim<-estim+rowMeans(Cf-Cg)
  }
  e<-exp(-r*T)*estim
  end<-Sys.time()
  return(c(mean(e),sd(e),end-start))
}

multiLevelQMC(N=10,epsilon = 0.05,func=torus)

########################################
#############  SYNTHESE  ###############
########################################

results<-data.frame(matrix(nrow=10,ncol=3),
                    row.names = c('mc','antithetic','ctrl','ctrl antithetic','qmc',
                                  'qmc antithetic','rqmc','rqmc antithetic',
                                  'multilevel','qmc multilevel'))

colnames(results)<-c('moyenne','ecart-type','temps de calcul')
results[1,]<-monteCarlo(antithet=F,N=10000)
results[2,]<-monteCarlo(antithet = T,N=10000)
results[3,]<-control(N=10000)
results[4,]<-control(antithet=T,N=10000)
results[5,]<-qmc(antithet = F,N=10000)
results[6,]<-qmc(antithet = T,N=10000)
results[7,]<-qmc(antithet = F,randomized = T,N=10000)
results[8,]<-qmc(antithet = T,randomized = T,N=10000)
results[9,]<-multiLevel(N=100,epsilon = 0.1)
results[10,]<-multiLevel(N=100,epsilon = 0.1)
results
stargazer(results)

######### PLOT #############

get_plot <- function(method, ...){
  x <- c(10,50,100,500,1000,5000,10000)
  mean <- vector("numeric", length(x))
  sdv <- vector("numeric", length(x))
  ci <- vector("numeric", length(x))
  for (i in 1:length(x)){
    res <- method(N=x[i], ...)
    mean[i] <- res[1]
    sdv[i] <- res[2]
    ci[i] <- 1.96 * sdv[i]/sqrt(x[i])
  }
  tib <- tibble(x, mean, sdv, ci)
  ggplot(data = tib, mapping = aes(x = x, y = mean)) +
    geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci),
                  colour="red", width=.05, size = .7) +
    geom_point(size = 2.2) +
    scale_x_continuous(trans = 'log10') +
    labs(x = 'nombre de simulations', 'mean') + 
    theme_bw()
}

get_plot(method = monteCarlo, antithet = F)


par(mfrow=c(1,2))
# Plot QMC
get_plot(method=qmc,antithet =F,randomized=F)
get_plot(method=qmc,antithet =T,randomized=F)
# Plot RQMC
get_plot(method=qmc,antithet =F,randomized=T)
get_plot(method=qmc,antithet =T,randomized=T)
# Plot Mutilevel
get_plot(method=multiLevel, epsilon = 0.05)



get_plot(method=qmc,antithet=T)

#get_plot(method=control,antithet=F)
