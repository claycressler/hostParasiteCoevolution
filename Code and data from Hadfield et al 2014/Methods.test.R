############################################################
## generates simulated data and tests permutation schemes ##
############################################################

library(MCMCglmm)
library(igraph)
library(asreml)
library(gdata)
source("~/Work/DoublePhylo/R/hommola.R")

results.filepath<-"~/Work/DoublePhylo/Data/Intermediate"  # file path where results are to be stored 

n.sim<-10            # number of simulated data sets

n<-10               # number of hosts
m<-10               # number of parasites

int<-0              # intercpet (0 = 50% incidence)
Vp<-1e-8            # parasite phylogenetic main effect variance
Vh<-1e-8            # host phylogenetic main effect variance
VpI<-1e-8           # parasite evolutionary interaction effect variance
VhI<-1e-8           # host evolutionary interaction effect variance
Vhp<-4              # coevolutionary interaction effect variance

incidence<-TRUE     # should the data be binary or gaussian

run.parafit<-TRUE   # should Legendre's statistic be calculated
run.hommola<-TRUE   # should Hommola's statistic be calculated
run.ives<-TRUE      # should Ives's statistic be calculated
run.module<-TRUE    # should a Krasnov-like module anlaysis be run

run.asreml<-TRUE    # should asreml be tested too
run.MCMC<-TRUE      # should MCMCglmm be tested too

nitt=13000*5        # number of MCMC iterations
thin=10*5           # thinning interval
burnin=3000*5       # MCMC burn-in 

nitt=130
thin=1
burnin=30


leg.sampling<-TRUE  # should Legendre permutation be used or Hommola permutation


HomS<-rep(NA, n.sim)
LegS<-rep(NA, n.sim)
IvesS<-rep(NA, n.sim)
ASRS<-matrix(NA, n.sim,3)
MCMS<-matrix(NA, n.sim,21)
MSP<-rep(NA, n.sim)  
MSH<-rep(NA, n.sim)               

# empty vectors for storing results

if(run.parafit==FALSE){
LegS<-NULL
}
if(run.hommola==FALSE){
HomS<-NULL
}
if(run.ives==FALSE){
IvesS<-NULL
}

c2 <- (16 * sqrt(3)/(15 * pi))^2  

# Term used when rescaling a variance to what would have been observed had the 
# (unidentified) residual variance been set to zero in a binary model

Vt<-Vp+Vh+VpI+VhI+Vhp  # Total variance

Ir<-Vt/(Vt+pi^2/3)     # ICC for all effects

if(incidence){
Ve<-0
}else{
Ve<-Vt*(1-Ir)/Ir       # Ve picked for Guassian data so ICC is the same as what it would be for binary data
}

if(incidence){

  # parameter expanded priors for the variances, with models varying in the number of variance parameters

  prior0<-list(B=list(V=(5+pi^2/3), mu=0), R=list(V=1, fix=1))

  prior1<-list(B=list(V=(5+pi^2/3), mu=0), R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

  prior2<-list(B=list(V=(5+pi^2/3), mu=0), R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

  prior3<-list(B=list(V=(5+pi^2/3), mu=0), R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

}else{

  prior0<-list(R=list(V=1, nu=0))

  prior1<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

  prior2<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

  prior3<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

}

for(i in 1:n.sim){

  host.tree<-rcoal(n)
  par.tree<-rcoal(m)
  # simulate trees

  ht<-vcv(host.tree, corr=T)
  pt<-vcv(par.tree, corr=T)
  # get A matrices

  Ch<-chol(ht)
  Cp<-chol(pt)
  # get cholesky decomposition of A matrices

  dat<-data.frame(y=rep(NA, n*m),In=rep(NA, n*m), host=rep(host.tree$tip.label, each=m), par=rep(par.tree$tip.label,n))
  dat$host.par<-paste(dat$host, dat$par, sep=".")
  dat$hostI.par<-paste(dat$host, dat$par, sep=".")
  dat$host.parI<-paste(dat$host, dat$par, sep=".")
  # get data-frame with empty data vectors

  h.e<-rbv(host.tree, Vh, nodes="TIPS")                       # simulate host phylogenetic main effects
  p.e<-rbv(par.tree, Vp, nodes="TIPS")                        # simulate parasite phylogenetic main effects

  hI.e<-c(t(sqrt(VhI)*t(Ch)%*%matrix(rnorm(n*m),n,m)))        # simulate host evolutionary interaction effects
  pI.e<-c(t(sqrt(VpI)*matrix(rnorm(n*m),n,m)%*%Cp))           # simulate parasite evolutionary interaction effects

  hp.e<-c(t(sqrt(Vhp)*t(Ch)%*%matrix(rnorm(n*m),n,m)%*%Cp))   # simulate coevolutionary interaction effects

  dat$y<-int+h.e[rep(1:n, each=m)]+p.e[rep(1:m, n)]+rnorm(m*n, 0, sqrt(Ve))+hp.e+hI.e+pI.e
  # combine effects to get latent variable

  if(incidence){
    dat$In<-rbinom(n*m,1,plogis(dat$y))
    # simulate binary data
  }

  In<-t(matrix(dat$In,m,n))
  # get incidence matrix

  rownames(In)<-rownames(ht)
  colnames(In)<-colnames(pt)
    
  hc<-hommola(In, ht,pt) 
  # get Hommola statistic

  htp<-pcoa(1-ht)$vectors
  ptp<-pcoa(1-pt)$vectors
  # get principal coordinates of phylogenetic distance matrices

  hte<-t(t(eigen(solve(ht))$vectors)*sqrt(eigen(solve(ht))$values))
  pte<-t(t(eigen(solve(pt))$vectors)*sqrt(eigen(solve(pt))$values))
  # get unnormalised eigenvectors of A matrices

  D<-t(htp)%*%In%*%ptp
  lc<-sum(diag(t(D)%*%D))
  # get Legendre's Parafit statistic

  iD<-t(hte)%*%(In-mean(In))%*%pte
  ic<-sum(diag(t(iD)%*%iD))
  # get Ives' MSEb statistic

  resL<-1:1000
  resH<-1:1000
  resI<-1:1000
  resM<-1:1000

  # vectors for storing statistics from permuted data
  if(run.parafit | run.hommola | run.ives){

    for(j in 1:1000){

      if(leg.sampling){
        In2<-apply(In, 2, sample)       # Legendre sampling
      }else{
        a<-sample(1:nrow(In))           # Hommola sampling
        b<-sample(1:ncol(In))
        In2<-In[a,b]
      }

      if(run.parafit){                  # Legendre statistic
        D<-t(htp)%*%In2%*%ptp          
        resL[j]<-sum(diag(t(D)%*%D))  
      }
      if(run.ives){                     # Ives statistic
        iD<-t(hte)%*%(In2-mean(In2))%*%pte
        resI[j]<-sum(diag(t(iD)%*%iD))
      }
      if(run.hommola){                  # Hommola statistic
        resH[j]<-hommola(In2, ht,pt)
      }
      if(j%%100==0){    
        print(paste(j/10, "%", sep="")) # print progress
      }
    }
    HomS[i]<-sum(resH>hc)/1000          
    LegS[i]<-sum(resL>lc)/1000
    IvesS[i]<-sum(resI>ic)/1000
    # store proportion of statistics from permuted data that are more extreme than atual statistic
  }

  if(run.MCMC | run.asreml){
    # Generate S^{-1} matrices
    hA<-inverseA(host.tree)$Ainv
    pA<-inverseA(par.tree)$Ainv
    hpA<-as(kronecker(hA,pA), "dgCMatrix")
    hpAS<-as(kronecker(hA, Diagonal(nrow(pA))), "dgCMatrix")  
    hpSA<-as(kronecker(Diagonal(nrow(hA)), pA), "dgCMatrix")  
    rownames(hpA)<-rownames(hpAS)<-rownames(hpSA)<-apply(expand.grid(rownames(pA), rownames(hA)), 1, function(x){paste(x[2],x[1], sep=".")})
  }

  if(run.MCMC){

    m1<-MCMCglmm(In~1, random=~host+par+host.par, ginverse=list(host=hA, par=pA, host.par=hpA), data=dat, prior=prior3, family=c("gaussian", "categorical")[incidence+1], slice=TRUE, nitt=nitt, thin=thin, burnin=burnin)
    # model with main effects and coevolutionary interaction

    m2<-MCMCglmm(In~1, random=~host.par, ginverse=list(host.par=hpA), data=dat, prior=prior1, family=c("gaussian", "categorical")[incidence+1], slice=T, nitt=nitt, thin=thin, burnin=burnin)
   # model with coevolutionary interaction
 
    m3<-MCMCglmm(In~1, random=~hostI.par+host.parI+host.par, ginverse=list(hostI.par=hpSA, host.parI=hpAS, host.par=hpA), data=dat, prior=prior3, family=c("gaussian", "categorical")[incidence+1], slice=TRUE, nitt=nitt, thin=thin, burnin=burnin)
   # model with evolutionary interactions and coevolutionary interaction

   MCMS[i,][1]<-posterior.mode(m1$VCV[,1]/(rowSums(m1$VCV)+(pi^2/3)*incidence))
   MCMS[i,][2:3]<-HPDinterval(m1$VCV[,1]/(rowSums(m1$VCV)+(pi^2/3)*incidence))
   MCMS[i,][4]<-posterior.mode(m1$VCV[,2]/(rowSums(m1$VCV)+(pi^2/3)*incidence))
   MCMS[i,][5:6]<-HPDinterval(m1$VCV[,2]/(rowSums(m1$VCV)+(pi^2/3)*incidence))
   MCMS[i,][7]<-posterior.mode(m1$VCV[,3]/(rowSums(m1$VCV)+(pi^2/3)*incidence))
   MCMS[i,][8:9]<-HPDinterval(m1$VCV[,3]/(rowSums(m1$VCV)+(pi^2/3)*incidence))
   MCMS[i,][10]<-posterior.mode(m2$VCV[,1]/(rowSums(m2$VCV)+(pi^2/3)*incidence))
   MCMS[i,][11:12]<-HPDinterval(m2$VCV[,1]/(rowSums(m2$VCV)+(pi^2/3)*incidence))
   MCMS[i,][13]<-posterior.mode(m3$VCV[,1]/(rowSums(m3$VCV)+(pi^2/3)*incidence))
   MCMS[i,][14:15]<-HPDinterval(m3$VCV[,1]/(rowSums(m3$VCV)+(pi^2/3)*incidence))
   MCMS[i,][16]<-posterior.mode(m3$VCV[,2]/(rowSums(m3$VCV)+(pi^2/3)*incidence))
   MCMS[i,][17:18]<-HPDinterval(m3$VCV[,2]/(rowSums(m3$VCV)+(pi^2/3)*incidence))
   MCMS[i,][19]<-posterior.mode(m3$VCV[,3]/(rowSums(m3$VCV)+(pi^2/3)*incidence))
   MCMS[i,][20:21]<-HPDinterval(m3$VCV[,3]/(rowSums(m3$VCV)+(pi^2/3)*incidence))
   # store ICC's
  } 

  if(run.asreml){
   m1.asreml<-asreml(In~1, random=~giv(host.par,TRUE)+giv(host,TRUE)+giv(par,TRUE), ginverse=list(host=sm2asreml(hA), par=sm2asreml(pA), host.par=sm2asreml(hpA)), data=dat, family=asreml.binomial())
   # model with main effects and coevolutionary interaction

   m2.asreml<-asreml(In~1, random=~giv(host.par,TRUE), ginverse=list(host.par=sm2asreml(hpA)), data=dat, family=asreml.binomial())
   # model with coevolutionary interaction

   m3.asreml<-asreml(In~1, random=~giv(host.par,TRUE)+giv(hostI.par,TRUE)+giv(host.parI,TRUE), ginverse=list(hostI.par=sm2asreml(hpSA), host.parI=sm2asreml(hpAS), host.par=sm2asreml(hpA)), data=dat, family=asreml.binomial())
   # model with evolutionary interactions and coevolutionary interaction

   ASRS[i,]<-c(summary(m1.asreml)$varcomp[,4][1],summary(m2.asreml)$varcomp[,4][1],summary(m3.asreml)$varcomp[,4][1])
   # store z-ratios for coevolutionary interaction term
  }

  if(run.module){

    keep.r<-which(rowSums(In)!=0)
    keep.c<-which(colSums(In)!=0)
    # retain species that have some observed interaction

    if(nrow(In[keep.r,keep.c, drop=FALSE])!=0){
      # if there are some observed interactions!

      G<-graph.incidence(In[keep.r,keep.c, drop=FALSE])
      MU<-walktrap.community(G)
      MU<-membership(MU)
      # assign to modules
 
      lth<-which(lower.tri(ht[keep.r,keep.r]))

      MSH[i]<-cor(outer(MU[1:length(keep.r)], MU[1:length(keep.r)], "==")[lth], c(1-ht[keep.r,keep.r])[lth])
      # correlation between comembership and phylogenetic distance of hosts

      ltp<-lower.tri(pt[keep.c,keep.c])

      MSP[i]<-cor(outer(MU[length(keep.r)+1:length(keep.c)], MU[length(keep.r)+1:length(keep.c)], "==")[ltp], c(1-pt[keep.c,keep.c])[ltp])
      # correlation between comembership and phylogenetic distance of parasites
    } 
  }
  print(i)
}

colnames(MCMS)<-c("h|h.p.hp", "h.lower|h.p.hp", "h.upper|h.p.hp", "p|h.p.hp", "p.lower|h.p.hp", "p.upper|h.p.hp", "hp|h.p.hp", "hp.lower|h.p.hp", "hp.upper|h.p.hp","hp|hp", "hp.lower|hp", "hp.upper|hp", "hIp|hIp.hpI.hp", "hIp.lower|hIp.hpI.hp", "hIp.upper|hIp.hpI.hp", "hpI|hIp.hpI.hp", "hpI.lower|hIp.hpI.hp", "hpI.upper|hIp.hpI.hp", "hp|hIp.hpI.hp", "hp.lower|hIp.hpI.hp", "hp.upper|hIp.hpI.hp")

simresults<-cbind(MSH,MSP,HomS, LegS, IvesS, MCMS)

save(simresults, file=paste(paste(data.filepath, "simresults", sep=.Platform$file.sep), c("H","L")[leg.sampling+1], "Rdata", sep=".")


