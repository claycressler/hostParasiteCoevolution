#####################################################
## generates simulated data and tests PQL adequacy ##
#####################################################

library(asreml)
library(MCMCglmm)
library(gdata)
library(igraph)

n.sim<-1         # number of simulations

data.filepath<-"~/Work/DoublePhylo/Data/Raw"              # file path to data files 
results.filepath<-"~/Work/DoublePhylo/Data/Intermediate"  # file path where results are to be stored 

ptree<-"krasnov"                 # which algorithm to use to make parasite tree ultrameteric
                                          # "krasnov" (that used in Krasnov 2012: Mesquite)
                                          # "grafen"  (Grafen's method - uses apes compute.brlen)
                                          # "chronpl" (Sanderson's method - uses chronopl)

#################
## Phylogenies ##
#################

  if(ptree=="krasnov"){
    parasite.tree<-read.tree(paste(data.filepath, "FleaTreeK.tre", sep=.Platform$file.sep))
  }else{
    parasite.tree<-read.tree(paste(data.filepath, "FleaTree.tre", sep=.Platform$file.sep))
    if(ptree=="grafen"){
      parasite.tree<-compute.brlen(parasite.tree)
    }
    if(ptree=="chronopl"){
      parasite.tree<-chronopl(parasite.tree)
    }
  }

  host.tree<-read.tree(paste(data.filepath, "HostTree.tre", sep=.Platform$file.sep))

  host.tree$node.label<-NULL
  parasite.tree$node.label<-NULL

###############################
## Form S^{-1} for each term ##
###############################

  parasiteA<-inverseA(parasite.tree)$Ainv                                        # parasite main effect
  hostA<-inverseA(host.tree)$Ainv                                                # host main effect
  host.parasiteA<-as(kronecker(hostA, parasiteA), "dgCMatrix")                   # coevolutionary effect
  host.parasiteAS<-as(kronecker(hostA, Diagonal(nrow(parasiteA))), "dgCMatrix")  # host evolutionary effect
  host.parasiteSA<-as(kronecker(Diagonal(nrow(hostA)), parasiteA), "dgCMatrix")  # parasite evolutionary effect

  rownames(host.parasiteA)<-apply(expand.grid(rownames(parasiteA), rownames(hostA)), 1, function(x){paste(x[2],x[1], sep=".")})
  rownames(host.parasiteAS)<-rownames(host.parasiteSA)<-rownames(host.parasiteA)

##########
## Data ##
##########

if(!"ndat"%in%ls()){   # read in and format data if not there

dat<-read.csv(paste(data.filepath, "DoublePhyloData.csv", sep=.Platform$file.sep))

################
## Clean Data ##
################

while(sum(duplicated(paste(dat$Host.species, dat$Region)))>0){
dat[,-(1:2)][which(duplicated(paste(dat$Host.species, dat$Region)))[1]-1,]<-dat[,-(1:2)][which(duplicated(paste(dat$Host.species, dat$Region)))[1]-1,]+dat[,-(1:2)][which(duplicated(paste(dat$Host.species, dat$Region)))[1],]
res<-any(dat[,1:2][which(duplicated(paste(dat$Host.species, dat$Region)))[1]-1,]!=dat[,1:2][which(duplicated(paste(dat$Host.species, dat$Region)))[1],])
if(res){stop("not consecuitive")}
dat<-dat[-which(duplicated(paste(dat$Host.species, dat$Region)))[1],]
}
# data contains 3 duplicate species/Region counts: amalgamate into single records

#################
## Format Data ##
#################

ndat<-reshape(dat, direction="long", varying=list(names(dat)[-c(1:3)]), times=names(dat)[-c(1:3)])
# turn into long format

# rename columns: 
names(ndat)[which(names(ndat)=="time")]<-"Parasite.species"
names(ndat)[which(names(ndat)=="Amalaraeus.arvicolae")]<-"counts"

# duplicate columns with new names: 
ndat$Parasite.species<-gsub("\\.", "_", ndat$Parasite.species)                    # phylogenetic main effect for parasites
ndat$Host.species<-gsub(" ", "_", ndat$Host.species)                              # phylogenetic main effect for hosts
ndat$Parasite.species.ide<-ndat$Parasite.species                                  # non-phylogenetic main effect for parasites
ndat$Host.species.ide<-ndat$Host.species                                          # non-phylogenetic main effect for hosts
ndat$Host.Parasite<-paste(ndat$Host.species, ndat$Parasite.species, sep=".")      # phylogentic coevolutionary effect
ndat$Host.Parasite.ide<-paste(ndat$Host.species, ndat$Parasite.species, sep=".")  # non-phylogentic interaction effect
ndat$Host.Parasite.ide2<-paste(ndat$Host.species, ndat$Parasite.species, sep=".") # phylogentic host evolutionary effect
ndat$Host.Parasite.ide3<-paste(ndat$Host.species, ndat$Parasite.species, sep=".") # phylogentic parasite evolutionary effect
ndat$present<-as.numeric(ndat$counts>0)                                           # Turn counts into incidence

# remove structural-zeros for parasites (i.e. zero counts for parasites on hosts in regions where the parasite has 
# not been found) the structural zeros for hosts are removed already (although there are some real zero's):

pr<-table(ndat$counts==0, paste(ndat$Parasite.species, ndat$Region))
mis<-colnames(pr)[which(pr[1,]==0)]
mis.p<-which(paste(ndat$Parasite.species, ndat$Region)%in%mis)
ndat<-ndat[-mis.p,]

# get number of host species per region 
nhosts<-tapply(ndat$Host.species, ndat$Region, function(x){length(unique(x))})
ndat$nhosts<-nhosts[match(ndat$Region, names(nhosts))]

# get number of parasites per host species/region
par.counts<-tapply(ndat$counts, paste(ndat$Parasite.species,ndat$Region), sum)
ndat$no.parasites.sampled<-par.counts[match(paste(ndat$Parasite.species,ndat$Region), names(par.counts))]

# create factors for unique region/host, region/parasite and region/host/parasite
# combinations for the anlayses ignoring between-region information:

ndat$Parasite.species.r<-paste(ndat$Region, ndat$Parasite.species, sep=".")
ndat$Host.species.r<-paste(ndat$Region, ndat$Host.species, sep=".") 
ndat$Host.Parasite.r<-paste(ndat$Region, ndat$Host.Parasite, sep=".")          
ndat$Host.Parasite.ide2.r<-paste(ndat$Region, ndat$Host.Parasite, sep=".")     
ndat$Host.Parasite.ide3.r<-paste(ndat$Region, ndat$Host.Parasite, sep=".")     
ndat$Host.species.ide.r<-paste(ndat$Region, ndat$Host.species.ide, sep=".") 
ndat$Parasite.species.ide.r<-paste(ndat$Region, ndat$Parasite.species.ide, sep=".")      
}

#####################
#### simulation #####
#####################

if(any(ls()=="mI.MCMCa")==FALSE){    # read in most recent main MCMC model

    mfiles<-list.files(results.filepath)
    mfiles<-mfiles[which(grepl("mI\\.MCMCa\\.", mfiles) & grepl(ptree, mfiles))]
    mfiles<-mfiles[which.max(as.Date(substr(mfiles, nchar(mfiles)-15, nchar(mfiles)-6), "%d-%m-%Y"))]

    load(paste(results.filepath, .Platform$file.sep, mfiles,".Rdata", sep=""))
}

if(any(ls()=="mI.ASReml")==FALSE){  # read in most recent main MCMC model
  load(paste(results.filepath, .Platform$file.sep, "mI.ASRemla.", ptree, ".Rdata",  sep=""))
}

host.name<-host.tree$tip.label
parasite.name<-parasite.tree$tip.label

n<-length(host.name)
m<-length(parasite.name)

host.par.name<-paste(rep(host.name, each=m), rep(parasite.name, n), sep=".")

ht<-vcv(host.tree, corr=T)
pt<-vcv(parasite.tree, corr=T)
# get A matrices

Ch<-chol(ht)
Cp<-chol(pt)
# get Cholesky decomposition of A matrices

res<-matrix(NA, n.sim,10)
resb<-matrix(NA, n.sim,11)


X<-mI.MCMC$X
Z<-mI.MCMC$Z[,which(grepl("Host.Parasite", colnames(mI.MCMC$Z)) & grepl("Host.Parasite.ide", colnames(mI.MCMC$Z))==FALSE & grepl("Node", colnames(mI.MCMC$Z))==FALSE)]
colnames(Z)<-gsub("Host.Parasite.", "", colnames(Z))
Z<-Z[,match(colnames(Z), host.par.name)]
# Fixed effect design matrix and Random effect design matrix for coevolutionary interaction

e.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"units"]))                  # posterior mode of residual variation (square rooted)
r.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Region"]))                 # posterior mode of regional variation (square rooted)
pi.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Parasite.species.ide"]))  # posterior mode of non-phylogenetic parasite variation (square rooted)
hi.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Host.species.ide"]))      # posterior mode of non-phylogenetic host variation (square rooted)
hipi.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Host.Parasite.ide"]))   # posterior mode of non-phylogenetic host-parasite intreraction variation (square rooted)
p.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Parasite.species"]))       # posterior mode of parasite phylogenetic main effect variation (square rooted)
h.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Host.species"]))           # posterior mode of host phylogenetic main effect variation (square rooted)
hp.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Host.Parasite"]))         # posterior mode of coeveolutioanry interaction effect variation (square rooted)
hpi.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Host.Parasite.ide2"]))   # posterior mode of host eveolutioanry interaction effect variation (square rooted)
hip.sd<-sqrt(posterior.mode(mI.MCMC$VCV[,"Host.Parasite.ide3"]))   # posterior mode of parasite eveolutioanry interaction effect variation (square rooted)

ndat$simy<-rep(NA, nrow(ndat))
ndat$simi<-rep(NA, nrow(ndat))
# empty vectors for the simulated latent variable (y) and the incidence data (i)

msimASReml.st<-asreml(simy~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+giv(Parasite.species, var=T)+giv(Host.species, var=T) +Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+Host.Parasite.ide+giv(Host.species, var=T):Parasite.species.ide+giv(Parasite.species, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species=sm2asreml(parasiteA), Host.species=sm2asreml(hostA), Host.Parasite=sm2asreml(host.parasiteA)), workspace = 8e+07, maxiter = 30, start.values=TRUE)

# getting starting values for Gaussian model

G.param1<-msimASReml.st$G.param
for(i in 1:9){
G.param1[[i]][[length(G.param1[[i]])]]$initial[1]<-posterior.mode(mI.MCMC$VCV)[i]
}

# give MCMC estimates from data-fitted model as starting values (reduces convergence time)

msimASReml2.st<-asreml(simi~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+Host.Parasite.ide+giv(Host.species, var=T):Parasite.species.ide+giv(Parasite.species, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species=sm2asreml(parasiteA), Host.species=sm2asreml(hostA), Host.Parasite=sm2asreml(host.parasiteA)), workspace = 8e+07, maxiter = 30, family = asreml.binomial(), start.values=TRUE)

# getting starting values for Binary model

G.param2<-msimASReml2.st$G.param
for(i in 1:9){
G.param2[[i]][[length(G.param2[[i]])]]$initial[1]<-summary(mI.ASRemla)$varcomp[,2][i]
}

# give ASReml estimates from data-fitted model as starting values (reduces convergence time)

for(i in 1:n.sim){

  r.e<-rnorm(length(unique(ndat$Region)), 0, r.sd)       # simulate regional effects
  pi.e<-rnorm(m, 0, pi.sd)                               # simulate non-phylogenetic parasite main effects
  hi.e<-rnorm(n, 0, hi.sd)                               # simulate non-phylogenetic host main effects
  hipi.e<-rnorm(n*m, 0, hipi.sd)                         # simulate non-phylogenetic host-parasite interaction effects

  p.e<-rbv(parasite.tree, p.sd^2, nodes="TIPS")          # simulate phylogenetic parasite main effects
  h.e<-rbv(host.tree, h.sd^2, nodes="TIPS")              # simulate phylogenetic parasite host effects

  hpi.e<-c(hpi.sd*matrix(rnorm(n*m), m,n)%*%Ch)          # simulate host evolutionary interaction effects
  hip.e<-c(hip.sd*t(Cp)%*%matrix(rnorm(n*m), m,n))       # simulate parasite evolutionary interaction effects
  hp.e<-c(hp.sd*t(Cp)%*%matrix(rnorm(n*m), m,n)%*%Ch)    # simulate coevolutionary interaction effects

  ndat$simy<-rnorm(nrow(ndat), 0, e.sd)                  # simulate residual effects

  ndat$simy<-ndat$simy+r.e[match(ndat$Region, unique(ndat$Region))]
  ndat$simy<-ndat$simy+pi.e[match(ndat$Parasite.species.ide, parasite.name)]
  ndat$simy<-ndat$simy+hi.e[match(ndat$Host.species.ide, host.name)]
  ndat$simy<-ndat$simy+hipi.e[match(ndat$Host.Parasite.ide, host.par.name)]
  ndat$simy<-ndat$simy+p.e[match(ndat$Parasite.species, parasite.name)]
  ndat$simy<-ndat$simy+h.e[match(ndat$Host.species, host.name)]
  ndat$simy<-ndat$simy+hp.e[match(ndat$Host.Parasite, host.par.name)]
  ndat$simy<-ndat$simy+hpi.e[match(ndat$Host.Parasite, host.par.name)]
  ndat$simy<-ndat$simy+hip.e[match(ndat$Host.Parasite, host.par.name)]
  # add all effects together in the appropriate way

  ndat$simy<-ndat$simy+posterior.mode(mI.MCMC$Sol[,1])
  ndat$simy<-ndat$simy+posterior.mode(mI.MCMC$Sol[,2])*log(ndat$no.hosts.sampled)
  ndat$simy<-ndat$simy+posterior.mode(mI.MCMC$Sol[,3])*log(ndat$no.parasites.sampled)
  # add on the fixed effects 

  ndat$simi<-rbinom(nrow(ndat), 1, plogis(ndat$simy))
  # simulate the binary data from the latent variable

  msimASReml<-asreml(simy~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+Host.Parasite.ide+giv(Host.species, var=T):Parasite.species.ide+giv(Parasite.species, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species=sm2asreml(parasiteA), Host.species=sm2asreml(hostA), Host.Parasite=sm2asreml(host.parasiteA)), workspace = 9e+07, maxiter = 30, G.param=G.param1)
  # fit Gaussian model

  res[i,]<-summary(msimASReml)$varcomp[,2]
  # store results 

  msimASReml2<-asreml(simi~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+Host.Parasite.ide+giv(Host.species, var=T):Parasite.species.ide+giv(Parasite.species, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species=sm2asreml(parasiteA), Host.species=sm2asreml(hostA), Host.Parasite=sm2asreml(host.parasiteA)), workspace = 9e+07, maxiter = 30, family = asreml.binomial(), G.param=G.param2)
  # fit Binary model

  resb[i,][1:10]<-summary(msimASReml2)$varcomp[,2]
  # store results 

  u0<-1/(1+exp(-(X%*%rev(msimASReml2$coefficients$fixed))))
  W0<-Diagonal(length(u0), (u0*(1-u0))@x)
  W1<-Diagonal(length(u0), (u0*(1-u0)*(1-2*u0))@x)
  W2<-Diagonal(length(u0), (u0*(1-u0)*(6*(u0^2)-6*u0+1))@x)

  Q<-diag(t(Z)%*%W0%*%Z)
  Qtilde<-matrix(Q, nrow(ht), nrow(pt))
  a<-t(Q)%*%c((ht*ht)%*%Qtilde%*%(pt*pt))/2

  q<-Z%*%Matrix(kronecker(diag(pt), diag(ht)))
  b<-(t(q)%*%W2%*%q-t(q)%*%W1%*%X%*%solve(t(X)%*%W0%*%X)%*%(t(X)%*%W1%*%q))/4

  # get correction factor (see supplemnetary material)
  resb[i,][11]<-(a/(a+b))@x

}

trait<-c(rep("gaussian", length(res)), rep("binary", length(resb)))
component<-c(rep(colnames(mI.MCMC$VCV), each=nrow(res)),rep(c(colnames(mI.MCMC$VCV), "Correction"), each=nrow(resb)))
asremlsim<-data.frame(trait=trait, component=component, estimate=c(c(res), c(resb)))
#save(asremlsim, file=paste(paste(data.filepath, "asremlsim.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))
# write results to file







