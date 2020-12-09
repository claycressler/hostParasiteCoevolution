library(asreml)
library(MCMCglmm)
library(gdata)
library(igraph)

incidence<-FALSE                 # run models for incidence of parasites on hosts
counts<-FALSE                    # run models for counts of parasites on hosts
abundance<-FALSE                 # run models for host abundance and parasite abundance 
run.MCMC<-FALSE                  # fit models using MCMC 
run.ASReml<-FALSE                # fit models using asreml
run.parafit<-FALSE               # perform parafit, Hommola and `Ives' tests for incidence data
run.module<-FALSE                # perform a krasnov-like test for incidence data
run.regional.ASReml<-FALSE       # fit incidence models for each region using asreml 
run.regional.parafit<-FALSE      # perform parafit, Hommola and `Ives' tests for incidence data for each region
run.regional.module<-FALSE       # perform a krasnov-like test for incidence data for each region
VKfigure<-FALSE                  # generate figure for the incidence matrix in the Volga-Kamma region
SIMfigure<-TRUE                  # generate figures of different covariance patterns for simulated data
eps<-TRUE                        # make eps (if FALSE pdf)
get.EB.correction<-FALSE         # calculate Engel/Buist Correction factor 
test.overdispersion<-FALSE       # test whether marginal counts are overdispersed 
ptree<-"krasnov"                 # which algorithm to use to make parasite tree ultrameteric
                                          # "krasnov" (that used in Krasnov 2012: Mesquite)
                                          # "grafen"  (Grafen's method - uses apes compute.brlen)
                                          # "chronpl" (Sanderson's method - uses chronopl)

data.filepath<-"~/ParasitePhylogeny/Coevolutionary_interactions/Code and data from Hadfield et al 2014/"              # file path to data files 
results.filepath<-"~/ParasitePhylogeny/Coevolutionary_interactions/Code and data from Hadfield et al 2014/"  # file path where results are to be stored 
graphs.filepath<-"~/ParasitePhylogeny/Coevolutionary_interactions/Code and data from Hadfield et al 2014/"                 # file path where graphs are to be stored 
# 
# source("~/Work/DoublePhylo/R/hommola.R")  # Function for calculating Hommola's statistic

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

  ## S^{-1} is a phylogenetic covariance matrix with ancestral nodes retained 
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

############################################
## S^{-1} ignoring between-region effects ## 
############################################

parasiteA.R<-bdiag()
hostA.R<-bdiag()
host.parasiteA.R<-bdiag()
host.parasiteAS.R<-bdiag()
host.parasiteSA.R<-bdiag()
pnames<-c()
hnames<-c()
hpnames<-c()

## going region-by-region
for(i in 1:length(unique(ndat$Region))){

  ## get all of the unique hosts found in this region
h.pres<-unique(ndat$Host.species[which(ndat$Region==unique(ndat$Region)[i])])
## get all of the unique parasites found in this region
p.pres<-unique(ndat$Parasite.species[which(ndat$Region==unique(ndat$Region)[i])])

## trim the host and parasite phylogenies to include only these species
host.tree.tmp<-drop.tip(host.tree, host.tree$tip.label[which(!host.tree$tip.label%in%h.pres)])
parasite.tree.tmp<-drop.tip(parasite.tree, parasite.tree$tip.label[which(!parasite.tree$tip.label%in%p.pres)])
# drop species not present in the region from the trees

# generate S^{-1}_{i} matrices for each term:

hostA.tmp<-inverseA(host.tree.tmp)$Ainv
parasiteA.tmp<-inverseA(parasite.tree.tmp)$Ainv
host.parasiteA.tmp<-as(kronecker(hostA.tmp, parasiteA.tmp), "dgCMatrix")
host.parasiteAS.tmp<-as(kronecker(hostA.tmp, Diagonal(nrow(parasiteA.tmp))), "dgCMatrix")  
host.parasiteSA.tmp<-as(kronecker(Diagonal(nrow(hostA.tmp)), parasiteA.tmp), "dgCMatrix")  

hnames<-c(hnames, paste(unique(ndat$Region)[i], rownames(hostA.tmp), sep="."))
pnames<-c(pnames, paste(unique(ndat$Region)[i], rownames(parasiteA.tmp), sep="."))
hpnames<-c(hpnames, apply(expand.grid(rownames(parasiteA.tmp), rownames(hostA.tmp)), 1, function(x){paste(unique(ndat$Region)[i], x[2],x[1], sep=".")}))

# create block diagonal Diag(S^{-1}_{1}, S^{-1}_{2} ... S^{-1}_{n}):

hostA.R<-bdiag(hostA.R,hostA.tmp)
parasiteA.R<-bdiag(parasiteA.R,parasiteA.tmp)
host.parasiteA.R<-bdiag(host.parasiteA.R,host.parasiteA.tmp )
host.parasiteAS.R<-bdiag(host.parasiteAS.R,host.parasiteAS.tmp )
host.parasiteSA.R<-bdiag(host.parasiteSA.R,host.parasiteSA.tmp )
}

colnames(hostA.R)<-rownames(hostA.R)<-hnames
colnames(parasiteA.R)<-rownames(parasiteA.R)<-pnames
colnames(host.parasiteA.R)<-rownames(host.parasiteA.R)<-hpnames
colnames(host.parasiteAS.R)<-rownames(host.parasiteAS.R)<-hpnames
colnames(host.parasiteSA.R)<-rownames(host.parasiteSA.R)<-hpnames

################
## Run models ## 
################

if(run.MCMC){

if(incidence){

##############################
##### MCMC Incidence Data #### 
##############################

priorI=list(R=list(V=1, fix=1))

priorI$G<-lapply(1:9, function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})
names(priorI$G)<-paste("G", 1:9, sep="")
# create parameter expanded priors for each variance component

mI.MCMCa<-MCMCglmm(present~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+Parasite.species+Host.species+Parasite.species.ide+Host.species.ide+Host.Parasite+Host.Parasite.ide+Host.Parasite.ide2+Host.Parasite.ide3, family="categorical", data=ndat, ginverse=list(Host.species=hostA, Host.Parasite=host.parasiteA, Host.Parasite.ide2=host.parasiteAS, Host.Parasite.ide3=host.parasiteSA), prior=priorI, slice=T, nitt=1000000, thin=400, burnin=200000)

# main model: includes between-region information and controls for sampling effort

save(mI.MCMCa, file=paste(paste(results.filepath, "mI.MCMCa", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

mI.MCMCb<-MCMCglmm(present~1, random=~Region+Parasite.species+Host.species+Parasite.species.ide+Host.species.ide+Host.Parasite+Host.Parasite.ide+Host.Parasite.ide2+Host.Parasite.ide3, family="categorical", data=ndat, ginverse=list(Host.species=hostA, Host.Parasite=host.parasiteA, Host.Parasite.ide2=host.parasiteAS, Host.Parasite.ide3=host.parasiteSA), prior=priorI, slice=T, nitt=400000, thin=400, burnin=200000)

# model b: includes between-region information but does not control for sampling effort

save(mI.MCMCb, file=paste(paste(results.filepath, "mI.MCMCb", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

priorI$G<-priorI$G[-9]

mI.MCMCc<-MCMCglmm(present~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+Parasite.species.r+Host.species.r+Parasite.species.ide.r+Host.species.ide.r+Host.Parasite.r+Host.Parasite.ide2.r+Host.Parasite.ide3.r, data=ndat, ginverse=list(Parasite.species.r=parasiteA.R, Host.species.r=hostA.R, Host.Parasite.r=host.parasiteA.R, Host.Parasite.ide2.r=host.parasiteAS.R, Host.Parasite.ide3.r=host.parasiteSA.R), prior=priorI, slice=T, nitt=1000000, thin=400, burnin=200000, family="categorical")

# model c: does not include between-region information but controls for sampling effort

save(mI.MCMCc, file=paste(paste(results.filepath, "mI.MCMCc", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

mI.MCMCd<-MCMCglmm(present~1, random=~Region+Parasite.species.r+Host.species.r+Parasite.species.ide.r+Host.species.ide.r+Host.Parasite.r+Host.Parasite.ide2.r+Host.Parasite.ide3.r, data=ndat, ginverse=list(Parasite.species.r=parasiteA.R, Host.species.r=hostA.R, Host.Parasite.r=host.parasiteA.R, Host.Parasite.ide2.r=host.parasiteAS.R, Host.Parasite.ide3.r=host.parasiteSA.R), prior=priorI, slice=T, nitt=1000000, thin=400, burnin=200000, family="categorical")

# model d: does not include between-region information and does not control for sampling effort

save(mI.MCMCd, file=paste(paste(results.filepath, "mI.MCMCd", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))
}

if(counts){

##########################
##### MCMC Count Data #### 
##########################

# These analyses not published - distribution of counts is assumed to be overdispersed zero-truncated Poisson

priorC=list(R=list(V=1, nu=0))

priorC$G<-lapply(1:9, function(x){list(V=1, nu=1, alpha.mu=0, alpha.V=1000)})
names(priorC$G)<-paste("G", 1:9, sep="")

mC.MCMCa<-MCMCglmm(counts~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+Parasite.species+Host.species+Parasite.species.ide+Host.species.ide+Host.Parasite+Host.Parasite.ide+Host.Parasite.ide2+Host.Parasite.ide3, family="ztpoisson", data=ndat, ginverse=list(Host.species=hostA, Host.Parasite=host.parasiteA, Host.Parasite.ide2=host.parasiteAS, Host.Parasite.ide3=host.parasiteSA), prior=priorC, nitt=1000000, thin=400, burnin=200000)

save(mC.MCMCa, file=paste(paste(results.filepath, "mC.MCMCa", sep=.Platform$file.sep),ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

# main count model: includes between-region information and controls for sampling effort

mC.MCMCb<-MCMCglmm(counts~1, random=~Region+Parasite.species+Host.species+Parasite.species.ide+Host.species.ide+Host.Parasite+Host.Parasite.ide+Host.Parasite.ide2+Host.Parasite.ide3, family="categorical", data=ndat, ginverse=list(Host.species=hostA, Host.Parasite=host.parasiteA, Host.Parasite.ide2=host.parasiteAS, Host.Parasite.ide3=host.parasiteSA), prior=priorC, nitt=1000000, thin=400, burnin=200000)

save(mC.MCMCb, file=paste(paste(results.filepath, "mC.MCMCb", sep=.Platform$file.sep),ptree, format(Sys.time(), "%d-%m-%Y"), "R", sep="."))

# count model b: includes between-region information but does not control for sampling effort

priorC$G<-priorC$G[-9]

mC.MCMCc<-MCMCglmm(counts~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+Parasite.species.r+Host.species.r+Parasite.species.ide.r+Host.species.ide.r+Host.Parasite.r+Host.Parasite.ide2.r+Host.Parasite.ide3.r, data=ndat, ginverse=list(Parasite.species.r=parasiteA.R, Host.species.r=hostA.R, Host.Parasite.r=host.parasiteA.R, Host.Parasite.ide2.r=host.parasiteAS.R, Host.Parasite.ide3.r=host.parasiteSA.R), prior=priorC, nitt=1000000, thin=400, burnin=200000)

save(mC.MCMCc, file=paste(paste(results.filepath, "mC.MCMCc", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "Rdata", sep="."))

# count model c: does not include between-region information but controls for sampling effort

mC.MCMCd<-MCMCglmm(counts~1, random=~Region+Parasite.species.r+Host.species.r+Parasite.species.ide.r+Host.species.ide.r+Host.Parasite.r+Host.Parasite.ide2.r+Host.Parasite.ide3.r, data=ndat, ginverse=list(Parasite.species.r=parasiteA.R, Host.species.r=hostA.R, Host.Parasite.r=host.parasiteA.R, Host.Parasite.ide2.r=host.parasiteAS.R, Host.Parasite.ide3.r=host.parasiteSA.R), prior=priorC, nitt=1000000, thin=400, burnin=200000)

save(mC.MCMCd, file=paste(paste(results.filepath, "mC.MCMCd", sep=.Platform$file.sep), ptree, format(Sys.time(), "%d-%m-%Y"), "Rdata", sep="."))

# count model d: does not include between-region information and does not control for sampling effort
}
}


if(run.ASReml){

if(incidence){

################################
##### ASReml Incidence Data #### 
################################

mI.ASRemla<-asreml(present~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+Host.Parasite.ide+giv(Host.species, var=T):Parasite.species.ide+giv(Parasite.species, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species=sm2asreml(parasiteA), Host.species=sm2asreml(hostA), Host.Parasite=sm2asreml(host.parasiteA)), workspace = 9e+07, family = asreml.binomial(), maxiter = 30)

save(mI.ASRemla, file=paste(paste(results.filepath, "mI.ASRemla.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))

# main model: includes between-region information and controls for sampling effort

mI.ASRemlb<-asreml(present~1, random=~Region+giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+Host.Parasite.ide+giv(Host.species, var=T):Parasite.species.ide+giv(Parasite.species, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species=sm2asreml(parasiteA), Host.species=sm2asreml(hostA), Host.Parasite=sm2asreml(host.parasiteA)), workspace = 9e+07, family = asreml.binomial(), maxiter = 30)

save(mI.ASRemlb, file=paste(paste(results.filepath, "mI.ASRemlb.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))

# model b: includes between-region information but does not control for sampling effort

mI.ASRemlc<-asreml(present~log(no.hosts.sampled)+log(no.parasites.sampled), random=~Region+giv(Parasite.species.r, var=T)+giv(Host.species.r, var=T)+Parasite.species.ide.r+Host.species.ide.r+giv(Host.Parasite.r, var=T)+giv(Host.species.r, var=T):Parasite.species.ide+giv(Parasite.species.r, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species.r=sm2asreml(parasiteA.R), Host.species.r=sm2asreml(hostA.R), Host.Parasite.r=sm2asreml(host.parasiteA.R)), workspace = 8e+07, family = asreml.binomial(), maxiter = 30)

save(mI.ASRemlc, file=paste(paste(results.filepath, "mI.ASRemlc.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))

# model c: does not include between-region information but controls for sampling effort

mI.ASRemld<-asreml(present~1, random=~Region+giv(Parasite.species.r, var=T)+giv(Host.species.r, var=T)+Parasite.species.ide.r+Host.species.ide.r+giv(Host.Parasite.r, var=T)+giv(Host.species.r, var=T):Parasite.species.ide+giv(Parasite.species.r, var=T):Host.species.ide, data=ndat, ginverse=list(Parasite.species.r=sm2asreml(parasiteA.R), Host.species.r=sm2asreml(hostA.R), Host.Parasite.r=sm2asreml(host.parasiteA.R)), workspace = 8e+07, family = asreml.binomial(), maxiter = 30)

save(mI.ASRemld, file=paste(paste(results.filepath, "mI.ASRemld.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))

# model d: does not include between-region information and does not control for sampling effort
}
if(counts){
stop("ASreml does not have a zero-truncated Poisson distribution for modelling the count data (could try Poisson?)")
}
}


######################
##### Abundanaces #### 
######################

if(abundance){

ndat.h<-ndat[which(duplicated(paste(ndat$Region, ndat$Host.species))==FALSE),]
ndat.p<-ndat[which(duplicated(paste(ndat$Region, ndat$Parasite.species))==FALSE),]

if(run.MCMC){

###############
##### MCMC #### 
###############

prior.mA.MCMC<-list(R=list(V=1, nu=0), G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G2=list(V=1, nu=1, alpha.mu=0, alpha.V=1000), G3=list(V=1, nu=1, alpha.mu=0, alpha.V=1000)))

mHA.MCMC<-MCMCglmm(log(no.hosts.sampled)~1, random=~Region+Host.species+Host.species.ide, data=ndat.h, ginverse=list(Host.species=hostA), prior=prior.mA.MCMC, nitt=13000*5, thin=10*2.5, burnin=3000*5)

save(mHA.MCMC, file=paste(results.filepath, "mHA.MCMC.Rdata", sep=.Platform$file.sep))

# host abundance model

mPA.MCMC<-MCMCglmm(log(no.parasites.sampled)~1, random=~Region+Parasite.species+Parasite.species.ide, ginverse=list(Parasite.species=parasiteA), prior=prior.mA.MCMC, nitt=13000*5, thin=10*2.5, burnin=3000*5,data=ndat.p)

save(mPA.MCMC, file=paste(paste(results.filepath, "mPA.MCMC.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))

# parasite abundance model
}

if(run.ASReml){

#################
##### ASReml #### 
#################

mHA.ASReml<-asreml(log(no.hosts.sampled)~1, random=~Region+giv(Host.species, var=T)+Host.species.ide, data=ndat.h, ginverse=list(Host.species=sm2asreml(hostA)), workspace = 8e+07, family = asreml.gaussian(), maxiter = 30)

save(mHA.ASReml, file=paste(results.filepath, "mHA.ASReml.Rdata", sep=.Platform$file.sep))

# host abundance model

mPA.ASReml<-asreml(log(no.parasites.sampled)~1, random=~Region+giv(Parasite.species, var=T)+Parasite.species.ide, data=ndat.p, ginverse=list(Parasite.species=sm2asreml(parasiteA)), workspace = 8e+07, family = asreml.gaussian(), maxiter = 30)

save(mPA.ASReml, file=paste(paste(results.filepath, "mPA.ASReml.", sep=""), ptree, ".Rdata", sep=.Platform$file.sep))

# parasite abundance model

}
}

######################################
## Run parafit on consolidated data ##
######################################

if(run.parafit){

  ht<-vcv(host.tree, corr=T)
  pt<-vcv(parasite.tree, corr=T)
  # A matrices for host and parasite

  Y<-table(ndat$Host.species, ndat$Parasite.species, ndat$present)[,,2]>0
  # incidence data

  hi<-match(rownames(Y), rownames(ht))
  pi<-match(colnames(Y), rownames(pt))

  ht<-ht[hi,hi]
  pt<-pt[pi,pi]
  # reorder A matrices so match row/columns of Y

  htp<-pcoa(1-ht)$vectors
  ptp<-pcoa(1-pt)$vectors
  # get principal coordinate of phylogenetic distance matrix

  hte<-t(t(eigen(solve(ht))$vectors)*sqrt(eigen(solve(ht))$values))
  pte<-t(t(eigen(solve(pt))$vectors)*sqrt(eigen(solve(pt))$values))

  # get unnormalised eigenvectors of A

  lD<-t(htp)%*%Y%*%ptp
  iD<-t(hte)%*%(Y-mean(Y))%*%pte
  # D matrices (see Appendix) for Legendre and Ives methods
 
  cl<-sum(diag(t(lD)%*%lD))
  ci<-sum(diag(t(iD)%*%iD))
  ch<-hommola(Y, ht,pt)
  # metrics (see Appendix) of Legendre, Ives & Hommola

  sl.l<-1:1000
  si.l<-1:1000
  sh.l<-1:1000
  # stoing metrics of Legendre, Ives & Hommola after Legendre permutations  
  sl.h<-1:1000
  si.h<-1:1000
  sh.h<-1:1000
  # stoing metrics of Legendre, Ives & Hommola after Hommola permutations  

  ##################
  ## Permutations ##
  ##################

  for(j in 1:1000){

    Y2<-apply(Y, 2, sample)                     # Legendre permtations

    lD<-t(htp)%*%Y2%*%ptp                       # Legendre D matrix
    iD<-t(hte)%*%(Y2-mean(Y2))%*%pte            # Ives D matrix

    sl.l[j]<-sum(diag(t(lD)%*%lD))              # Legendre metric (ParafitGlobal)
    si.l[j]<-sum(diag(t(iD)%*%iD))              # Ives metric (MSEb)
    sh.l[j]<-hommola(Y2, ht,pt)                 # Hommola metric

    Y2<-Y[sample(1:nrow(Y)),sample(1:ncol(Y))]  # Hommola sampling

    lD<-t(htp)%*%Y2%*%ptp
    iD<-t(hte)%*%(Y2-mean(Y2))%*%pte
 
    sl.h[j]<-sum(diag(t(lD)%*%lD)) 
    si.h[j]<-sum(diag(t(iD)%*%iD)) 
    sh.h[j]<-hommola(Y2, ht,pt)    
  }

    other.tests<-cbind(c(cl, ci, ch), c(sum(cl<sl.l)/1000, sum(ci>si.l)/1000, sum(ch<sh.l)/1000), c(sum(cl<sl.h)/1000, sum(ci>si.h)/1000, sum(ch<sh.h)/1000))
    # store metrics, and the proportion of times the metrics under permutation were greater 

    rownames(other.tests)<-c("L", "I", "H")
    colnames(other.tests)<-c("statistic", "L-pval", "H-pval")

    save(other.tests, file=paste(paste(results.filepath, "other.tests.", sep=.Platform$file.sep), ptree, ".Rdata",sep=""))
}

######################################
## Run module on consolidated data ##
######################################

if(run.module){

   ht<-vcv(host.tree, corr=T)
   pt<-vcv(parasite.tree, corr=T)
   # A matrices for host and parasite

   Y<-table(ndat$Host.species, ndat$Parasite.species, ndat$present)[,,2]>0
   # incidence data

   hi<-match(rownames(Y), rownames(ht))
   pi<-match(colnames(Y), rownames(pt))

   ht<-ht[hi,hi]
   pt<-pt[pi,pi]
   # reorder A matrices so match row/columns of Y

   G<-graph.incidence(Y)
   MU<-walktrap.community(G)
   MU<-membership(MU)
   # assign species to modules

   lth<-which(lower.tri(ht))

   chk<-cor(outer(MU[1:length(hi)], MU[1:length(hi)], "==")[lth], c(1-ht)[lth])
   # correlation between comembership and phylogenetic distance of hosts

   ltp<-lower.tri(pt)

   cpk<-cor(outer(MU[length(hi)+1:length(pi)], MU[length(hi)+1:length(pi)], "==")[ltp], c(1-pt)[ltp])
   # correlation between comembership and phylogenetic distance of parasites

   shk.l<-1:1000
   spk.l<-1:1000
   # storing host & parasite membership/phylogenetic-distance correlations after Legendre Permutations 
   shk.h<-1:1000
   spk.h<-1:1000
   # storing host & parasite membership/phylogenetic-distance correlations after Hommola Permutations 

   for(i in 1:1000){

     Y2<-apply(Y, 2, sample)                     # Legendre sampling

     G2<-graph.incidence(Y2)
     MU2<-walktrap.community(G2)
     MU2<-membership(MU2)

     shk.l[i]<-cor(outer(MU2[1:length(hi)], MU2[1:length(hi)], "==")[lth], c(1-ht)[lth])
     spk.l[i]<-cor(outer(MU[length(hi)+1:length(pi)], MU[length(hi)+1:length(pi)], "==")[ltp], c(1-pt)[ltp])

     Y2<-Y[sample(1:nrow(Y)),sample(1:ncol(Y))]  # Hommola sampling

     G2<-graph.incidence(Y2)
     MU2<-walktrap.community(G2)
     MU2<-membership(MU2)

     shk.h[i]<-cor(outer(MU2[1:length(hi)], MU2[1:length(hi)], "==")[lth], c(1-ht)[lth])
     spk.h[i]<-cor(outer(MU[length(hi)+1:length(pi)], MU[length(hi)+1:length(pi)], "==")[ltp], c(1-pt)[ltp])
     print(i)
   }

    k.tests<-cbind(c(chk, cpk), c(sum(chk>shk.l)/1000, sum(cpk>spk.l)/1000), c(sum(chk>shk.h)/1000, sum(cpk>spk.h)/1000))
    # store metrics, and the proportion of times the metrics under permutation were greater 

    rownames(k.tests)<-c("H", "P")
    colnames(k.tests)<-c("statistic", "L-pval", "H-pval")

    save(k.tests, file=paste(paste(results.filepath, "k.tests.", sep=.Platform$file.sep), ptree, ".Rdata", sep=""))

}

############################
## Site-specific analyses ##
############################

if(run.regional.module | run.regional.ASReml | run.regional.parafit){

resd<-matrix(NA, 51,6)
rese<-matrix(NA, 51,7*2)   # store results from ASReml models (controlling for sampling effort)
resf<-matrix(NA, 51,7*2)   # store results from ASReml models (not controlling for sampling effort)
resg<-matrix(NA, 51,9)

for(j in 1:length(unique(ndat$Region))){

ndat_tmp<-ndat[which(ndat$Region==unique(ndat$Region)[j]),]
# obtain site-specific data


ihosts<-unique(ndat_tmp$Host.species)
iparas<-unique(ndat_tmp$Parasite.species)

# obtain site-specific phylogenies

host.tree.tmp<-drop.tip(host.tree, host.tree$tip.label[which(!host.tree$tip.label%in%ihosts)])
parasite.tree.tmp<-drop.tip(parasite.tree, parasite.tree$tip.label[which(!parasite.tree$tip.label%in%iparas)])
# drop sepcies from the tree if not present at site i

ht<-vcv(host.tree.tmp, corr=T)     
pt<-vcv(parasite.tree.tmp, corr=T)  
# obtain site-specific A matrices

Y<-matrix(ndat_tmp$present,length(ihosts),length(iparas))
rownames(Y)<-ihosts
colnames(Y)<-iparas
# obtain site-specific incidence data

hi<-match(rownames(Y), rownames(ht))
pi<-match(colnames(Y), rownames(pt))

ht<-ht[hi,hi]
pt<-pt[pi,pi]
# reorder A matrices so they agree with incidence data


#######################################
## Site-specific ASReml anlayses ##
#######################################

if(run.regional.ASReml){

parasiteA.tmp<-inverseA(parasite.tree.tmp)$Ainv
hostA.tmp<-inverseA(host.tree.tmp)$Ainv
host.parasiteA.tmp<-as(kronecker(hostA.tmp, parasiteA.tmp), "dgCMatrix")
rownames(host.parasiteA.tmp)<-apply(expand.grid(rownames(parasiteA.tmp), rownames(hostA.tmp)), 1, function(x){paste(x[2],x[1], sep=".")})
# generate S^{-1}_{j} matrices
 
mI.ASRemle<-asreml(present~log(no.hosts.sampled)+log(no.parasites.sampled), random=~giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+giv(Parasite.species, var=T):Host.species.ide+giv(Host.species, var=T):Parasite.species.ide, data=ndat_tmp, ginverse=list(Parasite.species=sm2asreml(parasiteA.tmp), Host.species=sm2asreml(hostA.tmp), Host.Parasite=sm2asreml(host.parasiteA.tmp)), workspace = 8e+07, family = asreml.binomial(), maxiter = 30)

# model that controls for sampling effort

rese[j,]<-c(rbind(round(summary(mI.ASRemle)$varcomp[,1][-8],3),round(summary(mI.ASRemle)$varcomp[,3][-8],3)))
# store variance components and their standard errors

mI.ASRemlf<-asreml(present~1, random=~giv(Parasite.species, var=T)+giv(Host.species, var=T)+Parasite.species.ide+Host.species.ide+giv(Host.Parasite, var=T)+giv(Parasite.species, var=T):Host.species.ide+giv(Host.species, var=T):Parasite.species.ide, data=ndat_tmp, ginverse=list(Parasite.species=sm2asreml(parasiteA.tmp), Host.species=sm2asreml(hostA.tmp), Host.Parasite=sm2asreml(host.parasiteA.tmp)), workspace = 8e+07, family = asreml.binomial(), maxiter = 30)

# model that does not control for sampling effort

resf[j,]<-c(rbind(round(summary(mI.ASRemlf)$varcomp[,1][-8],3),round(summary(mI.ASRemlf)$varcomp[,3][-8],3)))
# store variance components and their standard errors

}

###################################
## Site-specific module anlayses ##
###################################

if(run.regional.module){

  r.keep<-which(rowSums(Y)>0)
  c.keep<-which(colSums(Y)>0)
  # only retain species with an observed interaction

  G<-graph.incidence(Y[r.keep,c.keep])
  MU<-walktrap.community(G)
  MU<-membership(MU)
  # assign module membership

  ht.m<-ht[r.keep,r.keep]
  pt.m<-pt[c.keep,c.keep]
  # remove species with no observed interaction from the A matrices

  hlt<-which(lower.tri(ht.m))
  plt<-which(lower.tri(pt.m))

  ch<-cor(outer(MU[1:length(ihosts[r.keep])], MU[1:length(ihosts[r.keep])], "==")[hlt], c(1-ht.m)[hlt])
  cp<-cor(outer(MU[length(ihosts[r.keep])+1:length(iparas[c.keep])], MU[length(ihosts[r.keep])+1:length(iparas[c.keep])], "==")[plt], c(1-pt.m)[plt])
  # calculate correlation between module membership and phylogenetic distance
 
  sh.l<-1:1000
  sp.l<-1:1000
  # storing host & parasite membership/phylogenetic-distance correlations after Legendre Permutations 
  sh.h<-1:1000
  sp.h<-1:1000
  # storing host & parasite membership/phylogenetic-distance correlations after Hommola Permutations 

  ##################
  ## Permute Data ##
  ##################

  for(k in 1:1000){

      Y2<-apply(Y, 2, sample)                    # Legendre sampling

      r.keep<-which(rowSums(Y2)>0)
      c.keep<-which(colSums(Y2)>0)

      G<-graph.incidence(Y2[r.keep,c.keep])
      MU<-walktrap.community(G)
      MU<-membership(MU)

      ht.m<-ht[r.keep,r.keep]
      pt.m<-pt[c.keep,c.keep]

      hlt<-which(lower.tri(ht.m))
      plt<-which(lower.tri(pt.m))

      sh.l[k]<-cor(outer(MU[1:length(ihosts[r.keep])], MU[1:length(ihosts[r.keep])], "==")[hlt], c(1-ht.m)[hlt])
      sp.l[k]<-cor(outer(MU[length(ihosts[r.keep])+1:length(iparas[c.keep])], MU[length(ihosts[r.keep])+1:length(iparas[c.keep])], "==")[plt], c(1-pt.m)[plt])

      Y2<-Y[sample(1:nrow(Y)),sample(1:ncol(Y))] # Hommola sampling

      r.keep<-which(rowSums(Y2)>0)
      c.keep<-which(colSums(Y2)>0)

      G<-graph.incidence(Y2[r.keep,c.keep])
      MU<-walktrap.community(G)
      MU<-membership(MU)

      ht.m<-ht[r.keep,r.keep]
      pt.m<-pt[c.keep,c.keep]

      hlt<-which(lower.tri(ht.m))
      plt<-which(lower.tri(pt.m))

      sh.h[k]<-cor(outer(MU[1:length(ihosts[r.keep])], MU[1:length(ihosts[r.keep])], "==")[hlt], c(1-ht.m)[hlt])

      sp.h[k]<-cor(outer(MU[length(ihosts[r.keep])+1:length(iparas[c.keep])], MU[length(ihosts[r.keep])+1:length(iparas[c.keep])], "==")[plt], c(1-pt.m)[plt])
    }

    resd[j,]<-c(ch,sum(ch>sh.l, na.rm=T)/sum(!is.na(sh.l)), sum(ch>sh.h, na.rm=T)/sum(!is.na(sh.h)), cp, sum(cp>sp.l, na.rm=T)/sum(!is.na(sp.l)), sum(cp>sp.h, na.rm=T)/sum(!is.na(sp.h)))
    # store metrics, and the proportion of times the metrics under permutation were greater. 
    # Note that under some permutations the correlation is not defined because there is only one module
}

###############################################
## Site-specific Legendre (Parafit) anlayses ##
###############################################

if(run.regional.parafit){

    htp<-pcoa(1-ht)$vectors
    ptp<-pcoa(1-pt)$vectors

    hte<-t(t(eigen(solve(ht))$vectors)*sqrt(eigen(solve(ht))$values))
    pte<-t(t(eigen(solve(pt))$vectors)*sqrt(eigen(solve(pt))$values))

    lD<-t(htp)%*%Y%*%ptp
    iD<-t(hte)%*%(Y-mean(Y))%*%pte


    cl<-sum(diag(t(lD)%*%lD))
    ci<-sum(diag(t(iD)%*%iD))
    ch<-hommola(Y, ht,pt)

    sl.l<-1:1000
    si.l<-1:1000
    sh.l<-1:1000
    sl.h<-1:1000
    si.h<-1:1000
    sh.h<-1:1000

    for(k in 1:1000){
      Y2<-apply(Y, 2, sample)  # Legendre sampling
      lD<-t(htp)%*%Y2%*%ptp
      iD<-t(hte)%*%(Y2-mean(Y2))%*%pte

      sl.l[k]<-sum(diag(t(lD)%*%lD))
      si.l[k]<-sum(diag(t(iD)%*%iD))
      sh.l[k]<-hommola(Y2, ht,pt)

      Y2<-Y[sample(1:nrow(Y)),sample(1:ncol(Y))] # Hommola sampling

      lD<-t(htp)%*%Y2%*%ptp
      iD<-t(hte)%*%(Y2-mean(Y2))%*%pte
 
      sl.h[k]<-sum(diag(t(lD)%*%lD))
      si.h[k]<-sum(diag(t(iD)%*%iD))
      sh.h[k]<-hommola(Y2, ht,pt)
    }

    resg[j,]<-t(cbind(c(cl, ci, ch), c(sum(cl<sl.l)/1000, sum(ci>si.l)/1000, sum(ch<sh.l)/1000), c(sum(cl<sl.h)/1000, sum(ci>si.h)/1000, sum(ch<sh.h)/1000)))
}


print(j)
}

##################################
## save region-specific results ##
##################################

if(run.regional.ASReml){
colnames(rese)<-c(rbind(rownames(summary(mI.ASRemle)$varcomp)[-8], paste(rownames(summary(mI.ASRemle)$varcomp)[-8], ".Z", sep="")))
rownames(rese)<-unique(ndat$Region)

colnames(resf)<-colnames(rese)
rownames(resf)<-rownames(resf)
save(rese, file=paste(paste(results.filepath, sitespecificA., sep=.Platform$file.sep), ptree, ".Rdata", sep=""))
save(resf, file=paste(paste(results.filepath, sitespecificB., sep=.Platform$file.sep), ptree, ".Rdata", sep=""))
}

if(run.regional.module){
rownames(resd)<-unique(ndat$Region)
resd[,2][which(is.na(resd[,1]))]<-NA
resd[,3][which(is.na(resd[,1]))]<-NA
resd[,5][which(is.na(resd[,4]))]<-NA
resd[,6][which(is.na(resd[,4]))]<-NA

colnames(resd)<-c("host-statistic", "host-pvalL", "host-pvalH", "para-statistic", "para-pvalL", "para-pvalH")  

save(resd, file=paste(paste(results.filepath, sitespecificK., sep=.Platform$file.sep), ptree, ".Rdata", sep=""))
}

if(run.regional.parafit){
rownames(resg)<-unique(ndat$Region)
colnames(resg)<-c("L-statistic", "L-pvalL", "L-pvalH", "I-statistic", "I-pvalL", "I-pvalH", "H-statistic", "H-pvalL", "H-pvalH") 

save(resg, file=paste(paste(results.filepath, sitespecificP., sep=.Platform$file.sep), ptree, ".Rdata", sep=""))
}

}


############################################
## Try out  Engel/Buist Correction factor ##
############################################

if(get.EB.correction){

if(any(ls()=="mI.MCMCa")==FALSE){    # read in most recent main MCMC model

    mfiles<-list.files(results.filepath)
    mfiles<-mfiles[which(grepl("mI\\.MCMCa\\.", mfiles) & grepl(ptree, mfiles))]
    mfiles<-mfiles[which.max(as.Date(substr(mfiles, nchar(mfiles)-15, nchar(mfiles)-6), "%d-%m-%Y"))]

    load(paste(results.filepath, .Platform$file.sep, mfiles,".Rdata", sep=""))
}

if(any(ls()=="mI.ASReml")==FALSE){  # read in most recent main MCMC model
  load(paste(results.filepath, .Platform$file.sep, "mI.ASRemla.", ptree, ".Rdata",  sep=""))
}

ht<-vcv(host.tree, corr=T)      
pt<-vcv(parasite.tree, corr=T)

# obtain A matrices

n<-dim(ht)[1]
m<-dim(pt)[1]

host.name<-host.tree$tip.label
parasite.name<-parasite.tree$tip.label
host.par.name<-paste(rep(host.name, each=m), rep(parasite.name, n), sep=".")

X<-mI.MCMC$X
# get fixed effect design matrix

Z<-mI.MCMC$Z[,which(grepl("Host.Parasite", colnames(mI.MCMC$Z)) & grepl("Host.Parasite.ide", colnames(mI.MCMC$Z))==FALSE & grepl("Node", colnames(mI.MCMC$Z))==FALSE)]
colnames(Z)<-gsub("Host.Parasite.", "", colnames(Z))
Z<-Z[,match(colnames(Z), host.par.name)]
# get random effect design matrix for the ceovolutionary interaction

u0<-1/(1+exp(-(X%*%rev(mI.ASRemla$coefficients$fixed))))

W0<-Diagonal(length(u0), (u0*(1-u0))@x)
W1<-Diagonal(length(u0), (u0*(1-u0)*(1-2*u0))@x)
W2<-Diagonal(length(u0), (u0*(1-u0)*(6*(u0^2)-6*u0+1))@x)

Q<-diag(t(Z)%*%W0%*%Z)
Qtilde<-matrix(Q, nrow(ht), nrow(pt))
a<-t(Q)%*%c((ht*ht)%*%Qtilde%*%(pt*pt))/2

q<-Z%*%Matrix(kronecker(diag(pt), diag(ht)))
b<-(t(q)%*%W2%*%q-t(q)%*%W1%*%X%*%solve(t(X)%*%W0%*%X)%*%(t(X)%*%W1%*%q))/4

f<-a/(a+b)

# calculate correction factor with terms defined in the Supplementary materials
}

####################################################
## Test whether marginal counts are overdispersed ##
####################################################

if(test.overdispersion){

PSR<-tapply(ndat$present, paste(ndat$Host.species, ndat$Region), sum)
# number of parasite species per host species in each region (parasite species richness)
HA<-tapply(ndat$no.hosts.sampled, paste(ndat$Host.species, ndat$Region), mean)
# number of individuals of each host species in each region (host abundance)
HR<-tapply(ndat$present, paste(ndat$Parasite.species, ndat$Region), sum)
# number of host species per parasite species in each region (host range)
PA<-tapply(ndat$no.parasites.sampled, paste(ndat$Parasite.species, ndat$Region), mean)
# number of individuals of each parasite species in each region (parasite abundance)

nH<-tapply(ndat$Host.species, ndat$Region, function(x){length(unique(x))})
# number of host species in each region
NH<-nH[match(tapply(ndat$Region, paste(ndat$Parasite.species, ndat$Region), function(x){as.character(x[1])}), names(nH))]
# number of available hosts for each parasite species in each region
nP<-tapply(ndat$Parasite.species, ndat$Region, function(x){length(unique(x))})
# number of parasite species in each region
NP<-nP[match(tapply(ndat$Region, paste(ndat$Host.species, ndat$Region), function(x){as.character(x[1])}),names(nP))]
# number of potential parasites for each host species in each region

pP<-tapply(PSR, names(NP), mean)/nP
# mean parasite species richness at each site divided by maximum possible parasite species richness
vP<-tapply(PSR, names(NP), var)
# variance in parasite species richness at each site
EvP<-nP*pP*(1-pP)
# expected variance in parasite species richness at each site in the absence of under/over-dispersion
summary(lm(vP~EvP))$coef[2,]
# variance in PSR greater than expected

pH<-tapply(HR, names(NH), mean)/nH
# mean host ramge at each site divided by maximum possible host range
vH<-tapply(HR, names(NH), var)
# variance in host-range at each site
EvH<-nH*pH*(1-pH)
# expected variance in host-range at each site in the absence of under/over-dispersion
summary(lm(vH~EvH))$coef[2,]
# variance in HR greater than expected
}

#######################
## get VolgaKama pdf ##
#######################

if(VKfigure){

  ndat_tmp<-ndat[which(ndat$Region=="Volga-Kama"),]

  ihosts<-unique(ndat_tmp$Host.species)
  iparas<-unique(ndat_tmp$Parasite.species)

  host.tree.tmp<-compute.brtime(multi2di(drop.tip(host.tree, host.tree$tip.label[which(!host.tree$tip.label%in%ihosts)])))
  parasite.tree.tmp<-compute.brtime(multi2di(drop.tip(parasite.tree, parasite.tree$tip.label[which(!parasite.tree$tip.label%in%iparas)])))

  # get trees for species found in the Volga-Kama region

  Y<-matrix(ndat_tmp$present,length(ihosts),length(iparas))

  Y<-Y[match(host.tree.tmp$tip.label, ihosts), match(parasite.tree.tmp$tip.label, iparas)]

  # get Incidence data for the Volga-Kama region
  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKama.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Y,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none", col=c("white", "black"),margins = c(11, 10), font.axis=3)
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKama.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()
}

######################################################################
## get heatmap representations of the different sources of variance ##
######################################################################

if(SIMfigure){

  set.seed(516)                        # this produces a pretty picture!

  n<-m<-10

  host.tree.tmp<-rcoal(10)
  parasite.tree.tmp<-rcoal(10)

  host.tree.tmp$tip.label<-gsub("t","h", host.tree.tmp$tip.label)
  parasite.tree.tmp$tip.label<-gsub("t","p", parasite.tree.tmp$tip.label)

  # generate host and parasite trees with 10 tips

  ht<-vcv(host.tree.tmp, corr=T)
  pt<-vcv(parasite.tree.tmp, corr=T)

  # get A matrices

  Ch<-chol(ht)
  Cp<-chol(pt)

  # get Cholesky decomposition of A

  Ytmp<-kronecker(pt, matrix(1,n,n))

  # get covariance structure for parasite phylogenetic main effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)
  # the covariance structure is 100x100 (100 = the number of interactions)
  # get the covariance between the 45th interaction (Parasite 3 and Host 8) and the rest, 
  # rank them by magnitude, and organise into a 10x10 matrix.

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim1.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim1.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(pt, ht)
  # get covariance structure for phylogenetic coevolutionary interaction effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim2.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim2.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(pt, diag(n))
  # get covariance structure for parasite phylogenetic evolutionary interaction effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim3.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim3.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(diag(m), ht)
  # get covariance structure for host phylogenetic evolutionary interaction effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim4.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim4.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(matrix(1,m,m), ht)
  # get covariance structure for host phylogenetic main effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim5.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim5.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(diag(m), matrix(1,n,n))

  # get covariance structure for paraiste non-phylogenetic main effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim6.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim6.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(matrix(1,m,m), diag(n))
  # get covariance structure for host non-phylogenetic main effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim7.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim7.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()

  Ytmp<-kronecker(diag(m), diag(n))
  # get covariance structure for non-phylogenetic interaction effects

  Ytmp<-matrix(match(Ytmp[,45], sort(unique(Ytmp[,45]))),n,m)

  rownames(Ytmp)<-rownames(ht)
  colnames(Ytmp)<-colnames(pt)

  if(!eps){
    pdf(paste(graphs.filepath, "VolgaKamaSim8.pdf", sep=.Platform$file.sep), width=8, height=8)
  }
  heatmap(Ytmp,  Rowv=as.dendrogram(as.hclust(host.tree.tmp)), Colv=as.dendrogram(as.hclust(parasite.tree.tmp)),labRow=gsub("_", " ", host.tree.tmp$tip.label),labCol=gsub("_", " ", parasite.tree.tmp$tip.label), scale="none",margins = c(5, 6), font.axis=3, col=rev(grey(seq(0,1,length=max(Ytmp)))), cexRow=2, cexCol=2) 
  if(eps){
    dev.copy2eps(file=paste(graphs.filepath, "VolgaKamaSim8.eps", sep=.Platform$file.sep), width=8, height=8)
  }
  dev.off()
}
