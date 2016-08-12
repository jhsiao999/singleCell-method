
## Squarem implementation of the Poisson Random Batch model

library(SQUAREM)
library(optimx)
library(gtools)
library(parallel)
library(lineprof)
library(lme4)
library(permute)
library(BioPhysConnectoR)

#######  Source the important files

source('poisson_loglik.R')
source('simplex_functions.R')
source('estimation_poisson_batch.R')
source('simulate_counts_poisson.R')

###  If you want to simulate from the batch model : try running the example_simulation.R file

source('example_simulation.R');

##### From now on we assume we only have the counts matrix available ##########################

source('poisson_topic_loglink.R')

### apply the main function for the modeling

out <- Poisson_topic.loglink(counts,n_clus=4,lab_batch,use_squarem = FALSE)

###
docweights=out$omega;
perm_set=rbind(1:K,allPerms(1:K));
diff=array(0,dim(perm_set)[1]);
for (p in 1:dim(perm_set)[1])
{
  temp=docweights[,perm_set[p,]];
  diff[p]=fnorm(temp,omega_true);
}

p_star=which(diff==min(diff));
docweights=docweights[,perm_set[p_star,]];

barplot(t(docweights),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

