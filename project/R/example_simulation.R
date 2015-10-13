
###  Building an example counts table

source('simulate_counts_poisson.R')

N_sim = 500;
N_genes =100;
N_clus =4;

alpha_true=matrix(rnorm(N_clus*N_genes,0.5,1),nrow=N_clus); ### the matrix of fixed effects

lab_batch=c(rep(1,N_sim/2),rep(2,N_sim/2));

B= max(lab_batch);

sigmab_true=2;

beta_true=matrix(0,B,N_genes);       ###  the matrix of the random effect

for(g in 1:N_genes)
{
  beta_true[,g]=rnorm(B,mean=0,sd=sigmab_true);
}

library(gtools)
T=N_sim/5;
omega_true=matrix(rbind(rdirichlet(T,c(3,4,2,6)),
                        rdirichlet(T,c(1,4,6,3)),
                        rdirichlet(T,c(4,1,2,2)),
                        rdirichlet(T,c(2,6,3,2)),
                        rdirichlet(T,c(3,3,5,4))), nrow=N_sim);


counts <- simulate_poisson_random_batch(N_sim, N_genes, N_clus, omega_true, alpha_true, beta_true, lab_batch)


library(maptpx)
N_clus=4
Topic_Clus=topics(counts,N_clus,kill=0,tol=0.001);
docweights_topics=Topic_Clus$omega;
barplot(t(docweights_topics),col=2:(N_clus+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",N_clus),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

barplot(t(omega_true),col=2:(N_clus+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",N_clus),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)



#barplot(t(omega_final),col=2:(K+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",k),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)


#omega_out <- as.matrix(read.table('omega_temp.txt'));
#barplot(t(omega_out),col=2:(n_clus+1),axisnames=F,space=0,border=NA,main=paste("No. of clusters=",N_clus),las=1,ylim=c(0,1),cex.axis=1.5,cex.main=1.4)

