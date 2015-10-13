
###  EM update iterative step:


library(parallel)
library(optimx)


update_param_poisson_randbatch <- function(param_vec,counts_vec,n_samples,n_genes,n_clus,lab_batch)
{
  counts=matrix(counts_vec,n_samples,n_genes);
  omega0=matrix(param_vec[1:(n_samples*n_clus)],n_samples,n_clus);
  B <- max(lab_batch);
  
##############  Estimating the alpha and beta ####################
  
  alpha=matrix(0,n_clus,n_genes);
  beta=matrix(0,B,n_genes);
  lab_batch_fac <- as.factor(lab_batch);
  system.time(out_parallel <- lapply(1:n_genes, 
                                       function(g) glmer (counts[,g] ~ omega0 + (1|lab_batch_fac)-1, family=poisson(),control=glmerControl(optCtrl=list(maxfun=300)))
                                       ));
  
# system.time(out_parallel <- lapply(1:n_genes, 
#                                      function(g) glm (counts[,g] ~ omega0 + as.factor(lab_batch_fac)-1, family=poisson(),contrasts = list(lab_batch_fac = "contr.sum"))
#                                          ));
  
# system.time(out_serial <- lapply(1:n_genes, 
#                              function(g) glmer (counts[,g] ~ omega0 + (1|lab_batch_fac)-1, family=poisson())
#                             ));
  
 

  alpha_out <- do.call(cbind,mclapply(1:n_genes, function(g) as.numeric(fixef(out_parallel[[g]]))));
  beta_out <-  do.call(cbind,mclapply(1:n_genes, function(g) as.numeric(as.matrix(ranef(out_parallel[[g]])$lab_batch_fac))))                          
  
 # alpha_out <- do.call(cbind,mclapply(1:n_genes, function(g) as.numeric(out_parallel[[g]]$coefficients[1:n_clus])));
#  beta_out <- do.call(cbind,mclapply(1:n_genes, function(g) as.numeric(out_parallel[[g]]$coefficients[-(1:n_clus)])));
  
  
  system.time(omega_out <- do.call(rbind,mclapply(1:n_samples,
                        FUN=function(n) 
                        {
                          res <- optim(reverse_transform(omega0[n,]),function(v) omega_loglik_poisson_onesample(transform(v),counts[n,],alpha_out,beta_out,lab_batch[n]))
                          #return(transform(as.numeric(res[,1:(n_clus-1)])))
                          return(transform(res$par))
                        },
                        mc.set.seed = TRUE,
                        mc.cores = detectCores())));
  
  param_vec_out=c(as.vector(omega_out),as.vector(alpha_out),as.vector(beta_out));
  return(param_vec_out)
}
