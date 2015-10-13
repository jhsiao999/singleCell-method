
Poisson_topic.loglink <- function(counts, n_clus, lab_batch, max_iter=200, prior_scale_omega=3, use_squarem=FALSE)
{
  n_samples <- dim(counts)[1];
  n_genes <- dim(counts)[2];
  omega_initial <- matrix(rdirichlet(n_samples,rep(prior_scale_omega/n_clus,n_clus)),nrow=n_samples);
  alpha_initial <- matrix(0,nrow=n_clus, ncol=n_genes);
  beta_initial <- matrix(0, nrow=max(lab_batch),ncol=n_genes);
  param_vec_in <- c(as.vector(omega_initial), as.vector(alpha_initial), as.vector(beta_initial));
  counts_vec <- as.vector(counts);
  
  if(use_squarem==FALSE)
  {
    options(warn=-1)
    param_vec <- param_vec_in
    for(iter in 1:max_iter)
    {
      param_vec <- update_param_poisson_randbatch(param_vec,counts_vec,n_samples,n_genes,n_clus,lab_batch);
      print(length(param_vec))
      cat("Iteration number: ", iter)
      temp <- param_vec
    }
    
    omega_final=matrix(param_vec[1:(n_samples*n_clus)],n_samples,n_clus);
    param_expr <- param_vec[-(1:(n_samples*n_clus))];
    alpha_final=matrix(param_expr[1:(n_clus*n_genes)],n_clus,n_genes)
    beta_final =matrix(param_expr[-(1:(n_clus*n_genes))],max(lab_batch),n_genes);
    
  }
  
  if(use_squarem==TRUE)
  {
    system.time(res <- squarem(par=as.numeric(param_vec_in),
                               fixptfn=update_param_poisson_randbatch,
                               objfn= omega_loglik_poisson_full,
                               counts_vec = counts_vec,
                               n_samples = n_samples,
                               n_genes = n_genes,
                               n_clus = n_clus,
                               lab_batch = lab_batch,
                               control=list(maxiter = max_iter, trace = FALSE)));
    
    omega_final=matrix(res$par[1:(n_samples*n_clus)],n_samples,n_clus);
    res_expr <- res$par[-(1:(n_samples*n_clus))];
    alpha_final=matrix(res_expr[1:(n_clus*n_genes)],n_clus,n_genes)
    beta_final =matrix(res_expr[-(1:(n_clus*n_genes))],max(lab_batch),n_genes);
    
  }
 
  out <- list('omega'=omega_final,'alpha'=alpha_final,'beta'=beta_final);
  
  return(out)
  
}
