
##  Log likelihood  of the Poisson Random Batch Model for each tissue/single cell sample

omega_loglik_poisson_onesample = function(omega_vec,counts_samp,alpha_mat,beta_mat,lab_batch_sample)
{
  out <- exp(omega_vec%*%alpha_mat +beta_mat[lab_batch_sample,]);
  sum_loglik <- sum(unlist(lapply(1:length(counts_samp),function(g) out[g] -counts_samp[g]*log(out[g]+1e-12))));
  return(sum_loglik);
}

omega_loglik_poisson_full = function(param_vec_in,counts_vec,n_samples, n_genes, n_clus, lab_batch)
{
  omega=matrix(param_vec_in[1:(n_samples*n_clus)],n_samples,n_clus);
  param_vec_expr <- param_vec_in[-(1:(n_samples*n_clus))];
  alpha= matrix(param_vec_expr[1:(n_clus*n_genes)],n_clus,n_genes);
  beta = matrix(param_vec_expr[-(1:(n_clus*n_genes))],max(lab_batch),n_genes);
  
  out <- sum(unlist(mclapply(1:n_samples, 
                  function(n) omega_loglik_poisson_onesample(omega[n,],counts[n,],alpha,beta,lab_batch[n]),
                  mc.set.seed = TRUE,
                  mc.cores = detectCores())));
  return(out)
    
}

## an example 

# omega_loglik_Poisson(counts[1,],omega_true[1,],alpha_true,beta_true,1)