
simulate_poisson_random_batch = function(n_samples,
                                          n_genes, 
                                          n_clus,
                                          topic_omega_sim,
                                          topic_alpha_sim,
                                          topic_beta_sim,
                                          label_batch_sim)
{
  
  # n_samples : the number of tissue or single cell samples
  # n_genes : the number of genes
  # n_clus: the number of topics/clusters
  # topic_omega_sim: omega (n_samples \times n_clus ) matrix of mixing proportions
  # topic_alpha_sim: alpha (n_clus \times n_genes) matrix of topic expressions
  # topic_beta_sim : beta (num_batches \times n_genes) matrix of batch effects
  # label_batch_sim: the labels of the batches startinf grom 1 to num_batches
  
  flag <- 0;
  if(n_clus !=dim(topic_omega_sim)[2]){
    print ("The number of clusters does not match with dimension of omega matrix");
    flag <- 1}
  if(n_clus !=dim(topic_alpha_sim)[1]){
    print ("The number of clusters does not match with dimension of alpha matrix");
    flag <- 1}
  if(max(label_batch_sim) !=dim(topic_beta_sim)[1]){
    print ("The number of batches does not match with the label vector");
    flag <- 1}
  if(flag==0) {
    out_temp=exp(topic_omega_sim%*%topic_alpha_sim +topic_beta_sim[label_batch_sim,]);
    out_temp_vec <- as.vector(out_temp);
    read_counts=matrix(unlist(lapply(1:length(out_temp_vec),function(x) rpois(1,out_temp_vec[x]))),nrow=n_samples);
    return(read_counts)}
  
}



