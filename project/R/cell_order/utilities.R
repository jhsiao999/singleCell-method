
## Bayesian cell re-ordering mechanism

## One iteration of cell reordering that takes the data as input and some
## current iterate of cell times and obtains new set of parameters of amplitude,
## phase and noise variation, and also refined estimates of cell times

atan3 <- function(beta2, beta1)
{
  if (beta1 > 0)
    v <- atan(beta2/beta1);
  if(beta2 >=0 & beta1 <0)
    v <- pi + atan(beta2/beta1);
  if(beta2 <0 & beta1 <0)
    v <- -pi + atan(beta2/beta1);
  if(beta2 >0 & beta1==0)
    v <- pi/2;
  if(beta2 <0 & beta1==0)
    v <- - (pi/2);
  if (v < 0)
    v <- v + 2*pi;
  return(v)
}

cell_reordering_iter <- function(cycle_data, cell_times_iter)
{
  if(length(unique(cell_times_iter))==1)
    stop("All the points have converged at same point on cycle");
  
  # cycle_data: a N \times G matrix, where N is number of cells, G number of genes
  # cell_times_iter:  the vector of cell times taken as input (a N \times 1)
  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];
  sigma <- array(0,G);
  amp <- array(0,G); phi <- array(0,G);
  
  for(g in 1:G)
  {
    fit <- lm(cycle_data[,g]  ~ sin(cell_times_iter) + cos(cell_times_iter) -1);
    sigma[g] <- sd(fit$residuals);
    beta1 <- fit$coefficients[1];
    beta2 <- fit$coefficients[2];
    amp[g] <- sqrt(beta1^2 + beta2^2);
    phi[g] <- atan3(as.numeric(beta2), as.numeric(beta1));
  }
  
  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  num_celltime_class <- length(cell_times_class);
  
  signal_intensity_per_class <- matrix(0, numcells, num_celltime_class)
  
  signal_intensity_per_class <- do.call(rbind,mclapply(1:numcells, function(cell) 
                        {
                            out <- array(0,length(cell_times_class));
                            for(times in 1:length(cell_times_class))
                            {
                                out[times] <- sum(mapply(dnorm, cycle_data[cell,], amp * sin(cell_times_class[times] + phi), sigma,log=TRUE));
                            }
                            return(out)
                        }, mc.cores=detectCores()));
  
  
  signal_intensity_class_exp <- do.call(rbind,lapply(1:dim(signal_intensity_per_class)[1], function(x) 
                                                                        {
                                                                            out <- exp(signal_intensity_per_class[x,]- max(signal_intensity_per_class[x,]));
                                                                            return(out)
                                                                        }));
  
  cell_times <- cell_times_class[unlist(lapply(1:dim(signal_intensity_class_exp)[1], function(x) 
                                                                                {
                                                                                  temp <- signal_intensity_class_exp[x,];
                                                                                  if(length(unique(signal_intensity_class_exp[x,]))==1)
                                                                                    out <- sample(1:dim(signal_intensity_class_exp)[2],1)
                                                                                  else
                                                                                    out <- which(rmultinom(1,1,signal_intensity_class_exp[x,])==1);
                                                                                  return(out)
                                                                                }))];

  out <- list("cell_times_iter"=cell_times, "amp_iter"=amp, "phi_iter"=phi, "sigma_iter"=sigma, "signal_intensity"=signal_intensity_per_class);
  return(out)
}

## calculate log likelihood under estimated cell times and amplitude
## phase and variation parameters 


loglik_cell_cycle <- function(cycle_data, cell_times, amp, phi, sigma)
{
  # cycle_data: a N \times G matrix, where N is number of cells, G number of genes
  # cell_times : a N \times 1 vector of cell times 
  # amp: the amplitude vector (G \times 1) over the genes 
  # phi: the G \times 1 vector of phase values over genes
  # sigma: the G \times 1 vector of gene variation
  
  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];
  sum <- 0;
  
  for(s in 1:numcells)
  {
    sum <- sum + sum(mapply(dnorm, cycle_data[s,],amp * sin(cell_times[s] + phi), sigma, log=TRUE));
  }
  
  return(sum)
}

## main workhorse function that takes in data and number of discrete levels 
## of cell times along with the number of iterations

cell_reordering_phase <- function(cycle_data, celltime_levels, num_iter)
{
  # cycle_data: a N \times G matrix, where N is number of cells, G number of genes
  # celltime_levels: number of discrete cell times used for estimation
  
  # We assume all the G genes are sinusoidal, if some are not, filter them
  
  G <- dim(cycle_data)[2];
  numcells <- dim(cycle_data)[1];
  
  celltimes_choice <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  cell_times_init <- sample(celltimes_choice, numcells, replace=TRUE);
  
  cell_times_iter <- cell_times_init;
  
  for(iter in 1:num_iter)
  {
    fun <- cell_reordering_iter(cycle_data, cell_times_iter);
    cell_times_iter <- fun$cell_times_iter;
    amp_iter <- fun$amp_iter;
    phi_iter <- fun$phi_iter;
    sigma_iter <- fun$sigma_iter;
    loglik_iter <- loglik_cell_cycle(cycle_data, cell_times_iter, amp_iter, phi_iter, sigma_iter);
    cat("The loglikelihood after iter", iter, "is:", loglik_iter,"\n")
  }
  
  out <- list("cell_times"=cell_times_iter, "amp"=amp_iter,"phi"=phi_iter, "sigma"=sigma_iter, "loglik"=loglik_iter)
  #save(out,file="../rdas/cell_order_ipsc_2.rda");
  return(out)
}

cell_reordering_full <- function(cycle_data, celltime_levels, cell_times, amp, phi, sigma)
{
  numcells <- dim(cycle_data)[1];
  G <- dim(cycle_data)[2];
  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1));
  sorted_cell_times_class <- sort(cell_times_class)
  order_class <- order(cell_times_class);
  
  signal_intensity <- do.call(rbind,mclapply(1:numcells, function(cell) 
  {
    out <- array(0,length(cell_times_class));
    for(times in 1:length(cell_times_class))
    {
      out[times] <- sum(mapply(dnorm, cycle_data[cell,], amp * sin(cell_times_class[times] + phi), sigma,log=TRUE));
    }
    return(out)
  }, mc.cores=detectCores()));
  
  
  signal_intensity <- signal_intensity[,order_class];
  
  cell_order_full <- array(0,numcells)
  
  for(cell in 1:numcells)
  {
    max_index <- which.max(signal_intensity[cell,]);
    if(max_index==1){
      denominator <- signal_intensity[cell,max_index] - signal_intensity[cell,(max_index+1)];
      numerator <- signal_intensity[cell,max_index] - signal_intensity[cell,celltime_levels];
      ratio <- numerator/(numerator+denominator);
      cell_order_full[cell] <- sorted_cell_times_class[celltime_levels] + ratio*4*pi/(celltime_levels-1);
    }else if(max_index==celltime_levels){
      denominator <- signal_intensity[cell,max_index] - signal_intensity[cell,1];
      numerator <- signal_intensity[cell,max_index] - signal_intensity[cell,(max_index-1)];
      ratio <- numerator/(numerator+denominator);
      cell_order_full[cell] <- sorted_cell_times_class[(max_index-1)] + ratio*4*pi/(celltime_levels-1);
    } else {
      denominator <- signal_intensity[cell,max_index] - signal_intensity[cell,(max_index+1)];
      numerator <- signal_intensity[cell,max_index] - signal_intensity[cell,(max_index-1)];
      ratio <- numerator/(numerator+denominator);
      cell_order_full[cell] <- sorted_cell_times_class[(max_index-1)] + ratio*4*pi/(celltime_levels-1);
    }
  }
  
  return(cell_order_full)
}

load_data <- get(load(file="../rdas/cell_order_ipsc.rda"));
cell_times <- load_data$cell_times;

phase_matching <- function(cell_times, cell_phase_vector)
{
  ## cell_times are the recovered times of the single cells from the single cell experiment
  ## cell_phase_vector is the vector of cell phases (G1, S, G2, M) etc for the cells
    # this can be obtained using the phase specific genes
  min_g1.s <- min(cell_times[which(cell_phase_vector=='G1.S')])
  cell_times_mod <- cell_times - min_g1.s;
  plot(cell_phase_vector, cell_times_mod)
  cos_cell_times <- tapply(cos(cell_times_mod), factor(cell_phase_vector), mean);
  sin_cell_times <- tapply(sin(cell_times_mod), factor(cell_phase_vector), mean);
  mean_angles <- array(0, length(sin_cell_times));
  for(l in 1:length(mean_angles)){
    mean_angles[l] <- atan3(sin_cell_times[l], cos_cell_times[l]);
  }
  plot(cell_phase_vector,cos(cell_times))
  plot(cell_phase_vector,sin(cell_times))
  plot(cell_phase_vector, cell_times)
  samp_recovered_order <- cbind.data.frame(rownames(cycle_data)[order(cell_times)])
  colnames(samp_recovered_order) = "recovered_order";
}



plot(amp_genes, amp_iter, col="red",xlab="true amplitudes", ylab="est amplitudes", main="amplitudes est, comparison")
plot(sigma_genes, sigma_iter, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")
plot(phi_genes, phi_iter, col="red",xlab="true phi", ylab="est phi", main="phase est, comparison");

library(plotrix)
radial.plot(lengths=1:length(cell_times_iter),radial.pos=cell_times_iter, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_iter)), lwd=2)
radial.plot(lengths=1:length(cell_times),radial.pos=cell_times, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times)), lwd=2)
radial.plot(lengths=1:length(sort(cell_order_full)),radial.pos=sort(cell_order_full), 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_order_full)), lwd=2)

