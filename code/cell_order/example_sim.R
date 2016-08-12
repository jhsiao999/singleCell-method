
## An example simulation run for the cell cycle ordering mechanism

G <- 500;
num_cells <- 400;
amp_genes <- rep(10, G);
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4);
cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);

cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);

celltime_levels <- 100;

out <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100)

  
plot(amp_genes, out$amp, col="red",xlab="true amplitudes", ylab="est amplitudes", main="amplitudes est, comparison")
plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")
plot(phi_genes, out$phi, col="red",xlab="true phi", ylab="est phi", main="phase est, comparison");

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(out$cell_times),radial.pos=out$cell_times[order(cell_times_sim)], 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(out$cell_times)), lwd=2)
radial.plot(lengths=1:length(cell_times_sim),radial.pos=sort(cell_times_sim), 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)


## Second example simulation run

G <- 600;
num_cells <- 400;
amp_genes <- rep(10, G);
sigma_genes <- rchisq(G, 4);
cell_times_sim1<- sample(seq(pi/3,2*(pi/3), 2*pi/(num_cells/2-1)), num_cells/2, replace=TRUE);
cell_times_sim2<- sample(seq(5*(pi/3),6*(pi/3), 2*pi/(num_cells/2-1)), num_cells/2, replace=TRUE);
cell_times_sim3<- sample(seq(3*(pi/3),4*(pi/3), 2*pi/(num_cells/2-1)), num_cells/2, replace=TRUE);

cell_times_sim <- c(cell_times_sim1, cell_times_sim2, cell_times_sim3);
radial.plot(lengths=1:length(cell_times_sim),radial.pos=cell_times_sim, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)
#phi_genes <- c(runif(G/3, pi/4, 3*(pi/4)), runif(G/3, 5*(pi/4), 7*(pi/4)));
phi_genes <- c(runif(G/3, pi/3, 2*(pi/3)), runif(G/3, 3*(pi/3), 4*(pi/3)),
               runif(G/3, 5*(pi/3), 2*pi));

source('cell_cycle_sim.R')
cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);

celltime_levels <- 100;

out <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100)

plot(amp_genes, out$amp, col="red",xlab="true amplitudes", ylab="est amplitudes", main="amplitudes est, comparison")
plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")
plot(phi_genes, out$phi, col="red",xlab="true phi", ylab="est phi", main="phase est, comparison");

library(plotrix)
library(RColorBrewer)

radial.plot(lengths=1:length(cell_times_sim),radial.pos=out$cell_times, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)

phi2_genes <- phi_genes;
phi2_genes <- c(rep(pi/2,G/3), rep(pi, G/3), rep(3 *pi/2, G/3));
phi2_genes <- c(runif(G/2, pi/4, 3*(pi/4)), runif(G/2, 5*(pi/4), 7*(pi/4)));
phi2_genes <- runif(G,0,2*pi);

out <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100,
                             fix.phase=TRUE, phase_in = phi2_genes)

plot(amp_genes, out$amp, col="red",xlab="true amplitudes", ylab="est amplitudes", main="amplitudes est, comparison")
plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")
plot(phi_genes, out$phi, col="red",xlab="true phi", ylab="est phi", main="phase est, comparison");

library(plotrix)
library(RColorBrewer)

radial.plot(lengths=1:length(cell_times_sim),radial.pos=out$cell_times, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)

#phi2_genes <- c(rep((pi/2),G/2),rep(3*(pi/2),G/2));
phi2_genes <- c(runif(G/2, (pi/3), 2*(pi/3)), runif(G/2, pi, 4*(pi/3)));

out <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100,
                             fix.phase=TRUE, phase_in = phi2_genes)
plot(amp_genes, out$amp, col="red",xlab="true amplitudes", ylab="est amplitudes", main="amplitudes est, comparison")
plot(sigma_genes, out$sigma, col="red",xlab="true sigma", ylab="est sigma", main="sigma(variation) est, comparison")
plot(phi_genes, out$phi, col="red",xlab="true phi", ylab="est phi", main="phase est, comparison");

library(plotrix)
library(RColorBrewer)

radial.plot(lengths=1:length(cell_times_sim),radial.pos=out$cell_times, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times_sim)), lwd=2)
