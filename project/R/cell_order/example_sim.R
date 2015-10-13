
## An example simulation run for the cell cycle ordering mechanism

G <- 500;
num_cells <- 400;
amp_genes <- rep(10, G);
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4);
cell_times <- seq(0,2*pi, 2*pi/(num_cells-1));

cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times);

celltime_levels <- 100;

