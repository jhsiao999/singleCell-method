
## An example simulation run for the cell cycle ordering mechanism

G <- 500;
num_cells <- 400;
amp_genes <- rep(10, G);
phi_genes <- runif(G, 0, 2*pi)
sigma_genes <- rchisq(G, 4);
cell_times <- seq(0,2*pi, 2*pi/(num_cells-1));

cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times);

celltime_levels <- 100;

molecules_single_cell_cycle <- read.table("../data/molecules_ipsc_single_cell_cycle.txt");
cycle_data <- t(molecules_single_cell_cycle);
cycle_data <- cycle_data[, which(colSums(cycle_data)!=0)]

## compare the SNR for Yoav's cell cycle data

load_data <- get(load(file="../rdas/cell_order_ipsc.rda"));
amp_genes <- load_data$amp;
sd_genes <- load_data$sigma;

ESS <- amp_genes^2; RSS <- sd_genes^2

SNR <- ESS/RSS;

plot(SNR, xlab="genes", ylab="SNR", main="Plot of SNR across all genes")

SNR_indices <- which(SNR > 10);

cycle_data <- cycle_data[,SNR_indices];
