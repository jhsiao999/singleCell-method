
## Yoav iPSC data

## compare the SNR for Yoav's cell cycle data

molecules_single_cell_cycle <- read.table("../data/molecules_ipsc_single_cell_cycle.txt");
cycle_counts_data <- t(molecules_single_cell_cycle);
cycle_counts_data <- cycle_counts_data[, which(colSums(cycle_counts_data)!=0)];

cycle_voom_data <- voom(cycle_counts_data)$E;
cycle_data <- apply(cycle_voom_data,2, function(x) x-mean(x));

cell_phase_vector <- as.vector(as.matrix(read.table("../data/cell_phase_vector_yoav.txt")));

cycle.pca <- princomp(cycle_data, scale.=TRUE, center=TRUE);
plot(cycle.pca, type="l", col="red", main="scree plot for cycle data")
library(rgl)
plot(cycle.pca$scores[,1],cycle.pca$scores[,2], col=as.factor(cell_phase_vector), 
     lwd=3, lty=1, pch=20, xlab="PC1", ylab="PC2", main="PC plot")
legend("topleft", legend=levels(as.factor(cell_phase_vector)), col=1:length(as.factor(cell_phase_vector)),
       pch=20)
plot(cycle.pca$scores[,2],cycle.pca$scores[,3], col=as.factor(cell_phase_vector), 
     lwd=3, lty=1, pch=20, xlab="PC2", ylab="PC3", main="PC plot")
legend("topleft", legend=levels(as.factor(cell_phase_vector)), col=1:length(as.factor(cell_phase_vector)),
       pch=20)
source('/Users/kushal/Documents/singleCell-method/project/R/cell_order/utilities.R')
out <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100, save_path="../rdas/cell_order_ipsc_2.rda")

rank_pca <- rank(cycle.pca$scores[,1]);

load_data <- get(load(file="../rdas/cell_order_ipsc.rda"));

cell_times <- load_data$cell_times;

rank_celltimes <- rank(cell_times);

plot(rank_celltimes, rank_pca, xlab="Rank (ordering)",
     ylab="Rank (PCA)", col="red")
amp_genes <- load_data$amp;
sd_genes <- load_data$sigma;
phi_genes <- load_data$phi;

ESS <- amp_genes^2; RSS <- sd_genes^2

SNR <- ESS/RSS;

plot(SNR, xlab="genes", ylab="SNR", main="Plot of SNR across all genes")


cell_cycle_genes <- read.table("../data/cellcyclegenes.txt", header = TRUE, sep="\t")

## create 5 lists of 5 phases (de-level and then remove "")
cell_cycle_genes_list <- lapply(1:5,function(x){
  temp <- as.character(cell_cycle_genes[,x])
  temp[temp!=""]
})


phi_vec <- array(0, 0);
factor_level <- array(0,0);
cell_phase <- colnames(cell_cycle_genes)
cell_phase <- factor(cell_phase, 
                     levels = c("G1.S", "S", "G2.M", "M", "M.G1"))

for(l in 1:5){
  labs <- match(cell_cycle_genes_list[[l]], colnames(cycle_data));
  labs <- labs[which(!is.na(labs))]
  phi_vec <- c(phi_vec,phi_genes[labs]);
  factor_level <- c(factor_level, rep(toString(cell_phase[l]), length(phi_genes[labs])));
}

phi_vec1 <- unique(phi_vec);
factor_vec <- factor_level[match(phi_vec1,phi_vec)]
par(mar=c(4,1,1,1))
plot(ordered(factor_vec, c("G1.S", "S", "G2.M", "M", "M.G1")), phi_vec1,
     xlab="factor_levels", ylab="phase (radians)", 
     main="Phase plot across cell phases", col="red")

amp_vec <- array(0, 0);
factor_level <- array(0,0);

for(l in 1:5){
  labs <- match(cell_cycle_genes_list[[l]], colnames(cycle_data));
  labs <- labs[!is.na(labs)]
  print(length(labs))
  amp_vec <- c(amp_vec,amp_genes[labs]);
  factor_level <- c(factor_level, rep(toString(cell_phase[l]), length(amp_genes[labs])));
}

amp_vec1 <- unique(amp_vec);
factor_vec <- factor_level[match(amp_vec1,amp_vec)]

par(mar=c(4,1,1,1))
plot(ordered(factor_vec, c("G1.S", "S", "G2.M", "M", "M.G1")), amp_vec1,
     xlab="factor_levels", ylab="amplitudes", 
     main="Amplitude plot across cell phases", col="red")

## SNR of genes across the cell phases

snr_vec <- array(0,0)
factor_level <- array(0,0)

for(l in 1:5){
  labs <- match(cell_cycle_genes_list[[l]], colnames(cycle_data));
  labs <- labs[!is.na(labs)]
  print(length(labs))
  snr_vec <- c(snr_vec,SNR[labs]);
  factor_level <- c(factor_level, rep(toString(cell_phase[l]), length(SNR[labs])));
}

snr_vec1 <- unique(snr_vec);
factor_vec <- factor_level[match(snr_vec1,snr_vec)]
par(mar=c(4,1,1,1))
plot(ordered(factor_vec, c("G1.S", "S", "G2.M", "M", "M.G1")), snr_vec1,
     xlab="factor_levels", ylab="SNR", 
     main="SNR plot across cell phases", col="red")


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

