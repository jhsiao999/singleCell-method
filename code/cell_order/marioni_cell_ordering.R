

## Marioni single cell re-ordering

## We apply our cell re-ordering algorithm on Marioni data where cells come from G1,
## G2 and S phases, the phases being validated by FACS sorting.

library(data.table)
setwd('/Users/kushal/Documents/singleCell-method/project/R/cell_order')
G1_single <- data.frame(fread('../../data/Marioni_data/G1_singlecells_counts.txt'), row.names=1);
G2M_single <- data.frame(fread('../../data/Marioni_data/G2M_singlecells_counts.txt'), row.names=1);
S_single <- data.frame(fread('../../data/Marioni_data/S_singlecells_counts.txt'), row.names=1);

ercc_start <- grep("ERCC", rownames(G1_single))[1]

G1_single <- G1_single[-(ercc_start:dim(G1_single)[1]),-(1:3)];
G2M_single <- G2M_single[-(ercc_start:dim(G2M_single)[1]),-(1:3)];
S_single <- S_single[-(ercc_start:dim(S_single)[1]),-(1:3)];

pooled_data <- t(cbind(G1_single, S_single, G2M_single));

cell_cycle_genes <- as.vector(as.matrix((read_excel('../../data/Marioni_data/cellcycle_genes_mouse.xlsx'))))

matched_indices <- match(cell_cycle_genes,colnames(pooled_data))
cycle_counts_data <- pooled_data[,matched_indices];
cycle_counts_data <- cycle_counts_data[, -which(colSums(cycle_counts_data)==0)]
library(limma)
cycle_voom_data <- voom(cycle_counts_data)$E;
library(scales)
#cycle_data <- apply(cycle_voom_data,2, function(x) rescale(x, to=c(-1,1)));
cycle_data <- apply(cycle_voom_data, 2, function(x) return((x-mean(x))))
cell_phase_vector <- c(rep("G1", 96), rep("S", 96), rep("G2M", 96))

cycle_data.frame <- cbind.data.frame(cycle_data, cell_phase_vector);
library(dplyr)
cycle_grp_mean <- cycle_data.frame %>% group_by(cell_phase_vector) %>% summarise_each(funs(mean))
cycle_grp_mean <- data.frame(cycle_grp_mean);
cycle_labs <- cycle_grp_mean[apply(cycle_grp_mean[,-1], 2, function(x) which.max(x)),1];
## cycle_data is the data we shall use for cell re-ordering scheme

phase_in <- array(0, length(cycle_labs));
seq <- seq(0, 2*pi, length.out=5); seq1 <- seq[2:4];
seq1 <- c(1.5, 4, 5)
phase_in[which(cycle_labs=="G1")]=runif(length(which(cycle_labs=="G1")),pi/3,2*(pi/3));
phase_in[which(cycle_labs=="S")]=runif(length(which(cycle_labs=="S")),pi,4*(pi/3));
phase_in[which(cycle_labs=="G2M")]=runif(length(which(cycle_labs=="G2M")),5*(pi/3),2*pi);


source('utilities.R')
out1 <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100, 
                             save_path="../../rdas/cell_order_marioni_phase_fixed.rda",
                             fix.phase=TRUE, phase_in=phase_in)

out2 <- cell_reordering_phase(cycle_data, celltime_levels = 100, num_iter=100, 
                             save_path="../../rdas/cell_order_marioni.rda"
)

cell_times_1 <- out1$cell_times;

library(vioplot)
vioplot(cell_times_1[which(cell_phase_vector=="G1")],
        cell_times_1[which(cell_phase_vector=="S")],
        cell_times_1[which(cell_phase_vector=="G2M")],
        names=c("G1","S","G2M"),
        col="red")

cell_times_2 <- out2$cell_times;

library(vioplot)
vioplot(cell_times_2[which(cell_phase_vector=="G1")],
        cell_times_2[which(cell_phase_vector=="S")],
        cell_times_2[which(cell_phase_vector=="G2M")],
        names=c("G1","S","G2M"),
        col="red")



amp_genes <- out$amp;
sd_genes <- out$sigma;
phi_genes <- out2$phi;


plot(density(phi_genes), col="red", main="Density plot of the phases")
plot(density(amp_genes), col="red")


##  qtlcharts representation of the Marioni data
library(qtlcharts)
iplotCurves(t(cycle_data[,1:100]))
