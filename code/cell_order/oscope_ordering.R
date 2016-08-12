
## Botstein Data analysis

setwd("/Users/kushal/Documents/singleCell-method/project/analysis/")
library(data.table)
data <- read.csv("../data/Oscope data/GSE64016_H1andFUCCI_normalized_EC.csv")
gene_names <- data[,1];
data <- data[,-1];
cell_phases <- sapply(colnames(data), function(x) strsplit(x,"_")[[1]][1]);
cell_data <- as.matrix(data[-1,which(cell_phases != "H1")]);
cycle_data <- t(cell_data);

cycle_data_norm <- apply(cycle_data,2,function(x)  return (x-mean(x))/sd(x))

library(cellcycleR)

celltime_levels <- 100;
library(parallel)
cycle_data_norm <- cycle_data_norm[, -which(colSums(cycle_data_norm)==0)]
out <- cell_ordering_class(cycle_data_norm, celltime_levels = 100, num_iter=50)

save(out,file="../rdas/Botstein_cell_cycle.rdata");

cell_order_full <- cell_ordering_full(out$signal_intensity, dim(cycle_data)[2])

cell_times <- out$cell_times;

library(plotrix)
library(RColorBrewer)
radial.plot(lengths=1:length(cell_times),radial.pos=cell_times, 
            line.col=colorRampPalette(brewer.pal(9,"Blues"))(length(cell_times)), lwd=2)

cell_phases <- cell_phases[which(cell_phases != "H1")];
library(vioplot)
vioplot((cell_times[which(cell_phases=="G1")]-1.5)%%6,
        (cell_times[which(cell_phases=="S")]-1.5)%%6,
        (cell_times[which(cell_phases=="G2")]-1.5)%%6,
        names=c("G1","S","G2"),
        col="red")

cell_times_ordered <- c(cell_times[which(cell_phases=="G1")], 
                        cell_times[which(cell_phases=="S")],
                        cell_times[which(cell_phases=="G2")]);

radial.plot(lengths=1:length(cell_times),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="G1"))),
                       rep("blue",length(which(cell_phases=="S"))),
                       rep("green",length(which(cell_phases=="G2")))),lwd=2)

cell_times_ordered <- c(cell_times[which(cell_phases=="G1")],
                        cell_times[which(cell_phases=="S")]);

radial.plot(lengths=1:length(cell_times_ordered),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="G1"))),
                       rep("green",length(which(cell_phases=="S")))),lwd=2)

cell_times_ordered <- c(cell_times[which(cell_phases=="G1")],
                        cell_times[which(cell_phases=="G2")]);

radial.plot(lengths=1:length(cell_times_ordered),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="G1"))),
                       rep("green",length(which(cell_phases=="G2")))),lwd=2)

cell_times_ordered <- c(cell_times[which(cell_phases=="S")],
                        cell_times[which(cell_phases=="G2")]);

radial.plot(lengths=1:length(cell_times_ordered),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="S"))),
                       rep("green",length(which(cell_phases=="G2")))),lwd=2)


amp_genes <- out$amp;
sd_genes <- out$sigma;
phi_genes <- out$phi;

plot(density(phi_genes), col="red", main="Density plot of the phases")
plot(density(amp_genes[-which.max(amp_genes)]), col="red", main="Density plot of the amplitudes")
plot(density(sd_genes[-which.max(sd_genes)]), col="red", main="Density plot of the standard deviation")

ESS <- amp_genes^2; RSS <- sd_genes^2

SNR <- ESS/RSS;

plot(SNR)

snr_high_indices <- which(SNR > 1);

cycle_data_norm_sinusoidal <- cycle_data_norm[,snr_high_indices];
out2 <- cell_ordering_class(cycle_data_norm_sinusoidal, celltime_levels = 100, num_iter=100)
cell_times <- out2$cell_times;

library(vioplot)
vioplot(cell_times[which(cell_phases=="G1")],
        cell_times[which(cell_phases=="S")],
        cell_times[which(cell_phases=="G2")],
        names=c("G1","S","G2"),
        col="red")

cell_times_ordered <- c(cell_times[which(cell_phases=="G1")],
                        cell_times[which(cell_phases=="S")]);

radial.plot(lengths=1:length(cell_times_ordered),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="G1"))),
                       rep("green",length(which(cell_phases=="S")))),lwd=2)

cell_times_ordered <- c(cell_times[which(cell_phases=="G1")],
                        cell_times[which(cell_phases=="G2")]);

radial.plot(lengths=1:length(cell_times_ordered),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="G1"))),
                       rep("green",length(which(cell_phases=="G2")))),lwd=2)

cell_times_ordered <- c(cell_times[which(cell_phases=="S")],
                        cell_times[which(cell_phases=="G2")]);

radial.plot(lengths=1:length(cell_times_ordered),radial.pos=cell_times_ordered, 
            line.col=c(rep("red",length(which(cell_phases=="S"))),
                       rep("green",length(which(cell_phases=="G2")))),lwd=2)
amp_genes <- out2$amp;
sd_genes <- out2$sigma;
phi_genes <- out2$phi;

plot(density(phi_genes), col="red", main="Density plot of the phases")

###  Top genes expression profile under ordered cell times 

top_genes <- which(SNR > 0.5);
library(qtlcharts)
iplotCurves(t(cycle_data_norm[order(cell_order_full),top_genes]))
iplotCurves(t(cycle_data_norm[c(which(cell_phases=="G1"),which(cell_phases=="S"), which(cell_phases=="G2")),top_genes]))
