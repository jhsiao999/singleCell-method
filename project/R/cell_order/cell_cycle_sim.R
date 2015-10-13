

## Simulate a model of the cell cycle patterns for oscillatory/sinusoidal genes 


sim_sinusoidal_cycle <- function(G, amp_genes, phi_genes, sigma_genes, cell_times)
{
  # G: the number of genes in the model
  # amp_genes : a G \times 1 vector of amplitudes of the genes
  # phi_genes : a G \times 1 vector of phases of the genes 
  # sigma_genes : a G \times 1 vector of variability of the genes 
  # cell_times : a N \times 1 vector (N being number of cells)
  ##             the time on (0, 2\pi) from where the cell was obtained in cell cycle
  
  signal <- matrix(0, length(cell_times), G);
  for(s in 1:length(cell_times))
  {
    signal[s,] <- mapply(rnorm, 1, amp_genes * sin(cell_times[s] + phi_genes), sigma_genes);
  }
  return(signal)
}