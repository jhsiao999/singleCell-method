#' @title Simulate sinusoidal gene expression along cell phases on cell cycle
#'
#' @param num_genes The number of sinusoidal genes to simulate
#' @param amp_genes The amplitude vector of the sinusoidal genes (of length equal to num_genes)
#' @param phi_genes The phase angle vector of the sinusoidal genes (of length equal to num_genes)
#' @param sigma_genes The noise variation vector of the sinusoidal genes (of length equal to num_genes)
#' @param cell_times The phases of the cells on the cell cycle (a vector of values between 0 to 360 degrees)
#'
#' @description The function simulates sinusoidal gene patterns over a cell cycle with the cell phases
#'              provided by the user.
#'
#'  @author  Kushal K Dey
#'
#'  @export
#'
#'
sim_sinusoidal_cycle <- function(num_genes, amp_genes, phi_genes, sigma_genes, cell_times)
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
