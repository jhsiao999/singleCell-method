
source("../../../classtpx/R/class_count.R")
source("../../../classtpx/R/class_topics.R")
source("../../../classtpx/R/class_tpx.R")
source("../../../classtpx/R/class_varselect.R")

quake_data <- get(load("../data/Quake_brain_cell/Quake_brain_data.rda"))
quake_metadata <- read.csv("../data/Quake_brain_cell/All_cell_info_brain.csv", sep=";")

samp_names_data <- colnames(quake_data)
samp_names_metadata <- as.character(quake_metadata[,1])

samp_names_data_filt <- unlist(lapply(samp_names_data, function(x) strsplit(x, "_")[[1]][2]))

metadata_indices <- match(samp_names_data_filt, samp_names_metadata);
quake_metadata <- quake_metadata[metadata_indices,];

cell_type <- quake_metadata$Cell_type

gene_names_quake <- rownames(quake_data);

### get Ensembl IDs from the gene names

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
listAttributes(mart)[1:20,]
results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 filters = "external_gene_name", values = as.vector(gene_names_quake),
                 mart = mart)

common_genes <- intersect(gene_names_quake, results[,2])
results_filtered <- results[match(common_genes, results[,2]),]

results_filtered <- results[complete.cases(results[match(gene_names_quake, results[,2]),]),]
match(results_filtered[,2], gene_names_quake)

library(data.table)
gtex_data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"));
matdata <- gtex_data[,-(1:2)];
gene_names <- as.vector(sapply(gtex_data[,2], function(x) substring(x,1,15)))
tissue_labels=read.table("../data/GTEX_V6/samples_id.txt")[,3];


common_gtex_genes <- intersect(results_filtered[,1], gene_names)

results_filtered_filtered <- results_filtered[match(common_gtex_genes, results_filtered[,1]),]

match(common_gtex_genes, results_filtered_filtered[,1])

matdata_filtered <- matdata[match(results_filtered_filtered[,1],gene_names),];
quake_data_filtered <- quake_data[match(results_filtered_filtered[,2],gene_names_quake),]

matdata_filtered_brain <- matdata_filtered[,grep("Brain",tissue_labels)];

library(dplyr)
cell_type_filtered <- cell_type[cell_type %in% c("Astrocytes", "Microglia", "OPC", "Oligodendrocytes", "Endothelial", "Neurons")]
quake_data_filtered_filtered <- quake_data_filtered[,cell_type %in% c("Astrocytes", "Microglia", "OPC", "Oligodendrocytes", "Endothelial", "Neurons")]

pooled_data <- t(cbind(as.matrix(quake_data_filtered_filtered), as.matrix(matdata_filtered_brain)));

library(maptpx)
library(classtpx)
source("../../../classtpx/R/class_count.R")
source("../../../classtpx/R/class_topics.R")
source("../../../classtpx/R/class_tpx.R")
source("../../../classtpx/R/class_varselect.R")

class_labs <- as.numeric(cell_type_filtered)
known_samples <- 1:length(class_labs)

Topic_clus <- class_topics(
  pooled_data, 
  K=length(unique(class_labs)), 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)

save(Topic_clus, file="../rdas/topic_clus_classtpx_quake_thetafix_shrink.rda")

omega <- Topic_clus$omega;
pooled_types_labels <- c(as.character(cell_type_filtered),as.character(tissue_labels[grep("Brain",tissue_labels)]))

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = pooled_types_labels
)

rownames(omega) <- annotation$sample_id;
CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(4, "Spectral")),
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


##################   Brain Cortex Data  ########################

cortex_indices <- tissue_labels %in% 
  c("Brain - Anterior cingulate cortex (BA24)","Brain - Frontal Cortex (BA9)","Brain - Cortex")

cortex_labels <- tissue_labels[cortex_indices]

cortex_data <- matdata_filtered[,cortex_indices]
quake_cortex_labels <- which(quake_metadata[,2]=="cortex")
cell_type_cortex <- cell_type[quake_cortex_labels]
cell_type_cortex_filtered <- cell_type_cortex[cell_type_cortex %in% c("Astrocytes", "Microglia", "OPC", "Oligodendrocytes", "Endothelial", "Neurons")]

quake_data_cortex_filtered <- quake_data_filtered[,quake_cortex_labels];
quake_data_cortex_filtered_filtered <- quake_data_cortex_filtered[,which(cell_type_cortex %in% c("Astrocytes", "Microglia", "OPC", "Oligodendrocytes", "Endothelial", "Neurons"))]

pooled_data_cortex_filtered_filtered <- cbind(as.matrix(quake_data_cortex_filtered_filtered), as.matrix(cortex_data))

pooled_labels <- c(as.character(cell_type_cortex_filtered), as.character(cortex_labels))

cell_type_cortex_filtered <- droplevels(cell_type_cortex_filtered)

class_labs <- as.numeric(cell_type_cortex_filtered);
known_samples <- 1:length(class_labs);

Topic_clus <- class_topics(
  t(pooled_data_cortex_filtered_filtered), 
  K=length(unique(cell_type_cortex_filtered)), 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)

save(Topic_clus, file="../rdas/topic_clus_classtpx_quake_thetafix_noshrunk.rda")
omega <- Topic_clus$omega;

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = pooled_labels
)

rownames(omega) <- annotation$sample_id;
CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(4, "Spectral")),
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))



###############   Brain Hippocampus Study  #######################

hippo_indices <- tissue_labels %in% 
  c("Brain - Hippocampus")

hippo_labels <- tissue_labels[hippo_indices]

hippo_data <- matdata_filtered[,hippo_indices]
quake_hippo_labels <- which(quake_metadata[,2]=="hippocampus")
cell_type_hippo <- cell_type[quake_hippo_labels]
cell_type_hippo_filtered <- cell_type_hippo[cell_type_hippo %in% c("Astrocytes", "Microglia", "OPC", "Oligodendrocytes", "Endothelial", "Neurons")]

quake_data_hippo_filtered <- quake_data_filtered[,quake_hippo_labels];
quake_data_hippo_filtered_filtered <- quake_data_hippo_filtered[,cell_type_hippo %in% c("Astrocytes", "Microglia", "OPC", "Oligodendrocytes", "Endothelial", "Neurons")]

pooled_data_hippo_filtered_filtered <- cbind(as.matrix(quake_data_hippo_filtered_filtered), as.matrix(hippo_data))

pooled_labels <- c(as.character(cell_type_hippo_filtered), as.character(hippo_labels))

cell_type_hippo_filtered <- droplevels(cell_type_hippo_filtered)

class_labs <- as.numeric(cell_type_hippo_filtered);
known_samples <- 1:length(class_labs);

Topic_clus <- class_topics(
  t(pooled_data_hippo_filtered_filtered), 
  K=length(unique(cell_type_hippo_filtered)), 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)

save(Topic_clus, file="../rdas/topic_clus_classtpx_quake_thetafix_noshrunk_hippo.rda")
omega <- Topic_clus$omega;

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = pooled_labels
)

rownames(omega) <- annotation$sample_id;
CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = c(RColorBrewer::brewer.pal(8, "Accent"),RColorBrewer::brewer.pal(4, "Spectral")),
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


