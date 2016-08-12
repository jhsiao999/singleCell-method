
## mash prep data

library(data.table)
library(e1071)
library(classtpx)
library(maptpx)

gtex_data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"));
matdata <- gtex_data[,-(1:2)];

tissue_labels=read.table("../data/GTEX_V6/samples_id.txt")[,3];

aorta_labels <- which(as.character(tissue_labels) == 'Artery - Aorta')
coronary_labels <- which(as.character(tissue_labels) == 'Artery - Coronary');
tibial_labels <- which(as.character(tissue_labels) == 'Artery - Tibial');
blood_labels <- which(as.character(tissue_labels) == 'Whole Blood')
uterus_labels <- which(as.character(tissue_labels) == 'Uterus')


pooled_data <- t(cbind(matdata[,aorta_labels],
                       matdata[,coronary_labels],
                       matdata[,tibial_labels],
                       matdata[,blood_labels],
                       matdata[,uterus_labels]))

class_labs <- rep(1:5, each=50);
length_vec <- c(length(aorta_labels), length(coronary_labels),
                length(tibial_labels), length(blood_labels),
                length(uterus_labels))
cum_length_vec <- cumsum(length_vec)

known_samples <- as.vector(outer(1:50, c(0, cum_length_vec[1:(length(cum_length_vec)-1)]), "+"))

voom2 <- function(counts){
  libsize.mat <- rep.col(rowSums(counts), dim(counts)[2]);
  voom.out <- log((counts+0.5), base=2) - log((libsize.mat+1), base=2)+ 6* log(10, base=2);
  return(voom.out)
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}



counts_class <- pooled_data[known_samples,];
voom_class <- voom2(counts_class);
model_mat <- model.matrix(~as.factor(class_labs)-1) 

beta_class <- matrix(0, dim(voom_class)[2], 5);
sebeta_class <- matrix(0, dim(voom_class)[2], 5)
for(k in 1:5){
  model_mat_temp <- cbind(rep(1,dim(voom_class)[1]), model_mat[,k]);
  limma.obj <- limma::lmFit(t(voom_class), model_mat_temp)
  limma.obj <- limma::eBayes(limma.obj)
  mean_genes_limma <- apply(limma.obj$coefficients, 1, mean)
  beta_class[,k] <- as.matrix(limma.obj$coefficients[,-1]);
  sebeta_class[,k] <- limma.obj$sigma*(as.matrix(limma.obj$stdev.unscaled[,-1]));
}

mash_beta <- read.table("../analysis/EE_biggridposterior.means.txt")[,-1];

lib_size <- rowSums(counts_class);
voom_class <- voom2(counts_class);
mean_voom_features <- apply(voom_class, 2, mean);

voom_shrunk_mean <- mean_voom_features + cbind.data.frame(mash_beta);
voom_shrunk_class_2 <- matrix(0, dim(counts_class)[1], dim(counts_class)[2])

for(i in 1:length(unique(class_labs))){
  voom_shrunk_class_2[which(class_labs==unique(class_labs)[i]),] <- rep.row(as.vector(voom_shrunk_mean[,i]), length(which(class_labs==unique(class_labs)[i])));
}


counts_shrunk_matrix <- (2^{voom_shrunk_class_2 - 6*log(10, base=2)})*(rep.col(lib_size+1, dim(voom_shrunk_class_2)[2])) - 0.5
counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;

mean_counts_shrunk_class <- do.call(rbind, lapply(1:dim(counts_shrunk_matrix)[2], function(l)
{
  mean_element <- tapply(counts_shrunk_matrix[,l], class_labs, mean);
  return(mean_element)
}))

ash_theta_class <- class.normalizetpx(mean_counts_shrunk_class+1e-20, byrow=FALSE);


class_labs <- rep(1:5, each=50);
length_vec <- c(length(aorta_labels), length(coronary_labels),
                length(tibial_labels), length(blood_labels),
                length(uterus_labels))
cum_length_vec <- cumsum(length_vec)

known_samples <- as.vector(outer(1:50, c(0, cum_length_vec[1:(length(cum_length_vec)-1)]), "+"))


Topic_clus <- classtpx::class_topics(
  pooled_data, 
  K=5, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=TRUE,
  shrink.method = 1,
  mash_user = ash_theta_class,
  tol=0.01,
  ord=FALSE)

omega <- Topic_clus$omega;

annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = factor(tissue_labels[c(aorta_labels, 
                                        coronary_labels,
                                        tibial_labels,
                                        blood_labels,
                                        uterus_labels)],
                        levels=c("Artery - Aorta", 
                                 "Artery - Coronary",
                                 "Artery - Tibial",
                                 "Whole Blood",
                                 "Uterus"))
)


CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = RColorBrewer::brewer.pal(8, "Accent"),
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


test_class_classtpx <- c("Aorta", "Coronary", "Tibial",
                         "Whole Blood", "Uterus")[apply(omega,1,function(x) which.max(x))];
tab_class_classtpx <- table(c(
  rep("Aorta", length(aorta_labels)), 
  rep("Coronary", length(coronary_labels)),
  rep("Tibial", length(tibial_labels)),
  rep("Whole Blood", length(blood_labels)),
  rep("Uterus", length(uterus_labels)))[-known_samples],
  test_class_classtpx[-known_samples]);

print(tab_class_classtpx)
misclass_classtpx <- sum(tab_class_classtpx[row(tab_class_classtpx)!=col(tab_class_classtpx)])/ sum(tab_class_classtpx)
misclass_classtpx

print(tab_class_classtpx)

##


Aorta Coronary Tibial Uterus Whole Blood
Aorta         153       21      0      0           0
Coronary        5       73      5      0           0
Tibial          2       38    240      2           0
Uterus          0        0      0     33           0
Whole Blood     0        6      0      1         336



ll <- list("betaclass"=beta_class, "sebetaclass"=sebeta_class,
           "tissue_names"=c("Artery-Aorta", "Artery-Coronary", "Artery-Tibial",
                            "Whole Blood", "Uterus"))
save(ll, file="../rdas/mash_prep_example_1.rda")

ll <- get(load("../rdas/mash_prep_example_1.rda"))
