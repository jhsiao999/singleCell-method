
###  classtpx gtex 


source("../../../classtpx/R/class_count.R")
source("../../../classtpx/R/class_topics.R")
source("../../../classtpx/R/class_tpx.R")
source("../../../classtpx/R/class_varselect.R")


library(data.table)
library(e1071)
library(maptpx)
library(limma)

gtex_data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"));
matdata <- gtex_data[,-(1:2)];

tissue_labels=read.table("../data/GTEX_V6/samples_id.txt")[,3];

gtex20_clus <- get(load("../rdas/gtex_all_tissues_classtpx_theta.fix_noshrink.rda"));
omega20 <- gtex20_clus$omega

tab_tissue_labels <- table(tissue_labels)
inds <- which(tab_tissue_labels > 100)
# gtex_data_filt <- matdata[,which(!is.na(match(tissue_labels, names(inds))))]
gtex_data_filt <- omega20[which(!is.na(match(tissue_labels, names(inds)))),]
tissue_labels_filt <- tissue_labels[which(!is.na(match(tissue_labels, names(inds))))]

tissue_labels_filt <- droplevels(tissue_labels_filt)
class_labs <- rep(1:length(levels(tissue_labels_filt)), each=50);
length_vec <- as.numeric(table(tissue_labels_filt))
tissue_labs <- names(table(tissue_labels_filt))
cum_length_vec <- cumsum(length_vec)
class_labs_full <- rep(1:length(levels(tissue_labels_filt)), length_vec)


known_samples <- as.vector(outer(1:50, 
                                 c(0, cum_length_vec[1:(length(cum_length_vec)-1)]), "+"))

#training_data <- t(gtex_data_filt[,known_samples]);

training_data <- gtex_data_filt[known_samples,];


mean_features <- apply(training_data, 2, function(x) return (tapply(x, as.factor(class_labs), function(y) return(round(mean(y))))));
mean_features <- class.normalizetpx(mean_features, byrow=TRUE)
top_features <- CountClust::ExtractTopFeatures(t(mean_features), top_features = 200,
                                               method="poisson", options="min")
top_features <- as.vector(top_features);
top_features <- top_features[!is.na(top_features)]

counts <- t(gtex_data_filt[top_features,]);

#counts <- gtex_data_filt;

#### SVM comparison ###

training_data <- out_pmd$u[known_samples,];
train.x <- training_data;
train.y <- class_labs;
test.x <- out_pmd$u[ -(known_samples), ]
test.y <- class_labs_full[-(known_samples)];




training.data.frame <- data.frame(cbind(train.y, train.x));
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)
colnames(test.x) <- colnames(training.data.frame)[-1]
test_class_svm <- predict(model, test.x)
tab_class_svm <- table(test_class_svm, test.y)
misclass_svm <- sum(tab_class_svm[row(tab_class_svm)!=col(tab_class_svm)])/ sum(tab_class_svm)
misclass_svm



shrink=TRUE
shrink.method=1
mash_user=NULL
shape=NULL 
initopics=NULL 
tol=0.1 
bf=FALSE 
kill=2 
ord=TRUE 
verb=1
tmax=10000 
wtol=10^(-4) 
qn=100
grp=NULL 
admix=TRUE 
nonzero=FALSE 
dcut=-10


Topic_clus <- class_topics(
  counts, 
  K=36, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="theta.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)


#save(Topic_clus, file="../rdas/gtex_all_tissues_classtpx_theta.fix_shrink.rda")

Topic_clus <- get(load("../rdas/gtex_all_tissues_classtpx_theta.fix_shrink.rda"))

par(mfrow=c(1,1))
par(mar=c(4,2,2,2))
omega <- Topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = tissue_labels_filt
)

cols1 <- c(rev(RColorBrewer::brewer.pal(12, "Set3")),
           RColorBrewer::brewer.pal(8, "Accent"),
           RColorBrewer::brewer.pal(9, "Set1"),
           RColorBrewer::brewer.pal(8, "Dark2"))


CountClust::StructureGGplot(omega = as.matrix(omega),
                            annotation = annotation,
                            palette = cols1,
                            yaxis_label = "Tissue type",
                            order_sample = TRUE,
                            axis_tick = list(axis_ticks_length = .1,
                                             axis_ticks_lwd_y = .1,
                                             axis_ticks_lwd_x = .1,
                                             axis_label_size = 7,
                                             axis_label_face = "bold"))


test_class_classtpx <- apply(omega,1,function(x) which.max(x));
tab_class_classtpx <- table(test_class_classtpx[-known_samples], test.y)
misclass_classtpx <- sum(tab_class_classtpx[row(tab_class_classtpx)!=col(tab_class_classtpx)])/ sum(tab_class_classtpx)
misclass_classtpx
