
###  classtpx gtex 


library(data.table)
library(e1071)
library(maptpx)
library(classtpx)
library(limma)

source("../../../classtpx/R/class_count.R")
source("../../../classtpx/R/class_topics.R")
source("../../../classtpx/R/class_tpx.R")
source("../../../classtpx/R/class_varselect.R")


gtex_data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"));
matdata <- gtex_data[,-(1:2)];

tissue_labels=read.table("../data/GTEX_V6/samples_id.txt")[,3];

tab_tissue_labels <- table(tissue_labels)
inds <- which(tab_tissue_labels > 100)
gtex_data_filt <- matdata[,which(!is.na(match(tissue_labels, names(inds))))]
tissue_labels_filt <- tissue_labels[which(!is.na(match(tissue_labels, names(inds))))]

tissue_labels_filt <- droplevels(tissue_labels_filt)
class_labs <- rep(1:length(levels(tissue_labels_filt)), each=50);
length_vec <- as.numeric(table(tissue_labels_filt))
tissue_labs <- names(table(tissue_labels_filt))
cum_length_vec <- cumsum(length_vec)
class_labs_full <- rep(1:length(levels(tissue_labels_filt)), length_vec)


known_samples <- as.vector(outer(1:50, 
                                 c(0, cum_length_vec[1:(length(cum_length_vec)-1)]), "+"))

training_data <- t(gtex_data_filt[,known_samples]);

mean_features <- apply(training_data, 2, function(x) return (tapply(x, as.factor(class_labs), function(y) return(round(mean(y))))));
mean_features <- class.normalizetpx(mean_features, byrow=TRUE)
top_features <- CountClust::ExtractTopFeatures(t(mean_features), top_features = 200,
                                               method="poisson", options="min")
top_features <- as.vector(top_features);
top_features <- top_features[!is.na(top_features)]

counts <- t(gtex_data_filt[top_features,]);

train.x <- training_data[,top_features]
train.y <- class_labs;
test.x <- t(gtex_data_filt[top_features, -(known_samples)])
test.y <- class_labs_full[-(known_samples)];


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
  method="omega.fix",
  shrink=FALSE,
  shrink.method = 1,
  tol=1,
  ord=FALSE)


#save(Topic_clus, file="../rdas/gtex_all_tissues_classtpx_omega.fix.rda")

Topic_clus <- get(load("../rdas/gtex_all_tissues_classtpx_theta.fix_shrink.rda"))
