
####  Deep Learning with the mxnet package ################################

library(mxnet)
data(Sonar, package="mlbench")
Sonar[,61] = as.numeric(Sonar[,61])-1
train.ind = c(1:50, 100:150)
train.x = data.matrix(Sonar[train.ind, 1:60])
train.y = Sonar[train.ind, 61]
test.x = data.matrix(Sonar[-train.ind, 1:60])
test.y = Sonar[-train.ind, 61]


mx.set.seed(0)
model <- mx.mlp(train.x, train.y, hidden_node=50, out_node=2, out_activation="softmax",
                num.round=50, array.batch.size=15, learning.rate=0.1, momentum=0.8, 
                eval.metric=mx.metric.accuracy)

graph.viz(model$symbol$as.json())
preds = predict(model, test.x)

pred.label = max.col(t(preds))-1
table(pred.label, test.y)

training.data.frame <- data.frame(cbind(train.y, train.x));
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)

test_class_svm <- predict(model, test.x)
tab_class_svm <- table(test_class_svm, test.y)

data <- mx.symbol.Variable("data")
fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=128)
act1 <- mx.symbol.Activation(fc1, name="relu1", act_type="relu")
fc2 <- mx.symbol.FullyConnected(act1, name="fc2", num_hidden=64)
act2 <- mx.symbol.Activation(fc2, name="relu2", act_type="relu")
fc3 <- mx.symbol.FullyConnected(act2, name="fc3", num_hidden=10)
softmax <- mx.symbol.SoftmaxOutput(fc3, name="sm")


data <- mx.symbol.Variable("data")
fc1 <- mx.symbol.FullyConnected(data, name="fc1", num_hidden=128)
act1 <- mx.symbol.Activation(fc1, name="tanh1", act_type="tanh")
fc2 <- mx.symbol.FullyConnected(act1, name="fc2", num_hidden=64)
act2 <- mx.symbol.Activation(fc2, name="tanh2", act_type="tanh")
fc3 <- mx.symbol.FullyConnected(act2, name="fc3", num_hidden=32)
act3 <- mx.symbol.Activation(fc3, name="tanh3", act_type="tanh")
fc4 <- mx.symbol.FullyConnected(act3, name="fc4", num_hidden=16)
act4 <- mx.symbol.Activation(fc4, name="tanh4", act_type="tanh")
fc5 <- mx.symbol.FullyConnected(act4, name="fc5", num_hidden=8)
act5 <- mx.symbol.Activation(fc5, name="tanh5", act_type="tanh")
fc6 <- mx.symbol.FullyConnected(act5, name="fc6", num_hidden=2)
softmax <- mx.symbol.SoftmaxOutput(fc6, name="sm")

devices <- mx.cpu()

mx.set.seed(0)
model <- mx.model.FeedForward.create(softmax, X=train.x, y=train.y,
                                     ctx=devices, num.round=10, array.batch.size=10,
                                     learning.rate=0.07, momentum=0.9,  eval.metric=mx.metric.accuracy,
                                     initializer=mx.init.uniform(0.07),
                                     epoch.end.callback=mx.callback.log.train.metric(100))

graph.viz(model$symbol$as.json())




library(mxnet)
data(BostonHousing, package="mlbench")
train.ind = seq(1, 506, 3)
train.x = data.matrix(BostonHousing[train.ind, -14])
train.y = BostonHousing[train.ind, 14]
test.x = data.matrix(BostonHousing[-train.ind, -14])
test.y = BostonHousing[-train.ind, 14]

data <- mx.symbol.Variable("data")

fc1 <- mx.symbol.FullyConnected(data, num_hidden=2)

lro <- mx.symbol.LinearRegressionOutput(fc1)
mx.set.seed(0)
model <- mx.model.FeedForward.create(lro, X=train.x, y=as.array(train.y),
                                     ctx=mx.cpu(), num.round=50, array.batch.size=20,
                                     learning.rate=2e-6, momentum=0.9, eval.metric=mx.metric.rmse)



########  Deep Learning with the h2o package ############################

library(h2o)
localH2O <- h2o.init(ip = "localhost", port = 54321, startH2O = TRUE)
data(Sonar, package="mlbench")
train.ind = c(1:50, 100:150)
train.x = data.matrix(Sonar[train.ind, 1:60])
train.y = Sonar[train.ind, 61]

training_dat <- data.frame(cbind(train.y, train.x));


h2o.init()
prosPath <- system.file("extdata", "prostate.csv", package="h2o")
prostate.hex <- h2o.uploadFile(path = prosPath)
as.data.frame(prostate.hex)

h2o.dat <- as.h2o(training_dat, destination_frame = "dat")

model <- 
  h2o.deeplearning(x = 2:60,  # column numbers for predictors
                   y = 1,   # column number for label
                   training_frame = as.h2o(training_dat, destination_frame = "dat"),
                  # data = dat, # data in H2O format
                   activation = "TanhWithDropout", # or 'Tanh'
                   input_dropout_ratio = 0.2, # % of inputs dropout
                   hidden_dropout_ratios = c(0.5,0.5,0.5), # % for nodes dropout
                   balance_classes = TRUE, 
                   hidden = c(50,50,50), # three layers of 50 nodes
                   epochs = 100) # max. no. of epochs

#################################################################################


########   Deep Learning model on Arteries data ########################

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

pooled_data <- t(cbind(matdata[,aorta_labels],
                       matdata[,coronary_labels],
                       matdata[,tibial_labels]))

class_labs <- c(rep(1,50), rep(2,50), rep(3,50));
known_samples <- c(1:50, 
                   length(aorta_labels) + (1:50), 
                   length(aorta_labels)+length(coronary_labels)+ (1:50));

training_data <- pooled_data[known_samples,];
mean_features <- apply(training_data, 2, function(x) return (tapply(x, class_labs, function(y) round(mean(y)))));
mean_features <- class.normalizetpx(mean_features, byrow=TRUE)
top_features <- CountClust::ExtractTopFeatures(t(mean_features), top_features = 1000,
                                               method="poisson", options="min")

voom_generate <- class.voom_generator(pooled_data, known_samples=known_samples, class_labs = class_labs, doshrink=TRUE)


train.x <- voom_generate$voom_shrunk_class;
lib_size_1 <- rowSums(pooled_data[known_samples, ]);
lib_size <- rep(mean(lib_size_1), length(lib_size_1));
counts_shrunk_matrix <- (2^{train.x - 6*log(10, base=2)})*(rep.col(lib_size_1+1, dim(train.x)[2])) - 0.5;
counts_shrunk_matrix[counts_shrunk_matrix < 0]=1e-08;
train.x.count <- counts_shrunk_matrix

#train.x <- voom2(pooled_data[known_samples, top_features])
train.y <- class_labs;
test.x <- pooled_data[-(known_samples),];
test.y <- c(rep(1, length(aorta_labels)), 
rep(2, length(coronary_labels)),
rep(3, length(tibial_labels)))[-known_samples]

mx.set.seed(0)
model <- mx.mlp(train.x, train.y, hidden_node=100, out_node=3, out_activation="softmax",
                num.round=500, array.batch.size=15, learning.rate=0.1, momentum=0.8, 
                eval.metric=mx.metric.accuracy)

preds = predict(model, test.x)

pred.label = max.col(t(preds))
tab_class_mlp <- table(pred.label, test.y)
misclass_mlp <- sum(tab_class_mlp[row(tab_class_mlp)!=col(tab_class_mlp)])/ sum(tab_class_mlp)
misclass_mlp

training.data.frame <- data.frame(cbind(train.y, train.x));
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)

test_class_svm <- predict(model, test.x)
tab_class_svm <- table(test_class_svm, test.y)
misclass_svm <- sum(tab_class_svm[row(tab_class_svm)!=col(tab_class_svm)])/ sum(tab_class_svm)
misclass_svm


out2 <- plsgenomics::pls.lda(train.x.count, train.y, test.x, ncomp=50)
test_class_pls <- out2$predclass
tab_class_pls <- table(test_class_pls, test.y)
misclass_pls <- sum(tab_class_pls[row(tab_class_pls)!=col(tab_class_pls)])/ sum(tab_class_pls)
misclass_pls

##########   Deep Learning on full gtex data  ##################################


library(data.table)
library(e1071)
library(classtpx)
library(maptpx)

gtex_data <- data.frame(fread("../data/GTEX_V6/cis_gene_expression.txt"));
matdata <- gtex_data[,-(1:2)];

tissue_labels=read.table("../data/GTEX_V6/samples_id.txt")[,3];

tab_tissue_labels <- table(tissue_labels)
inds <- which(tab_tissue_labels > 100)
#inds <- tab_tissue_labels[grep("Skin", names(tab_tissue_labels))]
#inds <- c("Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Breast - Mammary Tissue")
gtex_data_filt <- matdata[,which(!is.na(match(tissue_labels, names(inds))))]
tissue_labels_filt <- tissue_labels[which(!is.na(match(tissue_labels, names(inds))))]

out_pmd <- PMD(t(gtex_data_filt), K = 40, upos = TRUE,
                   vpos = TRUE, center = TRUE, sumabs = 1, niter = 2000,
                   sumabsu = sqrt(600))
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
#top_features <- CountClust::ExtractTopFeatures(t(mean_features), top_features = 1000,
#                                              method="poisson", options="min")
#top_features <- as.vector(top_features)
#top_features <- top_features[!is.na(top_features)]

top_features <- 1:16069

train.x <- training_data[,top_features]
train.y <- class_labs;
test.x <- t(gtex_data_filt[top_features, -(known_samples)])
test.y <- class_labs_full[-(known_samples)];

mx.set.seed(0)
model <- mx.mlp(train.x, train.y, hidden_node=500, out_node=36, out_activation="softmax",
                num.round=500, array.batch.size=15, learning.rate=0.1, momentum=0.8, 
                eval.metric=mx.metric.accuracy)

graph.viz(model$symbol$as.json())
#save(model, file="../rdas/mlp_model_gtex.rda")
preds = predict(model, test.x)

pred.label = max.col(t(preds))-1
table(pred.label, test.y)

tab_class_mlp <- table(pred.label, test.y)
misclass_mlp <- sum(tab_class_mlp[row(tab_class_mlp)!=col(tab_class_mlp)])/ sum(tab_class_mlp)
misclass_mlp
############################################################################


#######  PLS-LDA on GTEx V6 counts #########

out2 <- plsgenomics::pls.lda(train.x, train.y, test.x, ncomp=50)
test_class_pls <- out2$predclass
tab_class_pls <- table(test_class_pls, test.y)
misclass_pls <- sum(tab_class_pls[row(tab_class_pls)!=col(tab_class_pls)])/ sum(tab_class_pls)
misclass_pls

#######  PLS-LDA on GTEx V6 voom #########

out2 <- plsgenomics::pls.lda(voom2(train.x), train.y, voom2(test.x), ncomp=50)
test_class_pls <- out2$predclass
tab_class_pls <- table(test_class_pls, test.y)
misclass_pls <- sum(tab_class_pls[row(tab_class_pls)!=col(tab_class_pls)])/ sum(tab_class_pls)
misclass_pls

#######  PLS-LDA on GTEx V6 voom shrink #########

voom_generate <- class.voom_generator(train.x, class_labs = class_labs, doshrink=TRUE)

out2 <- plsgenomics::pls.lda(voom_generate$voom_shrunk_class, train.y, voom2(test.x), ncomp=50)
test_class_pls <- out2$predclass
tab_class_pls <- table(test_class_pls, test.y)
misclass_pls <- sum(tab_class_pls[row(tab_class_pls)!=col(tab_class_pls)])/ sum(tab_class_pls)
misclass_pls

#############   SVM on GTEx V6 counts data ###########################

training.data.frame <- data.frame(cbind(train.y, train.x));
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)
colnames(test.x) <- colnames(train.x)
test_class_svm <- predict(model, test.x)
tab_class_svm <- table(test_class_svm, test.y)
misclass_svm <- sum(tab_class_svm[row(tab_class_svm)!=col(tab_class_svm)])/ sum(tab_class_svm)
misclass_svm

#############  SVM on GTEx V6 voom data ##############################

voom.train.x <- voom2(train.x)
voom.test.x <- voom2(test.x)
training.data.frame <- data.frame(train.y, voom.train.x);
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)
colnames(voom.test.x) <- colnames(voom.train.x)
test_class_svm <- predict(model, voom.test.x)
tab_class_svm <- table(test_class_svm, test.y)
misclass_svm <- sum(tab_class_svm[row(tab_class_svm)!=col(tab_class_svm)])/ sum(tab_class_svm)
misclass_svm

############  SVM on GTEx V6 voom shrunk data ##########################

voom.train.x <- voom_generate$voom_shrunk_class;
voom.test.x <- voom2(test.x)
training.data.frame <- data.frame(train.y, voom.train.x);
model <- svm(as.factor(train.y) ~ ., data = training.data.frame)
colnames(voom.test.x) <- colnames(voom.train.x)
test_class_svm <- predict(model, voom.test.x)
tab_class_svm <- table(test_class_svm, test.y)
misclass_svm <- sum(tab_class_svm[row(tab_class_svm)!=col(tab_class_svm)])/ sum(tab_class_svm)
misclass_svm


#############   Model based framework #################################

#################  Normal Loglik model  #############################

counts <- t(gtex_data_filt[top_features,])
model1 <- class.model_clust(counts, known_samples = known_samples,
                            class_labs = class_labs, dist="normal")

test_class_normal <- apply(model1, 1, which.max);
tab_class_normal <- table(test_class_normal, test.y)
misclass_normal <- sum(tab_class_normal[row(tab_class_normal)!=col(tab_class_normal)])/ sum(tab_class_normal)
misclass_normal

#################   Poisson Loglik model  ############################

model2 <- class.model_clust(counts, known_samples = known_samples,
                            class_labs = class_labs, dist="poisson")

test_class_poisson <- apply(model2, 1, which.max);
tab_class_poisson <- table(test_class_poisson, test.y)
misclass_poisson <- sum(tab_class_poisson[row(tab_class_poisson)!=col(tab_class_poisson)])/ sum(tab_class_poisson)
misclass_poisson

#################  Neg Binomial model    #####################

model3 <- class.model_clust(counts, known_samples = known_samples,
                            class_labs = class_labs, dist="negbinom")

test_class_negbinom <- apply(model3, 1, which.max);
tab_class_negbinom <- table(test_class_negbinom, test.y)
misclass_negbinom <- sum(tab_class_negbinom[row(tab_class_negbinom)!=col(tab_class_negbinom)])/ sum(tab_class_negbinom)
misclass_negbinom

##################   classtpx (theta.fix)  ###############################


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
  K=2, 
  known_samples = known_samples,
  class_labs = class_labs,
  method="omega.fix",
  shrink=TRUE,
  shrink.method = 1,
  tol=0.001,
  ord=FALSE)


par(mfrow=c(1,1))
par(mar=c(4,2,2,2))
omega <- Topic_clus$omega;
annotation <- data.frame(
  sample_id = paste0("X", c(1:NROW(omega))),
  tissue_label = tissue_labels_filt
)

cols1 <- c(RColorBrewer::brewer.pal(12, "Set3"),
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

rownames(tab_class_classtpx) <- unique(tissue_labels_filt)
colnames(tab_class_classtpx) <-  rep("", 36)
myImagePlot(class.normalizetpx(tab_class_classtpx))