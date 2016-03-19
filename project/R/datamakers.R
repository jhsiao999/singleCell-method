# datamaker: a function to generate a G*(2N) count matrix, with G genes and 2N samples for 2 groups.
# First N samples for group/condition A, and the rest N samples for group B.

# Input: args is a list of arguments:
# - path: (Required!) Data path (I put the GTEx data files under the directory path/gtex/).
# - tissue: (Required!) String of one tissue name, or a vector of 2 tissue names. 
#           If inputting one tissue name, then all 2N samples are selected from this tissue's sample.
#           If inputting two tissue names, then first N samples are selected from tissue1 
#           and last N samples from tissue2 or from a mixture of tissue1 and tissue2 
#           (when breaksample==TRUE & 0<nullpi<1).
# - Nsamp: N, Number of samples in each group (so the total number of samples is 2*N). 
#          If Nsamp==NULL, select all samples of the tissue(s).
# - Ngene: Number of top high-expressed genes that will be selected. If Ngene==NULL, select all genes.
# - breaksample: Flag, whether select each gene's counts from different GTEx samples. 
#                This will break the possible within sample correlation. The default value is FALSE.
# - poisthin: Flag, whether use Poisson thinning to change some counts. 
#             Need specify 2 parameters for thinning: log2foldmean, log2foldsd,
#             and the proportion of genes for thinning: nullpi. The default value is FALSE.
# - log2foldmean: Only works when poisthin==TRUE. Mean parameter for log2-fold change.
# - log2foldsd: Only works when poisthin==TRUE. SD parameter for log2-fold change.
# - nullpi: Proportion of genes that are nulls. 
#           Only works when tissue is one tissue name and poisthin==TRUE, 
#           or tissue contains two tissue names and breaksample==TRUE.
# - pseudocounts: Add pseudocounts to the count matrix. The default value is 1.
# - RUV.k: Number of surrogate variables for RUV. The default value is round(log2(Nsamp)).

# Output: a list of input info and meta info for dscr.
# - input: a list of G*2N count matrix (count) and a 2N vector of group assignment (condition), 
#          and estimated effects (betahat.XXX), sd (sebetahat.XXX) and degree of freedom (df.XXX)
#          from several different methods.
# - meta: the true null/alternative
#         info for each gene (null), and the input arguments of datamaker function (after setting default).

# Example: 
# Compare 50 Lung samples against 50 Ovary samples (only select top 10000 highly-expressed genes):
# args=list(tissue=c("Lung","Ovary"), Nsamp=50, Ngene=10000, path="/mnt/lustre/home/mengyin")
# More examples in https://github.com/mengyin/dscr-gtex-total/blob/master/scenarios.R.

library(edgeR)
library(limma)
library(RUVSeq)
library(sva)
library(DESeq)
library(data.table)

datamaker = function(args){  
    dfargs = default_datamaker_args(args)
    
    # rawdata1 = readtissue(dfargs$path, dfargs$tissue[1])
    # rawdata1 = read.table(paste0(dfargs$path,"/gtex/tissues/",dfargs$tissue[1],".txt"),header=TRUE)
#    rawdata1 = as.matrix(fread(paste0(dfargs$path,"/gtex/tissues/",dfargs$tissue[1],".txt"),header=FALSE,drop=1,skip=1))
#    colnames(rawdata1) = as.character(read.table(paste0(dfargs$path,"/gtex/tissues/",dfargs$tissue[1],".txt"),
#                                                 nrows=1,stringsAsFactors=FALSE))
    rawdata1 <- dfargs$rawdata1

#    if (length(dfargs$tissue)>1){
     if (!is.null(dfargs$rawdata2)) {
        # rawdata2 = readtissue(dfargs$path, dfargs$tissue[2])
        # rawdata2 = read.table(paste0(dfargs$path,"/gtex/tissues/",dfargs$tissue[2],".txt"),header=TRUE)
        # rawdata2 = as.matrix(fread(paste0(dfargs$path,"/gtex/tissues/",dfargs$tissue[2],".txt"),header=FALSE,drop=1,skip=1))
        # colnames(rawdata1) = as.character(read.table(paste0(dfargs$path,"/gtex/tissues/",dfargs$tissue[2],".txt"),
        #                                              nrows=1,stringsAsFactors=FALSE))
        rawdata2 <- args$rawdata2
        
        if (is.null(dfargs$Nsamp)){
            dfargs$Nsamp = min(dim(rawdata1)[2],dim(rawdata2)[2])
        }
        if (dim(rawdata1)[2]<dfargs$Nsamp | dim(rawdata2)[2]<dfargs$Nsamp){
            stop("Not enough samples in the raw dataset!")
        }
        
        
        if (dfargs$nullpi==0){
            # All genes are alternatives
            counts1 = selectsample(rawdata1, dfargs$Nsamp, dfargs$breaksample)     
            counts2 = selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            
            counts = cbind(counts1, counts2)
        }else{
            # Some genes are nulls, some are alternatives
            counts1 = selectsample(rawdata1, 2*dfargs$Nsamp, dfargs$breaksample)      
            counts2 = selectsample(rawdata2, dfargs$Nsamp, dfargs$breaksample)
            counts = cbind(counts1, counts2)
        }
    }else{
        if (is.null(dfargs$Nsamp)){
            dfargs$Nsamp = floor(dim(rawdata1)[2]/2)
        }
        if (dim(rawdata1)[2] < 2*dfargs$Nsamp){
            stop("Not enough samples in the raw dataset!")
        }
        
        counts = selectsample(rawdata1, 2*dfargs$Nsamp, dfargs$breaksample)
    } 
    
    # Remove genes without any reads
    counts = counts[apply(counts,1,sum)>0,]
    
    # Take the top Ngene high-expressed genes
    if (!is.null(dfargs$Ngene)){
        dfargs$Ngene = min(dfargs$Ngene, dim(counts)[1])
        counts = counts[sort(order(rowSums(counts),decreasing=TRUE)[1:dfargs$Ngene]),]
    }
    dfargs$Ngene = dim(counts)[1]
    
    # Model's design: Nsamp samples for group A and Nsamp samples for group B
    condition = factor(rep(1:2,each=dfargs$Nsamp))
    design = model.matrix(~condition)
    
    # Ground truth of null hypotheses: beta_g=0
    null = rep(0,dfargs$Ngene)
    null[sample(dfargs$Ngene, round(dfargs$Ngene*dfargs$nullpi))] = 1  
    
    # Poisson thinning (optional)
    counts = pois_thinning(counts, dfargs, null)
    
    # Mix null and alternative genes from different samples (optional)
    counts = mix_sample(counts, dfargs, null)
    
    # Voom transformation
    voom = voom_transform(counts, condition)
    
    # Quasi-binomial glm
    qb = quasi_binom(counts, condition)
    
    # Myrna (Langmead et al. '10) & Quasi-binomial glm 
    # Use log(75th quantile of samples' counts) as covariate
    W.Myrna = apply(counts,2,function(x) log(quantile(x[x>0],0.75)))
    Myrnaqb = quasi_binom(counts, condition, W=W.Myrna)
    # Use log(75th quantile of samples' counts) as offset
    offset.Myrnaoff = apply(counts,2,function(x) quantile(x[x>0],0.75))
    Myrnaoffqb = quasi_binom(counts, condition, W=NULL, offset=offset.Myrnaoff)
    
    # RUV & voom
    halfnull = rep(0,length(null))  # Use half of the true nulls to do supervised RUV/SVA
    halfnull[which(null==1)[1:floor(length(which(null==1))/2)]] = 1
    W.RUV = RUV_factor(counts, dfargs, halfnull)
    RUVvoom = voom_transform(counts, condition, W=W.RUV)
    
    # SVA & voom
    W.SVA = SVA_factor(counts, condition, dfargs, halfnull)
    SVAvoom = voom_transform(counts, condition, W=W.SVA)
    
    # RUV & quasi-binomial glm
    RUVqb = quasi_binom(counts, condition, W=W.RUV)
    
    # SVA & quasi-binomial glm
    SVAqb = quasi_binom(counts, condition, W=W.SVA)
    
    # Get sebetahat from edgeR.glm (infer from betahat & pval)
    edgeRglm = edgeR_glmest(counts, condition, dfargs)
    
    # Get sebetahat from DESeq2 (infer from betahat & pval)
    DESeq2glm = DESeq2_glmest(counts, condition, dfargs)
    
    # meta data
    meta = list(null=null, dfargs=dfargs)
    
    # input data
    input = list(counts=counts, condition=condition,
                 v=voom$v, betahat.voom=voom$betahat, sebetahat.voom=voom$sebetahat, df.voom=voom$df,
                 betahat.RUVvoom=RUVvoom$betahat, sebetahat.RUVvoom=RUVvoom$sebetahat, df.RUVvoom=RUVvoom$df, W.RUV=W.RUV,
                 betahat.SVAvoom=SVAvoom$betahat, sebetahat.SVAvoom=SVAvoom$sebetahat, df.SVAvoom=SVAvoom$df, W.SVA=W.SVA,
                 betahat.qb=qb$betahat, sebetahat.qb=qb$sebetahat, df.qb=qb$df, dispersion.qb=qb$dispersion,
                 betahat.RUVqb=RUVqb$betahat, sebetahat.RUVqb=RUVqb$sebetahat, dispersion.RUVqb=RUVqb$dispersion, df.RUVqb=RUVqb$df, W.RUV=W.RUV,
                 betahat.SVAqb=SVAqb$betahat, sebetahat.SVAqb=SVAqb$sebetahat, dispersion.SVAqb=SVAqb$dispersion, df.SVAqb=SVAqb$df, W.SVA=W.SVA,
                 betahat.Myrnaqb=Myrnaqb$betahat, sebetahat.Myrnaqb=Myrnaqb$sebetahat, dispersion.Myrnaqb=Myrnaqb$dispersion, df.Myrnaqb=Myrnaqb$df, W.Myrna=W.Myrna,
                 betahat.Myrnaoffqb=Myrnaoffqb$betahat, sebetahat.Myrnaoffqb=Myrnaoffqb$sebetahat, dispersion.Myrnaoffqb=Myrnaoffqb$dispersion, df.Myrnaoffqb=Myrnaoffqb$df, offset.Myrnaoff=offset.Myrnaoff,
                 betahat.edgeRglm=edgeRglm$betahat, sebetahat.edgeRglm=edgeRglm$sebetahat, df.edgeRglm=edgeRglm$df,
                 betahat.DESeq2glm=DESeq2glm$betahat, sebetahat.DESeq2glm=DESeq2glm$sebetahat, df.DESeq2glm=DESeq2glm$df)
    
    data = list(meta=meta,input=input)
    return(data)
}

###<--- Set default arguments for datamaker function
default_datamaker_args = function(args){
    # poisthin: flag of Poisson thinning
    if (is.null(args$poisthin)){
        args$poisthin = FALSE
    }
    
    # log2foldmean, log2foldsd: Poisson thinning params
    if (args$poisthin==TRUE){
        if (is.null(args$log2foldmean)){
            args$log2foldmean = 0
        }
        if (is.null(args$log2foldsd)){
            args$log2foldsd = 1
        }
    }
    
    # breaksample: flag of each gene randomly select samples
    if (is.null(args$breaksample)){
        args$breaksample = FALSE
    }
    
    
    # nullpi: proportion of null genes
    # default if thinning: 90% null genes
    # default if not thinning and one tissue: 100% null genes
    # default if not thinning and two tissues: 0% null genes
    if (is.null(args$nullpi)){
        if (args$poisthin==TRUE){
            args$nullpi = 0.9
        }else if (length(args$tissue)==1){
            args$nullpi = 1
        }else if (length(args$tissue)>1){
            args$nullpi = 0
        }
    }
    
    # RUV.k: number of surrogate variables for RUV.
    #   if (is.null(args$RUV.k)){
    #     args$RUV.k = round(log2(args$Nsamp))
    #   }
    
    # pseudocounts: add pseudocounts to count matrix
    if (is.null(args$pseudocounts)){
        args$pseudocounts = 1
    }
    
    return(args)
}


###<--- Poisson thinning
pois_thinning = function(counts, args, null){
    if (args$poisthin==TRUE){ 
        log2foldchanges = rnorm(sum(!null), mean=args$log2foldmean, sd=args$log2foldsd)
        foldchanges = 2^log2foldchanges
        
        # thin group A
        counts[which(!null)[log2foldchanges>0],1:args$Nsamp]=matrix(rbinom(sum(log2foldchanges>0)*args$Nsamp, 
                                                                           size=c(as.matrix(counts[which(!null)[log2foldchanges>0],1:args$Nsamp])),
                                                                           prob=rep(1/foldchanges[log2foldchanges>0],args$Nsamp)),ncol=args$Nsamp)
        # thin group B
        counts[which(!null)[log2foldchanges<0],(args$Nsamp+1):(2*args$Nsamp)]=matrix(rbinom(sum(log2foldchanges<0)*args$Nsamp, 
                                                                                            size=c(as.matrix(counts[which(!null)[log2foldchanges<0],
                                                                                                                    (args$Nsamp+1):(2*args$Nsamp)])),
                                                                                            prob=rep(foldchanges[log2foldchanges<0],args$Nsamp)),
                                                                                     ncol=args$Nsamp)
        
    }
    return(counts)
}


###<--- Mix null and alternative genes from different samples
mix_sample = function(counts, args, null){
    if(args$nullpi<1 & args$nullpi>0 & args$breaksample==TRUE){
        newcounts = matrix(rep(0, args$Ngene*2*args$Nsamp),nrow=args$Ngene)
        newcounts[as.logical(null),] = counts[as.logical(null),1:(2*args$Nsamp)]
        newcounts[!null,] = counts[!null,c(1:args$Nsamp,(2*args$Nsamp+1):(3*args$Nsamp))]
        counts = newcounts
        rm(newcounts); 
    }
    return(counts=counts)
}


###<--- Voom transformation
voom_transform = function(counts, condition, W=NULL){
    dgecounts = calcNormFactors(DGEList(counts=counts,group=condition))
    
    if (is.null(W)){
        design = model.matrix(~condition)
    }else{
        design = model.matrix(~condition+W)
    }
    
    v = voom(dgecounts,design,plot=FALSE)
    lim = lmFit(v)
    #zdat.voom = apply(cbind(v$E,v$weights),1,wls.wrapper,g=condition)
    #betahat.voom = zdat.voom[1,]
    #sebetahat.voom = zdat.voom[2,]
    betahat.voom = lim$coefficients[,2]
    sebetahat.voom = lim$stdev.unscaled[,2]*lim$sigma
    df.voom = length(condition)-2-!is.null(W)
    
    return(list(betahat=betahat.voom, sebetahat=sebetahat.voom, df=df.voom, v=v))
}


###<--- Weighted least squares regression
# g: formula
# ynweights: matrix of response y and corresponding weights
wls.wrapper = function(ynweights,g,...){
    y = ynweights[1:(length(ynweights)/2)]
    weights = ynweights[(length(ynweights)/2+1):length(ynweights)]
    y.wls = lm(y~g,weights=weights,...)
    
    # slope estimate & standard error
    c = as.vector(t(summary(y.wls)$coeff[2,1:2]))
    return(c)
}


###<--- Quasi-binomial glm
quasi_binom = function(counts, condition, W=NULL, offset=NULL){
    zdat.qb = counts.associate(counts, condition, W=W, offset=offset)
    betahat = zdat.qb[3,]
    sebetahat = zdat.qb[4,]
    dispersion = zdat.qb[5,]
    df = length(condition)-2-!is.null(W)
    return(list(betahat = betahat, sebetahat = sebetahat,
                df = df, dispersion = dispersion))
}

# counts is a ngene (or nwindow) by nsample matrix of counts (eg RNAseq)
# g is an nsample vector of group memberships/covariate
# looks for association between rows and covariate
counts.associate=function(counts, g, W=NULL, offset=NULL, pseudocounts=1){
    y.counts=t(as.matrix(counts)) 
    col.sum = apply(y.counts, 2, sum)
    y.counts=y.counts[,col.sum>0] #remove 0 columns
    y.counts = y.counts+pseudocounts 
    if (is.null(offset)){
        offset = apply(y.counts,1,sum)
    }
    
    y.prop = y.counts/apply(y.counts,1,sum) # compute proportions
    zdat = apply(y.prop,2,glm.binomial.wrapper,g=g,W=W,weights=offset,epsilon=1e-6)  #estimate effect sizes and standard errors
    #zdat.ash = ash(zdat[3,],zdat[4,],df=2,method='fdr') #shrink the estimated effects
    #return(list(zdat=zdat,zdat.ash=zdat.ash))
    return(zdat)
}

glm.binomial.wrapper = function(y,g,W=NULL,...){
    if (is.null(W)){
        y.glm=safe.quasibinomial.glm(y~g,...)
    }else{
        y.glm=safe.quasibinomial.glm(y~g+W,...)
    }
    return(c(get.coeff(y.glm),summary(y.glm)$dispersion))
}


#fill NAs with 0s (or other value)
fill.nas=function(x,t=0){
    x[is.na(x)]=t
    return(x)
}

#get estimates and standard errors from a glm object
#return NAs if not converged
get.coeff=function(x.glm){
    c=as.vector(t(summary(x.glm)$coeff[,1:2]))  
    if(x.glm$conv){return(c)} else {return(rep(NA,length(c)))}
}

# use glm to fit quasibinomial, but don't allow for underdispersion!
safe.quasibinomial.glm=function(formula,forcebin=FALSE,...){
    if(forcebin){
        fm=glm(formula,family=binomial,...)
    } else{
        fm = glm(formula,family=quasibinomial,...)
        if(is.na(summary(fm)$dispersion) | summary(fm)$dispersion<1){
            fm=glm(formula,family=binomial,...)
        }
    }
    return(fm)
}

# randomly subsample data for each gene
# gene: a vector of reads for one gene
# Nsamp: # of samples wanted
sampleingene = function(gene, Nsamp){
    sample = sample(length(gene),Nsamp)
    return(c(gene[sample]e))
}

# Randomly select samples
# counts: full count matrix
# Nsamp: # of samples wanted
# breaksample: flag, if select different samples for each gene
selectsample = function(counts, Nsamp, breaksample){
    if (breaksample==FALSE){
        subsample = sample(1:dim(counts)[2],Nsamp)
        counts = counts[,subsample]
    }else{
        counts = t(apply(counts, 1, sampleingene, Nsamp=Nsamp))
    }
    return(counts)
}

# Use RUV to estimate confounding factor
RUV_factor = function(counts, args, null){
    W = NULL
    seq = newSeqExpressionSet(as.matrix(counts[as.logical(null),]))
    if (sum(null)>0){
        controls = rownames(seq)
        # differences = matrix(data=c(1:args$Nsamp, (args$Nsamp+1):(2*args$Nsamp)), byrow=TRUE, nrow=2)
        if (is.null(args$RUV.k)){
            k = num.sv(counts[as.logical(null),],rep(1,dim(counts)[2]))
        }else{
            k = args$RUV.k
        }
        
        if (k>0){
            #seqRUV = RUVs(seq, controls, k=k, differences)
            seqRUV = RUVg(seq, controls, k=k)
            W = as.matrix(pData(seqRUV))
        }
    }
    return(W = W) 
}

# Use SVA to estimate confounding factor
SVA_factor = function(counts, condition, args, null){
    mod1 = model.matrix(~condition)
    mod0 = cbind(mod1[,1])
    
    if (args$nullpi>0){
        svseq = svaseq(counts,mod1,mod0,control=null)
    }else{
        svseq = svaseq(counts,mod1,mod0)
    }
    
    if(svseq$n.sv>0){
        return(W = svseq$sv)
    }else{
        return(W = NULL)
    }
}

# Get sebetahat from edgeR.glm (infer from betahat & pval)
edgeR_glmest = function(counts, condition, args){
    design = model.matrix(~condition)
    y = DGEList(counts=counts+args$pseudocounts, group=condition)
    y = calcNormFactors(y)
    y = estimateGLMCommonDisp(y,design)
    y = estimateGLMTrendedDisp(y,design)
    y = estimateGLMTagwiseDisp(y,design)
    fit = glmFit(y,design)
    lrt = glmLRT(fit,coef=2)
    betahat.edgeRglm = fit$coef[,2]
    df.edgeRglm = length(condition)-2
    tscore = qt(1-lrt$table$PValue/2,df=df.edgeRglm)
    sebetahat.edgeRglm = abs(betahat.edgeRglm/tscore)
    
    return(list(betahat=betahat.edgeRglm, sebetahat=sebetahat.edgeRglm,
                df=df.edgeRglm)) 
}

# Get sebetahat from DESeq.glm (infer from betahat & pval)
DESeq_glmest = function(counts, condition, args){
    cds = newCountDataSet(counts+args$pseudocounts, condition )
    cds = estimateSizeFactors( cds )
    cds = try(estimateDispersions( cds ),silent=TRUE)
    if (class(cds)=="try-error"){
        betahat.DESeqglm = NA
        sebetahat.DESeqglm = NA
        df.DESeqglm = length(condition)-2
    }else{
        fit1 = fitNbinomGLMs( cds, count ~ condition )     
        fit0 = fitNbinomGLMs( cds, count ~ 1 )
        betahat.DESeqglm = fit1[,2]
        df.DESeqglm = length(condition)-2
        tscore = qt(1-nbinomGLMTest(fit1,fit0)/2,df=df.DESeqglm)
        sebetahat.DESeqglm = abs(betahat.DESeqglm/tscore)
    }
    return(list(betahat=betahat.DESeqglm, sebetahat=sebetahat.DESeqglm,
                df=df.DESeqglm)) 
}

# Get sebetahat from DESeq2 (infer from betahat & pval)
DESeq2_glmest = function(counts, condition, args){
    cond = condition
    dds = DESeqDataSetFromMatrix(counts+args$pseudocounts, DataFrame(cond), ~cond)
    dds = estimateSizeFactors(dds)
    dds = estimateDispersions(dds)
    dds = nbinomWaldTest(dds)
    res = results(dds,cooksCutoff=FALSE)
    pvalue = res$pvalue
    betahat.DESeq2 = res$log2FoldChange
    df.DESeq2 = length(condition)-2
    sebetahat.DESeq2 = res$lfcSE
    
    return(list(betahat=betahat.DESeq2, sebetahat=sebetahat.DESeq2,
                df=df.DESeq2)) 
}


# Extract dataset for a specific tissue from the GTEx reads txt file
# Note: use GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct_new.txt,
# which removes the first 3 lines from GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.txt
readtissue = function(path, tissue){
    tis = read.table(paste0(path,"/gtex/sample_tissue.txt"), header=TRUE)
    
    tissue.idx = grep(tissue,tis[,2],fixed=TRUE)
    cols = rep("NULL",dim(tis)[1])
    cols[tissue.idx]="numeric"
    
    data = read.table(paste0(path,"/gtex/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct_new.txt"),
                      colClasses = cols, header = TRUE)
    
    return(data)
}