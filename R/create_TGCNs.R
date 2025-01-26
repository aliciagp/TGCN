
#' getSamplingIndex - It splits samples into different folds for bootstrapped LASSO
#'
#' @param seed the seed number to ensure the reproducibility of the results
#' @param exprData a data frame representing a gene expression matrix with genes as rows and samples as columns
#' @param train.split the proportion of samples included in this fold
#' @param t the number of folds to create
#' @return a list with as many items as t and each one containing the samples index selected
#' @export
#' @examples

getSamplingIndex <- function(seed=1234,
                             exprData,
                             train.split=0.7,
                             t=10,
                             size=NULL) {

  set.seed(seed)
  samples <- list()

  for (i in 1:t) {
    samples[[i]] <- sample(1:ncol(exprData), train.split*ncol(exprData))
  }


  return(samples)
}


#' getMetrics - It returns the performance of the model created for both train and test sets
#'
#' @param model a model built with cv.glmnet function from glmnet package or train function from caret R package
#' @param data.train a data frame containing the data used to train the model with samples as rows and genes as columns
#' @param data.test a data frame containing the data reserved to test the model with samples as rows and genes as columns
#' @param target.train a vector with the response variable (target) used to train the model
#' @param target.test a vector with the response variable (target) to be predicted to test the model
#' @return a list containing 1) the metrics of the performance of the model with train and test sets; 2) the confusion matrix for the test set
#' @export
#' @examples

getMetrics <- function(model,
                       data.train,
                       data.test,
                       target.train=NULL,
                       target.test=NULL) {

  # Load library
  require(caret)

  # Define covariates
  # if(is.null(target.train)) {
  #   target.train <- data.train$covariate
  #   target.test <- data.test$covariate
  # }

  tab <- table(target.train)

  # Define parameters
  if(is.factor(target.train)) {
    type="class"
  } else {
    type="response"
  }

  # Get predictions
  if(is.null(model$lambda)) { # !cv.glmnet
    predict.train <- as.vector(predict(model, as.data.frame(data.train)))
    predict.test <- as.vector(predict(model, as.data.frame(data.test)))

  } else { # cv.glmnet
    predict.train <- predict(model, s = "lambda.min", newx = data.train, type=type)
    predict.test <- predict(model, s = "lambda.min", newx = data.test, type=type)
  }

  # Get the metrics
  if (type=="class") {
    train <- caret::confusionMatrix(as.factor(as.vector(predict.train)), as.factor(target.train), positive=levels(target.train)[2])
    test <- caret::confusionMatrix(as.factor(as.vector(predict.test)), as.factor(target.test), positive=levels(target.train)[2])
    results <- as.data.frame(rbind(c(as.vector(train$overall[c(1,2)]), as.vector(train$byClass[c(1,2)])),
                                   c(as.vector(test$overall[c(1,2)]), as.vector(test$byClass[c(1,2)]))))
    rownames(results) <- c("train", "test")
    colnames(results) <- c("Accuracy", "Kappa", "Sensitivity", "Specificity")
    cm_test <- test$table
  } else {
    train <- c(MLmetrics::R2_Score(predict.train, target.train), MLmetrics::RMSE(predict.train, target.train))
    test <- c(MLmetrics::R2_Score(predict.test, target.test), MLmetrics::RMSE(predict.test, target.test))
    results <- as.data.frame(rbind(train, test))
    colnames(results) <- c("R2", "RMSE")
    cm_test <- NULL

  }

  return(list(results=results, cm_test=cm_test))
}



#' getLASSOmodels - It runs bootstrapped LASSO t times with different train-test partitions and return the results
#'
#' @param exprData a data frame representing a gene expression matrix with genes as rows and samples as columns
#' @param target a vector representing the response variable to predict (numeric or categorical)
#' @param t the number of LASSO runs
#' @param nfolds the number of folds for the LASSO cross-validation
#' @param train.split the proportion of samples to use to train the model
#' @param seed the seed number to ensure the reproducibility of the results
#' @return a list containing for each run 1) the model created; 2) the coefficients of the model;
#' 3) the metrics of the performance of the model with train and test sets; 4) the confusion matrix for the test set;
#' 5) the samples used to train the model; 6) the samples used to test the model
#' @export
#' @examples

getLASSOmodels <- function(exprData,
                           target,
                           t = 10,
                           nfolds = 5,
                           train.split = 0.7,
                           seed=1234) {

  # Load libraries
  require(glmnet)
  require(doParallel)
  cores <- detectCores()/2
  registerDoParallel(cores)

  # Factor or numeric?
  if(is.factor(target)) {
    tab <- table(target)

    if(length(tab)==3) {
      family <- "multinomial"
      type <- "class"
    } else if (length(tab)==2) {
      family <- "binomial"
      type <- "class"
    }
  } else {
    family <- "gaussian"
    type <- "response"
  }

  # Data split
  ind.train.list <- getSamplingIndex(seed=seed,
                                     exprData=exprData,
                                     train.split=train.split,
                                     t=t)

  # Create an empty list
  models <- list()

  # For each sampling, create a model and test it
  for (iter in 1:t) {
    # cat("Iteration", iter, "\n")

    # Data split
    ind.train <- ind.train.list[[iter]]
    data.train <- t(exprData[, ind.train])
    target.train <- target[ind.train]
    data.test <- t(exprData[, -ind.train])
    target.test <- target[-ind.train]

    # Create the model
    cvfit <- glmnet::cv.glmnet(x=data.train,
                               y=target.train,
                               nfolds = nfolds,
                               alpha = 1,
                               family = family,
                               parallel=TRUE,
                               keep=TRUE)

    # coeffs <- getCoeffs(model=cvfit)

    coeffs <- data.table::setDT(as.data.frame(as.matrix(stats::coef(cvfit, s="lambda.min"))), keep.rownames=T)[-1, ]
    colnames(coeffs) <- c("genes", "coeff")
    coeffs <- coeffs[order(abs(coeffs$coef), decreasing=T), ]
    coeffs <- coeffs[coeffs$coeff!=0, ]

    metrics <- getMetrics(model=cvfit,
                          data.train=data.train,
                          data.test=data.test,
                          target.train=target.train,
                          target.test=target.test)


    models[[paste0("iter", iter)]] <- list(model=cvfit,
                                           coeffs=coeffs,
                                           metrics=metrics$results,
                                           confusionMatrixTest=metrics$cm_test,
                                           trainSamples=rownames(data.train),
                                           testSamples=rownames(data.test))
  }

  return(models)
}


#' getLMmodel - It builds a cross-validated linear regression model with the selected hub genes and return the results
#'
#' @param exprData a data frame representing a gene expression matrix with genes as rows and samples as columns
#' @param target a vector representing the response variable to predict (numeric or categorical)
#' @param nfolds the number of folds for the cross-validation
#' @param genes a vector with the genes to use as features in the final model (features were selected with bootstrapped LASSO)
#' @param train.split the proportion of samples to use to train the model
#' @param seed the seed number to ensure the reproducibility of the results
#' @return a list containing 1) the model created; 2) the coefficients of the model;
#' 3) the metrics of the performance of the model with train and test sets; 4) the confusion matrix for the test set;
#' 5) the samples used to train the model; 6) the samples used to test the model
#' @export
#' @examples

getLMmodel <- function(exprData,
                       target,
                       nfolds = 5,
                       train.split = 0.7,
                       genes,
                       seed=1234) {
  # Load libraries
  require(caret)

  # Factor or numeric?
  if(is.factor(target)) {
    tab <- table(target)

    if(length(tab)==3) {
      family <- "multinomial"
      metric <- "Accuracy"
      classProbs <- TRUE
    } else if (length(tab)==2) {
      family <- "binomial"
      metric <- "Accuracy"
      classProbs <- TRUE
    }
  } else {
    family <- "gaussian"
    classProbs <- FALSE
    metric <- "RMSE"
  }

  # Data split
  # genes <- gsub("\\.", "-", genes)
  mydata <- exprData[match(genes, rownames(exprData)), ]

  ind.train.list <- getSamplingIndex(seed=seed,
                                     exprData=mydata,
                                     train.split=train.split,
                                     t=1)[[1]]

  # Data split
  ind.train <- ind.train.list
  data.train <- t(mydata[, ind.train])
  target.train <- target[ind.train]
  data.test <- t(mydata[, -ind.train])
  target.test <- target[-ind.train]

  # Make names
  if(family=="binomial") {
    target.train <- factor(make.names(target.train), levels=c("X0", "X1"))
    target.test <- factor(make.names(target.test), levels=c("X0", "X1"))

  } else if(family=="multinomial") {
    target.train <- factor(make.names(target.train), levels=c("X0", "X1", "X2"))
    target.test <- factor(make.names(target.test), levels=c("X0", "X1", "X2"))
  }

  # fitControl
  fitControl <- trainControl(method = "repeatedcv",
                             number = nfolds,
                             repeats = 10,
                             allowParallel = TRUE,
                             savePredictions = T,
                             classProbs = classProbs,
                             summaryFunction = defaultSummary)

  # train the model
  model <- train(x=data.train,
                 y=target.train,
                 trControl = fitControl,
                 metric = metric,
                 method = "glmnet",
                 tuneGrid = expand.grid(.alpha=0,
                                        .lambda=10^seq(-1, -4, length = 10)),
                 preProcess=c("center", "scale"),
                 family=family)

  if(family!="multinomial") {
    coeffs <- as.data.frame(as.matrix(coef(model$finalModel, unlist(model$bestTune$lambda))))
  } else {
    coeffs <- as.matrix(coef(model$finalModel, unlist(model$bestTune$lambda))$X1)
  }

  coeffs <- data.frame(gene=rownames(coeffs), coeff=as.vector(unlist(coeffs)))
  coeffs <- coeffs[coeffs$gene!="(Intercept)", ]
  coeffs <- coeffs[order(abs(coeffs$coeff), decreasing=T), ]
  coeffs$coeff <- as.numeric(as.character(coeffs$coeff))

  ci <- model$results[model$results$lambda==model$bestTune$lambda, ]

  # Get metrics
  metrics <- getMetrics(model=model,
                        data.train=data.train,
                        data.test=data.test,
                        target.train=target.train,
                        target.test=target.test)

  return(list(model=model,
              coeffs=coeffs,
              crossValidation_ci=ci,
              metrics=metrics,
              confusionMatrixTest=metrics$cm_test,
              trainSamples=rownames(data.train),
              testSamples=rownames(data.test)))

}


#' getHubGenes - It applies bootstrapped LASSO t times and extract the most relevant features of each run.
#' Then, features are grouped based on their ratio of appearance.
#' For each ratio, a final model is built with the features that meet this ratio.
#'
#' @param exprData a data frame representing a gene expression matrix with genes as rows and samples as columns
#' @param target a vector representing the response variable to predict (numeric or categorical)
#' @param train.split the proportion of samples to use to train the model
#' @param nfolds the number of folds for the cross-validation
#' @param t the number of LASSO runs
#' @param seed the seed number to ensure the reproducibility of the results
#' @param cutoffs the ratios of appearance selected to create a final model with the corresponding genes
#' @param path the path where results are stored
#' @param tissueName the name of the tissue, cohort or dataset
#' @param targetName the name of the target
#' @param force if force=T, all the ratios of appearance with at least one gene will lead to the corresponding final linear regression model.
#' If force=F, only the ratios of appearance with at least 5 genes will lead to the corresponding final models. to create the final models.
#' @param save if save=T, results will be saved in the path indicated.
#' @return a list containing 1) LASSO models, one for each run; 2) Linear regression models, one for each ratio of appearance selected.
#' @export
#' @examples

getHubGenes <- function(exprData,
                        target,
                        train.split=0.7,
                        nfolds=5,
                        t=10,
                        seed=seed,
                        cutoffs=NULL,
                        path=getwd(),
                        tissueName="tissue1",
                        targetName="target1",
                        force=F,
                        save=T) {

  # Create directory if necessary
  if(!dir.exists(paste0(path, "/hubGenes/"))) {
    dir.create(paste0(path, "/hubGenes/"))
    file <- NA
  } else {# Check if the file is already created
    file <- list.files(path=paste0(path, "/hubGenes/"), full.names=T)
    file <- file[match(paste0(path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"), file)]
  }

  # If the file is not create or you want to overwrite it
  if(is.na(file)) {

    # Apply lasso t times for feature selection purposes
    cat("Apply LASSO", t, "times for feature selection \n")
    models_lasso <- getLASSOmodels(exprData=exprData,
                                   target=target,
                                   t = t,
                                   nfolds = nfolds,
                                   train.split = train.split,
                                   seed=seed)
    # Save lasso models
    result <- list()

    # Get hub genes enrichment for each cutoff
    sumLogPval <- getGeneSetEnrichment(lasso_models=models_lasso)
    enrichment <- sumLogPval$enrich
    sumLogPval <- sumLogPval$sumLogPval

    result[["LASSO"]] <- list(models=models_lasso, hubs_enrich=enrichment, hubs_enrich_stats=sumLogPval)

    # Select the cutoffs to evaluate
    hubs <- table(unlist(lapply(models_lasso, function(x) x$coeffs$genes)))
    hubGenesPerCutoff <- sapply(seq(1,t,1), function(x) names(hubs)[hubs>=x])
    names(hubGenesPerCutoff) <- paste0("cutoff", 1:t)

    # Cutoffs selected
    len <- rev(lengths(hubGenesPerCutoff))
    cat("Number of hubs per ratio of appearance\n")
    print(len)

    if(is.null(cutoffs)) {
      cutoffs <- sort(as.numeric(gsub("cutoff", "", names(hubGenesPerCutoff[lengths(hubGenesPerCutoff)>1]))), decreasing=T)
    } else {
      cutoffs <- intersect(cutoffs, as.numeric(gsub("cutoff", "", names(hubGenesPerCutoff[lengths(hubGenesPerCutoff)>1]))))
    }

    cat("Cutoffs selected are", cutoffs, "\n")
    lm_models <- list()

    # For each cutoff
    for(c in cutoffs) {

      # Get final model and cv error
      hubGenes <- hubGenesPerCutoff[[paste0("cutoff", c)]]
      cat("Get final model and the cv error for cutoff",  c, "\n")

      finalModel <- getLMmodel(exprData=exprData,
                              target=target,
                              nfolds = nfolds,
                              train.split = train.split,
                              genes=hubGenes,
                              seed=seed)

      lm_models[[paste0("cutoff", c)]] <- list(lm_models=finalModel, lm_genes_selected=finalModel$coeffs)

    }

    # Save the results
    result[["LM"]] <- lm_models

    if(save==T) {
      saveRDS(result, paste0(path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"))
    }

    } else {
      result <- readRDS(file)
    }

  return(result)

}



#' getGeneSetEnrichment - It applies a functional enrichment analysis for a gene set from both LASSO models or networks.
#'
#' @param lasso_models a list of LASSO models
#' @param net a data frame containing the structure of the network (hub gene, gene and correlation columns)
#' @param sources a vector with the names of the sources to be tested
#' @param hubs if the net parameter is not null and hubs=T, the enrichment of the hub genes will also be tested
#' @return a list containing 1) the annotations obtained for each query; 2) the summary of the enrichment for each query
#' @export
#' @examples

getGeneSetEnrichment <- function(lasso_models=NULL, net=NULL,
                                 sources = c("GO", "KEGG", "REAC", "HP"), hubs=T) {

  # Load libraries
  require(gprofiler2)
  all.genes <- list()

  if(!is.null(lasso_models)) { # Get the enrichment of the hub genes per cutoff

    hubGenes <- table(unlist(lapply(lasso_models, function(x) x$coeffs$genes)))
    hubGenesPerCutoff <- sapply(seq(1,length(lasso_models),1), function(x) names(hubGenes)[hubGenes>=x])
    names(hubGenesPerCutoff) <- paste0("cutoff", 1:length(lasso_models))
    all.genes <- hubGenesPerCutoff

  } else if(!is.null(net)) { # Get the enrichment for each module of the CS-GCN

    if(!is.null(net$hubGene)) {
      for(hub in unique(net$hubGene)) {
        all.genes[[hub]] <- net$genes[net$hubGene==hub]
      }
      if(hubs==T) {
        all.genes[["hubs"]] <- unique(net$hubGene)
      }

    } else { # Get the enrichment for each wgcna module

      for(hub in unique(as.vector(net$moduleColors))) {
        all.genes[[hub]] <- names(net$moduleColors)[net$moduleColors==hub]
      }
      if(hubs==T) {
        all.genes[["hubs"]] <- unique(as.vector(net$moduleColors))
      }
    }
  }

  # Get enrichment
  enrich <- gprofiler2::gost(all.genes,
                             correction_method="g_SCS",
                             sources=sources,
                             organism = "hsapiens",
                             exclude_iea = T)$result

  if(!is.null(enrich)) {
    enrich <- enrich[order(enrich$query, enrich$p_value, decreasing=F), ]
    enrich$query <- factor(enrich$query, levels=names(all.genes))

    # Sum the -log10(pval) per query (cutoff)
    enrich$logPval <- -log10(enrich$p_value)
    sumLogPval <- tapply(enrich$logPval, factor(enrich$query, levels=names(all.genes)), FUN="sum")
    sumLogPval[is.na(sumLogPval)] <- 0
    sumLogPval <- data.frame(query=names(sumLogPval), sum=round(as.numeric(sumLogPval), 3))

    # Divide the sum of log10(pval) by the number of hub genes
    numGenes <- lengths(all.genes)
    stopifnot(identical(names(numGenes), as.vector(sumLogPval$query)))
    sumLogPval$numGenes <- numGenes
    sumLogPval$sum_corrected <- sumLogPval$sum/sumLogPval$numGenes
  } else {
    enrich <- NULL
    sumLogPval <- NULL
  }


  return(list(enrichment=enrich, sumLogPval=sumLogPval))
}



#' getCorrelations - It estimates the correlations between the hub genes and the rest of genes in the expression matrix.
#'
#' @param exprData a data frame representing the gene expression matrix with genes as rows and samples as columns
#' @param hubGenes a vector with the name of the hub genes selected
#' @return a correlation matrix with as many rows as rows in the expression matrix and as many columns as the number of hub genes.
#' @export
#' @examples

getCorrelations <- function(exprData,
                            hubGenes) {

  # Get correlation matrix
  matrix_all <- t(exprData)
  matrix_hubs <- t(exprData[match(hubGenes, rownames(exprData)), ])
  df <- as.data.frame(t(cor(matrix_hubs, matrix_all)))

  return(df)
}

#' getModules - It adds the most correlated genes to each hub gene to get the modules. The number of genes added to each hub gene depends
#' on the approach selected: fixed, where all modules have the same size; coefficient, it depends on the hub gene coefficient in the final model;
#' enrichment, it looks for the smallest gene set with the maximum enrichment.
#'
#' @param hubs a data frame containing the genes selected for a specific ratio of appearance together with
#' their coefficients in the corresponding linear regression model
#' @param exprData a data frame representing the gene expression matrix with genes as rows and samples as columns
#' @param n the maximum size of a module for all the approaches
#' @param s the minimum size of a module for all the approaches
#' @param m the number of genes added in each iteration when approach="enrichment"
#' @param minCor the minimum correlation of a gene to be added to a module
#' @param maxTol the maximum number of tries to get a better enrichment without getting it
#' @param approach the approach to complete the modules (fixed, coefficient or enrichment)
#' @return a list containing 1) the network; 2) the plot of the correlation distribution for each module created.
#' If approach="enrichment", it also includes 3) the statistics of the process of creating the modules; 4) the plot of the statistics.
#' @export
#' @examples

getModules <- function(hubs,
                       exprData,
                       n=100,
                       s=10,
                       m=10,
                       minCor=0.3,
                       maxTol=3,
                       approach=c("fixed", "coefficient", "enrichment")) {

  # Load library
  require(psych)
  require(ggplot2)

  # Get correlation matrix
  M <- getCorrelations(exprData=exprData,
                       hubGenes=hubs$gene)

  # Get coeffs
  coeffs <- data.frame(gene=hubs$gene, estimate=hubs$coeff)
  net <- data.frame()
  stats <- data.frame()

  if(approach=="coefficient") {

    for (i in 1:nrow(coeffs)) {
      hub <- coeffs$gene[i]

      if(i==1) {
        nGenes <- n
      } else {
        nGenes <- round(abs((coeffs$estimate[i]*n)/coeffs$estimate[1]), 0)
        if(nGenes<s) {
          nGenes <- s
        }
      }

      df <- data.frame(genes=rownames(M), cor=M[, match(hub, colnames(M))])
      df <- df[order(abs(df$cor), decreasing=T), ]
      df <- df[df$cor>abs(minCor), ]
      df <- df[1:nGenes, ]
      df <- data.frame(hubGene=rep(hub, nrow(df)), df)
      df <- df[complete.cases(df), ]
      net <- rbind(net, df)
    }

    stats=NULL
    plot_stats=NULL

  } else if (approach=="fixed"){

    for (i in 1:nrow(coeffs)) {
      hub <- coeffs$gene[i]
      df <- data.frame(genes=rownames(M), cor=M[, match(hub, colnames(M))])
      df <- df[order(abs(df$cor), decreasing=T), ]
      df <- df[df$cor>abs(minCor), ]
      df <- df[1:n, ]
      df <- data.frame(hubGene=rep(hub, nrow(df)), df)
      df <- df[complete.cases(df), ]
      net <- rbind(net, df)
    }

    stats=NULL
    plot_stats=NULL

  } else if(approach=="enrichment") {

    # For each hub gene

    for(i in 1:nrow(coeffs)) {
      mystats <- data.frame()
      tol <- 0 # the number of times we try to get a better E(G) without getting it
      bestE <- 0 # best E(G*) until now
      bestG <- c() # gene set with best E(G*) until now

      h <- coeffs$gene[i]  # hub gene
      # print(h)

      C <- M[, match(h, colnames(M))]
      C <- data.frame(genes=rownames(M), cor=C)
      C$cor <- round(C$cor, 3)
      C <- C[C$cor>abs(minCor), ]
      C <- C[order(abs(C$cor), decreasing=T), ] # gene set ordered by decreasing correlation with h

      C$genes <- gsub("\\.", "-", C$genes)

      size <- s # initially, the size of the module (size) is the same as the minimum size (s)

      # If we found no genes correlated with h with a minCor of x, stop here
      if(nrow(C)<size & s>10) {
        stop("We found less than ", size, " genes with a minimum correlation of ", minCor, ", please, select other value for minCor or s\n")
      } else if(nrow(C)<size & s==10) {
        size <- nrow(C)
      }

      # While the size of the module (size) is lower than the maximum size (m)
      # and the number of times we try to get a better E(G) without getting it is lower than the maximum number of times allowed (maxTol)
      while(size<n & tol<maxTol & size<=nrow(C)) {
        h_expr <- as.vector(unlist(exprData[match(h, rownames(exprData)), ])) # hub gene expression
        C_expr <- as.data.frame(t(exprData[match(C$genes[1:size], rownames(exprData)), ])) # gene set expression
        c <- corr.test(C_expr, h_expr, adjust="BH") # correlation between the gene set and the hub
        G <- rownames(c$p.adj)[c$p.adj<0.05] # gene set with a minimum correlation > minCor and padj<0.05

        # If we find genes with p.adj<0.05
        if(length(G)>0) {
          mynet <- c()
          mynet$moduleColors <- rep(h, length(G))
          names(mynet$moduleColors) <- G

          E <- getGeneSetEnrichment(net=mynet, hubs=F)$sumLogPval # E(G)

          if(is.null(E)) { # if we found no enrichment, add this row to the table
            E <- data.frame(query=h, sum=0, numGenes=size, sum_corrected=0)
          }

          E <- E$sum_corrected # E(G)

          # If this is the first iteration, bestE is E and bestG is G
          # Or, if this is not the first iteration and bestE>E
          if(E>bestE || size==s || size==nrow(C)) {
            bestE <- E
            bestG <- G

            # If this is not the first iteration and E<bestE
          } else if(E<bestE) {
            tol <- tol + 1
          }

          mystats <- rbind(mystats, c(h, size, length(G), round(E, 3), round(min(abs(c$r)),3), format(max(c$p.adj[-1]), sci=T, digits=3), "NO"))
          size <- size+m

        } else {
          size <- n
          tol <- maxTol
          cat("We found no genes with padj<0.05 and corr>", minCor, "using", size, "size\n")
        }

      }
      # Which row of stats contain the bestG/bestE?
      index <- which(mystats[, 3]==length(bestG))
      mystats[index, ncol(mystats)] <- "YES"
      colnames(mystats) <- c("hubGene", "genes", "filteredGenes", "enrichment", "minCor", "maxpadj", "bestE")

      if(mystats$enrichment[mystats$bestE=="YES"]>0) {
        stats <- rbind(stats, mystats)
      }

      # Add the bestG as one module in the network
      mymodule <- cbind(hubGene=rep(h, length(bestG)), C[match(bestG, C$genes), ])
      net <- rbind(net, mymodule)
    }

    # Plot the enrichment for each gene set for each hub gene
    df <- stats
    levels <- sort(unique(df$hubGene))
    df$enrichment <- as.numeric(as.character(df$enrichment))
    df$hubGene <- factor(df$hubGene, levels=levels)

    plot_stats <- ggplot(df, aes(x=filteredGenes, y=enrichment, color=hubGene, group=hubGene)) +
                         geom_point(size = 1.8, aes(shape = bestE)) +
                         scale_shape_manual(values=c(16, 17)) +
                         geom_line() +
                         theme_classic() +
                         theme(text=element_text(size=14),
                              legend.key.size=unit(0.3, "cm"),
                              legend.text = element_text(size=8)) +
                         ylab(expression(paste("sum (-", log[10] ~ "Pval)/", "\n module size"))) +
                         xlab("Module size") +
                         labs(color="Seed") +
                         labs(shape="best enrichment")

    }

    # Plot correlation distribution for each module
    df <- net
    levels <- sort(unique(df$hubGene))
    df$hubGene <- factor(df$hubGene, levels=levels)
    df <- df[df$cor<1, ]

    plot_net <- ggplot(df, aes(x=hubGene, y=cor, color=hubGene)) +
                       geom_boxplot() +
                       theme_classic() +
                       theme(text=element_text(size=14),
                            legend.key.size=unit(0.3, "cm"),
                            legend.text = element_text(size=8),
                            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                       ylab("|Pearson correlation|") +
                       xlab("Seed") +
                       labs(color="Seed")

  net$hubGene <- factor(net$hubGene, levels=coeffs$gene)

  return(list(stats=stats, plot_stats=plot_stats,
              net=net, plot_net=plot_net))
}


#' saveResults - It saves LASSO models and linear regression models performance to files
#'
#' @param tissueName the name of the tissue, cohort or dataset
#' @param targetName the name of the target
#' @param hubGenes the getHubGenes output containing both LASSO models and linear regression models
#' @param path the path where the files will be saved
#' @return files saved in the path indicated by the user
#' @export
#' @examples

saveResults <- function(tissueName, targetName, hubGenes, path) {

  if(!is.null(hubGenes)) {
    lasso_models_stats <- do.call("rbind", lapply(hubGenes$LASSO$models, function(x) x$metrics))
    lasso_models_stats$iter <- gsub("\\..*.", "", rownames(lasso_models_stats))
    lasso_models_stats$set <- gsub("..*.\\.", "", rownames(lasso_models_stats))
    rownames(lasso_models_stats) <- NULL
    write.csv(lasso_models_stats, paste0(path, "/results/", targetName, "_", tissueName, "_", "lasso_stats.csv"), row.names = F)

    lm_models_stats <- do.call("rbind", lapply(hubGenes$LM, function(x) x$lm_models$metrics$results))
    lm_models_stats$cutoff <- gsub("\\..*.", "", rownames(lm_models_stats))
    lm_models_stats$set <- gsub("..*.\\.", "", rownames(lm_models_stats))
    rownames(lm_models_stats) <- NULL
    write.csv(lm_models_stats, paste0(path, "/results/", targetName, "_", tissueName, "_", "lm_stats.csv"), row.names = F)
  }

  return()
}


#' getModulesAnnotation - It saves LASSO models and linear regression models performance to files
#'
#' @param net the output of getModules function
#' @param exprData a data frame representing a gene expression matrix with genes as rows and samples as columns
#' @param covs a data frame with donors metadata where samples are the rows and features are the columns
#' @param path the path where the files will be saved
#' @param save if save=T, the results of the analysis will be saved in separate files
#' @param overwrite if overwrite=T, the analysis will be repeated and files will be overwritten
#' @param cutoff the name of the cutoff of this TGCN
#' @param tissueName the name of the tissue, cohort or dataset
#' @param targetName the name of the target
#' @return a list containing 1) the net, the eigengenes and the plot of the correlation per module; 2) the GO enrichment annotations and statistics; 3) the cell-type markers enrichment analysis results;
#' 4) the module-trait association analysis results; 5) the crossTabPlot analysis resylts between the modules of the TGCN.
#' @export
#' @examples

getModulesAnnotation <- function(net,
                                 exprData,
                                 covs,
                                 path,
                                 save=T,
                                 overwrite=T,
                                 cutoff,
                                 tissueName,
                                 targetName,
                                 reduced=F) {

  # TGCN already characterized?
  if(!dir.exists(paste0(path, "/results/"))) {
    dir.create(paste0(path, "/results/"))
    files <- NA
  } else {
    files <- list.files(path=paste0(path, "/results/"), pattern="csv", full.names=T)
    files <- files[grep(cutoff, files)]
  }

  # If not
  if(is.na(files) || overwrite==T) {

  # Get the enrichment per module
  mynet <- net$net
  enrich <- getGeneSetEnrichment(net=mynet)
  sumLogPval <- enrich$sumLogPval
  sumLogPval$name <- rep(cutoff, nrow(sumLogPval))

  terms <- enrich$enrichment
  terms$parents <- unlist(lapply(terms$parents, function(x) paste0(x, collapse=", ")))

  # Reduced GO terms
  go_results <- enrich$enrichment
  if(reduced==T) {
    go_reduced <- plotReducedGOterms(go_results[go_results$source=="GO:BP", ], module=F)
  } else {
    go_reduced <- NULL
  }

  # Module-trait corr
  colors <- mynet$hubGene
  names(colors) <- mynet$genes
  expr.data <- t(exprData)
  expr.data <- expr.data[, match(mynet$genes, colnames(expr.data))]
  MEs <- WGCNA::moduleEigengenes(expr.data, colors, excludeGrey=F)$eigengenes

  if(!is.null(covs) & save==T) {
    png(filename=paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_moduleTraitCorr.png"), width=12, height=15, unit="cm", res=500)
    traits <- plotModuleTraitCorr(MEs=MEs, covs=covs, ylab=paste0(targetName, " modules"), max=6,
                                    height=8, width=6)
    dev.off()
    colnames(traits$pval) <- paste0(colnames(traits$pval), "Pval")

    if(!is.null(traits$cor)) {
      colnames(traits$cor) <- paste0(colnames(traits$cor), "Cor")
    }

    traits_df <- as.data.frame(cbind(rownames(traits$pval), traits$cor, traits$pval))
    colnames(traits_df)[1] <- "trait"
    rownames(traits_df) <- NULL
    write.csv(traits_df, paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_moduleTraitCorr.csv"), row.names=F)

  } else {
    traits <- list()
    traits$cor <- NULL
    traits$pval <- NULL
    traits$htCor <- NULL
    traits$htPval <- NULL
  }

  # Get CT enrichment
  ct_enrich <- getCTenrich(mynet)

  if(!is.null(ct_enrich) & save==T) {
    png(filename=paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_CTenrich.png"), width=12, height=12, unit="cm", res=500)
    plot_ct <- plotCTenrich(ct=ct_enrich$matrix, target=targetName, height=10, width=10)
    dev.off()
  } else {
    plot_ct <- NULL
    ct_enrich$df <- NULL
    ct_enrich$matrix <- NULL
  }

  # Save files
  if(save==T) {
    write.csv(mynet, paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN.csv"), row.names=F)
    write.csv(terms, paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_GOenrich.csv"), row.names=F)
    write.csv(enrich$sumLogPval, paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_GOenrich_stats.csv"), row.names=F)
    write.csv(ct_enrich$df, paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_CTenrich.csv"))
  }

  # CrossTabPlot
  png(filename=paste0(path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_crossTabPlot.png"), width=15, height=15, unit="cm", res=500)
  ctp <- plotModulesOverlap(name1=targetName,
                            name2=targetName,
                            tgcn1=mynet,
                            tgcn2=mynet,
                            background1=rownames(exprData),
                            background2=rownames(exprData),
                            colsOrder=NULL,
                            significant=F,
                            moduleSize=F,
                            main="",
                            diag=T)

  dev.off()

  plots <- plotGOenrichSummary(enrich)

  # Save results in a list
  result <- list()
  result <- list(net=list(modules=mynet, MEs=MEs,
                          plotCorr=net$plot_net,
                          moduleSizeSelectionStats=net$stats,
                          moduleSizeSelectionPlot=net$plot_stats),
                 GOenrich=list(terms=enrich$enrichment[, -ncol(enrich$enrichment)],
                               stats=enrich$sumLogPval,
                               reducedTerms=go_reduced$TGCN$reducedTerms,
                               plotNterms=plots$nterms + labs(title=cutoff),
                               plotStats=plots$stats + labs(title=cutoff),
                               plotReduced=go_reduced$TGCN$scatterPlot),
                 CTenrich=list(sigTests=ct_enrich$df, allTests=ct_enrich$matrix, plot=plot_ct),
                 moduleTraitCorr=list(corr=traits$cor, plot_cor=traits$htCor,
                                      pval=traits$pval, plot_pval=traits$htPval),
                 crossTabPlot=list(pval=ctp$df, plot=paste0("see ", path, "/results/", targetName, "_", tissueName, "_", cutoff, "_TGCN_crossTabPlot.png figure")))
  } else {
    result <- readRDS(paste0(path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds"))
  }

  return(result)
}


#' testAllCutoffs - It applies bootstrapped LASSO t times and extract the most relevant genes of each iteration. Then, genes are
#' grouped based on their ratio of appearance. For each ratio, a linear regression model is built and then, each hub gene is completed
#' adding the most correlated genes to create a module.
#'
#' @param exprData a data frame representing a gene expression matrix with genes as rows and samples as columns
#' @param target a vector representing the response variable to predict (numeric or categorical)
#' @param covs a data frame with donors metadata where samples are the rows and features are the columns
#' @param targetName the name of the target
#' @param tissueName the name of the tissue, cohort or dataset
#' @param train.split train.split the proportion of samples used to train the models
#' @param nfolds the number of folds for the cross-validation
#' @param t the number of LASSO runs
#' @param cutoffs the ratios of appearance selected to create a final model and a TGCN with the corresponding genes
#' @param n the maximum size of a module for all the approaches
#' @param s the minimum size of a module for all the approaches
#' @param m the number of genes added in each iteration when approach="enrichment"
#' @param minCor the minimum correlation of a gene to be added to a module
#' @param maxTol the maximum number of tries to get a better enrichment without getting it
#' @param approach the approach to complete the modules (fixed, coefficient or enrichment)
#' @param seed the seed number to ensure the reproducibility of the results
#' @param save if save=T, the results of the analysis will be saved in separate files
#' @param overwrite if overwrite=T, the analysis will be repeated and files will be overwritten
#' @param path the path where results are stored
#' @return
#' @export
#' @examples

testAllCutoffs <- function(exprData,
                           target,
                           covs,
                           targetName="target",
                           tissueName="tissue",
                           train.split=0.7,
                           nfolds=5,
                           t=10,
                           cutoffs=10:1,
                           n=100,
                           m=10,
                           s=10,
                           minCor=0.3,
                           maxTol=3,
                           approach=c("fixed", "coefficient", "enrichment"),
                           seed=1234,
                           save=T,
                           overwrite=T,
                           path=getwd(),
                           reduced=F,
                           report=F) {

  # Load libraries
  require(WGCNA)
  require(stringr)
  require(RColorBrewer)
  require(ggplot2)
  require(caret)
  require(CoExpNets)
  require(ComplexHeatmap)
  require(circlize)
  require(glmnet)
  require(doParallel)
  require(gprofiler2)
  require(psych)
  require(rrvgo)

  # Hub genes already tested?
  if(!dir.exists(paste0(path, "/hubGenes/"))) {
    dir.create(paste0(path, "/hubGenes/"))
    file <- NA
  } else {# Check if the file is already created
    file <- list.files(path=paste0(path, "/hubGenes/"), full.names=T)
    file <- file[match(paste0(path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"), file)]
  }

  if(is.na(file)) {
    cat("- Step 1: Create a linear regression model for each gene set of hubs based on their ratio of appearance \n")
    # Step 1:
    # Apply LASSO 10 times
    # Create a linear regression model for each set of hubs based on the ratio of appearance
    hubGenes <- getHubGenes(exprData=exprData,
                            target=target,
                            train.split=train.split,
                            nfolds=nfolds,
                            t=t,
                            path=path,
                            tissueName=tissueName,
                            targetName=targetName,
                            seed=seed,
                            cutoffs=cutoffs,
                            save=F,
                            force=T)

    if(save==T) {
      saveRDS(hubGenes, paste0(path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"))
    }
  } else {
    hubGenes <- readRDS(paste0(path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"))
  }


  # All cutoffs already tested?
  if(!dir.exists(paste0(path, "/Net/"))) {
    dir.create(paste0(path, "/Net/"))
    file <- NA
  } else {# Check if the file is already created
    file <- list.files(path=paste0(path, "/Net/"), full.names=T)
    file <- file[match(paste0(path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds"), file)]
  }


  if(is.na(file) || overwrite==T) {
    cat("\n- Step 2: represent train and test error per ratio of appearance model\n")
    # Step 2:
    # per ratio of appearance, get the train and test error as well as the number of hubs used.
    results <- list()
    results[["selectRatio"]] <- selectAppearanceRatio(hubGenes)

    cat("\n- Step 3: creating a network for each ratio of appearance where the number of hubs is <=30\n")
    # Step 3:
    # Create a network for each cutoff where nHubs<=30
    nHubs <- lengths(lapply(hubGenes$LM, function(x) x$lm_genes_selected$gene))
    nHubs <- nHubs[nHubs<=30]
    names <- as.numeric(gsub("cutoff", "", names(nHubs)))

    for(cutoff in names) {
      cat("\n Step 3.1: TGCN creation for ratio of appearance", cutoff, "\n")

      # Complete the modules
      hubs <- hubGenes$LM[[paste0("cutoff", cutoff)]]
      net1 <- getModules(hubs=hubs$lm_genes_selected,
                         exprData=exprData,
                         n=n,
                         s=s,
                         m=m,
                         minCor=minCor,
                         maxTol=maxTol,
                         approach=approach)

      cat("\n Step 3.2: TGCN characterization for cutoff", cutoff, "\n")
      net2 <- getModulesAnnotation(net=net1,
                                   exprData=exprData,
                                   covs=covs,
                                   save=save,
                                   path=path,
                                   overwrite=T,
                                   cutoff=paste0("c", cutoff),
                                   tissueName=tissueName,
                                   targetName=targetName,
                                   reduced=reduced)

      results[["nets"]][[paste0("c", cutoff)]] <- net2

    }

    if(save==T) {
      saveRDS(results, paste0(path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds"))
    }

  } else {
    results <- readRDS(paste0(path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds"))
  }

  if(report==T) {
    template=list.files(system.file("report", "", package = "TGCN"), full.names=T)

    rmarkdown::render(input = template,
                      output_file = paste0(path, "/results/", targetName, "_", tissueName, "_TGCNs.html"),
                      params = list(target=targetName, tissue=tissueName, path=path))
  }

  return(results)

}




