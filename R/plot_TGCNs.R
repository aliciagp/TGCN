
#' selectAppearanceRatio - It represents, per ratio of appearance, the number of hubs included in the model together with the train and
#' test performance of the model (Accuracy/R2).
#'
#' @param hubGenes the output of getHubGenes function containing both LASSO models and linear regression models
#' @return a list containing 1) the plot of the number of hubs selected per ratio of appearance; 2) the plot of the train and test error
#' for the model created for each ratio of appearance
#' @export
#' @examples

selectAppearanceRatio <- function(hubGenes) {

  require(ggplot2)

  nHubs <- lengths(lapply(hubGenes$LM, function(x) x$lm_genes_selected$gene))
  nHubs <- rep(nHubs, each=2)
  error <- do.call("rbind", lapply(hubGenes$LM, function(x) x$lm_models$metrics$results))
  error$cutoff <- gsub("\\..*.", "", rownames(error))
  error$set <- gsub(".*.\\.", "", rownames(error))
  stopifnot(identical(names(nHubs), error$cutoff))

  df <- cbind(error, nHubs)
  df$cutoff <- gsub("cutoff", "", df$cutoff)
  df$cutoff <- factor(df$cutoff, levels=1:length(hubGenes$LM))
  df$set <- factor(df$set, levels=c("train", "test"))
  ylab <- colnames(df)[1]
  colnames(df)[1] <- "metric"

  p1 <- ggplot(df, aes(x=cutoff, y=nHubs)) +
               geom_point() +
               theme_classic() +
               theme(text=element_text(size=14)) +
               ylab("Number of hub genes") +
               xlab("Ratio of appearance")

  p2 <- ggplot(df, aes(x=cutoff, y=metric, color=set, group=set)) +
               geom_point() +
               geom_line() +
               theme_classic() +
               theme(text=element_text(size=14)) +
               ylab(ylab) +
               xlab("Ratio of appearance")

  return(list(nHubs=p1, stats=p2))
}

#' plotReducedGOterms - It plots the reduced GO:BP terms at TGCN level and module level.
#'
#' @param go_results a data frame containing the enrichment from gprofiler2::gost function.
#' @param module if TRUE, the reduced GO terms are plotted for each module separately.
#' @return a list containing 1) the similarity matrix; 2) the reduced GO terms; 3) the scatter plot of the reduced terms.
#' @export
#' @examples

plotReducedGOterms <- function(go_results, module=T) {

  # Load library
  require(rrvgo)

  # Remove hub genes enrichment
  go_results <- go_results[go_results$query!="hubs", ]

  mylist <- list()
  mylist[["TGCN"]] <- go_results
  names <- unique(go_results$query)

  if(module==T) {
    for(name in names) {
      mylist[[name]] <- go_results[go_results$query==name, ]
    }
  }

  results <- list()
  name <- names(mylist)

  for (name in names(mylist)) {
    simMatrix <- calculateSimMatrix(mylist[[name]]$term_id,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel")

    if(!is.null(dim(simMatrix))) {
      if(nrow(simMatrix)>5) {
        scores <- setNames(-log10(mylist[[name]]$p_value), mylist[[name]]$term_id)
        reducedTerms <- reduceSimMatrix(simMatrix,
                                        scores,
                                        threshold=0.7,
                                        orgdb="org.Hs.eg.db")

        results[[name]][["simMatrix"]] <- simMatrix
        results[[name]][["reducedTerms"]] <- reducedTerms
        results[[name]][["scatterPlot"]] <- scatterPlot(simMatrix, reducedTerms)
      }
    }
  }

  return(results)
}



#' plotModulesOverlap - It applies a Fisher exact test for each pair of modules and represent the results. Only significant
#' overlaps are colored. Each cell represents the number of genes that overlap and the color represents the -log10 p-value
#' adjusted.
#'
#' @param name1 the name of the TGCN1
#' @param name2 the name of the TGCN2 in case TGCN1!=TGCN2
#' @param tgcn1 the TGCN whose modules we want to compare with TGCN2
#' @param tgcn2 the TGCN whose modules we want to compare with TGCN1 in case TGCN1!=TGCN2
#' @param background1 the pool of genes used to create the TGCN1
#' @param background2 the pool of genes used to create the TGCN2 in case TGCN1!=TGCN2
#' @param colsOrder the order of the TGCN2 modules (columns) in the plot
#' @param significant if significant=T, the cells corresponding to significant overlaps are represented by *** instead
#'of showing the number of genes that overlaps.
#' @param moduleSize if moduleSize=T, the size of each module is included in the corresponding label in the plot
#' @param main the main of the plot
#' @param diag if diag=T, the diagonal of the plot will also be colored based on the significance, if not, only significant
#' overlaps different from the diagonal will be colored.
#' @return a list containing 1) the Fisher exact test results for each pair of modules compared;
#' 2) the plot of the Fisher exact test results.
#' @export
#' @examples

plotModulesOverlap <- function(name1="TGCN1",
                               name2="TGCN2",
                               tgcn1,
                               tgcn2=NULL,
                               background1=NULL,
                               background2=NULL,
                               colsOrder=NULL,
                               significant=F,
                               moduleSize=F,
                               main="",
                               diag=T) {

  # Load library
  require(WGCNA)

  # If we only have one TGCN
  if(is.null(tgcn2)) {
    tgcn2 <- tgcn1
    name2 <- name1
  }

  # If we are dealing with WGCNA networks, save module composition in the same way as our method does
  if(is.null(tgcn1$genes)) {
    mynet <- list()
    mynet <- stack(tgcn1$moduleColors)
    colnames(mynet) <- c("hubGene", "genes")
    tgcn1 <- mynet
    wgcna=T
  } else {
    tgcn1 <- tgcn1[, 1:2]
  }
  # The same thing for tgcn2
  if(is.null(tgcn2$genes)) {
    mynet <- list()
    mynet <- stack(tgcn2$moduleColors)
    colnames(mynet) <- c("hubGene", "genes")
    tgcn2 <- mynet
    wgcna=T
  } else {
    tgcn2 <- tgcn2[, 1:2]
  }

  # Create a fake module to include the non-common genes if necessary
  diff <- setdiff(tgcn1$genes, tgcn2$genes)
  index1 <- match(diff, tgcn1$genes)
  index2 <- match(diff, tgcn2$genes)
  newmodule <- data.frame(hubGene=rep("OTHER", length(diff)), genes=diff)

  if(sum(is.na(index2))>0) {
    tgcn2 <- rbind(tgcn2, newmodule)
  } else if(sum(is.na(index1))>0) {
    tgcn1 <- rbind(tgcn1, newmodule)
  }

  # Get overlap table (code from overlapTable function from WGCNA package adapted to this specific case)
  modules1 <- as.vector(unique(tgcn1$hubGene))
  modules2 <- as.vector(unique(tgcn2$hubGene))

  countMat <- matrix(0, length(modules1), length(modules2))
  pMat <- matrix(0, length(modules1), length(modules2))
  nAll <- length(intersect(background1, background2))

  for(m1 in 1:length(modules1)) {
    for (m2 in 1:length(modules2)) {
      genes1 <- tgcn1$genes[tgcn1$hubGene==modules1[m1]]
      genes2 <- tgcn2$genes[tgcn2$hubGene==modules2[m2]]
      .n1 = length(genes1)
      .n2 = length(genes2)
      .n12 = length(intersect(genes1, genes2))

      if (.n12 > 0) {
        pMat[m1, m2] <- phyper(.n12 - 1, .n1, nAll - .n1, .n2, lower.tail = FALSE, log.p = FALSE)
      } else {
        pMat[m1, m2] <- 1
      }
      countMat[m1, m2] = .n12
    }
  }

  XTbl <- list()
  XTbl$countTable <- countMat
  XTbl$pTable <- pMat

  rownames(XTbl$countTable) <- modules1
  colnames(XTbl$countTable) <- modules2

  rownames(XTbl$pTable) <- modules1
  colnames(XTbl$pTable) <- modules2

  # CrossTabPlot code
  XTbl$pTable[] = p.adjust(XTbl$pTable, method = "fdr")
  toreturn = XTbl$pTable
  XTbl$pTable <- -log10(XTbl$pTable)
  XTbl$pTable[XTbl$pTable > 50] = 50

  ModTotals.1 = table(factor(tgcn1$hubGene, levels=modules1))
  stopifnot(identical(names(ModTotals.1), modules1))

  ModTotals.2 = table(factor(tgcn2$hubGene, levels=modules2))
  stopifnot(identical(names(ModTotals.2), modules2))

  if(identical(modules1, modules2) & diag==F) {
    diag(XTbl$pTable) <- NA
  }

  XTbl$pTable[XTbl$pTable<1.3] <- 0

  if(!is.null(colsOrder)) {
    XTbl$pTable <- XTbl$pTable[, match(colsOrder, colnames(XTbl$pTable))]
    stopifnot(identical(colnames(XTbl$pTable), colsOrder))

    XTbl$countTable <- XTbl$countTable[, match(colsOrder, colnames(XTbl$countTable))]
    stopifnot(identical(colnames(XTbl$countTable), colsOrder))

    ModTotals.2 <- ModTotals.2[match(colsOrder, names(ModTotals.2))]
    stopifnot(identical(names(ModTotals.2), colsOrder))
  }

  if(significant==T) {
    XTbl$sig <- XTbl$pTable
    XTbl$sig[which(XTbl$sig>(-log10(0.001)))] <- "***"
    XTbl$sig[which(XTbl$sig>(-log10(0.01)))] <- "**"
    XTbl$sig[which(XTbl$sig>(-log10(0.05)))] <- "*"
    XTbl$sig[which(XTbl$sig==0)] <- ""
  } else {
    XTbl$sig <- XTbl$countTable
  }

  if(moduleSize==T) {
    ySymbols = paste0(names(ModTotals.1), ": ", as.vector(ModTotals.1))
    xSymbols = paste0(names(ModTotals.1), ": ", as.vector(ModTotals.1))
  } else {
    ySymbols = names(ModTotals.1)
    xSymbols = names(ModTotals.2)
  }

  # par(mar = c(5, 6, 4, 4) + 0.01)

  p <- labeledHeatmap(Matrix = XTbl$pTable,
                      yLabels = ySymbols,
                      xLabels = paste(" ", names(ModTotals.2)),
                      colorLabels = TRUE,
                      xColorWidth=0.03,
                      yColorWidth=0.02,
                      textMatrix = XTbl$sig,
                      colors = blueWhiteRed(100)[50:100],
                      legendLabel = expression(paste(-log[10], P[adj])),
                      naColor = "white",
                      xlab = paste0(name2, " modules"),
                      ylab = paste0(name1, " modules"),
                      cex.text = 0.6, cex.lab = 0.6, cex.main=0.8, setStdMargins = FALSE,
                      plotLegend = TRUE,
                      main=main)

  return(list(p=p, df=toreturn))
}



#' getLASSOinstability - It estimates the Jaccard index between two gene sets selected by LASSO in two different runs.
#'
#' @param lasso_models a list of LASSO models obtained from the getHubGenes or getLASSOmodels functions.
#' @return for each pair of LASSO runs, the Jaccard index between the genes selected. Results are returned both as a matrix
#' and as a data frame.
#' @export
#' @examples

getLASSOinstability <- function(lasso_models) {

  # Get the Jaccard index per pairs of iterations
  genes <- lapply(lasso_models, function(iter) iter$coeffs$genes) # Get genes selected at least once
  iters <- paste0("iter", 1:length(lasso_models))
  M <- data.frame() # to represent the corrplot
  df <- data.frame()

  # Estimate the instability
  for (iter1 in iters) {
    J_all <- c()
    for (iter2 in iters) {
      genesIter1 <- genes[[iter1]]
      genesIter2 <- genes[[iter2]]
      overlap <- length(intersect(genesIter1, genesIter2))
      all <- length(unique(c(genesIter1, genesIter2)))
      J_pair <- (overlap/all)
      J_all <- c(J_all, J_pair)

      df <- rbind(df, c(iter1, iter2, J_pair))
    }
    M <- rbind(M, J_all)
  }

  rownames(M) <- iters
  colnames(M) <- iters

  colnames(df) <- c("iter_n", "iter_m", "Jpair")
  df$Jpair <- as.numeric(as.character(df$Jpair))

  return(list(M=M, df=df))
}



#' getLASSOstats - It returns the main metrics of the bootstrapped LASSO.
#'
#' @param lasso_models a list of LASSO models obtained from the getHubGenes or getLASSOmodels functions.
#' @return it returns a data frame with the main metrics of the bootstrapped LASSO including: the mean number
#' of samples used to train the model (m), the number of genes used to train the model (p), the Jaccard index
#' between the genes selected for each pair of LASSO runs, the number of genes selected per run (nzero), the
#' train and test error as RMSE, the train and test R2 and the number of samples used to train and test the
#' model.
#' @export
#' @examples

getLASSOstats <- function(hubGenes) {

  require(stringr)
  results <- data.frame()

  # Get the short name
  name <- str_split(basename(hubGenes), "_")[[1]][1]

  # Read the file
  file <- readRDS(hubGenes)
  m <- mean(unlist(lapply(file$lasso_models, function(x) length(x$trainSamples))))
  p <- nrow(stats::coef(file$lasso_models$iter1$model, s="lambda.min"))

  results <- rbind(results, data.frame(values=m, metric="m", name=name))
  results <- rbind(results, data.frame(values=p, metric="p", name=name))

  # Get RMSE of the lasso models
  df <- stack(getLassoMetrics(file$lasso_models))
  df$name <- rep(name, nrow(df))
  df$set <- rep(c("Train", "Test"), nrow(df)/2)
  colnames(df)[2] <- "metric"
  df$metric <- paste0(df$metric, "_", df$set)
  df$set <- NULL
  results <- rbind(results, df)

  # Get instability of the lasso models
  inst <- getLASSOinstability(file$lasso_models)$M
  inst[upper.tri(inst, diag=T)] <- NA
  inst <- stack(inst)
  inst <- inst[!is.na(inst$values), ]
  inst$name <- rep(name, nrow(inst))
  colnames(inst)[2] <- "metric"
  inst$metric <- rep("Jaccard", nrow(inst))
  results <- rbind(results, inst)

  results$name <- gsub("\\s*\\([^\\)]+\\)","", results$name)
  results <- results[-which(results$metric %in% "nzero_Test"), ]
  results$metric <- gsub("nzero_Train", "nzero", results$metric)
  results$metric <- factor(results$metric, levels=c("m", "p", "Jaccard", "nzero", "RMSE_Train", "RMSE_Test", "R2_Train", "R2_Test",  "numSamples_Train", "numSamples_Test"))

  return(results)

}



#' getColors - It returns a color for each name of a vector.
#'
#' @param names a vector with the names to be colored
#' @param palette the name of the palette to be used to assign the colors
#' @return a vector with the names and the color assign to each one
#' @export
#' @examples

getColors <- function(names, palette="Set2") {

  require(RColorBrewer)
  mycolors <- rev(colorRampPalette(brewer.pal(5, palette))(length(names)))
  names(mycolors) <- names

  return(mycolors)
}



#' plotLASSOstats - It plots the main statistics of bootstrapped LASSO per tissue (and per group).
#'
#' @param stats the output of getLASSOstats function
#' @param colors a vector with the names of the tissue and the assigned colors to be used
#' in the plot. If colors=NULL, colors will be assigned using getColors function
#' @param order if order=T, tissues will be ordered based on the Jaccard index in descending
#' order
#' @return a plot with the main statistics of bootstrapped LASSO per tissue (and per group).
#' @export
#' @examples

plotLASSOstats <- function(stats, colors=NULL, order=F) {

  require(RColorBrewer)
  require(ggplot2)

  if(order==T) {
    mean <- aggregate(stats$values[stats$metric=="Jaccard"], list(stats$name[stats$metric=="Jaccard"]), FUN="mean")
    mean <- mean[order(mean$x, decreasing=F), ]
    stats$name <- factor(stats$name, levels=mean$Group.1)
  }

  if(is.null(colors)) {
    if(order==T) colors <- getColors(levels(stats$name))
    else colors <- getColors(unique(stats$name))
  }

  if(!is.null(stats$group)) {

    wrap_text <- function(x, chars = 100) {
      x <- gsub("_", " ", x)
      stringr::str_wrap(x, chars)
    }

    p <- ggplot(stats, aes(y=values, x=name, fill=group)) +
                facet_wrap(.~metric, scales="free_y", ncol=5, labeller=as_labeller(wrap_text)) +
                geom_bar(data=subset(stats, metric %in% "m"), stat="identity", position="dodge") +
                geom_boxplot(data = subset(stats, metric %in% c("Jaccard", "nzero", "RMSE_Train", "RMSE_Test")),
                             outlier.size=0.3, size=0.3) +
                theme_classic() +
                xlab("") +
                scale_fill_brewer(palette="Blues") +
                labs(fill="Group") +
                theme(text=element_text(size=14),
                      strip.text.x=element_text(size=15),
                      axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  } else {
    wrap_text <- function(x, chars = 5) {
      x <- gsub("_", " ", x)
      stringr::str_wrap(x, chars)
    }

  p <- ggplot(stats, aes(y=name, x=values, fill=name)) +
              facet_grid(.~metric, scales="free_x", labeller=as_labeller(wrap_text)) +
              geom_violin(data = subset(stats, metric %in% c("Jaccard", "Jaccard_corrected"))) +
              geom_boxplot(data = subset(stats, metric %in% c("nzero", "RMSE_Train", "RMSE_Test"))) +
              geom_bar(data=subset(stats, metric %in% c("m", "p")),stat="identity") +
              theme_classic() +
              theme(legend.position="none",
                    text=element_text(size=14),
                    panel.spacing = unit(0.7, "lines"),
                    strip.text.x = element_text(size=13),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=13)) +
              ylab("") +
              xlab("values") +
              scale_fill_manual(values=colors)

    }

  return(list(p=p, colors=colors))

}



#' selectNumberOfLASSOruns - It plots the new hub genes obtained after more than one run.
#'
#' @param hubs the output of getHubGenes function
#' @return a plot showing the number of new hubs non-selected in the first run but selected in successive LASSO
#' runs, grouped by ratio of appearance.
#' @export
#' @examples

selectNumberOfLASSOruns <- function(hubs) {

  if(is.character(hubs)) {
    hubs <- readRDS(hubs)
  }

  nzeroDF <- data.frame()
  nzeroList <- lapply(hubs$LASSO$models, function(x) x$coeffs$genes)
  initialPool <- unique(unlist(nzeroList[[1]]))
  max <- length(hubs$LASSO$models)

  for (iter in seq(1, max, 1)) {
    iterPool <- unique(unlist(nzeroList[c(1:iter)]))
    newPool <- setdiff(iterPool, initialPool)

    newPoolTab <- table(unlist(nzeroList[c(1:iter)]))
    newPoolTab <- newPoolTab[match(newPool, names(newPoolTab))]

    for(timesSelected in seq(1, max, 1)) {
      numOfGenes <- length(which(newPoolTab>=timesSelected))
      nzeroDF <- rbind(nzeroDF, c(iter, timesSelected, numOfGenes))
    }
  }

  colnames(nzeroDF) <- c("iter", "timesSelected", "numOfGenes")
  nzeroDF$numOfGenes <- as.numeric(as.character(nzeroDF$numOfGenes))

  nzeroDF <- nzeroDF[nzeroDF$timesSelected %in% c(1,3,5,7), ]
  nzeroDF$timesSelected <- as.factor(nzeroDF$timesSelected)
  mycolors <- getColors(unique(nzeroDF$timesSelected), palette="Spectral")
  names(mycolors) <- rev(names(mycolors))

  p <- ggplot(data=nzeroDF, aes(x=iter, y=numOfGenes, color=timesSelected)) +
              geom_line(size=1) +
              geom_point(size=2) +
              theme_classic() +
              ylab("Number of new relevant genes") +
              xlab("iteration") +
              scale_x_continuous(breaks=seq(1,max,1)) +
              scale_colour_manual(values = mycolors) +
              theme(text=element_text(size=12))

  return(p)

}



#' getLASSOcoeffs - For each gene selected at least once in all LASSO runs, it returns the corresponding
#' normalized coefficient for each run.
#'
#' @param hubs the output of getHubGenes function
#' @return it returns a data frame with as many columns as LASSO runs and as many rows as genes selected at
#' least once in any of these runs. Each cell represent the normalized coefficient for one gene in one
#' LASSO run.
#' @export
#' @examples

getLASSOcoeffs <- function(hubs) {

  # Load library
  require(caret)

  # Load files
  if(is.character(hubs)) {
    hubs <- readRDS(hubs)
  }

  # Get genes selected at least once
  genes <- unique(as.vector(unlist(lapply(hubs$LASSO$models, function(iter) iter$coeffs$genes))))

  # Create a data frame with the coefficient value for each gene in each iteration
  df <- data.frame()

  for(gene in genes) {
    coeffs <- lapply(hubs$LASSO$models, function(iter) iter$coeffs$coeff[match(gene, iter$coeffs$genes)])
    df <- rbind(df, coeffs)
  }

  rownames(df) <- genes

  # Get coefficients in absolute value
  df <- abs(df)

  # Normalize the coefficients between 0 and 1
  preProc <- preProcess(df, method = c("range"))
  dfn <- predict(preProc, df)

  return(dfn)

}



#' getLASSOmeanCoeff - For each gene selected at least once in all LASSO runs, it returns its mean coefficient
#' in absolute value.
#' @param hubs the output of getHubGenes function
#' @return For each gene selected at least once in all LASSO runs, it returns its mean coefficient in
#' absolute value.
#' @export
#' @examples

getLASSOmeanCoeffs <- function(hubs) {

  # Get coeffs per iter
  dfn <- getLASSOcoeffs(hubs)

  # Get the number of iterations where each gene has been selected as relevant
  numIters <- ncol(dfn)-rowSums(is.na(dfn))

  # Get the mean coefficient for each gene
  meanCoeffs <- rowMeans(dfn, na.rm=T)

  # Create a final data frame with this information
  finaldf <- as.data.frame(cbind(rownames(dfn), numIters, round(meanCoeffs,5)))
  rownames(finaldf) <- c()
  colnames(finaldf) <- c("genes", "numIters", "meanCoeff")
  finaldf$meanCoeff <- as.numeric(as.character(finaldf$meanCoeff))
  finaldf$numIters <- as.numeric(finaldf$numIters)
  finaldf <- finaldf[order(finaldf$numIters, finaldf$meanCoeff), ]
  finaldf$genes <- factor(finaldf$genes, levels=c(finaldf$genes))

  return(finaldf)
}



#' plotLMmodel - It plots the regression between the mean coefficient of a gene and the number of LASSO runs
#' where it was selected.
#' @param hubs the output of getHubGenes function
#' @return a list containing both 1) for each selected gene, the mean coefficient and the number of LASSO runs
#' where it was selected; 2) the regression plot between the mean coefficient of a gene and the number of LASSO
#' runs where it was selected.
#' @export
#' @examples

plotLMmodel <- function(hubs) {

  # Load library
  require(ggplot2)

  # Get mean coefficient
  df <- getLASSOmeanCoeffs(hubs)

  # Create a linear regression model
  colnames(df) <- c("genes", "predictor", "output")
  model <- lm(output~predictor, data=df)

  # Plot the results
  lm_eqn <- function(m){
    eq <- substitute(italic(y) == a + b * italic(x)*","~~italic(r)^2~"="~r2,
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
  }

  p <- ggplot(df, aes(x=predictor, y=output)) +
              geom_point(aes(color=output), size=2) +
              geom_smooth(aes(as.numeric(predictor), output), method = "lm")+
              theme_classic() +
              theme(text=element_text(size=15)) +
              scale_colour_gradient(low = "lightblue", high = "red") +
              geom_text(x = 5, y = max(df$output)-0.1, label = lm_eqn(model), parse = TRUE)

  return(list(data=df, plot=p))
}



#' plotGOenrichStats - It plots the GO enrichment for each tissue/condition.
#' @param hubs the name of the files containing the hub genes
#' @return a plot that represents the GO enrichment per file name as the sum(-log10Pval)/ module size
#' @export
#' @examples

plotGOenrichStats <- function(hubs) {

  require(ggplot2)
  require(stringr)

  superDF <- data.frame()

  for(h in files) {
    name <- str_split(basename(h), "_")[[1]][1]
    h <- readRDS(h)
    df <- h$lasso_hubGenes_enrichment_stats
    df$query <- gsub("cutoff", "", df$query)
    df$query <- factor(df$query, levels=seq(1,10,1))
    df$name <- as.factor(rep(name, nrow(df)))
    superDF <- rbind(superDF, df)
  }

  p <- ggplot(data=superDF, aes(x=query, y=sum_corrected, group=name, color=name)) +
              geom_point() +
              geom_line() +
              theme_classic() +
              xlab("cutoff") +
              ylab(expression(paste("sum (-", log[10] ~ "Pval)/", "\n module size"))) +
              theme(text=element_text(size=14),
                    legend.key.size = unit(0.5, 'cm')) +
              scale_y_continuous(limits=c(0,2), breaks=seq(0,2, 0.5)) +
              geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")


  return(p)
}



#' getCTenrich - It returns the cell type markers enrichment per TGCN module.
#' @param net a net file as an output of testAllCutoffs function
#' @return a list containing both 1) a matrix with cell types as rows and modules as columns and each cell
#' contains corrected p-values; 2) a data frame that only contains the significant tests;
#' @export
#' @examples

getCTenrich <- function(net) {

  require(CoExpNets)

  if(!is.null(net$hubGene)) {
    mc <- net$hubGene
    names(mc) <- net$genes

    mynet <- list()
    mynet$moduleColors <- mc
  } else {
    mynet <- net
  }

  ct <- as.data.frame(genAnnotationCellType(net.in=mynet,return.processed=F,which.one="new"))
  sigTests <- ct[,apply(ct,2,function(x){ any(x < 0.05)}),drop=FALSE]
  df <- data.frame()

  if(ncol(sigTests)>=1) {

    for(module in 1:ncol(sigTests)) {
      test <- sigTests[, module]
      index <- as.vector(which(test<0.05))

      for(i in index) {
        newrow <- c(rownames(ct)[i], format(as.vector(test[i]), sci=T, digits=4), colnames(sigTests)[module])
        names(newrow) <- c("ct", "pVal", "module")
        df <- rbind(df, newrow)
      }
    }
    colnames(df) <- c("ct", "correctedPvalue", "module")

    return(list(matrix=ct, df=df))
  } else {
    return()
  }

}



#' plotCTenrich - It plots the cell type markers enrichment per TGCN module.
#' @param ct the output of getCTenrich function
#' @param target the name of the target
#' @param height the height of the plot
#' @param width the width of the plot
#' @return a plot that shows the significant associations between cell type markers and TGCN modules.
#' @export
#' @examples

plotCTenrich <- function(ct, target="target", height=10, width=14) {

  require(circlize)
  require(complexHeatmap)

  col_fun = colorRamp2(c(0, 5, 10), c("white", "orange", "red"))

  if(!is.null(ct$matrix)) {
    matrix <- ct$matrix
  } else {
    matrix <- ct
  }

  matrix <- apply(matrix, 1, function(x) -log10(x))
  modulesNames <- names(which(rowSums(matrix)>0))
  matrix <- matrix[which(rowSums(matrix)>0), ]

  if(!is.null(dim(matrix))) {
    CTnames <- names(which(colSums(matrix)>0))
    matrix <- as.data.frame(matrix[, which(colSums(matrix)>0)])
    matrix <- t(matrix)
  } else {
    CTnames <- names(matrix)[which(matrix>0)]
    matrix <- matrix[which(matrix>0)]
    matrix <- as.data.frame(matrix)
  }
  rownames(matrix) <- CTnames
  rownames(matrix) <- gsub("_", " ", rownames(matrix))
  colnames(matrix) <- modulesNames

  if(nrow(matrix)>2) {
    ht <- Heatmap(matrix,
                  name = "-log10P",
                  column_title = target,
                  column_title_gp=grid::gpar(fontsize = 11),
                  row_title = "Cell type markers",
                  col = col_fun,
                  na_col = "black",
                  column_names_rot = 45,
                  border = TRUE,
                  row_gap = unit(1.5, "mm"),
                  column_gap = unit(1.5, "mm"),
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  column_names_gp = grid::gpar(fontsize = 8),
                  row_names_gp = grid::gpar(fontsize = 8),
                  heatmap_height = unit(height, "cm"),
                  heatmap_width = unit(width, "cm"),
                  cluster_rows = T,
                  cluster_columns = F,
                  show_row_dend = FALSE,
                  heatmap_legend_param = list(
                    legend_direction = "horizontal",
                    legend_width = unit(3, "cm")),
                  row_names_max_width = max_text_width(rownames(matrix),
                                                       gp = gpar(fontsize = 10)))

    draw(ht, heatmap_legend_side="bottom")
  } else {
    ht <- NULL
  }

  return(ht)
}



#' plotModuleTraitCorr - It plots the associations between traits and modules eigengenes.
#' @param MEs a data frame containing modules eigengenes with modules as columns and samples as rows
#' @param covs a data frame containing the donors traits
#' @param cex.lab.x the size of the x axis labels
#' @param cex.lab.y the size of the y axis labels
#' @param height the height of the plots
#' @param width the width of the plots
#' @param main the title of the plot
#' @param method for numeric covariates, the type of correlation (Pearson or Spearman)
#' @param max the -log10P over max are represented as max.
#' @param ylab the label of y axis
#' @return a plot that shows the associations between traits and modules eigengenes.
#' @export
#' @examples

plotModuleTraitCorr <- function(MEs, covs=NULL, cex.lab.x=0.9, cex.lab.y=0.7, height=10, width=8,
                           main="Module-trait relationships", method="spearman", max=10, ylab="Modules"){
  require(ComplexHeatmap)
  require(WGCNA)
  require(circlize)

  for(i in 1:ncol(covs)){
    if(typeof(covs[,i]) ==  "character")
      covs[,i] = as.factor(covs[,i])
  }
  factor.mask = unlist(lapply(covs,is.factor))
  cat("We will work with",sum(factor.mask),"factors\n")
  # stopifnot(sum(factor.mask) > 0)

  fcm = matrix(nrow=ncol(MEs),ncol=sum(factor.mask))
  index = 1
  for(i in which(factor.mask)){
    #cat("Factor",colnames(trait.data)[i],"\n")
    #cat("Levels",levels(trait.data[,i]))
    #print(trait.data[,i])
    for(j in 1:ncol(MEs)){
      #print(paste0(i,j)
      if(length(unique(covs[,i])) > 1){
        form = eg ~ cov
        data.in = data.frame(MEs[,j],covs[,i])
        colnames(data.in) = c("eg","cov")
        fcm[j,index] = anova(aov(form,data.in))$`Pr(>F)`[1]
      }else
        fcm[j,index] = 1

    }
    fcm[,index] = p.adjust(fcm[,index],method="BH")
    index = index + 1
  }

  if(sum(!factor.mask) > 0){
    moduleTraitCor = cor(MEs,covs[,!factor.mask,drop=FALSE],use="p", method=method)
    #Generate the p-values for significance of a given matrix of correlations, for all modules,
    #between traits data and eigengenes, both from samples
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor,nrow(MEs))
    moduleTraitPvalue = cbind(moduleTraitPvalue,fcm)
    colnames(moduleTraitPvalue) = c(colnames(covs)[!factor.mask],
                                    colnames(covs)[factor.mask])
  }else{
    moduleTraitPvalue = fcm
    colnames(moduleTraitPvalue) = colnames(covs)[factor.mask]
    rownames(moduleTraitPvalue) = gsub("ME","",names(MEs))
    moduleTraitCor <- NULL
  }

  pval <- moduleTraitPvalue
  cor <- moduleTraitCor

  moduleTraitPvalue = -log10(moduleTraitPvalue)
  moduleTraitPvalue[moduleTraitPvalue > 10] = 10
  moduleTraitPvalue[moduleTraitPvalue < -log10(0.05)] <- 0

  col_fun = colorRamp2(c(0, max/2, max), c("white", "orange", "red"))
  rownames(moduleTraitPvalue) <- gsub("^ME", "", names(MEs))

  htPval <- Heatmap(moduleTraitPvalue,
                    name = "-log10P",
                    column_title = "Traits",
                    row_title = ylab,
                    col = col_fun,
                    na_col = "black",
                    column_names_rot = 45,
                    border = TRUE,
                    row_gap = unit(1.5, "mm"),
                    column_gap = unit(1.5, "mm"),
                    clustering_method_rows = "ward.D2",
                    clustering_method_columns = "ward.D2",
                    column_names_gp = grid::gpar(fontsize = 8),
                    row_names_gp = grid::gpar(fontsize = 7),
                    heatmap_height = unit(height, "cm"),
                    heatmap_width = unit(width, "cm"),
                    cluster_rows = T,
                    cluster_columns = F,
                    heatmap_legend_param = list(
                      legend_direction = "horizontal",
                      legend_width = unit(2, "cm")),
  row_names_max_width = max_text_width(rownames(matrix),
                                       gp = gpar(fontsize = 9)))

  draw(htPval, heatmap_legend_side="bottom")

  if(!is.null(moduleTraitCor)) {

    col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    rownames(moduleTraitCor) <- gsub("ME","",rownames(moduleTraitCor))

    htCor <- Heatmap(moduleTraitCor,
                     name = "corr",
                     column_title = "Traits",
                     row_title = ylab,
                     col = col_fun,
                     na_col = "black",
                     column_names_rot = 45,
                     border = TRUE,
                     row_gap = unit(1.5, "mm"),
                     column_gap = unit(1.5, "mm"),
                     clustering_method_rows = "ward.D2",
                     clustering_method_columns = "ward.D2",
                     column_names_gp = grid::gpar(fontsize = 8),
                     row_names_gp = grid::gpar(fontsize = 7),
                     heatmap_height = unit(height, "cm"),
                     heatmap_width = unit(width, "cm"),
                     cluster_rows = T,
                     cluster_columns = F,
                     heatmap_legend_param = list(
                       legend_direction = "horizontal",
                       legend_width = unit(4, "cm")),
                     row_names_max_width = max_text_width(rownames(matrix),
                                                          gp = gpar(fontsize = 8)))

    draw(htCor, heatmap_legend_side="bottom")


  } else {
    htCor <- NULL
  }


  return(list(pval=pval, cor=cor, htPval=htPval, htCor=htCor))
}



#' plotGOenrichTerms - It plots the GO terms of a specific source per TGCN module.
#' @param enrich a data frame containing GO enrichment results from gost function
#' @param source the name of the source to plot the corresponding terms
#' @param name the name of the TGCN
#' @param height the height of the plots
#' @param width the width of the plots
#' @return a plot that shows the GO terms for each TGCN module
#' @export
#' @examples

plotGOenrichTerms <- function(enrich, source="GO:CC", name="APP", width=18, height=29) {

  # Load library
  require(ComplexHeatmap)

  # Get enrichment as a matrix
  enrich <- enrich[enrich$source==source, ]

  df <- data.frame()
  terms <- unique(enrich$term_name[enrich$source==source])
  modules <- unique(enrich$query)

  for (term in terms) {
    row <- rep(NA, length(modules))
    names(row) <- modules
    myenrich <- enrich[enrich$term_name==term, ]
    row[match(myenrich$query, modules)] <- myenrich$p_value
    df <- rbind(df, row)
  }

  colnames(df) <- modules
  rownames(df) <- terms
  df <- as.matrix(apply(df, 2, function(x) -log10(x)))

  if(!is.null(df)) {
    df[is.na(df)] <- 0 # Replace NA values by zeros
    dim(df)

    # How to filter the rows and columns
    if(nrow(df)<=20) {
      min_sum_row=1
    } else if (nrow(df)>20 & nrow(df)<=35) {
      min_sum_row=2
    } else if (nrow(df)>35) {
      min_sum_row=5
    }

    if(ncol(df)<=5) {
      min_sum_column=1
    } else if (ncol(df)>5 & ncol(df)<=10) {
      min_sum_column=2
    } else if (ncol(df)>10) {
      min_sum_column=4
    }

    # Filter rows
    termsToRemove <- which(rowSums(df)<min_sum_row)
    if(length(termsToRemove)>0) {
      df <- df[-termsToRemove, ]
      dim(df)
    }

    # Filter columns
    modulesToRemove <- which(colSums(df)<min_sum_column)
    if(length(modulesToRemove)>0) {
      df <- df[ ,-modulesToRemove]
      dim(df)
    }

    # Define colors
    col_fun = colorRamp2(c(0, 5, 10), c("white", "orange", "red"))

    if(ncol(df)<=2) {
      column_km=1
    } else if (ncol(df)>2 & ncol(df)<=10) {
      column_km=2
    } else if (ncol(df)>10) {
      column_km=3
    }

    # Create the heatmap
    ht <- Heatmap(df,
                  name = "-log10P",
                  column_title = paste0(name, " modules"),
                  row_title = paste0(source, " terms"),
                  col = col_fun,
                  na_col = "black",
                  column_names_rot = 45,
                  row_km = 4,
                  column_km = column_km,
                  border = TRUE,
                  row_gap = unit(1.5, "mm"),
                  column_gap = unit(1.5, "mm"),
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2",
                  column_names_gp = grid::gpar(fontsize = 10),
                  row_names_gp = grid::gpar(fontsize = 10),
                  heatmap_height = unit(height-2, "cm"),
                  heatmap_width = unit(width-0.5, "cm"),
                  cluster_rows = T,
                  cluster_columns = T,
                  heatmap_legend_param = list(
                    legend_direction = "horizontal",
                    legend_width = unit(4, "cm")),
                  row_names_max_width = max_text_width(
                    rownames(df),
                    gp = gpar(fontsize = 12)))
    draw(ht, heatmap_legend_side="bottom")

  } else {
    ht <- NULL
  }

  return(ht)
}



#' plotGOenrichSummary - It plots the GO terms of a specific source per TGCN module.
#' @param enrich a data frame containing GO enrichment results from gost function
#' @return a list containing both 1) a plot that shows the sum(-log10P)/module size per module;
#' 2) a plot that shows the number of terms of each source per module.
#' @export
#' @examples

plotGOenrichSummary <- function(enrich) {

  stats <- enrich$sumLogPval
  stats <- stats[order(stats$sum_corrected, decreasing=T), ]
  stats$query <- factor(stats$query, levels=unique(stats$query))

  stats$group <- rep("module", nrow(stats))
  stats$group[stats$query=="hubs"] <- "hubs"

  p1 <- ggplot(data=stats, aes(x=query, y=sum_corrected, fill=group)) +
               geom_bar(stat="identity") +
               theme_classic() +
               theme(text=element_text(size=12),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
               xlab("Gene sets") +
               ylab("sum(-log10P) /\n module size") +
               labs(title="Enrichment per module") +
               geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red")

    terms <- enrich$enrichment
    terms <- as.data.frame(table(factor(terms$query, levels=levels(terms$query)), factor(terms$source)))
    colnames(terms) <- c("module", "source", "nTerms")
    levels <- aggregate(terms$nTerms, list(terms$module), FUN="sum")
    levels <- levels[order(levels$x, decreasing=T), ]
    terms$module <- factor(terms$module, levels=levels$Group.1)


  p2 <- ggplot(data=terms, aes(x=module, y=nTerms, fill=source)) +
               geom_bar(stat="identity") +
               theme_classic() +
               theme(text=element_text(size=12),
                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
               xlab("Modules") +
               ylab("Number of terms") +
               labs(title="Number of terms\n per module")

  return(list(stats=p1, nterms=p2))

}


