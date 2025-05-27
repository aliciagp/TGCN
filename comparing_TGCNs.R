#' selectBasicFeatures - It represents, given the transcriptomic and proteomic data, the dimensionality ratio: molecules / samples, and the Venn diagram showing the percentage of shared molecules between both datasets.
#'
#' @param transcriptomicData the dataset used to create transcriptomic TGCNs.
#' @param proteomicData the dataset used to create proteomic TGCNs.
#' @param targetName name of the target.
#' @param tissueName name of the dataset or cohort.
#' @param cr_path path where the results will be saved.
#' @param overwrite
#' @return a list containing 1) the plot of the number of ratio dimensionality; 2) the plot with the Venn diagram
#' @export
#' @examples

selectBasicFeatures <- function(transcriptomic_data,
                                proteomic_data,
                                targetName,
                                tissueName,
                                cr_path,
                                overwrite=F) {
  
  # Load libraries
  require(ggVennDiagram) 
  require(ggplot2)
  
  # file names for the resulting plots:
  dimratio_file <- paste0(cr_path, "/EarlyReport/", targetName, "_", tissueName, "_dimRatio.png")
  venn_file <- paste0(cr_path, "/EarlyReport/", targetName, "_", tissueName, "_sharedMolecules.png")
  
  
  # Create the directory if it does not exist:
  if(!dir.exists(paste0(cr_path, "/EarlyReport/"))) {
    dir.create(paste0(cr_path, "/EarlyReport/"))
  } 
  
  if(!file.exists(dimratio_file) || !file.exists(venn_file) || overwrite==T) {
    # compute dimensionality ratio 
    dim_proteomic <- dim(proteomic_data)  
    dim_transcript <- dim(transcriptomic_data)  
    
    ratio_proteomic <- dim_proteomic[1] / dim_proteomic[2]  # Proteins / Samples
    ratio_transcript <- dim_transcript[1] / dim_transcript[2]  # Genes / Samples
    
    # df for the first plot
    df_ratio <- data.frame(
      Dataset = c("Proteomic", "Transcriptomic"),
      Ratio = c(ratio_proteomic, ratio_transcript)
    )
    
    # dimensionality ratio plot
    p1 <- ggplot(df_ratio, aes(x = Dataset, y = Ratio, fill = Dataset, alpha=0.5)) +
      scale_y_continuous(expand = expansion(mult = 0.2)) +
      geom_bar(stat = "identity", width = 0.6) +  
      geom_text(aes(label = round(Ratio, 2)), vjust = -0.5, size = 6, color = "black") +  
      scale_fill_manual(values = c("Proteomic" = "red", "Transcriptomic" = "blue")) +
      theme_minimal() + 
      ggtitle("Dimensionality ratio comparison") +
      labs(y = "Ratio") + 
      theme(
        axis.title.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray", linetype = "dashed"), 
        panel.grid.major.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none"
      )
    if (!file.exists(dimratio_file) || overwrite==T) {
      ggsave(dimratio_file, bg="White", plot = p1, width = 8, height = 8, dpi=300)
    }  
    
    # get genes and proteins names (uniques)
    genes <- unique(rownames(transcriptomic_data))  
    proteins <- unique(rownames(proteomic_data)) 
    
    # Intersect and get the percentage
    shared_molecules <- intersect(genes, proteins)
    percentage_shared <- (length(shared_molecules) / min(length(genes), length(proteins))) * 100
    
    # Crear lista de conjuntos
    venn_data <- list(
      Transcriptomic = genes,
      Proteomic = proteins
    )
    
    # Crear el diagrama de Venn con ggplot2
    p2 <- ggVennDiagram(venn_data, label_alpha = 0, set_size = 6) + 
      scale_fill_distiller(palette = "RdBu") +
      guides(fill = guide_legend(title = "Number of molecules")) +
      labs(title = "Shared molecules") +
      theme(
        legend.position = "bottom"
      ) +
      scale_x_continuous(expand = expansion(mult = 0.5)) +  
      scale_y_continuous(expand = expansion(mult = 0.2))
    
    
    # Save the plot if necessary
    if (!file.exists(venn_file) || overwrite==T) {
      ggsave(venn_file, bg="White", plot = p2, width = 8, height = 8, dpi=300)
    }
    # Return the list with both plots.
    return(list(dimratio=p1, venn=p2))
  }
}

#'
#' selectSeedMolecules - It plots the R^2 for both transcriptomic and proteomic
#' TGCNs to compare them over training and test sets.
#'
#' @param: transcriptomic_path: location for the transcriptomic TGCNs output.
#' @param: proteomic_path: location for the proteomic TGCNs output.
#' @param: cr_path: location where the report will be saved.
#' @param: targetName: the name of the targeted variable.
#' @param: tissueName: the name of the tissue, cohort or dataset.
#' @param: overwrite: If overwrite=T, it overwrites the previous results.
#' @param: transcriptomic_hubs: TGCn hub generated using transcriptomic data.
#' @param: proteomic_hubs: TGCn hub generated using proteomic data.
#' @param: selectedCutoffs: Array with the selected cutoff for at least one of the TGCNs.
#' @return: a list contained all the generated plots for each cutoff.
#'
selectSeedMolecules <- function(transcriptomic_path,
                                proteomic_path,
                                cr_path,
                                targetName,
                                tissueName,
                                overwrite=F,
                                selectedCutoffs) {
  
  # load libraries
  require(ggplot2)
  require(fs)
  
  
  # create the directory if it does not exist:
  output_dir <- file.path(cr_path, paste0(targetName, "_", tissueName, "_mTGCN_seedsPerformance"))
  dir_create(output_dir)
  
  # path for the output plot
  output_file <- file.path(output_dir, paste0(targetName, "_", tissueName, "_", "_mTGCN_seedR2.png"))
  
  transcriptomic_hubs <- readRDS(file=paste0(transcriptomic_path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"))
  
  proteomic_hubs <- readRDS(file = paste0(proteomic_path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"))
  
  # if not tested yet, the function generates the plots
  if (!file.exists(output_file) || overwrite==T) {
    # dataframe for storing the data
    results_df <- data.frame()
    
    for (cutoff in selectedCutoffs) {
      cutoffn <- paste0("cutoff", cutoff)
      
      if (!is.null(transcriptomic_hubs$LM[[cutoffn]]) & !is.null(proteomic_hubs$LM[[cutoffn]])) {
        
        # Is Accuracy or R2?
        if ("Accuracy" %in% names(transcriptomic_hubs$LM[[cutoffn]]$lm_models$metrics$results) ||
            "Accuracy" %in% names(proteomic_hubs$LM[[cutoffn]]$lm_models$metrics$results)) {
          
          metric <- "Accuracy"
          
        } else {
          
          metric <- "R2"
          
        }
        
        # select the results for train and test
        transcriptomic_test  <- transcriptomic_hubs$LM[[cutoffn]]$lm_models$metrics$results[[metric]][2]
        transcriptomic_train <- transcriptomic_hubs$LM[[cutoffn]]$lm_models$metrics$results[[metric]][1]
        proteomic_test       <- proteomic_hubs$LM[[cutoffn]]$lm_models$metrics$results[[metric]][2]
        proteomic_train      <- proteomic_hubs$LM[[cutoffn]]$lm_models$metrics$results[[metric]][1]
        
        # Aggregate
        results_df <- rbind(results_df,
                            data.frame(Cutoff = cutoff, Dataset = "Transcriptomic", Set = "Train", Metric = transcriptomic_train),
                            data.frame(Cutoff = cutoff, Dataset = "Transcriptomic", Set = "Test", Metric = transcriptomic_test),
                            data.frame(Cutoff = cutoff, Dataset = "Proteomic", Set = "Train", Metric = proteomic_train),
                            data.frame(Cutoff = cutoff, Dataset = "Proteomic", Set = "Test", Metric = proteomic_test))
      }
    }
    
    colors <- c("Transcriptomic" = "blue", "Proteomic" = "red")
    linetypes <- c("Train" = "dashed", "Test" = "solid")
    
    # Label de la métrica
    metric_label <- ifelse(metric == "Accuracy", "Accuracy", expression(R^2))
    print(tissueName)
    # Plot
    p <- ggplot(results_df, aes(x = Cutoff, y = Metric, color = Dataset, linetype = Set, group = interaction(Dataset, Set))) +
      geom_line(size = 1) +         
      geom_point(size = 3) +
      scale_color_manual(values = colors, name = "Dataset") +
      scale_linetype_manual(values = linetypes, name = "Set") +
      theme_minimal() +
      labs(title="Seed Molecules Performance", x = "Cutoff", y = metric_label) +
      scale_x_continuous(breaks = selectedCutoffs) +
      theme(
        legend.position = "top",
        text = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 18, margin = margin(r = 15)),  
        legend.key.width = unit(1.5, "cm"),  
        legend.key.height = unit(1, "cm"),   
        
      ) +
      guides(
        color = guide_legend(override.aes = list(shape = NA, size = 1.5)),
        linetype = guide_legend(override.aes = list(size = 1.5, alpha = 1, linetype = c("solid","dashed")))
      )
    
    ggsave(output_file, bg="White",  p, width = 8, height = 8, dpi = 300)
    return(p)
    
  }
}

#' mgetModulesOverlap-It creates the plots representing the modules overlap for the multi-omic data
#' in each selected cutoff.
#'
#' @param: cr_path: path where the comparison results are saved.
#' @param: transcriptomic_hubs: TGCn hubs generated using transcriptomic data.
#' @param: proteomic_hubs: TGCn hubs generated using proteomic data.
#' @param: targetName: name of the target or cohort.
#' @param: tissueName: name of the tissue.
#' @param: transcriptomicData: expData of the transcriptomic dataset.
#' @param: proteomicData: expData of the proteomic dataset.
#' @param: n: 
#' @param: m:
#' @param: s:
#' @param: minCor:
#' @param: maxTol:
#' @param: approach:
#' @param: cutoff: selected cutoff for generating the plot.
#' @param: overwrite:
#' @return: It return the plot with the modules overlap for the selected cutoff.
mgetModulesOverlap <- function(main_path, 
                               targetName, 
                               tissueName, 
                               transcriptomicData, 
                               proteomicData,
                               cutoff){
  require(ggplot2)
  require(fs)
  
  # Create output directory if it doesn't exist
  output_dir <- file.path(main_path, paste0("comparison_report/", targetName, "_", tissueName, "_mTGCN_modulesOverlap"))
  if(!dir.exists(output_dir)) {
    dir_create(output_dir)
  }
  
  # Define output file name
  plot_filename <- file.path(output_dir, paste0(targetName, "_", tissueName, "_", cutoff, "_TGCN_crossTabPlot.png"))
  
  # Load networks
  tnet <- readRDS(file.path(main_path, "transcriptomics/Net", paste0(targetName, "_", tissueName, "_TGCNs.rds")))
  pnet <- readRDS(file.path(main_path, "proteomics/Net", paste0(targetName, "_", tissueName, "_TGCNs.rds")))
  
  # Extract modules for the specified cutoff
  tnetc <- tnet$nets[[paste0("c", cutoff)]]$net$modules
  pnetc <- pnet$nets[[paste0("c", cutoff)]]$net$modules
  
  # Get background genes
  tback <- rownames(transcriptomicData)
  pback <- rownames(proteomicData)
  
  # Set up plotting parameters
  par(mgp=c(5,0,0))
  par(mar = c(6, 6, 2.5, 1) + 0.5)
  
  # Call plotModulesOverlap function
  overlap_result <- plotModulesOverlap(
    name1 = "Transcriptomics",
    name2 = "Proteomics",
    tgcn1 = tnetc,
    tgcn2 = pnetc,
    background1 = tback,
    background2 = pback,
    significant = FALSE,
    moduleSize = TRUE,
    main = paste("Module Overlap (cutoff =", cutoff, ")")
  )
  
  # Get the raw p-values matrix before FDR adjustment
  # This is what comes back from overlap_result$df
  pvalues_matrix <- overlap_result$df
  
  # Find significant p-values (p < 0.05)
  significant_pvalues <- which(pvalues_matrix < 0.05, arr.ind = TRUE)
  
  # Create a dataframe with significant results
  if(length(significant_pvalues) > 0) {
    significant_results <- data.frame(
      Transcriptomic_Module = rownames(pvalues_matrix)[significant_pvalues[, 1]],
      Proteomic_Module = colnames(pvalues_matrix)[significant_pvalues[, 2]],
      p_value = pvalues_matrix[significant_pvalues],
      p_value_adjusted = p.adjust(pvalues_matrix[significant_pvalues], method = "fdr"),
      minuslog10_pvalue = -log10(pvalues_matrix[significant_pvalues])
    )
    # Sort by p-value (most significant first)
    significant_results <- significant_results[order(significant_results$p_value), ]

  } else {
    cat("No significant module overlaps found (p < 0.05).\n")
  }
  
  # Save the plot
  png(
    filename = plot_filename,
    width = 8,
    height = 6.5,
    units = "in",
    res = 150
  )
  # recall the function to save the plot.
  plotModulesOverlap(
    name1 = "Transcriptomics",
    name2 = "Proteomics",
    tgcn1 = tnetc,
    tgcn2 = pnetc,
    background1 = tback,
    background2 = pback,
    significant = FALSE,
    moduleSize = TRUE,
    main = paste("Module Overlap (cutoff =", cutoff, ")")
    )
  dev.off()
  
  # Return both the plot and significant p-values
  return(list(
    plot = overlap_result$p,
    pvalues = pvalues_matrix,
    significant_results = if(exists("significant_results")) significant_results else NULL
  ))
}



#' plotEnrichmentDistribution - It plots enrichment distribution for each TGCN using 
#' boxplots. Each point is a module with its enrichment.
#' @param: cr_path: location where the report will be saved.
#' @param: transcriptomic_path: location for the transcriptomic TGCNs output.
#' @param: transcriptomic_net: TGCN net generated using transcriptomic data.
#' @param: proteomic_path: location for the proteomic TGCNs output.
#' @param: proteomic_net: TGCN net generated using proteomic data.
#' @param: targetName: the name of the targeted variable.
#' @param: tissueName: the name of the tissue, cohort or dataset.
#' @param: overwrite: If overwrite=T, it overwrites the previous results.
#' @param: cutoff: selected cutoff for generating the plot.
#' @return: it returns the plot with both boxplots (transcriptomic and proteomic).
#'
plotEnrichmentDistribution <- function (cr_path, 
                                        transcriptomic_path,
                                        proteomic_path,
                                        targetName,
                                        tissueName,
                                        overwrite=F,
                                        cutoff) {
  
  library(ggplot2)
  library(dplyr)
  
  # create output directory if it does not exist
  output_dir <- file.path(cr_path, paste0(targetName, "_", tissueName, "_mTGCN_modulesEnrichment"))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  cutoffn <- paste0("c", cutoff)
  output_file <- file.path(output_dir, paste0(targetName, "_", tissueName, "_", cutoffn, "_mTGCN_ModulesEnrichmentBoxplot.png"))
  
  transcriptomic_net <- readRDS(file = paste0(transcriptomic_path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds")) 
  
  proteomic_net <- readRDS(file = paste0(proteomic_path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds")) 
  
  # if not tested or overwrite==T
  if (!file.exists(output_file) || overwrite==T) {
    # select the stats for the transcriptomic data and include the source 
    # new column that specifies the source of each row, transcriptomic or proteomic.
    if (!is.null(transcriptomic_net$nets[[cutoffn]])) {
      tdata <- transcriptomic_net$nets[[cutoffn]]$GOenrich$stats %>% 
        mutate(Source = "Transcriptomic")
      pdata <- proteomic_net$nets[[cutoffn]]$GOenrich$stats %>% 
        mutate(Source = "Proteomic")
      
      # join the data in one dataframe    
      all_data <- bind_rows(tdata, pdata)
      # generate boxplot
      p <- ggplot(all_data, aes(x = Source, y = sum_corrected, fill = Source)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, alpha = 0.5) +
        labs(title = paste("Cutoff:", cutoff),
             y = "Functional Abundance Distribution",
             x = "") +
        theme_minimal() +
        theme(legend.position = "none", text = element_text(size=18))
      ggsave(output_file, plot = p, bg="White", width = 8, height = 6, dpi = 300)
      return(p)
    }
  }
}

#' mtestAllCutoffs - It generalizes the function testAllCutoffs. First, it creates 3 paths: one for transcriptomic TGCN, one for proteomic TGCN and one for the comparison report, where the results will be saved. Then it recalls the function testAllCutoff to create both transcriptomic and proteomic TGCNs. Finally, (if not done yet or overwrite==T), it generates the plots and results comparing the performance and features of the generated networks.
#' @param proteomicData a dataframe representing a protein expression matrix with proteins as rows and samples as columns.
#' @param proteomicTarget a vector representing the response variable to predict in the proteomic dataset(numeric or categorical).
#' @param proteomicCovs a data frame with donors metadata where samples are the rows and features are the columns (for proteomic data).
#' @param transcriptomicData a dataframe representing a gene expression matrix with genes as rows and samples as columns.
#' @param transcriptomicTarget a vector representing the response variable to predict in the transcriptomic dataset(numeric or categorical).
#' @param trannscriptomicCovs a data frame with donors metadata where samples are the rows and features are the columns (for proteomic data).
#' @param: targetName the name of the targeted variable.
#' @param: tissueName the name of the tissue, cohort or dataset.
#' @param: train.split the proportion of samples used to train the models in TGCNs.
#' @param nfolds the number of folds for the cross-validation in TGCNs.
#' @param t the number of LASSO runs in both TGCNs.
#' @param: overwrite: If overwrite=T, it overwrites the previous results.
#' @param: cutoff: selected cutoff for generating the plot.
#' @param cutoffs the ratios of appearance selected to create a final model and the TGCNs with the corresponding genes and proteins.
#' @param n the maximum size of a module for all the approaches in TGCNs.
#' @param s the minimum size of a module for all the approaches.
#' @param m the number of genes (or proteins) added in each iteration when approach="enrichment".
#' @param minCor the minimum correlation of a gene or protein to be added to a module in TGCNs.
#' @param maxTol the maximum number of tries to get a better enrichment without getting it.
#' @param approach the approach to complete the modules in both TGCNs (fixed, coefficient or enrichment).
#' @param seed the seed number to ensure the reproducibility of the results.
#' @param save if save=T, the results of the analysis and the comparative results will be saved in separate files.
#' @param overwrite if overwrite=T, the analysis will be repeated and files will be overwritten
#' @param useNonCodingOnly if useNonCodingOnly=T, TCGN for transcriptomic will use only non-coding transcripts.
#' @param core_path the path used as root to create the aforementioned subpaths.
#' @param report if true, an automated report for each TGCN will be created (¡add the comparison report too!)
#' @return
#' 
mtestAllCutoffs <- function(proteomicData,
                            proteomicTarget,
                            proteomicCovs,
                            transcriptomicData,
                            transcriptomicTarget,
                            transcriptomicCovs,
                            targetName="target",
                            tissueName="tissue",
                            train.split=0.7,
                            nfolds=5,
                            t=10,
                            cutoffs=10:1,
                            n=100,
                            s=10,
                            m=10,
                            minCor=0.3,
                            maxTol=3,
                            approach=c("fixed", "coefficient", "enrichment"),
                            seed=1234,
                            save=T,
                            overwrite=F,
                            core_path=getwd(),
                            reduced=F,
                            report=F,
                            useNonCodingOnly=F) {
  
  results <- list()
  
  # Directory is created for transcriptomics and proteomics TGCNs.
  transcriptomic_path <- paste0(core_path, "/transcriptomics/")
  if(!dir.exists(transcriptomic_path)) dir.create(transcriptomic_path, recursive = TRUE)
  
  proteomic_path <- paste0(core_path, "/proteomics/")
  if(!dir.exists(proteomic_path)) dir.create(proteomic_path, recursive = TRUE)
  
  # testAllCutoffs from package TGCN.
  
  # transcriptomic TGCN
  transcriptomic_r <- testAllCutoffs(exprData=transcriptomicData,
                                     target=transcriptomicTarget,
                                     covs=transcriptomicCovs,
                                     train.split=train.split,
                                     nfolds=nfolds,
                                     t=t,
                                     path=transcriptomic_path,
                                     targetName=targetName,
                                     tissueName=tissueName,
                                     seed=seed,
                                     cutoffs=cutoffs,
                                     n=n,
                                     m=m,
                                     s=s,
                                     minCor=minCor,
                                     maxTol=maxTol,
                                     save=save,
                                     overwrite=overwrite,
                                     approach=approach,
                                     report = report, 
                                     useNonCodingOnly = useNonCodingOnly)
  
  # proteomic TGCN
  proteomics_r <- testAllCutoffs(exprData=proteomicData,
                                 target=proteomicTarget,
                                 covs=proteomicCovs,
                                 train.split=train.split,
                                 nfolds=nfolds,
                                 t=t,
                                 path=proteomic_path,
                                 targetName=targetName,
                                 tissueName=tissueName,
                                 seed=seed,
                                 cutoffs=cutoffs,
                                 n=n,
                                 m=m,
                                 s=s,
                                 minCor=minCor,
                                 maxTol=maxTol,
                                 save=save,
                                 overwrite=overwrite,
                                 approach=approach,
                                 report=report,
                                 useNonCodingOnly = F)
  
  # create the directory for the comparison report
  cr_path <- paste0(core_path, "/comparison_report")
  
  if(!dir.exists(cr_path)) {
    dir.create(cr_path, recursive = TRUE)
    files_exist <- FALSE
  } else { 
    subfolders <- list.dirs(cr_path, recursive = FALSE)
    files_exist <- any(sapply(subfolders, function(folder) {
      length(list.files(folder, full.names = TRUE)) > 0
    })) }
  
  if (files_exist==F || report==T) {
    # Generate results and plots, saving them in results.
    cat("\n- Generando gráficas de selectBasicFeatures...\n")
    basic_plots <- selectBasicFeatures(transcriptomic_data=transcriptomicData,
                                       proteomic_data=proteomicData,
                                       cr_path=cr_path,
                                       tissueName = tissueName,
                                       targetName = targetName,
                                       overwrite=overwrite)
    
    # save the plots
    results[["basicPlots"]][["dimRatio"]] <- basic_plots$dimratio
    results[["basicPlots"]][["vennPlot"]] <- basic_plots$venn
    
    # readRDS with the hubs and nets generated in TGCNs for both transcriptomic and proteomic.
    transcriptomic_df <- readRDS(file = paste0(transcriptomic_path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds")) 
    transcriptomic_net <- readRDS(file = paste0(transcriptomic_path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds")) 
    proteomic_df <- readRDS(file = paste0(proteomic_path, "/hubGenes/", targetName, "_", tissueName, "_hubGenes.rds"))
    proteomic_net <- readRDS(file = paste0(proteomic_path, "/Net/", targetName, "_", tissueName, "_TGCNs.rds")) 
    
    # select the relevant cutoff for both one of the TGCNs
    transcriptomic_nHubs <- lengths(lapply(transcriptomic_df$LM, function(x) x$lm_genes_selected$gene))
    transcriptomic_nHubs <- transcriptomic_nHubs[transcriptomic_nHubs <= 30]
    transcriptomic_names <- as.numeric(gsub("cutoff", "", names(transcriptomic_nHubs)))
    
    proteomic_nHubs <- lengths(lapply(proteomic_df$LM, function(x) x$lm_genes_selected$gene))
    proteomic_nHubs <- proteomic_nHubs[proteomic_nHubs <= 30]
    proteomic_names <- as.numeric(gsub("cutoff", "", names(proteomic_nHubs)))
    
    if (length(transcriptomic_names) <= length(proteomic_names)) {
      selected_cutoffs <- transcriptomic_names
    } else {
      selected_cutoffs <- proteomic_names
    }
    
    # performance plot
    molPLot <- selectSeedMolecules(transcriptomic_path = transcriptomic_path,
                                   proteomic_path = proteomic_path,
                                   tissueName = tissueName,
                                   targetName = targetName,
                                   cr_path = cr_path,
                                   overwrite = overwrite,
                                   selectedCutoffs = selected_cutoffs)
    
    results[["R2_plot"]] <- molPLot
    
    for (cutoff in selected_cutoffs) {
      cat("\n Creating comparison plots for cutoff", cutoff, "\n")
      
      # modules overlap betweeen hubs.
      modOverlap <- mgetModulesOverlap(main_path = core_path,
                                           targetName = targetName,
                                           tissueName = tissueName,
                                           transcriptomicData = transcriptomicData,
                                           proteomicData = proteomicData,
                                           cutoff=cutoff
      )
      
      results[["modulesOverlap"]][[paste0("c", cutoff)]] <- list(
        plot = modOverlap$plot,
        significant_results = modOverlap$significant_results
      )
      
      # enrichment distribution for both nets.
      enrichmentDistrib <- plotEnrichmentDistribution(cr_path=cr_path,
                                                      transcriptomic_path=transcriptomic_path,
                                                      proteomic_path=proteomic_path,
                                                      targetName=targetName,
                                                      tissueName=tissueName,
                                                      overwrite=overwrite,
                                                      cutoff=cutoff)
      
      results[["enrichmentDistribution"]][[paste0("c", cutoff)]] <- enrichmentDistrib
      
    }
    
    if (save==T) {
      dir.create(paste0(cr_path, "/results/"), recursive = TRUE, showWarnings = FALSE)
      saveRDS(results, paste0(cr_path, "/results/", targetName, "_", tissueName, "_mTGCNs.rds"))
    }
  } else {
    results <- readRDS(paste0(cr_path, "/results/", targetName, "_", tissueName, "_mTGCNs.rds"))
  }
  # if(report==T) {
  #   template=list.files() ...
  # }
  return(results) 
}