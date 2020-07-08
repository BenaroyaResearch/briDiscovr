
#' Perform metaclustering for a DISCOV-R experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment},
#' clustered using \code{clusterDiscovrExperiment}, and metaclustered with
#' \code{metaclusterDiscovrExperiment}
#' @param parentSubset A string indicating which subset represents the parent population (eg: all CD8+ cells).
#' If this is left NA, then it's assumed the first subset listed in the experiment's fcsInfo file is the parent
#' population and all other subsets are rare populations. (default: NA)
#' @param childSubsets A vector of strings indicating which subsets represent rare populations (eg: Tmr+ cells).
#' If this is left NA, then it's assumed the first subset listed in the experiment's fcsInfo file is the parent
#' population and all other subsets are rare populations. (default: NA)
#' @param filenamePrefix A string indicating the prefix to use for each file that gets saved.
#' This can include a directory path. (default: paste0(format(Sys.Date(), "%y%m%d"), "-Zscores-perSubj"))
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{setupDiscovrExperiment}} \code{\link{clusterDiscovrExperiment}} \code{\link{metaclusterDiscovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import dplyr
#' @importFrom flowCore exprs
#' @importFrom igraph graph.data.frame cluster_louvain membership modularity
#' @export
makeMetaclusterHeatmaps <- function(
  experiment,
  parentSubset = NA,
  childSubsets = NA,
  filenamePrefix = paste0(format(Sys.Date(), "%y%m%d"), "-Zscores-perSubj"),
  verbose = TRUE
){
  # check for experiment data object
  if(!is.discovrExperiment(experiment)){
    stop(
      "The object passed to this function is not a valid DISCOV-R experiment object. ",
      "Please create your experiment using the 'setupDiscovrExperiment' function ",
      "and perform the initial clustering using the 'clusterDiscovrExperiment' function. "
    )
  }
  if(experiment$status != "metaclustered"){
    stop(
      "The experiment must have already been metaclustered in order to print heatmaps using this function. ",
      "The current experiment has a status of ", experiment$status, ". ",
      "Please see the function 'metaclusterDiscovrExperiment' for details on metaclustering an experiment."
    )
  }
  if(!all(c("hmapDfAllSubsets", "allSubsetAllSubjectZscores", "allSubsetAllSubjectArcsinh") %in% names(experiment))){
    stop(
      "The object does not have the data expected in a metaclustered experiment. ",
      "Please re-run the clustering and metaclustering steps and try again."
    )
  }

  # Check requested subsets
  subsets <- unique(experiment$fcsInfo$cellSubset)

  if(is.na(parentSubset)){
    message("No parent population specified. Assumuning parent population is: ", subsets[1])
    parentSubset <- subsets[1]
  }

  if(is.na(childSubsets)){
    message(
      "No child populations specified. Assuming child populations are: ",
      paste0(subsets[2:length(subsets)], collapse = "; ")
    )
    childSubsets <- subsets[2:length(subsets)]
  }

  if(!parentSubset %in% subsets){
    stop(
      "The requested parent subset '", parentSubset, "' was not found in the list of cell subsets: ",
      paste0(subsets, collapse = "; ")
    )
  }
  if(!all(childSubsets %in% subsets)){
    stop(
      "At least one of the requested child subsets (",
      paste0(childSubsets, collapse = ";"),
      ") was not found in the list of all cell subsets: ",
      paste0(subsets, collapse = "; ")
    )
  }

  # TODO: check color palettes
  # TODO: make arguments to set new pals easily
  ### ln 225-250

  # blank out the legend titles by setting params to an empty list
  legendParams = list()
  exportWidth = 900
  exportHeight = 900
  titleFontParam = grid::gpar(fontface = "bold", fontsize = 15)
  marker_label_gp = grid::gpar(fontsize = 13)

  # Set up a palette for the heatmap
  my_zscore_pal <- circlize::colorRamp2(c(-2,0,2), c("#00ffffFF", "#000000FF", "#FDE725FF"))
  my_arcsinh_pal <- circlize::colorRamp2(c(0,5), c("#440154FF", "#FDE725FF"))

  # get subject ids, assign colors to them
  subject_ids <- sort(unique(experiment$clusterRarePopCts$samp))
  subject_id_colors <- data.frame(subject = subject_ids,
                                  #color = rainbow(length(subject_ids)))
                                  color = c("#000000FF",
                                            viridisLite::plasma(length(subject_ids)+2)[4:(length(subject_ids)+2)]))

  # pass_threshold_colors <- data.frame(pass = c(T,F), color = c("#FF0000", "#888888"))
  # pass_threshold_colors <- data.frame(pass = c("no", "yes", "enriched"), color = c("#888888", "#FF0000", "#0000FF"))
  # pass_threshold_colors <- data.frame(pass = c(0:8), color = c("#FFFFFF", "#3F054E", "#409E9B", "#CC0101", "#3F054E", "#409E9B", "#CC0101", "#5DC762", "#FFFFFF"))
  #
  # duration_colors <- data.frame(status = c("LS", "RO"), color = c("#2D297A", "#BB9C3E"))
  # progression_colors <- data.frame(status = c("Pre-slow", "Slow", "Unclassified", "Rapid", "Post-rapid"),
  #                                  color = c("#184349", "#368BA1", "#D3D3D3", "#5E9130", "#21632E"))

  #########################################################################
  # Section 3.d from original SOP - create heatmaps
  #########################################################################

  allSubsetZscoreAnnoDf <-
    experiment$hmapDfAllSubsets %>%
    dplyr::select(sample, !!subsets) %>%
    unique %>%
    dplyr::slice(match(colnames(experiment$allSubsetAllSubjectZscores), .$sample)) %>%
    dplyr::left_join(
      data.frame(group = experiment$colIndices) %>% rownames_to_column("sample"),
      by = "sample"
    ) %>%
    dplyr::mutate(subject = str_replace(.$sample, "_[0-9]+$", "")) %>%
    dplyr::select(subject, group, !!childSubsets)


  # TODO: add an argument to facilitate more annotations
  # TODO: add an argument for a palette for each additional anno
  # %>%
  # left_join(select(clinical_vars, ID, progression = `Progression rate`, duration = Duration), by = c("subject" = "ID"))

  # set up colors as a named list
  allSubsetZscoreAnnoColors <-
    list(
      subject = setNames(
        subject_id_colors$color,
        as.character(subject_id_colors$subject)
      ),
      # for 1-12MCs
      group = setNames(
        RColorBrewer::brewer.pal(experiment$kGroups, "Set3"),
        as.character(unique(experiment$colIndices))
      )
    )
  #group = setNames(c(brewer.pal(12, "Set3"), brewer.pal(kGroups-12, "Set2")), #for 12-20MCs
  #         as.character(unique(colIndices)))) #,
  # group = setNames(viridis(kGroups),
  #          as.character(unique(colIndices)))) # original viridis scheme

  # add metacluster information to tmr counting information
  subsetEventCounting <- left_join(
    experiment$subsetEventCounting,
    data.frame(
      sample = colnames(experiment$allSubsetAllSubjectZscores),
      group = experiment$colIndices
    )
  )

  ## Use for pct rare subset color scheme
  for (t in childSubsets){
    allSubsetZscoreAnnoColors[[t]] = c(setNames("white", 0),
                                        setNames(colorRampPalette(c("grey95", "grey0"))(100), seq(0.01, 1.00, 0.01)))
  }

  allSubsetZscoreAnno <-
    ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(allSubsetZscoreAnnoDf %>% mutate_all(funs(ifelse(is.na(.), 0, .)))),
      col = allSubsetZscoreAnnoColors,
      show_annotation_name = T,
      show_legend = F)
  # annotation_legend_param = list(Islet = list(ncol = 4, title_position = "topcenter"),
  #                                group = list(ncol = 2, title_position = "topcenter"),
  #                                subject = list(ncol = 4, title_position = "topcenter")),
  # show_legend = c(F, T, T, T, T))

  zscore_hmap <- ComplexHeatmap::Heatmap(
    experiment$allSubsetAllSubjectZscores,
    col = my_zscore_pal,
    name = "z-score",
    column_title = paste0("Cluster Phenotypes from All Subjects"),
    column_title_gp = titleFontParam,
    #column_dend_height=unit(10,"cm"),
    clustering_method_columns = experiment$linkage,
    cluster_rows = T,
    clustering_method_rows = experiment$linkage,
    show_column_names = F,
    row_names_gp = marker_label_gp,
    top_annotation = allSubsetZscoreAnno,
    heatmap_legend_param = legendParams
  )

  png(
    filename = paste0(
      filenamePrefix, "_", experiment$linkage, "_", experiment$distance, "_",
      "_AllClusters_cutree", experiment$nMetaclusters, "_pct_Islet.png"
    ),
    width = exportWidth,
    height = exportHeight
  )
  print(zscore_hmap)
  dev.off()

  marker_order = ComplexHeatmap::row_order(zscore_hmap)#[[1]]


  #MFI plot did not work
  png(
    filename = paste0(
      filenamePrefix, "_", experiment$linkage, "_", experiment$distance, "_",
      "_AllClusters_cutree", experiment$nMetaclusters, "_pct_tmrs_MFI.png"
    ),
    width = exportWidth,
    height = exportHeight
  )
  print(ComplexHeatmap::Heatmap(
    experiment$allSubsetAllSubjectArcsinh,
    col = my_arcsinh_pal,
    name = "MFI",
    column_title = paste0("Cluster Phenotypes from All Subjects"),
    column_title_gp = titleFontParam,
    cluster_columns = F,
    column_order = ComplexHeatmap::column_order(zscore_hmap),
    cluster_rows = F,
    row_order = experiment$metaclusterMarkers[marker_order],
    #split = rep(c("", " "), c(length(marker_order), 2)),
    gap = unit(5, "mm"),
    show_column_names = F,
    row_names_gp = marker_label_gp,
    top_annotation = allSubsetZscoreAnno,
    heatmap_legend_param = legendParams)
  )
  dev.off()


  # Weighted Average phenotype hmap - average of arcsinch values in each metacluster
  per_metax_avg <-
    data.frame(t(experiment$allSubsetAllSubjectArcsinh)) %>%
    dplyr::mutate(sample = rownames(.)) %>%
    merge(subsetEventCounting, by = "sample") %>%
    dplyr::rowwise() %>%
    mutate(
      n_event = sum(c_across(subsets)),
      subj = str_remove(sample, "_.*"),
      RP_clust = str_remove(sample, ".*_")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(subj) %>%
    dplyr::mutate(subj_event = sum(n_event)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = n_event/subj_event) %>%
    dplyr::mutate_at(experiment$metaclustergMarkers, funs(.*proportion)) %>%
    dplyr::select(-sample, -!!subsets, -n_event, -subj, -RP_clust, -subj_event) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate_at(experiment$metaclusterMarkers, funs(./proportion)) %>%
    dplyr::select(-group, -proportion) %>%
    t

  colnames(per_metax_avg) = 1:experiment$nMetaclusters

  # set up data frame for annotation bar for average hmap
  tmr_avg_anno_df <-
    data.frame(group = colnames(per_metax_avg))

  # set up colors as a named list
  tmr_avg_anno_colors <- list(
    #1-12MCs
    group = setNames(
      RColorBrewer::brewer.pal(experiment$kGroups, "Set3"),
      as.character(colnames(per_metax_avg))
    )
  )
  #group = setNames(c(brewer.pal(12, "Set3"), brewer.pal(kGroups-12, "Set2")),        #12-20MCs
  #               as.character(colnames(per_metax_avg))))

  tmr_avg_anno <- ComplexHeatmap::HeatmapAnnotation(df = tmr_avg_anno_df,
                                    col = tmr_avg_anno_colors,
                                    show_annotation_name = T,
                                    show_legend = T)


  #Did not work Error in .local(object, ...) : Number of rows in the matrix are not the same as the length of the cluster or the row orders.
  tmr_hm <- ComplexHeatmap::Heatmap(
    per_metax_avg,
    col = my_arcsinh_pal,
    name = "MFI",
    column_title = paste0("Clusters Enriched for ", experiment$nMetaclusters, " Tmrs"),
    column_title_gp = titleFontParam,
    cluster_columns = F,
    cluster_rows = F,
    row_order = experiment$metaclusterMarkers[marker_order],
    #split = rep(c("", " "), c(length(marker_order), 2)),
    gap = unit(5, "mm"),
    #combined_name_fun = NULL,
    show_column_names = F,
    row_names_gp = marker_label_gp,
    top_annotation = tmr_avg_anno,
    heatmap_legend_param = legendParams
  )


  png(
    filename = paste0(
      filenamePrefix, "_", experiment$linkage, "_rpheno_", experiment$nMetaclusters, "_arcsinh_weightedAvgByClust.png"
    ),
    width = 2500,
    height = 2500,
    res = 300
  )
  print(tmr_hm)
  dev.off()

  ## Weighted Average phenotype for total CD8s
  # weighted by proportion of cells in a subject (not by number of cells) in order to allow
  # equivalent contributions from subjects with disparate numbers of collected events
  cd8_avg <-
    data.frame(t(experiment$allSubsetAllSubjectArcsinh)) %>%
    dplyr::mutate(sample = rownames(.)) %>%
    merge(experiment$subsetEventCounting, by = "sample") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_event = sum(c_across(subsets)),
      subj = str_remove(sample, "_.*"),
      RP_clust = str_remove(sample, ".*_")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(subj) %>%
    dplyr::mutate(subj_event = sum(n_event)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = n_event/subj_event) %>%
    dplyr::mutate_at(experiment$metaclusterMarkers, funs(.*proportion)) %>%
    dplyr::select(-sample, -!!subsets, -n_event, -subj, -RP_clust, -subj_event) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate_at(experiment$metaclusterMarkers, funs(./proportion)) %>%
    dplyr::select(-proportion) %>%
    t

  colnames(cd8_avg) = "CD8"

  #Did not work
  cd8_hm <- ComplexHeatmap::Heatmap(
    cd8_avg,
    col = my_arcsinh_pal,
    name = "MFI",
    column_title = paste0("Total CD8+"),
    column_title_gp = titleFontParam,
    cluster_columns = F,
    cluster_rows = F,
    row_order = experiment$metaclusterMarkers[marker_order],
    #split = rep(c("", " "), c(length(marker_order), 2)),
    gap = unit(5, "mm"),
    #combined_name_fun = NULL,
    show_column_names = F,
    row_names_gp = marker_label_gp,
    heatmap_legend_param = legendParams
  )

  png(filename = paste0(filenamePrefix, "_", experiment$linkage, "_CD8_arcsinh_weightedAvg.png"),
      width = 700,
      height = 2500,
      res = 300)
  print(cd8_hm)
  dev.off()

  #########################################################################
  # Section 3.e from original SOP - export data
  #########################################################################
}
