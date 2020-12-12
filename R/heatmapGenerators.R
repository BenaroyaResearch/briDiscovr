## Copyright (C) 2020  Mario Rosasco and Benaroya Research Institute
##
## This file is part of the briDiscovr package

#' Generate heatmaps from a metaclustered experiment
#'
#' @param experiment A discovrExperiment created using \code{setupDiscovrExperiment},
#' clustered using \code{clusterDiscovrExperiment}, and metaclustered with
#' \code{metaclusterDiscovrExperiment}
#' @param childSubsets A vector of strings indicating which subsets represent rare populations (eg: Tmr+ cells).
#' If this is left NA, then it's assumed all non-parent subsets are rare populations. (default: NA)
#' @param dropMarkers A string or vector of strings indicating names of markers to exclude
#' from display. Note that printing the experiment object will display a list of all
#' markers. By default all markers will be included. Markers not used in metaclustering will be shown in
#' a separate heatmap. Also note that you cannot drop markers that were used in metaclustering, and should instead
#' re-run the metclustering step, dropping the markers at that stage. (default: NA)
#' @param filenamePrefix A string indicating the prefix to use for each file that gets saved.
#' This can include a directory path. If left as NA, will use the current date as the prefix
#' with the format YYMMDD. (default: NA)
#' @param parentTitle A string to use to label the parent population heatmaps (default: "Parent")
#' @param metaclusterColors An optional vector of color values to use to label metaclusters in heatmaps.
#' If left as the default 'NA', a random list of visually distinct colors will be generated.
#' Note that if a color vector is supplied, it must contain at least as many colors as there are metaclusters. (default: NA)
#' @param zScoreBreaks A vector of numbers indicating the values used to calibrate the Z score color scale. This vector
#' must have the same number of values as \code{zScoreColors}. (default: c(-2,0,2))
#' @param zScoreColors A vector of strings indicating the colors used to represent the Z score value range. This vector
#' must have the same number of values as \code{zScoreBreaks}. (default: c("#00FFFF", "#000000", "#FDE725"))
#' @param intensityBreaks A vector of numbers indicating the values used to calibrate the transformed marker intensity
#' scale. This vector must have the same number of values as \code{intensityColors} (default: c(0,2,4,6))
#' @param intensityColors A vector of strings indicating the colors used to represent the transformed marker intensity
#' value range. This vector must have the same number of values as \code{intensityBreaks}
#' (default: c("#4A0E53", "#925F54", "#D6B343", "#FDE727"))
#' @param columnDendHeight A numeric indicating the column dendrogram height (in mm)
#' for the Z score cluster phenotypes heatmap (default: 10)
#' @param rowDendWidth A numeric indicating the row dendrogram width (in mm)
#' for the Z score cluster phenotypes heatmap (default: 10)
#' @param verbose A boolean specifying whether to display processing messages (default: TRUE)
#' @return An S3 object of class \code{discovrExperiment}
#'
#' @seealso \code{\link{setupDiscovrExperiment}} \code{\link{clusterDiscovrExperiment}} \code{\link{metaclusterDiscovrExperiment}}
#' @author Mario G Rosasco, \email{mrosasco@@benaroyaresearch.org}, Virginia Muir
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_remove str_extract
#' @importFrom rlang .data
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom grid unit
#' @importFrom circlize colorRamp2
#' @export
makeMetaclusterHeatmaps <- function(
  experiment,
  childSubsets = NA,
  dropMarkers = NA,
  filenamePrefix = NA,
  parentTitle = "Parent",
  metaclusterColors = NA,
  zScoreBreaks = c(-2,0,2),
  zScoreColors = c("#00FFFF", "#000000", "#FDE725"),
  intensityBreaks = c(0,2,4,6),
  intensityColors = c("#4A0E53", "#925F54", "#D6B343", "#FDE727"),
  columnDendHeight = 10,
  rowDendWidth = 10,
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

  # Check marker names
  if(!all(is.na(dropMarkers))){
    if (!all(dropMarkers %in% experiment$markerInfo$commonMarkerName)){
      stop(
        "The following markers were not found in the dataset: ",
        paste0(dropMarkers[!dropMarkers %in% experiment$markerInfo$commonMarkerName], collapse = ", "),
        ". Please remove these markers from the 'dropMarkers' argument."
      )
    }
    if (any(dropMarkers %in% experiment$metaclusterMarkers)){
      stop(
        "The markers ", dropMarkers[dropMarkers %in% experiment$metaclusterMarkers],
        " were listed in 'dropMarkers' but were used for metaclustering and cannot be dropped. ",
        "Please re-run metaclustering without these markers to remove them."
      )
    }
  }

  # Check requested subsets
  subsets <- unique(experiment$fcsInfo$cellSubset)

  parentSubset <- experiment$parentPopulation

  if(is.na(childSubsets)){
    childSubsets <- subsets[subsets != parentSubset]
    message(
      "No child populations specified. Assuming child populations are: ",
      paste0(childSubsets, collapse = "; ")
    )
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

  # Set up the file name prefix
  if(is.na(filenamePrefix)){
    filenamePrefix <- format(Sys.Date(), "%y%m%d-")
  }
  filenamePrefix = paste0(
    filenamePrefix,
    experiment$linkage, "_",
    experiment$distance, "_",
    experiment$pctInClusterThreshold, "pctThresh_",
    experiment$nMetaclusters, "metaclusters"
  )

  # Set up color palettes
  if(all(is.na(metaclusterColors))){
    metaclusterColors = getColorList(experiment$kGroups)
  }
  if (length(metaclusterColors) != experiment$kGroups){
    stop("The length of 'metaclusterColors' must be the same as the number of metaclusters. ",
         length(metaclusterColors), " colors were provided for ", experiment$kGroups, " metaclusters.")
  }

  ### ln 225-250

  # Display Parameters
  exportWidth = 900
  exportHeight = 900
  titleFontParam = grid::gpar(fontface = "bold", fontsize = 15)
  marker_label_gp = grid::gpar(fontsize = 13)

  # Set up palettes for the heatmap
  if(length(zScoreBreaks) != length(zScoreColors)){
    stop("zScoreBreaks and zScoreColors must be the same length.")
  }
  if(length(intensityBreaks) != length(intensityColors)){
    stop("intensityBreaks and intensityColors must be the same length.")
  }
  my_zscore_pal <- circlize::colorRamp2(
    breaks = zScoreBreaks,
    colors = zScoreColors
  )
  my_arcsinh_pal <- circlize::colorRamp2(
    breaks = intensityBreaks,
    colors = intensityColors
  )

  # Set up legend labels to accurately represent the data
  zScoreMaxLabel = ifelse(
    max(zScoreBreaks) < max(experiment$allSubsetAllSubjectZscores),
    paste0("> ", max(zScoreBreaks)),
    max(zScoreBreaks)
  )
  zScoreMinLabel = ifelse(
    min(zScoreBreaks) > min(experiment$allSubsetAllSubjectZscores),
    paste0("< ", min(zScoreBreaks)),
    min(zScoreBreaks)
  )
  intensityMaxLabel = ifelse(
    max(intensityBreaks) < max(experiment$allSubsetAllSubjectArcsinh),
    paste0("> ", max(intensityBreaks)),
    max(intensityBreaks)
  )
  intensityMinLabel = ifelse(
    min(intensityBreaks) > min(experiment$allSubsetAllSubjectArcsinh),
    paste0("< ", min(intensityBreaks)),
    min(intensityBreaks)
  )
  intensityBreaks = sort(intensityBreaks)
  zScoreBreaks = sort(zScoreBreaks)

  intensityBreakLabels = intensityBreaks
  intensityBreakLabels[1] = intensityMinLabel
  intensityBreakLabels[length(intensityBreakLabels)] = intensityMaxLabel
  zScoreBreakLabels = zScoreBreaks
  zScoreBreakLabels[1] = zScoreMinLabel
  zScoreBreakLabels[length(zScoreBreakLabels)] = zScoreMaxLabel

  zScoreLegendParam = list(
    at = zScoreBreaks,
    labels = zScoreBreakLabels
  )
  intensityLegendParam = list(
    at = intensityBreaks,
    labels = intensityBreakLabels
  )

  # get subject ids, assign colors to them
  subject_ids <- sort(unique(experiment$clusterRarePopCts$samp))
  subject_id_colors <- data.frame(
    subject = subject_ids,
    color = viridisLite::plasma(length(subject_ids), begin = 0.2)
  )

  #########################################################################
  # Section 3.d from original SOP - create heatmaps
  #########################################################################

  allSubsetZscoreAnnoDf <-
    experiment$hmapDfAllSubsets %>%
    dplyr::select(.data$sample, !!subsets) %>%
    unique %>%
    dplyr::slice(match(colnames(experiment$allSubsetAllSubjectZscores), .data$sample)) %>%
    dplyr::left_join(
      data.frame(group = experiment$colIndices) %>% rownames_to_column("sample"),
      by = "sample"
    ) %>%
    dplyr::mutate(subject = str_replace(.data$sample, "_[0-9]+$", "")) %>%
    dplyr::select(.data$subject, .data$group,  parentSubset, !!childSubsets)

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
      # assign metacluster names to colors
      group = setNames(
        metaclusterColors,
        as.character(unique(experiment$colIndices))
      )
    )

  # add metacluster information to tmr counting information
  subsetEventCounting <- left_join(
    experiment$subsetEventCounting,
    data.frame(
      sample = colnames(experiment$allSubsetAllSubjectZscores),
      group = experiment$colIndices
    )
  )

  ## Use for pct rare subset color scheme
  for (t in c(parentSubset, childSubsets)){
    allSubsetZscoreAnnoColors[[t]] = c(
      setNames("white", 0),
      setNames(colorRampPalette(c("grey95", "grey0"))(100), seq(0.01, 1.00, 0.01))
    )
  }

  allSubsetZscoreAnno <-
    ComplexHeatmap::HeatmapAnnotation(
      df = data.frame(allSubsetZscoreAnnoDf %>% mutate_all(function(val){ifelse(is.na(val), 0, val)})),
      col = allSubsetZscoreAnnoColors,
      show_annotation_name = T,
      show_legend = F)
  # annotation_legend_param = list(Islet = list(ncol = 4, title_position = "topcenter"),
  #                                group = list(ncol = 2, title_position = "topcenter"),
  #                                subject = list(ncol = 4, title_position = "topcenter")),
  # show_legend = c(F, T, T, T, T))

  zscore_hmap <- ComplexHeatmap::Heatmap(
    experiment$allSubsetAllSubjectZscores[experiment$metaclusterMarkers,],
    col = my_zscore_pal,
    name = "z-score",
    # column styling
    column_title = paste0("Cluster Phenotypes from All Samples"),
    column_title_gp = titleFontParam,
    column_dend_height=unit(columnDendHeight, "mm"),
    clustering_method_columns = experiment$linkage,
    show_column_names = FALSE,
    # row styling
    cluster_rows = TRUE,
    clustering_method_rows = experiment$linkage,
    row_dend_width=unit(rowDendWidth, "mm"),
    row_names_gp = marker_label_gp,
    # overall styling
    top_annotation = allSubsetZscoreAnno,
    heatmap_legend_param = zScoreLegendParam
  )

  # write the primary z-score heatmap to a file
  png(
    filename = paste0(filenamePrefix, "_allClusters_zScore.png"),
    width = exportWidth,
    height = exportHeight
  )
  print(zscore_hmap)
  dev.off()

  #################################################################################
  # set marker order for other plots based on the z-score heatmap,
  # and identify any non-clustering markers to show in detached heatmaps
  markerOrderIndices = ComplexHeatmap::row_order(zscore_hmap)
  zscoreRownames = rownames(experiment$allSubsetAllSubjectZscores[experiment$metaclusterMarkers,])
  markerOrder = zscoreRownames[markerOrderIndices]

  additionalMarkers = setdiff(experiment$hmapDfAllSubsets$marker, experiment$metaclusterMarkers)
  additionalMarkers = setdiff(additionalMarkers, dropMarkers)

  # MFI plot - metaclustered markers
  png(
    filename = paste0(filenamePrefix,"_allClusters_metacluster_MFI.png"),
    width = exportWidth,
    height = exportHeight
  )
  print(ComplexHeatmap::Heatmap(
    experiment$allSubsetAllSubjectArcsinh[markerOrder,],
    col = my_arcsinh_pal,
    name = "MFI",
    column_title = paste0("Cluster Phenotypes from All Samples"),
    column_title_gp = titleFontParam,
    cluster_columns = F,
    column_order = ComplexHeatmap::column_order(zscore_hmap),
    cluster_rows = F,
    row_order = markerOrder,
    gap = unit(5, "mm"),
    show_column_names = F,
    row_names_gp = marker_label_gp,
    top_annotation = allSubsetZscoreAnno,
    heatmap_legend_param = intensityLegendParam)
  )
  dev.off()

  # Weighted Average phenotype hmap - average of arcsinch values in each metacluster
  perMetaclusterAvg <-
    data.frame(t(experiment$allSubsetAllSubjectArcsinh)) %>%
    tibble::rownames_to_column("sample") %>%
    merge(subsetEventCounting, by = "sample") %>%
    dplyr::rowwise() %>%
    mutate(
      n_event = sum(c_across(subsets)),
      subj = stringr::str_remove(.data$sample, "_[0-9]+$"),
      RP_clust = stringr::str_extract(.data$sample, "[0-9]+$")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$subj) %>%
    dplyr::mutate(subj_event = sum(.data$n_event)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = .data$n_event/.data$subj_event) %>%
    dplyr::mutate_at(experiment$markerInfo$commonMarkerName, list(~.*.data$proportion)) %>%
    dplyr::select(
      -.data$sample,
      -!!subsets,
      -.data$n_event,
      -.data$subj,
      -.data$RP_clust,
      -.data$subj_event) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate_at(experiment$markerInfo$commonMarkerName, list(~./.data$proportion)) %>%
    dplyr::select(-.data$group, -.data$proportion) %>%
    t

  colnames(perMetaclusterAvg) = 1:experiment$nMetaclusters

  # set up data frame for annotation bar for average hmap
  tmr_avg_anno_df <-
    data.frame(group = colnames(perMetaclusterAvg))

  # set up colors as a named list
  tmr_avg_anno_colors <- list(
    #1-12MCs
    group = setNames(
      metaclusterColors,
      as.character(colnames(perMetaclusterAvg))
    )
  )

  tmr_avg_anno <- ComplexHeatmap::HeatmapAnnotation(df = tmr_avg_anno_df,
                                    col = tmr_avg_anno_colors,
                                    show_annotation_name = T,
                                    show_legend = T)


  # Metacluster average heatmap - metaclustering markers
  tmr_hm <- ComplexHeatmap::Heatmap(
    perMetaclusterAvg[markerOrder,],
    col = my_arcsinh_pal,
    name = "MFI",
    column_title = paste(experiment$nMetaclusters, "Metaclusters - Weighted Average MFI"),
    column_title_gp = titleFontParam,
    cluster_columns = F,
    cluster_rows = F,
    row_order = markerOrder,
    gap = unit(5, "mm"),
    show_column_names = F,
    row_names_gp = marker_label_gp,
    top_annotation = tmr_avg_anno,
    heatmap_legend_param = intensityLegendParam
  )


  png(
    filename = paste0(filenamePrefix, "_metaclusters_arcsinh_weightedAvgByClust.png"),
    width = 2500,
    height = 2500,
    res = 300
  )
  print(tmr_hm)
  dev.off()



  ## Weighted Average phenotype for total CD8s
  # weighted by proportion of cells in a subject (not by number of cells) in order to allow
  # equivalent contributions from subjects with disparate numbers of collected events
  parentPopulationAvg <-
    data.frame(t(experiment$allSubsetAllSubjectArcsinh)) %>%
    tibble::rownames_to_column("sample") %>%
    merge(experiment$subsetEventCounting, by = "sample") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      n_event = sum(c_across(subsets)),
      subj = stringr::str_remove(.data$sample, "_[0-9]+$"),
      RP_clust = stringr::str_extract(.data$sample, "[0-9]+$")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$subj) %>%
    dplyr::mutate(subj_event = sum(.data$n_event)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = .data$n_event/.data$subj_event) %>%
    dplyr::mutate_at(experiment$markerInfo$commonMarkerName, list(~.*.data$proportion)) %>%
    dplyr::select(
      -.data$sample,
      -!!subsets,
      -.data$n_event,
      -.data$subj,
      -.data$RP_clust,
      -.data$subj_event
    ) %>%
    dplyr::summarise_all(sum) %>%
    dplyr::mutate_at(experiment$markerInfo$commonMarkerName, list(~./.data$proportion)) %>%
    dplyr::select(-.data$proportion) %>%
    t

  colnames(parentPopulationAvg) = parentTitle

  # Total Marker Average - metaclustering markers
  parentWeightedAvgHeatmap <- ComplexHeatmap::Heatmap(
    parentPopulationAvg[markerOrder,],
    col = my_arcsinh_pal,
    name = "MFI",
    column_title = paste0("Total ", parentTitle),
    column_title_gp = titleFontParam,
    cluster_columns = F,
    cluster_rows = F,
    row_order = markerOrder,
    gap = unit(5, "mm"),
    show_column_names = F,
    row_names_gp = marker_label_gp,
    heatmap_legend_param = intensityLegendParam
  )

  png(filename = paste0(filenamePrefix, "_total", parentTitle, "_arcsinh_weightedAvg.png"),
      width = 700,
      height = 2500,
      res = 300)
  print(parentWeightedAvgHeatmap)
  dev.off()


  # if there are additional markers to make plots for, do so
  if(!all(is.na(additionalMarkers)) && length(additionalMarkers) > 0){
    # MFI plot - additional markers
    png(
      filename = paste0(filenamePrefix, "_allClusters_nonClusteringMarkers_MFI.png"),
      width = exportWidth,
      height = exportHeight * (length(additionalMarkers)/length(markerOrder))
    )
    print(ComplexHeatmap::Heatmap(
      experiment$allSubsetAllSubjectArcsinh[additionalMarkers,, drop=FALSE],
      col = my_arcsinh_pal,
      name = "MFI",
      column_title = paste0("Cluster Phenotypes from All Samples (Non-Clustering Markers)"),
      column_title_gp = titleFontParam,
      cluster_columns = F,
      column_order = ComplexHeatmap::column_order(zscore_hmap),
      cluster_rows = F,
      row_order = additionalMarkers,
      gap = unit(5, "mm"),
      show_column_names = F,
      row_names_gp = marker_label_gp,
      top_annotation = allSubsetZscoreAnno,
      heatmap_legend_param = intensityLegendParam)
    )
    dev.off()

    # Total Marker Average - additional markers
    parentWeightedAvgHeatmap <- ComplexHeatmap::Heatmap(
      parentPopulationAvg[additionalMarkers,, drop=FALSE],
      col = my_arcsinh_pal,
      name = "MFI",
      column_title = paste0("Total ", parentTitle),
      column_title_gp = titleFontParam,
      cluster_columns = F,
      cluster_rows = F,
      row_order = additionalMarkers,
      gap = unit(5, "mm"),
      show_column_names = F,
      row_names_gp = marker_label_gp,
      heatmap_legend_param = intensityLegendParam
    )

    png(
      filename = paste0(filenamePrefix, "_total", parentTitle, "_nonClusteringMarkers_arcsinh_weightedAvg.png"),
      width = 700,
      height = 2500 * (length(additionalMarkers)/length(markerOrder)),
      res = 300
    )
    print(parentWeightedAvgHeatmap)
    dev.off()

    # Metacluster average heatmap - additional markers
    tmr_hm <- ComplexHeatmap::Heatmap(
      perMetaclusterAvg[additionalMarkers,, drop=FALSE],
      col = my_arcsinh_pal,
      name = "MFI",
      column_title = paste(experiment$nMetaclusters, "Metaclusters - Weighted Average MFI (Non-Clustering Markers)"),
      column_title_gp = titleFontParam,
      cluster_columns = F,
      cluster_rows = F,
      row_order = additionalMarkers,
      gap = unit(5, "mm"),
      show_column_names = F,
      row_names_gp = marker_label_gp,
      top_annotation = tmr_avg_anno,
      heatmap_legend_param = intensityLegendParam
    )

    png(
      filename = paste0(filenamePrefix, "_nonClusteringMarkers_arcsinh_weightedAvgByClust.png"),
      width = 2500,
      height = 2500 * (length(additionalMarkers)/length(markerOrder)),
      res = 300
    )
    print(tmr_hm)
    dev.off()
  }


  #########################################################################
  # Section 3.e from original SOP - export data
  #########################################################################
}
