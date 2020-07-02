
#   ### ln 225-250
#
#   # Set up a palette for the heatmap
#   my_zscore_pal <- colorRamp2(c(-2,0,2), c("#00ffffFF", "#000000FF", "#FDE725FF"))
#   my_arcsinh_pal <- colorRamp2(c(0,5), c("#440154FF", "#FDE725FF"))
#
#   # get subject ids, assign colors to them
#   subject_ids <- sort(unique(RPsubsetEventCounting$samp))
#   subject_id_colors <- data.frame(subject = subject_ids,
#                                   #color = rainbow(length(subject_ids)))
#                                   color = c("#000000FF",
#                                             plasma(length(subject_ids)+2)[4:(length(subject_ids)+2)]))
#
#   # pass_threshold_colors <- data.frame(pass = c(T,F), color = c("#FF0000", "#888888"))
#   # pass_threshold_colors <- data.frame(pass = c("no", "yes", "enriched"), color = c("#888888", "#FF0000", "#0000FF"))
#   # pass_threshold_colors <- data.frame(pass = c(0:8), color = c("#FFFFFF", "#3F054E", "#409E9B", "#CC0101", "#3F054E", "#409E9B", "#CC0101", "#5DC762", "#FFFFFF"))
#   #
#   # duration_colors <- data.frame(status = c("LS", "RO"), color = c("#2D297A", "#BB9C3E"))
#   # progression_colors <- data.frame(status = c("Pre-slow", "Slow", "Unclassified", "Rapid", "Post-rapid"),
#   #                                  color = c("#184349", "#368BA1", "#D3D3D3", "#5E9130", "#21632E"))
#
#   # Use Hannah/Virginia's color palette for matching to cluster colors
#   cb_pal <- colorblind_pal()(8)
#   vm_pal <- cb_pal[-1] # This is Hannah's favorite color scheme
#
#   #########################################################################
#   # Section 3.d from original SOP - create heatmaps
#   #########################################################################
#
#   # TODO: take n_metaclusters as argument
#   # TODO: check color palette
#
#
#   all_tmr_zscore_anno_df <-
#     hmapDfAllSubsets %>%
#     select(sample, !!subsets) %>%
#     unique %>%
#     slice(match(colnames(allSubsetAllSubjectZscores), .$sample)) %>%
#     mutate(sample = str_replace(.$sample, "_[0-9]+$", "")) %>%
#     rename(subject = sample) # %>%
#   # left_join(select(clinical_vars, ID, progression = `Progression rate`, duration = Duration), by = c("subject" = "ID"))
#
#   my_hm <- Heatmap(allSubsetAllSubjectZscores,
#                    clustering_method_columns = linkage,
#                    clustering_method_rows = linkage)
#
#   # cut the dendrogram to get phenotypic metaclusters
#   n_metaclusters = 12
#   col_clust <- as.hclust(column_dend(my_hm))
#   col_indices <- cutree(col_clust, k = n_metaclusters)
#   ordered_indices <- col_indices[column_order(my_hm)]
#   k_groups <- length(unique(col_indices))
#
#   plot_tmrs = c(#"A1_Islet", "A1_Insulin", "A1_Virus",
#     "A2_Islet", "A2_Insulin", "A2_Virus") # , "Insulin", "Virus"
#
#   # set up colors as a named list
#   all_tmr_zscore_anno_colors <- list(subject = setNames(subject_id_colors$color,
#                                                         as.character(subject_id_colors$subject)),
#                                      group = setNames(brewer.pal(k_groups, "Set3"),  #for 1-12MCs
#                                                       as.character(unique(col_indices))))
#   #group = setNames(c(brewer.pal(12, "Set3"), brewer.pal(k_groups-12, "Set2")), #for 12-20MCs
#   #         as.character(unique(col_indices)))) #,
#   # group = setNames(viridis(k_groups),
#   #          as.character(unique(col_indices)))) # original viridis scheme
#
#   # add metacluster information to tmr counting information
#   subsetEventCounting = subsetEventCounting %>%
#     left_join(data.frame(sample = colnames(allSubsetAllSubjectZscores), group = col_indices))
#
#   ## Use for presence/enrichment of tmrs color scheme
#   # for (t in plot_tmrs){
#   #   all_tmr_zscore_anno_colors[[t]] = setNames(pass_threshold_colors$color, pass_threshold_colors$pass)
#   # }
#
#   ## Use for number tmrs color scheme
#   # for (t in plot_tmrs){
#   #   tmr_vect = all_tmr_zscore_anno_df %>% select(t) %>% unlist(use.names = FALSE)
#   #   all_tmr_zscore_anno_colors[[t]] = c(setNames("#888888", 0), setNames(colorRampPalette(brewer.pal(9, "Greens")[2:9])(max(tmr_vect)),
#   #                             1:(max(tmr_vect))))
#   # }
#   # all_tmr_zscore_anno_df_counting = data.frame(subject = all_tmr_zscore_anno_df[,1], group = col_indices,
#   #                                     all_tmr_zscore_anno_df[,2:ncol(all_tmr_zscore_anno_df)])
#
#   ## Use for pct tmrs color scheme
#   for (t in plot_tmrs){
#     tmr_vect = all_tmr_zscore_anno_df %>% select(t) %>% unlist(use.names = FALSE)
#     # all_tmr_zscore_anno_colors[[t]] = c(setNames("white", 0),
#     #                                     setNames(colorRampPalette(c("grey90", "grey30"))(19), seq(0.01, 0.19, 0.01)),
#     #                                     setNames(rep("black", 81), seq(0.20, 1.00, 0.01)))
#     all_tmr_zscore_anno_colors[[t]] = c(setNames("white", 0),
#                                         setNames(colorRampPalette(c("grey95", "grey0"))(100), seq(0.01, 1.00, 0.01)))
#   }
#
#
#   all_tmr_zscore_anno_df %<>% mutate(group = col_indices) %>%
#     select(subject, group, #A1_Islet, A1_Insulin, A1_Virus,
#            A2_Islet, A2_Insulin, A2_Virus)#, Islet
#
#   all_tmr_zscore_anno <-
#     HeatmapAnnotation(
#       df = data.frame(all_tmr_zscore_anno_df %>% mutate_all(funs(ifelse(is.na(.), 0, .)))),
#       col = all_tmr_zscore_anno_colors,
#       show_annotation_name = T,
#       show_legend = F)
#   # annotation_legend_param = list(Islet = list(ncol = 4, title_position = "topcenter"),
#   #                                group = list(ncol = 2, title_position = "topcenter"),
#   #                                subject = list(ncol = 4, title_position = "topcenter")),
#   # show_legend = c(F, T, T, T, T))
#
#   zscore_hmap <- Heatmap(allSubsetAllSubjectZscores,
#                          col = my_zscore_pal,
#                          name = "z-score",
#                          column_title = paste0("Cluster Phenotypes from All Subjects"),
#                          column_title_gp = title_font_gp,
#                          #column_dend_height=unit(10,"cm"),
#                          clustering_method_columns = linkage,
#                          cluster_rows = T,
#                          clustering_method_rows = linkage,
#                          show_column_names = F,
#                          row_names_gp = marker_label_gp,
#                          top_annotation = all_tmr_zscore_anno,
#                          heatmap_legend_param = legend_params)
#
#   png(filename = paste0(fname_prefix_hmap, "_", linkage, "_", distance, "_",
#                         "_AllClusters_cutree", n_metaclusters, "_pct_Islet.png"),
#       width = export_width,
#       height = export_height)
#   print(zscore_hmap)
#   dev.off()
#
#   marker_order = row_order(zscore_hmap)#[[1]]
#
#
#   #MFI plot did not work
#   png(filename = paste0(fname_prefix_hmap, "_", linkage, "_", distance, "_",
#                         "_AllClusters_cutree", n_metaclusters, "_pct_tmrs_MFI.png"),
#       width = export_width,
#       height = export_height)
#   print(Heatmap(allSubsetAllSubjectArcsinh,
#                 col = my_arcsinh_pal,
#                 name = "MFI",
#                 column_title = paste0("Cluster Phenotypes from All Subjects"),
#                 column_title_gp = title_font_gp,
#                 cluster_columns = F,
#                 column_order = column_order(zscore_hmap),
#                 cluster_rows = F,
#                 row_order = c(markers[marker_order], "GRZMB", "CD160"),
#                 split = rep(c("", " "), c(length(marker_order), 2)),
#                 gap = unit(5, "mm"),
#                 show_column_names = F,
#                 row_names_gp = marker_label_gp,
#                 top_annotation = all_tmr_zscore_anno,
#                 heatmap_legend_param = legend_params))
#   dev.off()
#
#
#   # Weighted Average phenotype hmap - average of arcsinch values in each metacluster
#   per_metax_avg <-
#     data.frame(t(allSubsetAllSubjectArcsinh)) %>%
#     mutate(sample = rownames(.)) %>%
#     merge(subsetEventCounting, by = "sample") %>%
#     mutate(n_event = #A1_CD8 + A1_Islet + A1_Insulin + A1_Virus +
#              A2_CD8 + A2_Islet + A2_Insulin + A2_Virus,
#            subj = str_remove(sample, "_.*"),
#            RP_clust = str_remove(sample, ".*_")) %>%
#     group_by(subj) %>%
#     mutate(subj_event = sum(n_event)) %>%
#     ungroup() %>%
#     mutate(proportion = n_event/subj_event) %>%
#     mutate_at(markers, funs(.*proportion)) %>%
#     select(-sample, -!!subsets, -n_event, -subj, -RP_clust, -subj_event) %>%
#     group_by(group) %>%
#     summarise_all(sum) %>%
#     mutate_at(markers, funs(./proportion)) %>%
#     select(-group, -proportion) %>%
#     t
#
#   colnames(per_metax_avg) = 1:n_metaclusters
#
#   # set up data frame for annotation bar for average hmap
#   tmr_avg_anno_df <-
#     data.frame(group = colnames(per_metax_avg))
#
#   # set up colors as a named list
#   tmr_avg_anno_colors <- list(group = setNames(brewer.pal(k_groups, "Set3"), as.character(colnames(per_metax_avg)))) #1-12MCs
#   #group = setNames(c(brewer.pal(12, "Set3"), brewer.pal(k_groups-12, "Set2")),        #12-20MCs
#   #               as.character(colnames(per_metax_avg))))
#
#   tmr_avg_anno <- HeatmapAnnotation(df = tmr_avg_anno_df,
#                                     col = tmr_avg_anno_colors,
#                                     show_annotation_name = T,
#                                     show_legend = T)
#
#
#   #Did not work Error in .local(object, ...) : Number of rows in the matrix are not the same as the length of the cluster or the row orders.
#   tmr_hm <- Heatmap(per_metax_avg,
#                     col = my_arcsinh_pal,
#                     name = "MFI",
#                     column_title = paste0("Clusters Enriched for ", n_metaclusters, " Tmrs"),
#                     column_title_gp = title_font_gp,
#                     cluster_columns = F,
#                     cluster_rows = F,
#                     row_order = c(markers[marker_order], "GRZMB", "CD160"),
#                     split = rep(c("", " "), c(length(marker_order), 2)),
#                     gap = unit(5, "mm"),
#                     #combined_name_fun = NULL,
#                     show_column_names = F,
#                     row_names_gp = marker_label_gp,
#                     top_annotation = tmr_avg_anno,
#                     heatmap_legend_param = legend_params)
#
#
#   png(filename = paste0(fname_prefix_hmap, "_", linkage, "_rpheno_", n_metaclusters, "_arcsinh_weightedAvgByClust.png"),
#       width = 2500,
#       height = 2500,
#       res = 300)
#   print(tmr_hm)
#   dev.off()
#
#   ## Weighted Average phenotype for total CD8s
#   # weighted by proportion of cells in a subject (not by number of cells) in order to allow
#   # equivalent contributions from subjects with disparate numbers of collected events
#   cd8_avg <-
#     data.frame(t(allSubsetAllSubjectArcsinh)) %>%
#     mutate(sample = rownames(.)) %>%
#     merge(subsetEventCounting, by = "sample") %>%
#     mutate(n_event = #A1_CD8 + A1_Islet + A1_Insulin + A1_Virus +
#              A2_CD8 + A2_Islet + A2_Insulin + A2_Virus,
#            subj = str_remove(sample, "_.*"),
#            RP_clust = str_remove(sample, ".*_")) %>%
#     group_by(subj) %>%
#     mutate(subj_event = sum(n_event)) %>%
#     ungroup() %>%
#     mutate(proportion = n_event/subj_event) %>%
#     mutate_at(markers, funs(.*proportion)) %>%
#     select(-sample, -!!subsets, -n_event, -subj, -RP_clust, -subj_event) %>%
#     summarise_all(sum) %>%
#     mutate_at(markers, funs(./proportion)) %>%
#     select(-group, -proportion) %>%
#     t
#
#   colnames(cd8_avg) = "CD8"
#
#   #Did not work
#   cd8_hm <- Heatmap(cd8_avg,
#                     col = my_arcsinh_pal,
#                     name = "MFI",
#                     column_title = paste0("Total CD8+"),
#                     column_title_gp = title_font_gp,
#                     cluster_columns = F,
#                     cluster_rows = F,
#                     row_order = c(markers[marker_order], "GRZMB", "CD160"),
#                     split = rep(c("", " "), c(length(marker_order), 2)),
#                     gap = unit(5, "mm"),
#                     #combined_name_fun = NULL,
#                     show_column_names = F,
#                     row_names_gp = marker_label_gp,
#                     heatmap_legend_param = legend_params)
#
#   png(filename = paste0(fname_prefix_hmap, "_", linkage, "_CD8_arcsinh_weightedAvg.png"),
#       width = 700,
#       height = 2500,
#       res = 300)
#   print(cd8_hm)
#   dev.off()
#
#   #########################################################################
#   # Section 3.e from original SOP - export data
#   #########################################################################
# }
