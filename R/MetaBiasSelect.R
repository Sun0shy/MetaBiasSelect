#' MetaBiasSelect: helper to extract metadata
#'
#' Internal helper. Accepts either a data.frame-like metadata table or a Seurat-like
#' object containing a @meta.data slot.
#'
#' @noRd
.extract_meta <- function(object, sample_col, cancer_col, cluster_col) {
  if (is.data.frame(object)) {
    meta <- object
  } else {
    meta <- tryCatch(object@meta.data, error = function(e) NULL)
    if (is.null(meta)) {
      stop(
        "`object` must be either a metadata data.frame or a Seurat-like object with @meta.data."
      )
    }
  }

  meta <- tibble::as_tibble(meta)

  if (!cluster_col %in% colnames(meta)) {
    stop("`cluster_col` not found in metadata: ", cluster_col)
  }
  if (!cancer_col %in% colnames(meta)) {
    stop("`cancer_col` not found in metadata: ", cancer_col)
  }

  if (is.null(sample_col) || !sample_col %in% colnames(meta)) {
    sample_col <- cancer_col
  }

  meta_std <- meta %>%
    dplyr::transmute(
      sample  = as.character(.data[[sample_col]]),
      cancer  = as.character(.data[[cancer_col]]),
      cluster = as.character(.data[[cluster_col]])
    ) %>%
    dplyr::filter(
      !is.na(.data$sample),
      !is.na(.data$cancer),
      !is.na(.data$cluster)
    )

  list(meta = meta_std, sample_col_used = sample_col)
}

#' Build sample-by-cluster count matrix
#'
#' @noRd
.build_count_table <- function(meta) {
  tab <- table(meta$sample, meta$cluster)
  count_df <- as.data.frame.matrix(tab, stringsAsFactors = FALSE)
  count_df$sample <- rownames(count_df)
  rownames(count_df) <- NULL

  sample_meta <- meta %>%
    dplyr::select(.data$sample, .data$cancer) %>%
    dplyr::distinct()

  list(count_df = count_df, sample_meta = sample_meta)
}

#' Compute sample-level cluster proportions
#'
#' @noRd
.compute_sample_proportions <- function(count_df, sample_meta) {
  count_df %>%
    tidyr::pivot_longer(
      cols = -.data$sample,
      names_to = "cluster",
      values_to = "count"
    ) %>%
    dplyr::left_join(sample_meta, by = "sample") %>%
    dplyr::group_by(.data$sample) %>%
    dplyr::mutate(
      sample_total = sum(.data$count),
      prop = dplyr::case_when(
        .data$sample_total > 0 ~ .data$count / .data$sample_total,
        TRUE ~ 0
      )
    ) %>%
    dplyr::ungroup()
}

#' Compute cancer-level enrichment score
#'
#' @noRd
.compute_enrichment <- function(prop_df, pseudocount = 1e-6) {
  cancer_mean <- prop_df %>%
    dplyr::group_by(.data$cancer, .data$cluster) %>%
    dplyr::summarise(
      Observed = mean(.data$prop),
      .groups = "drop"
    )

  overall_prop <- prop_df %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::summarise(
      Expected = mean(.data$prop),
      .groups = "drop"
    )

  cancer_mean %>%
    dplyr::left_join(overall_prop, by = "cluster") %>%
    dplyr::mutate(
      log2OE = log2((.data$Observed + pseudocount) / (.data$Expected + pseudocount))
    )
}

#' Create an empty pairwise DA table
#'
#' @noRd
.empty_pairwise_da <- function() {
  tibble::tibble(
    cluster = character(),
    contrast = character(),
    cancer_high = character(),
    cancer_low = character(),
    logFC = numeric(),
    logCPM = numeric(),
    F = numeric(),
    PValue = numeric(),
    FDR = numeric()
  )
}

#' Run pairwise pseudobulk abundance testing across cancers
#'
#' This function treats each cluster as a feature and each sample as a replicate,
#' then runs pairwise edgeR quasi-likelihood tests across eligible cancer groups.
#' Cancers with fewer than `min_samples_per_group` samples are excluded from DA and
#' can still be retained later as candidate-biased clusters.
#'
#' @param count_df A sample-by-cluster count data.frame generated inside the pipeline.
#' @param sample_meta A data.frame with columns `sample` and `cancer`.
#' @param min_samples_per_group Minimum number of samples required in each cancer group.
#'
#' @return A data.frame with pairwise abundance statistics.
#' @export
run_pairwise_cluster_da <- function(count_df,
                                    sample_meta,
                                    min_samples_per_group = 3) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required. Please install it with BiocManager::install('edgeR').")
  }

  count_mat <- count_df[, !(colnames(count_df) %in% "sample"), drop = FALSE]
  rownames(count_mat) <- count_df$sample
  count_mat <- t(as.matrix(count_mat))

  sample_meta <- sample_meta[match(colnames(count_mat), sample_meta$sample), , drop = FALSE]
  sample_meta$cancer <- as.character(sample_meta$cancer)

  group_sizes <- table(sample_meta$cancer)
  eligible_cancers <- names(group_sizes)[group_sizes >= min_samples_per_group]

  if (length(eligible_cancers) < 2) {
    out <- .empty_pairwise_da()
    attr(out, "group_sizes") <- group_sizes
    attr(out, "eligible_cancers") <- eligible_cancers
    attr(out, "min_samples_per_group") <- min_samples_per_group
    return(out)
  }

  keep_samples <- sample_meta$cancer %in% eligible_cancers
  sample_meta_use <- sample_meta[keep_samples, , drop = FALSE]
  count_mat_use <- count_mat[, keep_samples, drop = FALSE]
  sample_meta_use$cancer <- factor(sample_meta_use$cancer, levels = eligible_cancers)

  dge <- edgeR::DGEList(counts = count_mat_use)
  dge$samples$group <- sample_meta_use$cancer
  dge <- edgeR::calcNormFactors(dge)

  design <- stats::model.matrix(~0 + group, data = dge$samples)
  colnames(design) <- sub("^group", "", colnames(design))

  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)

  cancers <- colnames(design)
  pair_mat <- utils::combn(cancers, 2)

  out <- purrr::map_dfr(seq_len(ncol(pair_mat)), function(i) {
    high <- pair_mat[1, i]
    low  <- pair_mat[2, i]

    contrast_vec <- rep(0, ncol(design))
    names(contrast_vec) <- colnames(design)
    contrast_vec[high] <- 1
    contrast_vec[low]  <- -1

    qlf <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
    tab <- edgeR::topTags(qlf, n = Inf)$table

    tibble::rownames_to_column(tab, var = "cluster") %>%
      dplyr::transmute(
        cluster = .data$cluster,
        contrast = paste0(high, "_vs_", low),
        cancer_high = high,
        cancer_low = low,
        logFC = .data$logFC,
        logCPM = .data$logCPM,
        F = .data$F,
        PValue = .data$PValue,
        FDR = .data$FDR
      )
  })

  attr(out, "group_sizes") <- group_sizes
  attr(out, "eligible_cancers") <- eligible_cancers
  attr(out, "min_samples_per_group") <- min_samples_per_group
  out
}

#' Summarise pairwise DA support for each cluster-cancer pair
#'
#' @noRd
.summarise_da_support <- function(pairwise_da,
                                  cancers,
                                  logfc_cutoff = 1,
                                  fdr_cutoff = 0.05) {
  if (nrow(pairwise_da) == 0) {
    return(tibble::tibble(
      cluster = character(),
      cancer = character(),
      da_support_n = integer(),
      da_tested_n = integer(),
      da_support_prop = numeric()
    ))
  }

  purrr::map_dfr(cancers, function(target) {
    tmp_high <- pairwise_da %>%
      dplyr::filter(.data$cancer_high == target) %>%
      dplyr::transmute(
        cluster = .data$cluster,
        cancer = target,
        support = .data$logFC >= logfc_cutoff & .data$FDR <= fdr_cutoff,
        comparison = .data$contrast
      )

    tmp_low <- pairwise_da %>%
      dplyr::filter(.data$cancer_low == target) %>%
      dplyr::transmute(
        cluster = .data$cluster,
        cancer = target,
        support = (-.data$logFC) >= logfc_cutoff & .data$FDR <= fdr_cutoff,
        comparison = .data$contrast
      )

    dplyr::bind_rows(tmp_high, tmp_low) %>%
      dplyr::group_by(.data$cluster, .data$cancer) %>%
      dplyr::summarise(
        da_support_n = sum(.data$support, na.rm = TRUE),
        da_tested_n = dplyr::n(),
        da_support_prop = dplyr::case_when(
          .data$da_tested_n > 0 ~ .data$da_support_n / .data$da_tested_n,
          TRUE ~ 0
        ),
        .groups = "drop"
      )
  })
}

#' Classify cancer-biased and shared clusters
#'
#' This function integrates three pieces of evidence for each cluster-cancer pair:
#' enrichment, cancer-level mean proportion, and pairwise pseudobulk abundance support.
#' For cancers with insufficient sample replication for DA, clusters can still be
#' labeled as candidate-biased if enrichment and proportion criteria are both met.
#'
#' @param evidence_table A cluster-cancer evidence table generated inside the pipeline.
#' @param enrich_cutoff Minimum `log2(O/E)` required to support cancer bias.
#' @param prop_cutoff Minimum cancer-level mean proportion required to support cancer bias.
#' @param min_rules Minimum number of criteria that must be satisfied.
#' @param da_min_support Minimum number of supporting pairwise DA comparisons.
#'   If `NULL`, all tested other cancers must be outperformed.
#' @param shared_label Label assigned to clusters not selected as biased.
#'
#' @return A list with `evidence_table` and `cluster_summary`.
#' @export
classify_metabias_clusters <- function(evidence_table,
                                       enrich_cutoff = 0.3,
                                       prop_cutoff = 0.1,
                                       min_rules = 2,
                                       da_min_support = NULL,
                                       shared_label = "Shared") {

  evidence_out <- evidence_table

  evidence_out <- evidence_out %>%
    dplyr::mutate(
      cond_enrich = .data$log2OE >= enrich_cutoff,
      cond_prop   = .data$prop >= prop_cutoff
    )

  if (is.null(da_min_support)) {
    evidence_out <- evidence_out %>%
      dplyr::mutate(
        cond_da = dplyr::if_else(
          .data$da_available,
          .data$da_support_n >= .data$da_tested_n,
          FALSE
        )
      )
  } else {
    evidence_out <- evidence_out %>%
      dplyr::mutate(
        cond_da = dplyr::if_else(
          .data$da_available,
          .data$da_support_n >= pmin(da_min_support, .data$da_tested_n),
          FALSE
        )
      )
  }

  evidence_out <- evidence_out %>%
    dplyr::mutate(
      cond_sum = (.data$cond_enrich * 1L) +
        (.data$cond_prop * 1L) +
        (.data$cond_da * 1L),

      selected = .data$da_available & (.data$cond_sum >= min_rules),

      candidate_selected = (!.data$da_available) &
        .data$cond_enrich &
        .data$cond_prop,

      selection_status = dplyr::case_when(
        .data$selected ~ "Biased",
        .data$candidate_selected ~ "Candidate_biased",
        TRUE ~ shared_label
      ),

      plot_symbol = dplyr::case_when(
        .data$selected ~ "*",
        .data$candidate_selected ~ "†",
        TRUE ~ ""
      )
    )

  cluster_summary <- evidence_out %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::arrange(
      dplyr::desc(.data$selected),
      dplyr::desc(.data$candidate_selected),
      dplyr::desc(.data$cond_sum),
      dplyr::desc(.data$log2OE),
      dplyr::desc(.data$prop),
      dplyr::desc(.data$da_support_n),
      .by_group = TRUE
    ) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      classification = dplyr::case_when(
        .data$selected ~ .data$cancer,
        .data$candidate_selected ~ .data$cancer,
        TRUE ~ shared_label
      ),
      decision_reason = dplyr::case_when(
        .data$selected ~ paste0(
          "Selected by ", .data$cancer,
          " (", .data$cond_sum, "/3 rules passed; DA available)"
        ),
        .data$candidate_selected ~ paste0(
          "Candidate-biased in ", .data$cancer,
          " (enrichment + proportion passed; DA unavailable due to insufficient replicates)"
        ),
        TRUE ~ paste0(shared_label, ": no cancer met the minimum rule threshold")
      )
    ) %>%
    dplyr::select(
      .data$cluster,
      assigned_cancer = .data$cancer,
      .data$classification,
      .data$selection_status,
      .data$plot_symbol,
      .data$log2OE,
      .data$prop,
      .data$da_available,
      .data$da_support_n,
      .data$da_tested_n,
      .data$da_support_prop,
      .data$cond_enrich,
      .data$cond_prop,
      .data$cond_da,
      .data$cond_sum,
      .data$selected,
      .data$candidate_selected,
      .data$decision_reason
    )

  list(
    evidence_table = evidence_out,
    cluster_summary = cluster_summary,
    da_min_support = da_min_support
  )
}

#' Build score matrix for heatmap
#'
#' @noRd
.build_score_matrix <- function(enrich_df, cancers) {
  score_matrix <- enrich_df %>%
    dplyr::select(.data$cancer, .data$cluster, .data$log2OE) %>%
    tidyr::pivot_wider(
      names_from = .data$cancer,
      values_from = .data$log2OE,
      names_prefix = "score_"
    )

  for (cn in cancers) {
    nm <- paste0("score_", cn)
    if (!nm %in% colnames(score_matrix)) {
      score_matrix[[nm]] <- 0
    }
  }

  score_matrix <- score_matrix %>%
    dplyr::select(.data$cluster, dplyr::all_of(paste0("score_", cancers)))

  score_cols <- paste0("score_", cancers)

  score_matrix$shared_score <- apply(
    score_matrix[, score_cols, drop = FALSE],
    1,
    function(x) -stats::sd(x, na.rm = TRUE)
  )

  score_matrix
}

#' Annotate high-confidence shared clusters
#'
#' @noRd
.annotate_shared_clusters <- function(cluster_summary,
                                      evidence_table,
                                      score_matrix,
                                      cancers,
                                      shared_score_cutoff = -1.0,
                                      shared_prop_cutoff = 0.03,
                                      shared_min_cancers = 2,
                                      shared_symbol = "•") {

  prop_wide <- evidence_table %>%
    dplyr::select(.data$cluster, .data$cancer, .data$prop) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = .data$cancer,
      values_from = .data$prop,
      names_prefix = "prop_"
    )

  for (cn in cancers) {
    nm <- paste0("prop_", cn)
    if (!nm %in% colnames(prop_wide)) {
      prop_wide[[nm]] <- 0
    }
  }

  prop_cols <- paste0("prop_", cancers)

  prop_wide <- prop_wide %>%
    dplyr::mutate(
      mean_prop = rowMeans(dplyr::across(dplyr::all_of(prop_cols)), na.rm = TRUE),
      n_cancers_prop_ge = rowSums(
        as.data.frame(dplyr::across(dplyr::all_of(prop_cols), ~ .x >= shared_prop_cutoff)),
        na.rm = TRUE
      )
    ) %>%
    dplyr::select(.data$cluster, .data$mean_prop, .data$n_cancers_prop_ge)

  out <- cluster_summary %>%
    dplyr::select(-dplyr::any_of(c(
      "shared_score", "mean_prop", "n_cancers_prop_ge",
      "shared_selected", "plot_symbol_shared"
    ))) %>%
    dplyr::left_join(
      score_matrix %>% dplyr::select(.data$cluster, .data$shared_score),
      by = "cluster"
    ) %>%
    dplyr::left_join(prop_wide, by = "cluster") %>%
    dplyr::mutate(
      mean_prop = tidyr::replace_na(.data$mean_prop, 0),
      n_cancers_prop_ge = tidyr::replace_na(.data$n_cancers_prop_ge, 0L),
      shared_selected = (.data$classification == "Shared") &
        (.data$shared_score >= shared_score_cutoff) &
        (.data$mean_prop >= shared_prop_cutoff) &
        (.data$n_cancers_prop_ge >= shared_min_cancers),
      plot_symbol_shared = dplyr::if_else(.data$shared_selected, shared_symbol, ""),
      plot_symbol = dplyr::case_when(
        .data$plot_symbol != "" ~ .data$plot_symbol,
        .data$shared_selected ~ shared_symbol,
        TRUE ~ ""
      ),
      decision_reason = dplyr::case_when(
        .data$shared_selected ~ paste0(
          "High-confidence Shared (shared_score >= ", shared_score_cutoff,
          ", mean_prop >= ", shared_prop_cutoff,
          ", present in >= ", shared_min_cancers, " cancers)"
        ),
        TRUE ~ .data$decision_reason
      )
    )

  out
}

#' Plot integrated heatmap for cancer bias and shared stability
#'
#' @param result Output from [run_metabias_pipeline()].
#' @param cluster_order Optional manual cluster order.
#' @param class_order Optional order of classes on y-axis grouping logic.
#' @param shared_low Low color for shared stability.
#' @param shared_high High color for shared stability.
#' @param bias_low Low color for cancer bias.
#' @param bias_mid Mid color for cancer bias.
#' @param bias_high High color for cancer bias.
#' @param bias_limits Optional numeric length-2 vector controlling display range of bias score.
#' @param show_symbols Whether to overlay selection symbols on tiles.
#' @param selected_symbol Symbol used for fully supported biased clusters.
#' @param candidate_symbol Symbol used for candidate-biased clusters lacking DA.
#' @param shared_symbol Symbol used for high-confidence shared clusters.
#' @param symbol_size Text size for symbols.
#' @param symbol_color Text color for symbols.
#' @param shared_score_cutoff Cutoff for high-confidence shared clusters.
#' @param shared_prop_cutoff Minimum mean proportion for shared clusters.
#' @param shared_min_cancers Minimum number of cancers passing proportion cutoff.
#' @param axis_title_size Axis title size.
#' @param axis_text_x_size X-axis text size.
#' @param axis_text_y_size Y-axis text size.
#' @param legend_title_size Legend title size.
#' @param legend_text_size Legend text size.
#' @param axis_line_width Axis line width.
#' @param tile_line_width Tile border width.
#'
#' @return A ggplot object.
#' @export
plot_metabias_heatmap <- function(result,
                                  cluster_order = NULL,
                                  class_order = NULL,
                                  shared_low = "#f7fcf5",
                                  shared_high = "#31a354",
                                  bias_low = "#3B4CC0",
                                  bias_mid = "#F7F7F7",
                                  bias_high = "#B40426",
                                  bias_limits = NULL,
                                  show_symbols = TRUE,
                                  selected_symbol = "*",
                                  candidate_symbol = "†",
                                  shared_symbol = "•",
                                  symbol_size = 5,
                                  symbol_color = "black",
                                  shared_score_cutoff = -1.0,
                                  shared_prop_cutoff = 0.03,
                                  shared_min_cancers = 2,
                                  axis_title_size = 16,
                                  axis_text_x_size = 13,
                                  axis_text_y_size = 12,
                                  legend_title_size = 12,
                                  legend_text_size = 11,
                                  axis_line_width = 0.8,
                                  tile_line_width = 0.8) {

  score_matrix <- result$score_matrix
  cluster_summary <- result$cluster_summary
  evidence_table <- result$evidence_table
  cancers <- result$params$cancers

  cancer_df <- score_matrix %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("score_"),
      names_to = "cancer",
      values_to = "score"
    ) %>%
    dplyr::mutate(cancer = gsub("^score_", "", .data$cancer))

  cluster_summary2 <- .annotate_shared_clusters(
    cluster_summary = cluster_summary,
    evidence_table = evidence_table,
    score_matrix = score_matrix,
    cancers = cancers,
    shared_score_cutoff = shared_score_cutoff,
    shared_prop_cutoff = shared_prop_cutoff,
    shared_min_cancers = shared_min_cancers,
    shared_symbol = shared_symbol
  )

  shared_df <- score_matrix %>%
    dplyr::select(.data$cluster, .data$shared_score) %>%
    dplyr::left_join(
      cluster_summary2 %>%
        dplyr::select(.data$cluster, .data$plot_symbol_shared),
      by = "cluster"
    ) %>%
    dplyr::mutate(
      plot_symbol_shared = tidyr::replace_na(.data$plot_symbol_shared, "")
    )

  ord_df <- cluster_summary2 %>%
    dplyr::select(-dplyr::any_of("shared_score")) %>%
    dplyr::left_join(
      score_matrix %>% dplyr::select(.data$cluster, .data$shared_score),
      by = "cluster"
    ) %>%
    dplyr::mutate(
      class_for_order = .data$classification,
      primary_score = dplyr::case_when(
        .data$classification == "Shared" ~ .data$shared_score,
        TRUE ~ .data$log2OE
      )
    )

  if (is.null(class_order)) {
    class_order <- c("Shared", cancers)
  }

  if (is.null(cluster_order)) {
    cluster_order <- ord_df %>%
      dplyr::mutate(class_for_order = factor(.data$class_for_order, levels = class_order)) %>%
      dplyr::arrange(
        .data$class_for_order,
        dplyr::desc(.data$primary_score),
        dplyr::desc(.data$cond_sum),
        .data$cluster
      ) %>%
      dplyr::pull(.data$cluster)
  }

  cluster_levels <- rev(unique(cluster_order))

  cancer_df$cluster <- factor(cancer_df$cluster, levels = cluster_levels)
  shared_df$cluster <- factor(shared_df$cluster, levels = cluster_levels)

  symbol_df <- evidence_table %>%
    dplyr::mutate(
      plot_symbol = dplyr::case_when(
        .data$selected ~ selected_symbol,
        .data$candidate_selected ~ candidate_symbol,
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(.data$cluster, .data$cancer, .data$plot_symbol)

  cancer_df <- cancer_df %>%
    dplyr::left_join(symbol_df, by = c("cluster", "cancer")) %>%
    dplyr::mutate(plot_symbol = tidyr::replace_na(.data$plot_symbol, ""))

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = cancer_df,
      ggplot2::aes(x = .data$cancer, y = .data$cluster, fill = .data$score),
      color = "white",
      linewidth = tile_line_width
    )

  if (isTRUE(show_symbols)) {
    p <- p +
      ggplot2::geom_text(
        data = cancer_df %>% dplyr::filter(.data$plot_symbol != ""),
        ggplot2::aes(x = .data$cancer, y = .data$cluster, label = .data$plot_symbol),
        size = symbol_size,
        color = symbol_color,
        fontface = "bold"
      )
  }

  if (is.null(bias_limits)) {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = bias_low,
        mid = bias_mid,
        high = bias_high,
        midpoint = 0,
        name = "Cancer bias score"
      )
  } else {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = bias_low,
        mid = bias_mid,
        high = bias_high,
        midpoint = 0,
        limits = bias_limits,
        oob = scales::squish,
        name = "Cancer bias score"
      )
  }

  p <- p +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      data = shared_df,
      ggplot2::aes(x = "Shared", y = .data$cluster, fill = .data$shared_score),
      color = "white",
      linewidth = tile_line_width
    )

  if (isTRUE(show_symbols)) {
    p <- p +
      ggplot2::geom_text(
        data = shared_df %>% dplyr::filter(.data$plot_symbol_shared != ""),
        ggplot2::aes(x = "Shared", y = .data$cluster, label = .data$plot_symbol_shared),
        size = symbol_size,
        color = symbol_color,
        fontface = "bold"
      )
  }

  p <- p +
    ggplot2::scale_fill_gradient(
      low = shared_low,
      high = shared_high,
      name = "Shared stability"
    ) +
    ggplot2::scale_x_discrete(
      limits = c("Shared", cancers),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_discrete(
      expand = c(0, 0)
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.title.y = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.text.x  = ggplot2::element_text(size = axis_text_x_size, angle = 0, hjust = 0.5, color = "black"),
      axis.text.y  = ggplot2::element_text(size = axis_text_y_size, color = "black"),
      axis.line.x  = ggplot2::element_line(color = "black", linewidth = axis_line_width),
      axis.line.y  = ggplot2::element_line(color = "black", linewidth = axis_line_width),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = axis_line_width),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = axis_line_width),
      axis.ticks.length = grid::unit(0.2, "cm"),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.text = ggplot2::element_text(size = legend_text_size),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = "Cancer type",
      y = "Cell cluster"
    )

  p
}

#' Run the standard workflow for multi-cancer metastatic bias cluster selection
#'
#' This is the main user-facing function of the package. It accepts a Seurat-like
#' object or a metadata table, computes sample-level cluster abundance, enrichment,
#' pairwise edgeR abundance tests (when possible), integrated cluster classification,
#' and the final heatmap. Cancers lacking sufficient biological replicates are
#' automatically downgraded to candidate-only status rather than stopping the workflow.
#'
#' @param object A Seurat-like object with `@meta.data`, or a metadata data.frame.
#' @param cluster_col Column name containing cluster annotation.
#' @param cancer_col Column name containing cancer type.
#' @param sample_col Column name containing sample ID. If `NULL` or missing, `cancer_col`
#'   will be used, but sample-level replication is strongly recommended.
#' @param cancers Optional vector specifying which cancer groups to include and their order.
#' @param enrich_cutoff Enrichment threshold for `log2(O/E)`.
#' @param prop_cutoff Mean cancer-level proportion threshold.
#' @param logfc_cutoff Pairwise logFC threshold for DA support.
#' @param fdr_cutoff Pairwise FDR threshold for DA support.
#' @param min_rules Minimum number of rules passed among enrichment, proportion, and DA.
#' @param da_min_support Minimum number of pairwise DA comparisons supporting the target cancer.
#'   If `NULL`, all tested other cancers must be outperformed.
#' @param pseudocount Small constant used in enrichment calculation.
#' @param min_samples_per_group Minimum number of samples per cancer group for edgeR.
#' @param cluster_order Optional manual cluster order for the final heatmap.
#' @param return_plot Whether to generate a ggplot heatmap.
#' @param shared_score_cutoff Cutoff for high-confidence shared clusters.
#' @param shared_prop_cutoff Minimum mean proportion for shared clusters.
#' @param shared_min_cancers Minimum number of cancers passing proportion cutoff.
#' @param shared_symbol Symbol used for high-confidence shared clusters.
#'
#' @return A list containing parameters, intermediate results, summary tables, and plot.
#' @export
run_metabias_pipeline <- function(object,
                                  cluster_col,
                                  cancer_col,
                                  sample_col = NULL,
                                  cancers = NULL,
                                  enrich_cutoff = 0.3,
                                  prop_cutoff = 0.1,
                                  logfc_cutoff = 1,
                                  fdr_cutoff = 0.05,
                                  min_rules = 2,
                                  da_min_support = NULL,
                                  pseudocount = 1e-6,
                                  min_samples_per_group = 3,
                                  cluster_order = NULL,
                                  return_plot = TRUE,
                                  shared_score_cutoff = -1.0,
                                  shared_prop_cutoff = 0.03,
                                  shared_min_cancers = 2,
                                  shared_symbol = "•") {

  meta_obj <- .extract_meta(
    object = object,
    sample_col = sample_col,
    cancer_col = cancer_col,
    cluster_col = cluster_col
  )
  meta <- meta_obj$meta

  if (!is.null(cancers)) {
    meta <- meta %>% dplyr::filter(.data$cancer %in% cancers)
    meta$cancer <- factor(meta$cancer, levels = cancers) %>% as.character()
  } else {
    cancers <- unique(meta$cancer)
  }

  if (length(unique(meta$cancer)) < 2) {
    stop("At least two cancer groups are required.")
  }

  count_obj <- .build_count_table(meta)
  count_df <- count_obj$count_df
  sample_meta <- count_obj$sample_meta

  sample_counts <- sample_meta %>%
    dplyr::count(.data$cancer, name = "n_samples") %>%
    dplyr::mutate(da_eligible = .data$n_samples >= min_samples_per_group)

  eligible_cancers <- sample_counts %>%
    dplyr::filter(.data$da_eligible) %>%
    dplyr::pull(.data$cancer)

  prop_df <- .compute_sample_proportions(count_df, sample_meta)
  enrich_df <- .compute_enrichment(prop_df, pseudocount = pseudocount)

  prop_cancer <- prop_df %>%
    dplyr::group_by(.data$cancer, .data$cluster) %>%
    dplyr::summarise(
      prop = mean(.data$prop),
      .groups = "drop"
    )

  pairwise_da <- run_pairwise_cluster_da(
    count_df = count_df,
    sample_meta = sample_meta,
    min_samples_per_group = min_samples_per_group
  )

  da_support <- .summarise_da_support(
    pairwise_da = pairwise_da,
    cancers = cancers,
    logfc_cutoff = logfc_cutoff,
    fdr_cutoff = fdr_cutoff
  )

  evidence_table <- enrich_df %>%
    dplyr::left_join(prop_cancer, by = c("cancer", "cluster")) %>%
    dplyr::left_join(da_support, by = c("cluster", "cancer")) %>%
    dplyr::left_join(sample_counts, by = "cancer") %>%
    dplyr::mutate(
      prop = tidyr::replace_na(.data$prop, 0),
      da_support_n = tidyr::replace_na(.data$da_support_n, 0L),
      da_tested_n = tidyr::replace_na(.data$da_tested_n, 0L),
      da_support_prop = tidyr::replace_na(.data$da_support_prop, 0),
      n_samples = tidyr::replace_na(.data$n_samples, 0L),
      da_eligible = tidyr::replace_na(.data$da_eligible, FALSE),
      da_available = .data$da_eligible & (.data$da_tested_n > 0)
    )

  class_obj <- classify_metabias_clusters(
    evidence_table = evidence_table,
    enrich_cutoff = enrich_cutoff,
    prop_cutoff = prop_cutoff,
    min_rules = min_rules,
    da_min_support = da_min_support,
    shared_label = "Shared"
  )

  score_matrix <- .build_score_matrix(enrich_df, cancers = cancers)

  cluster_summary <- class_obj$cluster_summary %>%
    dplyr::left_join(
      score_matrix %>% dplyr::select(.data$cluster, .data$shared_score),
      by = "cluster"
    )

  cluster_summary <- .annotate_shared_clusters(
    cluster_summary = cluster_summary,
    evidence_table = class_obj$evidence_table,
    score_matrix = score_matrix,
    cancers = cancers,
    shared_score_cutoff = shared_score_cutoff,
    shared_prop_cutoff = shared_prop_cutoff,
    shared_min_cancers = shared_min_cancers,
    shared_symbol = shared_symbol
  )

  biased_clusters_table <- cluster_summary %>%
    dplyr::filter(.data$selection_status %in% c("Biased", "Candidate_biased")) %>%
    dplyr::arrange(
      dplyr::desc(.data$selection_status == "Biased"),
      .data$assigned_cancer,
      dplyr::desc(.data$log2OE)
    )

  shared_clusters_table <- cluster_summary %>%
    dplyr::filter(.data$shared_selected) %>%
    dplyr::arrange(dplyr::desc(.data$shared_score))

  params <- list(
    cluster_col = cluster_col,
    cancer_col = cancer_col,
    sample_col = meta_obj$sample_col_used,
    cancers = cancers,
    enrich_cutoff = enrich_cutoff,
    prop_cutoff = prop_cutoff,
    logfc_cutoff = logfc_cutoff,
    fdr_cutoff = fdr_cutoff,
    min_rules = min_rules,
    da_min_support = class_obj$da_min_support,
    pseudocount = pseudocount,
    min_samples_per_group = min_samples_per_group,
    eligible_cancers = eligible_cancers,
    shared_score_cutoff = shared_score_cutoff,
    shared_prop_cutoff = shared_prop_cutoff,
    shared_min_cancers = shared_min_cancers,
    shared_symbol = shared_symbol
  )

  result <- list(
    params = params,
    meta = meta,
    sample_meta = sample_meta,
    cancer_sample_counts = sample_counts,
    count_df = count_df,
    sample_proportions = prop_df,
    enrichment_table = enrich_df,
    pairwise_da = pairwise_da,
    evidence_table = class_obj$evidence_table,
    cluster_summary = cluster_summary,
    biased_clusters_table = biased_clusters_table,
    shared_clusters_table = shared_clusters_table,
    score_matrix = score_matrix,
    plot = NULL
  )

  if (isTRUE(return_plot)) {
    result$plot <- plot_metabias_heatmap(
      result = result,
      cluster_order = cluster_order,
      shared_score_cutoff = shared_score_cutoff,
      shared_prop_cutoff = shared_prop_cutoff,
      shared_min_cancers = shared_min_cancers,
      shared_symbol = shared_symbol
    )
  }

  class(result) <- c("metabias_result", class(result))
  result
}

#' Save standard outputs from the MetaBiasSelect workflow
#'
#' @param result Output from [run_metabias_pipeline()].
#' @param outdir Output directory.
#' @param prefix Prefix for exported files.
#' @param width Plot width.
#' @param height Plot height.
#'
#' @return Invisibly returns the output directory path.
#' @export
save_metabias_outputs <- function(result,
                                  outdir = "MetaBiasSelect_output",
                                  prefix = "metabias",
                                  width = 8,
                                  height = 6) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  param_table <- tibble::tibble(
    parameter = names(result$params),
    value = vapply(result$params, function(x) paste(x, collapse = ","), character(1))
  )

  utils::write.csv(
    param_table,
    file = file.path(outdir, paste0(prefix, "_parameters.csv")),
    row.names = FALSE
  )

  if (!is.null(result$cancer_sample_counts)) {
    utils::write.csv(
      result$cancer_sample_counts,
      file = file.path(outdir, paste0(prefix, "_cancer_sample_counts.csv")),
      row.names = FALSE
    )
  }

  utils::write.csv(
    result$cluster_summary,
    file = file.path(outdir, paste0(prefix, "_cluster_summary.csv")),
    row.names = FALSE
  )

  if (!is.null(result$biased_clusters_table)) {
    utils::write.csv(
      result$biased_clusters_table,
      file = file.path(outdir, paste0(prefix, "_biased_clusters_table.csv")),
      row.names = FALSE
    )
  }

  if (!is.null(result$shared_clusters_table)) {
    utils::write.csv(
      result$shared_clusters_table,
      file = file.path(outdir, paste0(prefix, "_shared_clusters_table.csv")),
      row.names = FALSE
    )
  }

  utils::write.csv(
    result$evidence_table,
    file = file.path(outdir, paste0(prefix, "_cluster_cancer_evidence.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    result$pairwise_da,
    file = file.path(outdir, paste0(prefix, "_pairwise_DA.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    result$score_matrix,
    file = file.path(outdir, paste0(prefix, "_score_matrix.csv")),
    row.names = FALSE
  )

  if (!is.null(result$plot)) {
    ggplot2::ggsave(
      filename = file.path(outdir, paste0(prefix, "_heatmap.pdf")),
      plot = result$plot,
      width = width,
      height = height
    )
    ggplot2::ggsave(
      filename = file.path(outdir, paste0(prefix, "_heatmap.png")),
      plot = result$plot,
      width = width,
      height = height,
      dpi = 300
    )
  }

  saveRDS(result, file = file.path(outdir, paste0(prefix, "_full_result.rds")))

  invisible(outdir)
}