#' Class markers
#'
#' This function identifies both on and off markers genes for classes for use with patchSeqQC library.
#' 'On markers, are genes that are highly expressed in the cell type of interest with enriched expression
#' relative to other cell types. The second class, Off markers, are expected to be expressed at low levels
#' in a given patch-seq cell type.'  Note that these markers are based on a relevant reference data set.
#'
#' @param datRef a matrix (rows=genes x columns=samples) of gene expression data (e.g., scRNA-seq)
#' @param onClasses a character (or factor) vector indicated the on class for each sample in datRef
#' @param offClasses a character (or factor) vector indicated the off class for each sample in datRef
#' @param numMarkers number of markers per class to return (default = 50)
#'
#' @return a 3 x count matrix of the top confused pairs of clusters with the three columns corresponding
#'   to mapped cluster, assigned cluster, and fraction of cells incorrectly mapped, respectively.
#'
defineClassMarkers <- function(datRef, onClasses, offClasses, numMarkers = 50) {
  # Data prep and errors
  if (is.null(colnames(datRef))) 
    colnames(datRef) = as.character(1:length(colnames(datRef)))
  samples = colnames(datRef)
  
  if (length(samples) != length(onClasses)) 
    return("Error: onClasses is the wrong length.")
  if (length(samples) != length(offClasses)) 
    return("Error: onClasses is the wrong length.")
  
  offClasses <- factor(offClasses)
  onClasses <- factor(onClasses)
  names(onClasses) <- names(offClasses) <- samples
  
  # Caclulate proportionsa and medians
  propExpr <- do.call("cbind", tapply(names(onClasses), onClasses, 
    function(x) rowMeans(datRef[, x] > 1)))
  propExpr <- propExpr[, levels(onClasses)]
  medianExpr <- do.call("cbind", tapply(names(onClasses), onClasses, 
    function(x) rowMeans(datRef[, x])))
  medianExpr <- log2(medianExpr[, levels(onClasses)] + 1)
  rownames(propExpr) <- rownames(medianExpr) <- rownames(datRef)
  
  propExprC <- do.call("cbind", tapply(names(offClasses), offClasses, 
    function(x) rowMeans(datRef[, x] > 1)))
  propExprC <- propExprC[, levels(offClasses)]
  medianExprC <- do.call("cbind", tapply(names(offClasses), offClasses, 
    function(x) rowMeans(datRef[, x])))
  medianExprC <- log2(medianExprC[, levels(offClasses)] + 1)
  rownames(propExprC) <- rownames(medianExprC) <- rownames(datRef)
  
  # Define and return markers
  markers = list()
  
  for (cn in colnames(propExpr)) {
    a = (propExpr[, cn] - apply(propExpr[, colnames(propExpr) != 
      cn], 1, mean))
    b = ((medianExpr[, cn] - rowMeans(medianExpr[, colnames(medianExpr) != 
      cn]))/(medianExpr[, cn] + 1))
    kp = a * b * (a > 0) * (b > 0) * propExpr[, cn] * medianExpr[, 
      cn] * (medianExpr[, cn] >= 5) * (propExpr[, cn] >= 0.5)
    markers[[paste0(cn, "_on")]] <- make.names(names(head(-sort(-kp), 
      numMarkers)))
  }
  
  for (cn in colnames(propExprC)) {
    a = (propExprC[, cn] - apply(propExprC[, colnames(propExprC) != 
      cn], 1, max))
    b = ((medianExprC[, cn] - apply(medianExprC[, colnames(medianExprC) != 
      cn], 1, max))/(medianExprC[, cn] + 1))
    kp = a * b * (a > 0) * (b > 0) * sqrt(medianExprC[, cn])
    markers[[cn]] <- make.names(names(head(-sort(-kp), numMarkers)))
  }
  
  return(markers)
  
}




#' Calculate PatchSeq QC Metrics
#'
#' This function identifies is the same as calculatePatchSeqQCMetrics from patchSeqQC, except that it
#' allows for any user-inputted comparison data set, and fixes some other errors.  Importantly, it
#' outputs the same quality score, marker sum, and contamination score.
#'
#' @param pat_df a matrix (rows=samples x columns=genes + meta-data) of gene expression data and
#'   meta-data for patch-seq data (e.g., the data data set for QCing)
#' @param facs_df an equivalent matrix of reference data
#' @param markers a list of marker genes (calculated using defineClassMarkers)
#'
#' @return a table containing all of the qc metrics:
#' sample_id: name of the samples.
#' major_type: cell type identities (provided by Cadwell2016)
#' contam_type: cell type identities (normalized to cell type names in markers)
#' marker_sum: Summed expression of 'On' cell type marker genes (with cell type defined by contam_type)
#' marker_sum_norm: Normalized summed expression of 'on'-type marker genes, normalized to median expression of same cell type in dissociated-cell reference data
#' contam_sum: Contamination score, defined as the sum of normalized expression across all 'off' cell types defined in compare_cell_types_inh
#' quality_score: Quality score, defined as the Spearman correlation of marker expression between markers expressed in single cell with mean expression of markers in dissociated cell reference data
#' This function also outputs normalized expression of each 'off'-cell type (defined in compare_cell_types_inh) and we can use the function plotContamHeatmap to show these (each column is one single cell)
#'
calculatePatchSeqQCMetrics2 <- function(pat_df, facs_df, markers) {
  
  # This calculates some comparison matrix between each pair of
  # types
  facs_df$contam_type = facs_df$major_type
  aibs_contam_all_broad = calcContamAllTypes(facs_df, markers)
  aibs_contam_all_broad$contam_type = factor(facs_df$contam_type)
  aibs_med_exprs_broad = aibs_contam_all_broad %>% dplyr::group_by(contam_type) %>% 
    dplyr::summarize_all(median) %>% as.data.frame()
  rownames(aibs_med_exprs_broad) = aibs_med_exprs_broad$contam_type
  
  facs_df$contam_type = paste0(facs_df$contam_type, "_on")
  aibs_contam_all_sub = calcContamAllTypes(facs_df, markers)
  aibs_contam_all_sub$contam_type = factor(facs_df$contam_type)
  aibs_med_exprs_sub = aibs_contam_all_sub %>% dplyr::group_by(contam_type) %>% 
    dplyr::summarize_all(median) %>% as.data.frame()
  rownames(aibs_med_exprs_sub) = aibs_med_exprs_sub$contam_type
  
  aibs_med_exprs = rbind(aibs_med_exprs_broad, aibs_med_exprs_sub)
  
  
  dataset_cell_types = unique(pat_df$major_type) %>% as.character()
  marker_sums = lapply(dataset_cell_types, function(cell_type) {
    curr_marker_type = pat_df[pat_df$major_type == cell_type, 
      "major_type"][1] %>% as.character
    # curr_marker_type = paste0(curr_marker_type, '_on')
    curr_marker_list = markers[[curr_marker_type]]
    cell_inds = which(pat_df$major_type == curr_marker_type)
    df = pat_df[cell_inds, ]
    rownames(df) = df$sample_id
    marker_expr = sumExpression(df, curr_marker_list)
    compare_cell_types = setdiff(c("Astro", "Macrophage", "Endo", 
      "Oligo", "Peri", "SMC", "VLMC", "Glutamatergic", "GABAergic"), 
      cell_type)
    
    mks <- markers[c(curr_marker_type, compare_cell_types)] %>% 
      unlist %>% unique
    compare_expr_profile = facs_df[facs_df$major_type == cell_type, 
      mks] %>% log2 %>% colMeans
    contam_values = calcContamAllTypes(df, markers)
    contam_values$contam_type = cell_type
    
    
    expected_expr = aibs_med_exprs[contam_values$contam_type[1], 
      compare_cell_types] %>% unlist()  # Expression of all cell types in cluster of interest
    compare_cell_type_exprs = aibs_med_exprs[compare_cell_types, 
      compare_cell_types] %>% as.matrix %>% diag
    normalizedContamValues = apply(contam_values[, compare_cell_types], 
      1, function(x) (x - expected_expr)/(compare_cell_type_exprs - 
        expected_expr))
    
    contam_sum = normalizedContamValues %>% repWZero() %>% colSums
    
    marker_sum_norm = normalizeContam(contam_values, aibs_med_exprs, 
      c(curr_marker_type, compare_cell_types))
    marker_sum_norm_vec = contam_values[, curr_marker_type]/aibs_med_exprs[curr_marker_type, 
      curr_marker_type]
    
    
    quality_score = rbind(compare_expr_profile %>% t() %>% as.data.frame, 
      df[, mks] %>% log2) %>% t() %>% cor(method = "spearman")
    quality_score = quality_score[1, -1]
    
    out_df = data.frame(sample_id = rownames(df), marker_sum = marker_expr, 
      marker_sum_norm = marker_sum_norm_vec, contam_sum = contam_sum, 
      quality_score = quality_score)
    out_df = cbind(out_df, marker_sum_norm %>% t())
    return(out_df)
  })
  marker_sums = dplyr::bind_rows(marker_sums) %>% dplyr::select(sample_id, 
    dplyr::everything())
  marker_sums = merge(pat_df %>% dplyr::select(dplyr::one_of("sample_id", 
    "major_type", "contam_type")), marker_sums, by = "sample_id")
  marker_sum_df = marker_sums
  marker_sum_df = marker_sum_df[order(match(marker_sum_df$sample_id, 
    pat_df$sample_id)), ]
  
  return(marker_sum_df)
}
