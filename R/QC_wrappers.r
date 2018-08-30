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
defineClassMarkers <- function(datRef, onClasses, 
  offClasses, numMarkers = 50) {
  # Data prep and errors
  if (is.null(colnames(datRef))) 
    colnames(datRef) = as.character(1:length(colnames(datRef)))
  samples = colnames(datRef)
  
  if (length(samples != length(onClasses))) 
    return("Error: onClasses is the wrong length.")
  if (length(samples != length(offClasses))) 
    return("Error: onClasses is the wrong length.")
  
  offClasses <- factor(offClasses)
  onClasses <- factor(anumMarkersoFACs2$subclass_label)
  names(onClasses) <- names(offClasses) <- samples
  
  # Caclulate proportionsa and medians
  propExpr <- do.call("cbind", tapply(names(onClasses), 
    onClasses, function(x) rowMeans(datRef[, x] > 
      1)))
  propExpr <- propExpr[, levels(onClasses)]
  medianExpr <- do.call("cbind", tapply(names(onClasses), 
    onClasses, function(x) rowMeans(datRef[, x])))
  medianExpr <- log2(medianExpr[, levels(onClasses)] + 
    1)
  rownames(propExpr) <- rownames(medianExpr) <- rownames(datRef)
  
  propExprC <- do.call("cbind", tapply(names(offClasses), 
    offClasses, function(x) rowMeans(datRef[, 
      x] > 1)))
  propExprC <- propExprC[, levels(offClasses)]
  medianExprC <- do.call("cbind", tapply(names(offClasses), 
    offClasses, function(x) rowMeans(datRef[, 
      x])))
  medianExprC <- log2(medianExprC[, levels(offClasses)] + 
    1)
  rownames(propExprC) <- rownames(medianExprC) <- rownames(datRef)
  
  # Define and return markers
  markers = list()
  
  for (cn in colnames(propExpr)) {
    a = (propExpr[, cn] - apply(propExpr[, colnames(propExpr) != 
      cn], 1, mean))
    b = ((medianExpr[, cn] - rowMeans(medianExpr[, 
      colnames(medianExpr) != cn]))/(medianExpr[, 
      cn] + 1))
    kp = a * b * (a > 0) * (b > 0) * propExpr[, 
      cn] * medianExpr[, cn] * (medianExpr[, 
      cn] >= 5) * (propExpr[, cn] >= 0.5)
    markers[[paste0(cn, "_on")]] <- make.names(names(head(-sort(-kp), 
      numMarkers)))
  }
  
  for (cn in colnames(propExprC)) {
    a = (propExprC[, cn] - apply(propExprC[, colnames(propExprC) != 
      cn], 1, max))
    b = ((medianExprC[, cn] - apply(medianExprC[, 
      colnames(medianExprC) != cn], 1, max))/(medianExprC[, 
      cn] + 1))
    kp = a * b * (a > 0) * (b > 0) * sqrt(medianExprC[, 
      cn])
    markers[[cn]] <- make.names(names(head(-sort(-kp), 
      numMarkers)))
  }
  
  return(markers)
  
}
