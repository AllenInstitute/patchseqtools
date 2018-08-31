#' Reorder a matrix
#'
#' This function reorders a matrix by rows of columns
#'
#' @param matrix1 a matrix (rows=genes x columns=samples) of gene expression data (e.g., scRNA-seq)
#' @param by.rows By rows (TRUE; default) or by columns
#'
#' @return a reordered matrix
#'
reorder_matrix <- function(matrix1, by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
  return(matrix1.reordered)
}


#' Compare two cluster sets matched to CCA
#'
#' This function takes cluster calls defined in two different data sets and then determines to what
#' extent these cluster calls match up with cluster calls from CCA.
#'
#' @param cl a matrix (rows=genes x columns=samples) of gene expression data (e.g., scRNA-seq)
#' @param by.rows By rows (TRUE; default) or by columns
#'
#' @return a reordered matrix
#'
compareClusterCalls <- function(FACs_Patch_cca,
                              plot.title = NA, plot.silent = TRUE,
                              heat.colors = colorRampPalette(c("grey99", "orange", "red"))(100),
                              row.cl.num = min(length(unique(cl)),
                                               length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)

  # Set up the variables
  combined.cl <- FACs_Patch_cca@meta.data %>%
    mutate(ref.cl = paste(protocol, ref.cl))
  colnames(combined.cl) <- gsub("res.0.8","seurat.cl",colnames(combined.cl))

  cl <- setNames(combined.cl$ref.cl, rownames(FACs_Patch_cca@meta.data))
  ref.cl <- setNames(combined.cl$seurat.cl, rownames(FACs_Patch_cca@meta.data))

  combined.cl.anno <- combined.cl %>%
    group_by(ref.cl) %>%
    summarise(protocol = dplyr::first(protocol),
              mean.nGene = mean(nGene),
              #            mean.nUMI = mean(nUMI),
              n.cells = n()) %>%
    as.data.frame()

  rownames(combined.cl.anno)  <- combined.cl.anno$ref.cl
  combined.cl.anno <- combined.cl.anno[,-1]

  # Set up more variables
  conf1 <- table(cl, ref.cl)
  conf1 <- sweep(conf1, 1, rowSums(conf1), "/")
  conf2 <- reorder_matrix(conf1)

  # Cluster co-occurence
  cl.prop.cocl <- apply(conf1, 2, function(x) {
    grid1 <- expand.grid(x, x)
    min.prop <- apply(grid1, 1, min)
  })
  cl.prop.cocl.total <- apply(cl.prop.cocl, 1, sum)
  cl.prop.cocl.m <- matrix(cl.prop.cocl.total, nrow(conf1), nrow(conf1),
                           dimnames = list(rownames(conf1), rownames(conf1)))

  # Heatmap
  ph1 <- pheatmap(conf2, cutree_rows = row.cl.num, clustering_method = "ward.D2",
                  annotation_row = combined.cl.anno,
                  color = heat.colors, fontsize = 6,
                  main = plot.title, silent = plot.silent)
  return(list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1))
}
