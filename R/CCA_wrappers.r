#' Reorder a matrix
#'
#' This function reorders a matrix by rows of columns
#'
#' @param matrix1 a matrix (rows=genes x columns=samples) of gene expression data
#'   (e.g., scRNA-seq)
#' @param by.rows By rows (TRUE; default) or by columns
#'
#' @return a reordered matrix
#'
#' @export
reorder_matrix <- function(matrix1,
                           by.rows = TRUE) {
  if (by.rows == TRUE) {
    conf.order <- order(apply(matrix1, 1, which.max))
    matrix1.reordered <- matrix1[conf.order, ]
  } else {
    conf.order <- order(apply(matrix1, 2, which.max))
    matrix1.reordered <- matrix1[, conf.order]
  }
  matrix1.reordered
}


#' Compare two cluster sets matched to CCA
#'
#' This function takes cluster calls defined in two different data sets and then
#' determines to what extent these cluster calls match up with cluster calls from CCA.
#'
#' @param cl a matrix (rows=genes x columns=samples) of gene expression data
#'   (e.g., scRNA-seq)
#' @param by.rows By rows (TRUE; default) or by columns
#'
#' @return a reordered matrix
#'
#' @export
compareClusterCalls <- function(cl,
                                ref.cl,
                                cl.anno,
                                plot.title = NA, plot.silent = TRUE,
                                heat.colors = colorRampPalette(c("grey99", "orange", "red"))(100),
                                row.cl.num = min(length(unique(cl)), length(unique(ref.cl)))) {
  library(grid)
  library(pheatmap)
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
    dimnames = list(rownames(conf1), rownames(conf1))
  )

  ph1 <- pheatmap(conf2,
    cutree_rows = row.cl.num, clustering_method = "ward.D2",
    annotation_row = cl.anno, color = heat.colors, fontsize = 6,
    main = plot.title, silent = plot.silent
  )
  list(conf = conf2, cocl = cl.prop.cocl.m, ph = ph1)
}


#' Get some summary statistics
#'
#' This does the summary. For each group return a vector with N, mean, and sd.
#' This is called by qcPlot.
#'
#' @param data Annotation data frame
#' @param measurevar what variables to measure
#' @param groupvars what to group by
#' @param na.rm how to treat NA
#' @param conf.interval confidence interval (default = .95)
#' @param .drop drop something?
#' @param roundall should everything be rounded
#'
#' @return a data frame of statistics
#'
#' @export
summarySE <- function(data = NULL,
                      measurevar,
                      groupvars = NULL,
                      na.rm = FALSE,
                      conf.interval = .95,
                      .drop = TRUE,
                      roundall = F) {
  require(dplyr)
  #
  names(data)[names(data) == measurevar] <- "measurevar"

  datac <- data %>%
    select(one_of(groupvars, "measurevar")) %>%
    filter(ifelse(na.rm == T, !is.na(measurevar), T)) %>%
    mutate(measurevar = as.numeric(measurevar)) %>%
    group_by_(c(groupvars)) %>%
    summarise(
      N = n(),
      median = median(measurevar),
      mean = mean(measurevar),
      max = max(measurevar),
      sd = ifelse(N == 1, 0, sd(measurevar)),
      q25 = as.numeric(quantile(measurevar, 0.25)),
      q75 = as.numeric(quantile(measurevar, 0.75))
    ) %>%
    mutate(se = sd / sqrt(N))
  # %>% mutate(ci =  se * qt(conf.interval/2 + 0.5, N-1))

  if (roundall) {
    roundcols <- c("median", "mean", "max", "sd", "q25", "q75", "se", "ci")
    datac[roundcols] <- round(datac[roundcols], 3)
  }
  # datac <- datac %>% mutate(xpos = 1:n())

  datac
}


#' Make QC plots
#'
#' Makes QC plots
#'
#' @param anno Annotation data frame for patch-seq
#' @param dendcluster_anno what to cluster by
#' @param groupvars what to group by
#' @param scaleLimits scaleLimits
#' @param scaleBreaks scaleBreaks
#' @param scaleLabels scaleLabels
#' @param ylab ylab
#' @param fileName fileName
#'
#' @return the plot is returned
#'
#' @export
qcPlot <- function(anno,
                   dendcluster_anno,
                   name,
                   scaleLimits = c(-5000, 12000),
                   scaleBreaks = seq(0, 12000, 2000),
                   scaleLabels = seq(0, 12, 2),
                   ylab = "value",
                   fileName = gsub("\\.", "_", gsub("_label", "", name)),
                   outputFolder) {

  # dendcluster_id is the annotation for cluster ordering based on the current, bootstrapped dendrogram
  stats <- summarySE(data = anno, measurevar = name, groupvars = "dendcluster_id")

  genes_plot <- ggplot() +
    # geom_quasirandom from the ggbeeswarm package
    # makes violin-shaped jittered point plots
    geom_quasirandom(
      data = anno,
      aes(
        x = dendcluster_id,
        y = eval(parse(text = name))
      ),
      color = "skyblue",
      # Need to set position_jitter height = 0 to prevent
      # jitter on the y-axis, which changes data representation
      position = position_jitter(width = .3, height = 0), size = 0.1
    ) +
    # Errorbars built using stats values
    geom_errorbar(
      data = stats,
      aes(x = dendcluster_id, ymin = q25, ymax = q75),
      size = 0.2
    ) +
    # Median points from stats
    geom_point(
      data = stats,
      aes(x = dendcluster_id, y = median),
      color = "red",
      size = 0.5
    ) +
    # Cluster labels as text objects
    geom_text(
      data = dendcluster_anno,
      aes(
        x = dendcluster_id, y = 0, label = dendcluster_label,
        color = dendcluster_color
      ),
      angle = 90,
      hjust = 2,
      vjust = 0.3,
      size = 2 * 5 / 6
    ) +
    scale_color_identity() +
    # Expand the y scale so that the labels are visible
    scale_y_continuous(ylab,
      limits = scaleLimits,
      breaks = scaleBreaks,
      labels = scaleLabels
    ) +
    # Remove X-axis title
    scale_x_continuous("") +
    theme_bw() +
    # Theme tuning
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

  ggsave(paste0(outputFolder, fileName, "_QC.pdf"), genes_plot, width = 8, height = 4, useDingbats = F)
  genes_plot
}
