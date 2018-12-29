


#' make plots for publication
#'
#' @param x Path to the csv file of DEseq2 output, eg, function
#'   DESeq2_for_featureCounts, result "transcripts_deseq2.csv"
#'
#' @param path_pdf Directory to save the pdf file, contains
#'   scatter, MA, Volcano plots
#'
#' @param save2df A logical value, whether or not save the plot to PDF file,
#'   default: TRUE
#'
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
DESeq2_publish_plot <- function(x, path_pdf,
                                gene_labels = c("nxf2", "piwi", "CG9754"),
                                save2pdf = TRUE) {
  stopifnot(file.exists(x))

  # read file
  df <- DESeq2_csv2df(x)
  names(df) <- gsub("[^A-Za-z0-9]", "_", names(df)) # remove unsupport characters
  x.name <- colnames(df)[2]
  y.name <- colnames(df)[3]

  # make plots
  p1 <- DESeq2_de_scatter(df, x.name, y.name,
                          show.sig.genes = FALSE,
                          labels = gene_labels,
                          pval_cutoff = 0.05, max_labels = 3)

  p2 <- DESeq2_de_ma(df, show.sig.genes = FALSE,
                     labels = gene_labels,
                     pval_cutoff = 0.05, max_labels = 3,
                     show_volcano = FALSE)

  p3 <- DESeq2_de_ma(df, show.sig.genes = TRUE,
                     labels = gene_labels,
                     pval_cutoff = 0.05, max_labels = 5,
                     show_volcano = TRUE)
  # add legend
  figure_legend <- "Figure. Differentially expression analysis. (a-b) Comparasion of
  gene expression is shown as rpm from two conditions. Dashed lines indicate
  twofold change. a. scatter plot, b. MA plot. (c) Volcano plot showing
  enrichment values and corresponding significance levels."

  pg1 <- cowplot::plot_grid(p1, p2, p3, NULL, align = "hv",
                            ncol = 2, labels = "auto")
  pg2 <- cowplot::add_sub(pg1, figure_legend, x = 0, hjust = 0, size = 10)
  # cowplot::ggdraw(pg2)

  # save plot
  if (isTRUE(save2pdf)) {
    if (! dir.exists(path_pdf)) {
      dir.create(path_pdf, showWarnings = FALSE, recursive = TRUE, mode = "0740")
    }
    rpt_fname <- paste0("DESeq2.", x.name, "_vs_", y.name, ".publish.pdf")
    rpt_file  <- file.path(path_pdf, rpt_fname)
    pdf(rpt_file, width = 6, height = 6, paper = "a4")
    print(cowplot::ggdraw(pg2))
    dev.off()
  } else {
    return(pg2)
  }
}




#' generate scatter plot for DESeq2 analysis
#'
#' @param data A data frame with gene expression,
#'
#' @param x.name Name in column names of data frame, present on x axis
#'
#' @param y.name Same as x.name, but present on y axis
#'
#' @param show.sig.genes A logical or int value, whether or not to show
#'   significantly changes genes, default: FALSE
#'
#' @param labels A set of gene ids to show in plot, default: NULL
#'
#' @param pval_cutoff A value to select significantly changes genes, default: 0.05
#'
#' @param max_labels A integer value, how many genes to be labeled on plot,
#'   default: 3
#'
#' @import ggplot2
#' @import ggrepeal
#'
#' @export
#'
DESeq2_de_scatter <- function(data, x.name, y.name, show.sig.genes = FALSE,
                              labels = NULL, pval_cutoff = 0.05, max_labels = 3) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("id", x.name, y.name, "padj") %in% names(data)))

  # rename column
  data <- dplyr::rename(data, pval = padj)

  # re-arrange data.frame
  data_sig <- dplyr::filter(data, pval <= pval_cutoff) %>% dplyr::arrange(pval)

  # pick labels
  data_label <- dplyr::filter(data, id %in% labels)
  if (isTRUE(show.sig.genes)) {
    data_label <- rbind(data_label, data_sig)
  }
  if (nrow(data_label) > 0) {
    max_labels <- ifelse(nrow(data_label) > max_labels, max_labels, nrow(data_label))
    data_label <- data_label[1:max_labels, ]
  }

  p <- plot_scatter(data = data,
                    data_sig = data_sig,
                    x.name = x.name,
                    y.name = y.name,
                    data_label = data_label)
  return(p)
}



#' create MA or Volcano plot
#'
#' @param data A data frame with gene expression,
#'
#' @param show.sig.genes A logical or int value, whether or not to show
#'   significantly changes genes, default: FALSE
#'
#' @param labels A set of gene ids to show in plot, default: NULL
#'
#' @param pval_cutoff A value to select significantly changes genes, default: 0.05
#'
#' @param max_labels A integer value, how many genes to be labeled on plot,
#'   default: 3
#'
#' @import ggplot2
#' @import ggrepeal
#'
#' @export
#'
DESeq2_de_ma <- function(data, show.sig.genes = FALSE, labels = NULL,
                         pval_cutoff = 0.05, max_labels = 3,
                         show_volcano = FALSE) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("id", "baseMean", "log2FoldChange", "padj") %in% names(data)))

  # rename columns
  data <- dplyr::rename(data, logFC = log2FoldChange, pval = padj)

  # re-arrange data.frame
  data_sig <- dplyr::filter(data, pval <= pval_cutoff) %>% dplyr::arrange(pval)

  data_non_sig <- dplyr::filter(data, pval <= pval_cutoff)

  # pick labels
  data_label <- dplyr::filter(data, id %in% labels)
  if (isTRUE(show.sig.genes)) {
    data_label <- rbind(data_label, data_sig)
  }
  if (nrow(data_label) > 0) {
    max_labels <- ifelse(nrow(data_label) > max_labels, max_labels, nrow(data_label))
    data_label <- data_label[1:max_labels, ]
  }

  # volcano plot
  if (isTRUE(show_volcano)) {
    p <- plot_volcano(data, data_sig, data_label)
  } else {
    # change log2foldchange
    data$logFC <- ifelse(data$logFC > 4, 4,
                         ifelse(data$logFC < -4, -4, data$logFC))

    p <- plot_ma(data, data_sig, data_label)
  }
  return(p)
}



#' generate scatter plot for DE analysis
#'
#' @param data A data frame for the full version of data, required
#'
#' @param data_sig A data frame for significantly changes genes, required
#'
#' @param data_label A data frame for points to be labeled, default: NULL
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_scatter <- function(data, data_sig, x.name, y.name, x.label = NULL,
                         y.label = NULL, data_label = NULL) {
  # data_sig is required
  stopifnot(is.data.frame(data))
  stopifnot(is.data.frame(data_sig))
  stopifnot(is.character(x.name))
  stopifnot(is.character(y.name))
  stopifnot(all(c("id", x.name, y.name) %in% names(data)))
  stopifnot(all(c("id", x.name, y.name) %in% names(data_sig)))

  # labels on axis
  x.label <- ifelse(is.null(x.label), x.name, x.label)
  y.label <- ifelse(is.null(y.label), y.name, y.label)
  xn <- bquote(.(x.label)~"["*log["10"]~"rpm]")
  yn <- bquote(.(y.label)~"["*log["10"]~"rpm]")

  # pre-process
  if (nrow(data) == 0) {
    return(NULL)
  }

  # function
  to_log10 <- function(data, x, y) {
    stopifnot(is.data.frame(data))
    stopifnot(all(c("id", x, y) %in% names(data)))
    # select columns
    df <- dplyr::select_(data, .dots = c("id", x, y)) %>%
      as.data.frame()
    # remove old row names
    rownames(df) <- NULL
    df <- tibble::column_to_rownames(df, "id") %>% log10()
    df <- df[complete.cases(df), ]
    return(df)
  }

  df <- to_log10(data, x.name, y.name)
  if (nrow(df) == 0) {
    stop("empty dataframe")
  }

  p <- ggplot(df, aes_string(x.name, y.name)) +
    geom_abline(slope = 1, intercept = 0, color = "grey10") +
    geom_abline(slope = 1, intercept = log10(2),
                color = "grey30", linetype = 2) +
    geom_abline(slope = 1, intercept = -log10(2),
                color = "grey30", linetype = 2) +
    geom_point(color = "grey60", size = .4) +
    xlab(xn) + ylab(yn) +
    scale_x_continuous(limits = c(0, 6),
                       breaks = seq(0, 6, 1),
                       labels = seq(0, 6, 1),
                       expand = c(0, 0, 0, 0)) +
    scale_y_continuous(limits = c(0, 6),
                       breaks = seq(0, 6, 1),
                       labels = seq(0, 6, 1),
                       expand = c(0, 0, 0, 0)) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
          plot.title   = element_text(color = "black", hjust = .5, size = 14),
          panel.grid   = element_blank(),
          axis.line    = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          axis.ticks   = element_line(color = "black", size = .5),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12))

  if (nrow(data_sig) > 0) {
    df_sig <- to_log10(data_sig, x.name, y.name)
    p <- p +
      geom_point(mapping = aes_string(x.name, y.name),
                 data    = df_sig,
                 color   = "red",
                 size    = .6)
  }

  if (is.data.frame(data_label) & nrow(data_label) > 0) {
    df_label <- to_log10(data_label, x.name, y.name) %>%
      tibble::rownames_to_column("id")
    stopifnot(all(c("id", x.name, y.name) %in% names(data)))
    p <- p +
      geom_text_repel(
        mapping       = aes_string(x.name, y.name),
        data          = df_label,
        label         = df_label$id,
        # size          = 10,
        box.padding   = 1,
        segment.size  = 0.4,
        segment.color = "grey50",
        direction     = "both") +
      geom_point(mapping = aes_string(x.name, y.name),
                 data    = data_label,
                 color   = "red",
                 size    = 1.2)
  }
  return(p)
}



#' generate volcano plot for DE analysis
#'
#' @param data A data frame for the full version of data, required
#'
#' @param data_sig A data frame for significantly changes genes, required
#'
#' @param data_label A data frame for points to be labeled, default: NULL
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_ma <- function(data, data_sig, data_label = NULL) {
  # data_sig is required
  stopifnot(is.data.frame(data))
  stopifnot(is.data.frame(data_sig))
  stopifnot(all(c("id", "baseMean", "logFC", "pval") %in% names(data)))
  stopifnot(all(c("id", "baseMean", "logFC", "pval") %in% names(data_sig)))

  xn <- expression(paste(log["10"], "(baseMean)"))
  yn <- expression(paste(log["2"], "(Fold-change)"))

  data_non_sig <- dplyr::filter(data, ! id %in% data_sig$id)

  # plot
  p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey10", size = .7) +
    geom_point(mapping = aes(baseMean, logFC),
               data    = data_non_sig,
               color   = "grey60",
               size    = .6) +
    geom_point(mapping = aes(baseMean, logFC),
               data    = data_sig,
               color   = "red",
               size    = .6) +
    scale_x_log10() +
    scale_y_continuous(limits = c(-4, 4)) +
    xlab(xn) + ylab(yn) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
          plot.title   = element_text(color = "black", hjust = .5, size = 14),
          panel.grid   = element_blank(),
          axis.line    = element_blank(),
          axis.ticks.length = unit(.2, "cm"),
          axis.ticks   = element_line(color = "black", size = .5),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12))

  if (is.data.frame(data_label) & nrow(data_label) > 0) {
    stopifnot(all(c("id", "baseMean", "logFC", "pval") %in% names(data_label)))
    p <- p +
      geom_text_repel(
        mapping       = aes(baseMean, logFC),
        data          = data_label,
        label         = data_label$id,
        # size          = 10,
        box.padding   = 1,
        segment.size  = 0.4,
        segment.color = "grey50",
        direction     = "both"
      ) +
      geom_point(mapping = aes(baseMean, logFC),
                 data    = data_label,
                 color   = "red",
                 size    = 1.2)
  }
  return(p)
}



#' generate volcano plot for DE analysis
#'
#' @param data A data frame for the full version of data
#'
#' @param data_sig A data frame for significantly changed genes
#'
#' @param data_label A data frame for points to be labeled
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
plot_volcano <- function(data, data_sig = NULL, data_label = NULL) {
  # data_sig not required
  stopifnot(is.data.frame(data))
  stopifnot(all(c("id", "logFC", "pval") %in% names(data)))

  xn <- expression(paste(log["2"], "(mean ratio of Treatment/Control)"))
  yn <- expression(paste(-log["10"], "(", italic("P"), " value)"))

  data$logPval <- -log10(data$pval)

  p <- ggplot() +
    geom_vline(xintercept = 0, size = .6, color = "black") +
    geom_point(mapping = aes(logFC, logPval),
               data    = data,
               color   = "grey70",
               size    = .6) +
    xlab(xn) + ylab(yn) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid   = element_blank(),
          axis.line    = element_line(color = "black", size = .7),
          axis.ticks.length = unit(.2, "cm"),
          axis.ticks   = element_line(color = "black", size = .7),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12))

  if (is.data.frame(data_sig)) {
    stopifnot(all(c("id", "logFC", "pval") %in% names(data)))
    data_sig$logPval <- -log10(data_sig$pval)
    p <- p + geom_point(mapping = aes(logFC, logPval),
                        data    = data_sig,
                        # color   = "chocolate2",
                        color   = "red",
                        size    = .6)
  }

  if (is.data.frame(data_label) & nrow(data_label) > 0) {
    stopifnot(all(c("id", "logFC", "pval") %in% names(data_label)))
    data_label$logPval <- -log10(data_label$pval)
    p <- p +
      geom_text_repel(
        mapping       = aes(logFC, logPval),
        data          = data_label,
        label         = data_label$id,
        # size          = 10,
        nudge_x       = 1.2,
        nudge_y       = .02,
        point.padding = .4,
        box.padding   = .4,
        segment.size  = .6,
        segment.color = "grey60",
        direction     = "both") +
      geom_point(mapping = aes(logFC, logPval),
                 data    = data_label,
                 # color   = "chocolate2",
                 color   = "red",
                 size    = 1.2)
  }
  return(p)
}
