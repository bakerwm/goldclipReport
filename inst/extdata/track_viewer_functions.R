
#' make overlay plot
#'
#' make overlay plot for two or more
#'
tkplot2 <- function(data1, data2, nameA = "A", nameB = "B") {
  # for 2 samples, overlay
  stopifnot(is.data.frame(data1))
  stopifnot(is.data.frame(data2))
  stopifnot(is.character(nameA) & length(nameA) == 1)
  stopifnot(is.character(nameB) & length(nameB) == 1)
  # tags, TE
  tagA <- unique(data1$name)
  tagB <- unique(data2$name)
  stopifnot(tagA == tagB & length(tagA) == 1)
  col_sel <- c("name","sample", "position", "score", "y")
  stopifnot(all(col_sel %in% names(data1)))
  stopifnot(all(col_sel %in% names(data2)))
  # add sample column
  data1$sample <- nameA
  data2$sample <- nameB
  # merge data.frame
  data <- dplyr::bind_rows(data1, data2)
  data$sample <- factor(data$sample, levels = c(nameA, nameB))
  # plot
  mytheme <- pick_theme(2)
  p2 <- ggplot(data, aes(position, y, height = score,
                         color = sample, fill = sample)) +
    geom_ridgeline(show.legend = TRUE) +
    guides(color = guide_legend(title = NULL),
           fill = guide_legend(title     = NULL,
                               keywidth  = .4,
                               keyheight = .2,
                               label.position = "right")) +
    ylab("base coverage [rpm]") +
    ggtitle(tagA) +
    scale_fill_manual(values = c("orange2", "red3")) +
    scale_color_manual(values = c("orange2", "red3")) +
    mytheme
  return(p2)
}



#' create plots
#'
#'
te_plot_overlay <- function(grA, grB, id, nameA = "A", nameB = "B") {
  stopifnot(is(grA, "GRanges"))
  stopifnot(is(grB, "GRanges"))
  stopifnot(is.character(id) & length(id) == 1)
  grAs <- as.character(unique(seqnames(grA)))
  grBs <- as.character(unique(seqnames(grB)))
  stopifnot(id %in% grAs & id %in% grBs)

  # choose te
  x <- grA[seqnames(grA) == id]
  y <- grB[seqnames(grB) == id]

  # convert to data.frame
  dfA <- gr2df(x)
  dfB <- gr2df(y)

  # fix TE names
  dfA <- tidyr::separate(dfA, "chr", c("id", "name"), "_")
  dfB <- tidyr::separate(dfB, "chr", c("id", "name"), "_")
  #
  #   # single-plot
  #   tkplot1(df_het)
  #   tkplot1(df_mut)
  # return(list(dfA, dfB))
  p <- tkplot2(dfA, dfB, nameA, nameB)
  return(p)
}


#' choose theme
#'
#'
pick_theme <- function(x) {
  # single track
  t1 <- theme_linedraw() +
    theme(panel.border = element_rect(fill = NULL, color = "black", size = .7),
          panel.grid   = element_blank(),
          plot.title   = element_text(color = "black", hjust = .5, size = 12),
          axis.line    = element_blank(),
          # axis.line    = element_line(color = "black", size = .7),
          axis.ticks.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_line(color = "black", size = .7),
          axis.title.y = element_text(color = "black", size = 10, vjust = .5))
  # overlay tracks
  t2 <- t1 + theme(legend.position = c(0.9, 0.8))
  # choose theme
  tm <- list(t1, t2)
  names(tm) <- seq_len(length(tm))
  if (x %in% names(tm)) {
    return(tm[[x]])
  } else {
    return(tm[[1]])
  }
}


#' make coverage plot for one sample
#'
#'
tkplot1 <- function(data) {
  # for 1 sample
  # no legend
  stopifnot(is.data.frame(data))
  col_sel <- c("name", "position", "score", "y")
  stopifnot(all(col_sel %in% names(data1)))
  # tags/TE
  tag  <- unique(data$name)
  stopifnot(length(tag) == 1)
  # plot
  mytheme <- pick_theme(1)
  p1 <- ggplot(data, aes(position, y, height = score,
                         color = name, fill = name)) +
    geom_ridgeline(show.legend = FALSE) +
    guides(fill = guide_legend(title = NULL, label.position = "right")) +
    ylab("base coverage [rpm]") +
    ggtitle(tag) +
    scale_fill_manual(values = "orange") +
    scale_color_manual(values = "orange") +
    mytheme
  return(p1)
}



#' convert GRanges to data.frame
#'
#' expand ranges to 1-base contnet
#'
gr2df <- function(x, chr = NULL, start = NULL, end = NULL, gr = NULL) {
  # checkpoint, much too large dataset
  stopifnot(is(x, "GRanges"))
  df <- as.data.frame(x)
  # column score
  col_req <- c("seqnames", "start", "end", "width", "strand", "score")
  stopifnot(all(col_req %in% names(df)))
  # expand data.frame
  dd <- lapply(seq_len(nrow(df)), function(i){
    dx <- df[i, ]
    dy <- data.frame(chr      = dx$seqnames,
                     position = dx$start:dx$end,
                     score    = dx$score,
                     strand   = dx$strand,
                     stringsAsFactors = FALSE)
    return(dy)
  })
  # merge data.frame
  df_bp <- dplyr::bind_rows(dd)
  df_bp$y <- 1 # for ridges plots
  return(df_bp)
}


merge_list_of_dataframes <- function(list, by) {
  stopifnot(is.list(list))
  # print(by)
  # check
  if (length(list) == 1) {
    # print("AAA")
    return(list[[1]])
  } else if (length(list) == 2) {
    # print("BBB")
    df <- merge(list[[1]], list[[2]], by = by, all = TRUE)
    df[is.na(df)] <- 0
    return(df)
  } else {
    # print("CCC")
    df <- merge(list[[1]], list[[2]], by = by, all = TRUE)
    df[is.na(df)] <- 0
    # drop 1, 2
    list <- list.remove(list, c(1, 2))
    list_new <- list.append(list, df)
    # list[[1]] <- NULL
    # list[[2]] <- NULL
    # list
    # list_new <- c(list, df)
    # print(length(list_new))
    return(merge_list_of_dataframes(list_new, by = by))
  }
}
