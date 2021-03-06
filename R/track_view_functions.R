
# Functions for bigWig viewer

#' examples for bigWig viewer
#'
#' ctl_ip_bw <- "het_ip.bigWig"
#' ctl_input_bw <- "het_input.bigWig"
#' tre_ip_bw <- "mut_ip.bigWig"
#' tre_input_bw <- "mut_input.bigWig"
#' ctl_label <- "sampleA"
#' tre_label <- "sampleB"
#' fasize <- "te.sizes"
#' # n <- "HeT-A"
#' # n <- TRUE
#' # n <- "Penelope"
#' n <- c("HeT-A", "blood", "mdg1", "Burdock")
#'
#'
#' # make plots for single sample
#' df_single <- chipseq_bw_parser(ip_bw = ctl_ip_bw, fasize = fasize,
#'                                input_bw = ctl_input_bw, n = n)
#'
#' plist_single <- lapply(seq_len(length(df_single)), function(i) {
#'   p <- coverage_plot_single(df_single[[i]])
#'   return(p)
#' })
#'
#' plot_n_pages(plist_single, nrow = 5, ncol = 2, pdf_out = "foo-single.pdf")
#'
#' # make plots for dual samples
#' df_dual <- chipseq_bw_parser2(ctl_ip_bw, tre_ip_bw, fasize, ctl_input_bw,
#'                               tre_input_bw, ctl_label = "CG9754 het",
#'                               tre_label = "CG9754 mut", n = n)
#'
#' plist_dual <- lapply(seq_len(length(df_dual)), function(i) {
#'   p <- coverage_plot_dual(df_dual[[i]])
#'   return(p)
#' })
#'
#' plot_n_pages(plist_dual, nrow = 5, ncol = 2, pdf_out = "foo-dual.pdf")
#'


#' te_rename
#'
#' rename the TE ids, from FBgn0000004_17.6 to 17.6
#'
te_rename <- function(x){
  x_new <- x
  if(inherits(x, "character")) {
    if(all(grepl("_", x))) {
      ma <- stringr::str_split(x, "_", simplify = TRUE)
      x_new <- ma[, 2]
    }
  }
  return(x_new)
}


#' read bigWig file as GRanges object
#'
#' 1. normalize by input
#' 2. filter by GRrange object
#'
#' @param x path to the bigWig file, generally, the IP of ChIP-seq results
#'
#' @param y path to a bigWig file, to normalize x, generally, the Input
#'   sample of ChIP-seq, default NULL
#'
#' @param which GRanges object, subset of bigWig records
#'
#' @param normalize.to.y logical, whether normalize IP to Input bigWig
#'
#' @import rtracklayer
#' @import trackViewer
#' @import dplyr
#'
#' @export
bw_reader <- function(x, y = NULL, which = NULL, normalize.to.y = FALSE) {
  stopifnot(file.exists(x))
  stopifnot(endsWith(tolower(x), "bw") | endsWith(tolower(x), "bigwig"))

  if (is.null(which)) {
    # read entire bigWig file
    # warning, large bigwig file
    gr_x <- rtracklayer::import.bw(x)
    bw   <- list(x = gr_x)

    # normalize to y
    # operate "-"
    if(inherits(y, "character") & isTRUE(normalize.to.y)) {
      stopifnot(file.exists(y))
      stopifnot(endsWith(tolower(x), "bw") | endsWith(tolower(x), "bigwig"))

      gr_y <- rtracklayer::import.bw(y)
      bw   <- list(x = gr_x, y = gr_y)

      gr_norm <- trackViewer::GRoperator(gr_x, gr_y,
                                         col = "score",
                                         operator = "-")
      bw   <- list(x = gr_x, y = gr_y, norm = gr_norm)
    }
  } else if(is(which, "GRanges")) {
    # read subset of bigWig
    # ip bigWig
    gr_x <- rtracklayer::import.bw(x, which = which)
    bw   <- list(gr_x)

    # normalize to y
    # operate GRoperator() "-"
    if(inherits(y, "character") & isTRUE(normalize.to.y)) {
      stopifnot(file.exists(y))
      stopifnot(endsWith(tolower(x), "bw") | endsWith(tolower(x), "bigwig"))

      gr_y <- rtracklayer::import.bw(y, which = which)
      bw   <- list(x = gr_x, y = gr_y)

      if(length(gr_x) > 0) {
        if (length(gr_y) > 0) {
          gr_norm <- trackViewer::GRoperator(gr_x, gr_y,
                                             col = "score",
                                             operator = "-")
        } else {
          gr_norm <- gr_x
        }
        bw <- list(x = gr_x, y = gr_y, norm = gr_norm)
      } else {
        bw <- NULL
      }
    }
  }
  return(bw)
}


#' read fa_size
#'
#' read fa_size as data.frame
#' or convert to GRanges
#'
#' @param x path or file to the fasize, at least two columns
#'
#' @param to_GRanges logical, whether convert to GRanges object, defult: FALSE
#'
#' @import dplyr
#' @import readr
#' @import tidyr
#' @import GenomicRanges
#'
#' @export
#'
fasize_reader <- function(x, to_GRanges = FALSE) {
  stopifnot(file.exists(x))

  # read file
  df <- readr::read_delim(x, "\t", col_names = c("chr", "length"),
                          col_types = readr::cols()) %>%
    dplyr::mutate(start = 1, end = length, strand = "+") %>%
    dplyr::select(chr, start, end, length, strand) %>%
    dplyr::mutate(name = te_rename(chr)) # split name

  # # for dm3 TEs, FBgn0000004_17.6, split
  # if (all(grepl("_", df2$chr))) {
  #   df2 <- df2 %>%
  #     tidyr::separate("chr", c("id", "name"), sep = "_", remove = FALSE) %>%
  #     dplyr::select(chr, start, end, length, strand, name)
  # } else {
  #   df2$name <- df2$chr
  # }

  # convert to GRanges
  if (isTRUE(to_GRanges)) {
    out <- GenomicRanges::GRanges(seqnames = df$chr,
                                  ranges = IRanges::IRanges(start = df$start,
                                                            end   = df$end,
                                                            width = df$length,
                                                            names = df$name),
                                  strand = df$strand)
  } else {
    out <- df
  }
  return(out)
}


#' pick chromosome from GRanges
#'
#' @param x a GRanges object
#'
#' @param g a list of chromosome names to extract
#'
#' @param exclude logical, whether extract or exclude the specific chromosome
#'   names from bigWig file. default: FALSE
#'
#' @import GenomicRanges
#' @import dplyr
#'
#' @export
#'
GRanges_seqname_picker <- function(x, g, exclude = FALSE) {
  stopifnot(is(x, "GRanges"))
  stopifnot(all(is.character(g)))
  stopifnot(is.logical(exclude))

  # all seqnames in x
  n <- as.character(droplevels(GenomicRanges::seqnames(x)))

  # checkpoint, if x is empty
  if (length(n) == 0) {
    return(NULL)
  }

  # subset seqnames
  n0 <- n[! n %in% g]
  n1 <- n[n %in% g]

  # checkpoint, if g is not in x
  if (identical(n1, character(0))) {
    return(NULL)
  }

  # subset GRanges
  if (isTRUE(exclude)) {
    hit = n0
    hit_null = n1
  } else {
    hit = n1
    hit_null = n0
  }

  # GRanges object, drop levels
  # x1 <- x[seqnames(x) == hit]
  x1 <- x[as.character(x@seqnames) %in% hit] # !!! %in% now working
  x1 <- dropSeqlevels(x1, hit_null, pruning.mode = "coarse")

  return(x1)
}


#' convert GRanges to data.frame
#'
#' expand ranges to 1-base contnet
#'
#' @param x GRanges object, require extra column "score"
#'
#' @import GenomicRanges
#' @import dplyr
#'
#' @export
#'
GRanges_to_data.frame <- function(x) {
  # checkpoint, much too large dataset
  stopifnot(is(x, "GRanges"))
  df <- as.data.frame(x)

  # column score
  col_required <- c("seqnames", "start", "end", "width", "strand", "score")
  stopifnot(all(col_required %in% names(df)))

  # expand data.frame
  if (nrow(df) == 0) {
    df_bp <- data.frame(chr      = factor(),
                        position = integer(),
                        score    = double(),
                        strand   = factor(),
                        name     = factor(),
                        y        = integer(),
                        sample   = character(),
                        stringsAsFactors = FALSE)
  } else {
    dd <- lapply(seq_len(nrow(df)), function(i){
      dx <- df[i, ]
      dy <- data.frame(chr      = dx$seqnames,
                       position = dx$start:dx$end,
                       score    = dx$score,
                       strand   = dx$strand,
                       name     = dx$seqnames,
                       stringsAsFactors = FALSE)
      return(dy)
    })
    # merge data.frame
    df_bp <- dplyr::bind_rows(dd)
    df_bp$y <- 1 # for ridges plots
  }
  return(df_bp)
}


#' read bigWig files and return in data.frame format
#'
#' @param ip_bw the bigWig file of IP sampels of ChIP-seq experiment
#'
#' @param fasize the chromosome size file of the bigWig, tab separated
#'
#' @param input_bw the bigWig file of Input samples of ChIP-seq experiment
#'
#' @param n array or logical, a list of chromosomes from the bigWig files, or
#'   all chromosomes. default: TRUE
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import dplyr
#'
#' @example  n <- c("FBgn0000199_blood")
#' @example  n <- c("297", "FBgn0000199_blood")
#'
#' @export
#'
chipseq_bw_parser <- function(ip_bw, fasize, input_bw = NULL, n = TRUE) {
  stopifnot(file.exists(ip_bw))
  stopifnot(file.exists(fasize))

  # determine how many seqnames to read from bigWig files
  chr_all <- fasize_reader(fasize, to_GRanges = TRUE)
  if (isTRUE(n)) {
    hits <- chr_all
  } else if (all(is.character(n))) {
    # filt by chr names, for TE only
    dfx <- BiocGenerics::as.data.frame(chr_all)
    dfx_names <- rownames(dfx)
    dfx_chrs  <- as.character(dfx$seqnames)
    # overlap
    h1 <- dfx[dfx_names %in% n, ]
    h2 <- dfx[dfx_chrs %in% n, ]
    chr_hits <- rbind.data.frame(h1, h2)$seqnames
    chr_hits <- as.character(droplevels(chr_hits))
    hits <- GRanges_seqname_picker(chr_all, chr_hits)
  } else {
    stop("unknown type of argument, n=")
  }

  if (is.null(hits)) {
    stop("unable to extract hits")
  }

  # the seqnames of hits
  hits_seqnames <- as.character(GenomicRanges::seqnames(hits))

  # number of seqnames records, hits
  if (length(hits) == 0) {
    stop("no chromosome records selected, check, fasize and n options")
  }

  # read bigWig file
  if (is.null(input_bw)) {
    gr_list <- bw_reader(x = ip_bw, which = hits)
    gr_bw <- gr_list[[1]]  # the first one is ip_bw
  } else {
    gr_list <- bw_reader(x = ip_bw, y = input_bw, which = hits,
                         normalize.to.y = TRUE)
    gr_bw <- gr_list$norm # normalized bigWig
  }

  hits_n <- length(hits_seqnames)
  bb <- lapply(seq_len(hits_n), function(i){
    i_name <- hits_seqnames[i]
    print(glue::glue("{i} of {hits_n}, {i_name}"))

    # parse bigWig records
    empty_df <- data.frame(chr      = factor(),
                           position = integer(),
                           score    = double(),
                           strand   = factor(),
                           name     = factor(),
                           y        = integer(),
                           sample   = character(),
                           stringsAsFactors = FALSE)

    # ctl sample
    if (is(gr_bw, "GRanges")) {
      gr_hit <- GRanges_seqname_picker(gr_bw, i_name)
      if (is.null(gr_hit)) {
        gr_df <- empty_df
      } else {
        gr_df  <- GRanges_to_data.frame(gr_hit)
      }
    } else {
      gr_df  <- empty_df
    }

    # report
    return(gr_df)
  })

  # rename list of GRanges objects
  names(bb) <- hits_seqnames

  # a list of data.frames
  return(bb)
}


#' process pair-experiment files
#' read bigWig files and return in data.frame format
#'
#' @param ctl_ip_bw the bigWig file of IP samples in control experiment
#'
#' @param ctl_input_bw the bigWig file of Input samples in control experiment
#'
#' @param tre_ip_bw the bigWig file of IP samples in treatment experiment
#'
#' @param tre_input_bw the bigWig file of Input samples in treatment experiment
#'
#' @param fasize the chromosome size file of the bigWig, tab separated
#'
#' @param n a list of chromosomes names to extract from the bigWig files,
#'   default, TRUE, use all seqnames
#'
#' @return a list of data.frame
#'
#' @import dplyr
#'
#' @export
chipseq_bw_parser2 <- function(ctl_ip_bw, tre_ip_bw, fasize,
                               ctl_input_bw = NULL, tre_input_bw = NULL,
                               ctl_label = NULL, tre_label = NULL,
                               n = TRUE) {
  stopifnot(file.exists(ctl_ip_bw))
  stopifnot(file.exists(tre_ip_bw))
  stopifnot(file.exists(fasize))

  # determine the labels in plot
  ctl_prefix <- gsub(".bigwig|.bw", "", basename(ctl_ip_bw), ignore.case = TRUE)
  tre_prefix <- gsub(".bigwig|.bw", "", basename(tre_ip_bw), ignore.case = TRUE)
  if (is.null(ctl_label)) {ctl_label <- ctl_prefix}
  if (is.null(tre_label)) {tre_label <- tre_prefix}

  # determine how many seqnames to read from bigWig files
  chr_all <- fasize_reader(fasize, to_GRanges = TRUE)
  hits <- NULL
  if (isTRUE(n)) {
    hits <- chr_all
  } else if (all(is.character(n))) {
    hitA <- names(chr_all) %in% n
    hitB <- as.character(chr_all@seqnames) %in% n
    hits <- chr_all[hitA | hitB, ]
  } else {
    stop("unknown type of argument, n=")
  }
  # checkpoint
  if (is.null(hits)) {
    stop("records not found in fasize, check n= option")
  }

  # the seqnames of hits
  hits_seqnames <- as.character(seqnames(hits))

  # number of seqnames records, hits
  if (length(hits) == 0) {
    stop("no chromosome records selected, check, fasize and n options")
  }

  # read bigWig file
  if (is.null(ctl_input_bw) & is.null(tre_input_bw)) {
    ctl_gr_list <- bw_reader(x = ctl_ip_bw, which = hits)
    tre_gr_list <- bw_reader(x = tre_ip_bw, which = hits)
    ctl_gr_bw <- ctl_gr_list[[1]] # the first one is ip_bw
    tre_gr_bw <- tre_gr_list[[1]] # the first one is ip_bw
  } else {
    stopifnot(file.exists(ctl_input_bw) & file.exists(tre_input_bw))
    ctl_gr_list <- bw_reader(x = ctl_ip_bw, y = ctl_input_bw, which = hits,
                             normalize.to.y = TRUE)
    tre_gr_list <- bw_reader(x = tre_ip_bw, y = tre_input_bw, which = hits,
                             normalize.to.y = TRUE)
    ctl_gr_bw <- ctl_gr_list$norm # normalized bigWig, NULL or GRanges
    tre_gr_bw <- tre_gr_list$norm # normalized bigWig, NULL or GRanges
  }

  # iterate all the chr_names
  hits_n <- length(hits_seqnames)
  bb <- lapply(seq_len(hits_n), function(i){
    i_name <- hits_seqnames[i]
    print(glue::glue("{i} of {hits_n}, {i_name}"))

    # parse bigWig records
    empty_df <- data.frame(chr      = factor(),
                           position = integer(),
                           score    = double(),
                           strand   = factor(),
                           name     = factor(),
                           y        = integer(),
                           sample   = factor(),
                           stringsAsFactors = FALSE)
    # ctl sample
    if (is(ctl_gr_bw, "GRanges")) {
      ctl_gr_hit <- GRanges_seqname_picker(ctl_gr_bw, i_name)
      if (is.null(ctl_gr_hit)) {
        ctl_gr_df <- empty_df
      } else {
        ctl_gr_df  <- GRanges_to_data.frame(ctl_gr_hit)
      }
      # add label to smaple
      if (nrow(ctl_gr_df) > 0) {
        ctl_gr_df$sample <- ctl_label
      }
    } else {
      ctl_gr_df  <- empty_df
    }

    # tre sample
    if (is(tre_gr_bw, "GRanges")) {
      tre_gr_hit <- GRanges_seqname_picker(tre_gr_bw, i_name)
      if (is.null(tre_gr_hit)) {
        tre_gr_df <- empty_df
      } else {
        tre_gr_df  <- GRanges_to_data.frame(tre_gr_hit)
      }
      # add label to smaple
      if (nrow(tre_gr_df) > 0) {
        tre_gr_df$sample <- tre_label
      }
    } else {
      tre_gr_df <- empty_df
    }

    # merge data.frame
    gr_df <- rbind.data.frame(ctl_gr_df, tre_gr_df)

    # set levels
    gr_df$sample <- factor(gr_df$sample, levels = c(ctl_label, tre_label))

    # remove minus values
    gr_df <- dplyr::filter(gr_df, score > 0)

    # report data.frame
    return(gr_df)
  })

  # rename list of GRanges objects
  names(bb) <- hits_seqnames

  # a list of data.frames
  return(bb)
}



#' process pair-experiment files
#' read bigWig files and return in data.frame format
#'
#' @param ctl_ip_bw the bigWig file of IP samples in control experiment
#'
#' @param ctl_input_bw the bigWig file of Input samples in control experiment
#'
#' @param tre_ip_bw the bigWig file of IP samples in treatment experiment
#'
#' @param tre_input_bw the bigWig file of Input samples in treatment experiment
#'
#' @param fasize the chromosome size file of the bigWig, tab separated
#'
#' @param n a list of chromosomes names to extract from the bigWig files,
#'   default, TRUE, use all seqnames
#'
#' @param pdf_out logical, save the plots in pdf file
#'
#' @import dplyr
#'
#' @export
#'
chipseq_bw_parser3 <- function(ctl_ip_bw, tre_ip_bw, fasize,
                               ctl_input_bw = NULL, tre_input_bw = NULL,
                               ctl_label = NULL, tre_label = NULL,
                               n = TRUE, pdf_out = NULL) {
  if (is.null(pdf_out)) {
    pdf_dual <- as.character(glue::glue("ChIPseq.{ctl_label}.vs.{tre_label}.5TEs.trackview.pdf"))
  } else {
    pdf_dual <- pdf_out
  }

  # make plots for dual samples
  df_dual <- chipseq_bw_parser2(ctl_ip_bw, tre_ip_bw, fasize,
                                ctl_input_bw, tre_input_bw,
                                ctl_label = ctl_label,
                                tre_label = tre_label,
                                n = n)

  plist_dual <- lapply(seq_len(length(df_dual)), function(i) {
    df_dual[[i]] <- dplyr::filter(df_dual[[i]], score > 0)
    p <- coverage_plot_dual(df_dual[[i]])
    return(p)
  })
  plot_n_pages(plist_dual, nrow = 5, ncol = 2, pdf_out = pdf_dual)
}








#' process pair-experiment files
#' one TE per page
#' read bigWig files and return in data.frame format
#'
#' @param ctl_ip_bw list of bigWig files, control ip
#'
#' @param ctl_input_bw list of bigWig files, control input
#'
#' @param tre_ip_bw list of bigWig files, treatment ip
#'
#' @param tre_input_bw list of bigWig files, treatment input
#'
#' @param fasize the chromosome size file of the bigWig, tab separated
#'
#' @param n a list of chromosomes names to extract from the bigWig files,
#'   default, TRUE, use all seqnames
#'
#' @param pdf_out logical, save the plots in pdf file
#'
#' @import dplyr
#'
#' @export
#'
chipseq_bw_parser4 <- function(ctl_ip_bw, tre_ip_bw, fasize,
                               ctl_input_bw = NULL, tre_input_bw = NULL,
                               ctl_label = NULL, tre_label = NULL,
                               n = TRUE, pdf_out = NULL) {
  stopifnot(length(ctl_ip_bw) == length(tre_ip_bw))

  # pick chr
  chr_gr  <- fasize_reader(fasize, to_GRanges = TRUE)
  chr_all <- as.character(seqnames(chr_gr))
  if(isTRUE(n)) {
    chr_hit <- chr_all
  } else {
    chr_hit <- n
  }

  if(is.null(pdf_out)) {
    # pdf_out <- as.character(glue::glue("ChIPseq.{ctl_label}.vs.{tre_label}.{n_chr}.{n_stage}.trackview.pdf"))
    pdf_out <- getwd()
  }

  if(! dir.exists(pdf_out)) {
    dir.create(pdf_out)
  }

  for(chr_n in chr_hit){
    pdf_n_name <- paste0("trackview-", chr_n, ".pdf")
    pdf_n      <- file.path(pdf_out, pdf_n_name)

    # make plots for dual samples
    plist <- lapply(seq_len(length(ctl_ip_bw)), function(i){
      # stage
      t1 <- stringr::str_match(ctl_ip_bw[i], "_(\\w+)_NSR")[1, 2]
      t2 <- stringr::str_match(tre_ip_bw[i], "_(\\w+)_NSR")[1, 2]
      stopifnot(t1 == t2)

      # strand
      direction <- basename(dirname(ctl_ip_bw[i]))

      # read bw
      df_dual <- chipseq_bw_parser2(ctl_ip_bw = ctl_ip_bw[i],
                                    tre_ip_bw = tre_ip_bw[i],
                                    fasize = fasize,
                                    ctl_input_bw = ctl_input_bw[i],
                                    tre_input_bw = tre_input_bw[i],
                                    ctl_label = ctl_label,
                                    tre_label = tre_label,
                                    n = chr_n)
      df <- df_dual[[1]]

      # add stage
      df$stage  <- t1
      df$direction <- direction

      # make plot
      p <- coverage_plot_dual(df, title = "stage")

      # return
      return(p)
    })

    plot_n_pages(plist, nrow = 5, ncol = 2, pdf_out = pdf_n)
  }
}












#' choose theme for coverage plot
#'
#' @param x a integer, 1 for single coverage plot, 2 for overlay of
#'   dual-coverage plot
#'
#' @import ggplot2
#'
theme_picker <- function(x) {
  # single track
  t1 <- theme_linedraw() +
    theme(panel.border = element_rect(fill = NULL, color = "black", size = .7),
          panel.grid   = element_blank(),
          plot.title   = element_text(color = "black", hjust = .5, size = 12),
          axis.line    = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_line(color = "black", size = .7),
          axis.title.y = element_text(color = "black", size = 10, vjust = .5))
  # overlay tracks
  # t2 <- t1 + theme(legend.position = c(0.8, 0.8))
  t2 <- t1 + theme(legend.position = "top")
  # choose theme
  tm <- list(t1, t2)
  names(tm) <- seq_len(length(tm))
  if (x %in% names(tm)) {
    return(tm[[x]])
  } else {
    return(tm[[1]])
  }
}


#' generate coverage plot for single chromosome
#'
#' @param data a data.frame contains position and coverage values, require
#'   "sample" to group the tracks
#'
#' @param fill.color the color for fill and line of coverage plot,
#'   default: orange
#'
#' @param exclude.minus.scores logical, whether exclude minus values in plot
#'
#' @import ggplot2
#' @import ggridges
#'
#' @export
#'
coverage_plot_single <- function(data, fill.color = "orange",
                                 y.scale.max = "auto",
                                 exclude.minus.scores = TRUE) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("name", "position", "score") %in% names(data)))

  # number of records of data
  if (nrow(data) == 0) {
    stop("empty data.frame in input")
  }

  # exclude minus values
  if (isTRUE(exclude.minus.scores)) {
    data <- dplyr::filter(data, score >= 0)
  }

  # number of names
  n_names <- as.character(unique(data$name))
  if (length(n_names) > 1) {
    stop("multiple names detected in data.frame, exiting...")
  }

  # convert id_name to name
  if (grepl("_", n_names)) {
    n_names <- unlist(strsplit(n_names, "_"))[2]
  }

  # add y scale for ridgeline
  data$y <- 1

  # chooose theme
  mytheme <- theme_picker(1)

  # determine y-max
  y.max <- ceiling(max(data$score / 10)) * 10
  if (is.numeric(y.scale.max) & y.scale.max > y.max) {
    y.max <- y.scale.max
  }

  # generate plot
  p1 <- ggplot(data, aes(position, y, height = score,
                         color = name, fill = name)) +
    ggridges::geom_ridgeline(show.legend = FALSE) +
    guides(fill = guide_legend(title = NULL, label.position = "right")) +
    scale_y_continuous(limits = c(0, y.max)) +
    ylab("base coverage [rpm]") +
    ggtitle(n_names) +
    scale_fill_manual(values = fill.color) +
    scale_color_manual(values = fill.color) +
    mytheme

  return(p1)
}


#' generate coverage plot for dual samples
#'
#' @param data a data.frame contains position and coverage values, require
#'   "sample" to group the tracks
#'
#' @param fill.color the color for fill and line of coverage plot,
#'   default: orange
#'
#' @param exclude.minus.scores logical, whether exclude minus values in plot
#'
#' @param y.scale.max numerical or "auto", specify the y-max in plot.
#'
#' @import ggplot2
#' @import ggridges
#'
#' @export
#'
coverage_plot_dual <- function(data,
                               fill.color = c("orange", "red2"),
                               exclude.minus.scores = TRUE,
                               y.scale.max = "auto",
                               title = "name") {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("sample", "name", "position", "score") %in% names(data)))
  stopifnot(length(fill.color) == 2)

  # number of records of data
  if (nrow(data) == 0) {
    warning("empty data.frame in input")
    return(NULL)
  }

  # exclude minus values
  if (isTRUE(exclude.minus.scores)) {
    data <- dplyr::filter(data, score >= 0)
  }

  # number of names
  n_names <- as.character(unique(data$name))
  if (length(n_names) > 1) {
    stop("multiple names detected in data.frame, exiting...")
  }

  # convert id_name to name
  if (grepl("_", n_names)) {
    n_names <- unlist(strsplit(n_names, "_"))[2]
  }

  if("stage" %in% colnames(data)) {
    n_stage <- as.character(data$stage)[1]
    if("direction" %in% colnames(data)) {
      direction <- as.character(data$direction)[1]
      n_stage <- paste0(n_stage, "_", direction)
    }
  } else {
    n_stage <- "stage"
  }

  if(title == "name") {
    plot_title = n_names
  } else if(title == "stage") {
    plot_title = n_stage
  } else {
    plot_title = n_chr
  }

  # number of samples
  n_samples <- unique(data$sample)
  if (! length(n_samples == 2)) {
    stop("2 samples in data.frame are expected, failed, exiting...")
  }

  # add y scale for ridgeline
  data$y <- 1

  # chooose theme
  mytheme <- theme_picker(2)

  # determine y-max
  y.max <- ceiling(max(data$score / 10)) * 10
  if (is.numeric(y.scale.max) & y.scale.max > y.max) {
    y.max <- y.scale.max
  }

  # generate plot
  p2 <- ggplot(data, aes(position, y, height = score,
                         color = sample, fill = sample)) +
    ggridges::geom_ridgeline(show.legend  = TRUE, alpha = .5) +
    guides(color = guide_legend(title     = NULL),
           fill  = guide_legend(title     = NULL,
                                label.vjust = 1,
                                keywidth  = .6,
                                keyheight = .2,
                                label.position = "right")) +
    scale_y_continuous(limits = c(0, y.max)) +
    ylab("base coverage [rpm]") +
    ggtitle(plot_title) +
    scale_fill_manual(values = fill.color) +
    scale_color_manual(values = fill.color) +
    mytheme

  return(p2)
}


#' save multiple plots in one page
#'
#' generate PDF file for multiple plots
#' arrange multiple n_row x n_col plots in one page
#'
#' @param plist a list of ggplot objects
#'
#' @param nrow Number of rows in plot grid, default: 1
#'
#' @param ncol Number of cols in plot grid, default: 1
#'
#' @param page.width int, default 8
#'
#' @param page.height int, default 10
#'
#' @param pdf_out The pdf file to save the plots
#'
#'
#' @import cowplot
#' @import ggplot2
#'
#' @export
#'
plot_n_pages <- function(plotlist = NULL, nrow = 2, ncol = 5,
                         page.width = 8, page.height = 10,
                         pdf_out = NULL) {
  # make a list of plots
  plots <- plotlist
  num_plots <- length(plots)

  # determine number of pages
  num_pages <- ceiling(length(plots) / (nrow * ncol))

  # determine the indexes of plots for each page
  page_index <- lapply(seq_len(num_pages), function(i){
    s <- nrow * ncol * (i - 1) + 1
    e <- nrow * ncol * i
    e <- ifelse(e > num_plots, num_plots, e)
    return(s:e)
  })

  # check pdf file
  if (is.null(pdf_out)){
    pdf_out <- tempfile(pattern = "plot_", tmpdir = ".", fileext = ".pdf")
  } else {
    dir.create(dirname(pdf_out), showWarnings = FALSE, mode = "0755")
  }

  # generate plots
  pdf(pdf_out, width = page.width, height = page.height, paper = "a4")
  for (n in seq_len(num_pages)){
    idx <- page_index[[n]]
    idx_plist <- plots[idx]
    pg <- cowplot::plot_grid(plotlist = idx_plist, align = "hv",
                             nrow = nrow, ncol = ncol, labels = "auto")
    print(cowplot::ggdraw(pg))
  }
  dev.off()

  # return the pdf filename
  return(pdf_out)
}


# EOF
