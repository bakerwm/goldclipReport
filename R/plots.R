





#' map_stat_plot
#'
#' @param x data.frame of mapping stat
#' @param stat identity or percentage
#'
#' @import ggplot2
#' @import ggridges
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @export
map_stat_plot <- function(x, stat = "percentage") {

  if(inherits(x, "character")){
    df <- map_stat_read(x)
  } else if(inherits(x, "data.frame")) {
    df <- x
  } else {
    stop("unknown stat file")
  }

  stopifnot(all(c("id", "type", "count") %in% names(df)))

  if(stat == "percentage") {
    if(! "pct" %in% names(df)) {
      df <- dplyr::group_by(df, id) %>%
        dplyr::mutate(pct = count / sum(count) * 100)
    }
    p1 <- ggplot(df, aes(x = id, y = pct, fill = type)) +
      geom_bar(stat = "identity", color = "black", size = .5)
  } else if(stat == "identity") {
    p1 <- ggplot(df, aes(x = id, y = count, fill = type)) +
      geom_bar(stat = "identity", color = "black", size = .5)
  } else {
    stop("unknown stat=")
  }

  p2 <- p1 +
    coord_flip() +
    xlab("") + ylab("Percentage %") +
    scale_x_discrete(limits = rev(levels(df$id))) +
    scale_y_continuous(position = "right") +
    expand_limits(x = 0, y = 0) +
    scale_fill_discrete(guide = guide_legend(
      nrow = 1,
      reverse = TRUE,
      direction = "horizontal",
      title = NULL,
      title.position = "top",
      label.position = "bottom",
      label.hjust = 0.5,
      label.vjust = 1,
      label.theme = element_text(angle = 0))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = .5),
          legend.key.width = unit(0.4, "cm"),
          legend.key.height = unit(0.4, "cm"),
          # panel.border = element_blank(),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(color = "black", size = 12),
          axis.line = element_line(color = "black", size = .7),
          axis.ticks = element_line(color = "black", size = 1),
          panel.grid = element_blank(), legend.position = "top")

  return(p2)
}


#' map_stat_read
#'
#' @param x file or list of files of mapping stat
#' @param origin where keep the original format
#'
#' @importFrom readr cols
#' @importFrom readr read_delim
#' @importFrom magrittr %>%
#'
#'
#'
#' @export
map_stat_read <- function(x, origin = FALSE) {

  if(inherits(x, "character")) {
    dx <- lapply(x, function(f){
      readr::read_delim(f, "\t", col_types = readr::cols())
    })
    df1 <- dplyr::bind_rows(dx)
  } else {
    stop("unknown x=")
  }

  # default header
  t1 <- c("genome_rRNA", "genome_unique", "genome_multiple", "spikein_rRNA",
          "spikein_unique", "spikein_multi", "unmap")
  t2 <- c("rRNA", "unique", "multiple", "rRNA.sp", "unique.sp", "multiple.sp", "unmap")
  names(df1) <- plyr::mapvalues(names(df1), from = t1, to = t2, warn_missing = FALSE)
  colnames(df1)[1] <- "id" # rename

  if(! isTRUE(origin)) {
    dfx <- df1 %>%
      dplyr::select(-total) %>%
      tidyr::gather(key = "type", value = "count", -1) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(pct = count / sum(count) * 100)
    if(all(t2 %in% dfx$type)) {
      dfx <- dplyr::mutate(dfx, type = factor(type, levels = rev(t2)))
    } else {
      dfx <- dplyr::mutate(dfx, type = factor(type, levels = sort(unique(type))))
    }
    return(dfx)
  } else {
    return(df1)
  }
}



#' BedLenBarplot
#'
#' @param data data.frame of BED file
#'
#' @import ggplot2
#' @import ggridges
#'
#' @export
BedLenBarplot <- function(data) {
  # length distribution of length
  stopifnot(any(class(data) %in% "data.frame"))
  stopifnot(all(c("id", "length", "count") %in% names(data)))
  # convert formats
  df <- dplyr::group_by(data, id) %>%
    mutate(density = count / sum(count)) %>%
    ungroup() %>%
    mutate(id = as.character(id),
           length = as.numeric(length))
  # plot
  p <- ggplot(dfy, aes(x = length, y = id, fill = id)) +
    geom_density_ridges(scale = 2) + theme_ridges() +
    xlab("Length of peaks (nt)") + ylab(NULL) +
    theme(panel.grid = element_blank())
  return(p)
}


#' BedNumBarplot
#'
#' barplot for numbers
#'
#' @param data data.frame of BED, contains length, count, id columns
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
BedNumBarplot <- function(data) {
  # barplot for number of peaks
  # input: data.frame, <length> <count> <id>
  stopifnot(any(class(data) %in% "data.frame"))
  stopifnot(all(c("id", "length", "count") %in% names(data)))
  # convert formats
  df <- dplyr::group_by(data, id) %>%
    dplyr::summarise(total = sum(count))
  nmax <- 10^nchar(max(df$total) - 1)
  ymax <- ceiling(max(df$total) / nmax) * nmax
  # plot
  p <- ggplot(df, aes(id, total)) +
    geom_bar(stat = "identity", fill = "grey80", color = "black", size = .5) +
    geom_text(aes(label = total), vjust = .5, hjust = -0.1) +
    scale_y_continuous(limits = c(0, ymax),
                       breaks = seq(0, ymax, length.out = 5),
                       labels = seq(0, ymax / nmax, length.out = 5)) +
    coord_flip() +
    xlab(NULL) + ylab(paste0("Number of peaks (x", nmax, ")")) +
    theme_bw() +
    theme(panel.grid   = element_blank(),
          panel.border = element_blank(),
          plot.title   = element_text(color = "black", size = rel(1.5),
                                      face = "bold", hjust = .5),
          axis.line    = element_line(color = "black", size = .5),
          axis.title   = element_text(color = "black", size = rel(1.2)),
          axis.text.x  = element_text(color = "black", size = rel(1.2)),
          axis.text.y  = element_text(color = "black", size = rel(1.2)))
  return(p)
}


#' BedAnnoBarplot
#'
#' @param data data.frame of the mapping stat
#' @param type percentile or count; identity, percentage
#' @category reads, peaks
#'
#' @import ggplot2
#' @import dplyr
#'
#' @return a ggplot object for the plot
#'
#' @export
BedAnnoBarplot <- function(data, type = "identity", genome = "hg19",
                           category = "reads") {
  stopifnot(any(class(data) %in% "data.frame"))
  stopifnot(all(c("id", "type", "count") %in% names(data)))
  # remove na rows
  data <- data[complete.cases(data$type), ]
  #mytheme <- pickThemes(type = "bar1")
  if(genome == "dm3") {
    mycols  <- ColorPicker("anno2", 10)
  } else {
    mycols  <- ColorPicker("anno1", 9)
    if(category == "peaks") {
      mycols <- ColorPicker("anno2", 10)
    }
  }
  mytheme <- theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.border = element_blank(),
          plot.title   = element_text(color = "black", size = rel(1.5),
                                      face = "bold", hjust = .5),
          axis.line    = element_line(color = "black", size = .5),
          axis.title   = element_text(color = "black", size = rel(1.2)),
          axis.text.x  = element_text(color = "black", size = rel(1.2), angle = 90,
                                      vjust = .5, hjust = 1),
          axis.text.y  = element_text(color = "black", size = rel(1.2)),
          strip.text.x = element_text(size = rel(1.2)),
          strip.text.y = element_text(angle = 0, size = rel(1.2)),
          strip.background = element_rect(fill = NA, colour = NA))
  ## determine ymax
  data <- dplyr::group_by(data, id) %>%
    mutate(total = sum(count)) %>%
    mutate(pct = count / total * 100)
  nmax <- 10^nchar(max(data$total) - 1)
  ymax <- ceiling(max(data$total) / nmax) * nmax
  p1 <- ggplot(data, aes(x = id, y = count, fill = type)) +
    geom_bar(stat = "identity", colour = "grey10", size = .3) +
    scale_y_continuous(limits = c(0, ymax),
                       breaks = seq(0, ymax, length.out = 5),
                       labels = seq(0, ymax / nmax, length.out = 5)) +
    xlab(NULL) + ylab(paste0("Number of count (x", nmax, ")")) +
    scale_fill_manual(values = mycols) +
    mytheme
  p2 <- ggplot(data, aes(x = id, y = pct, fill = type)) +
    geom_bar(stat = "identity", colour = "grey10", size = .3) +
    scale_fill_manual(values = mycols) +
    scale_y_continuous(breaks = seq(0, 100, by = 20),
                       labels = seq(0, 100, by = 20)) +
    xlab(NULL) + ylab("Percentage (%)") +
    mytheme
  if(type == "identity") {
    return(p1)
  } else if(type == "percentage") {
    return(p2)
  } else {
    stop("unknown option for type:")
  }
}



#' align_report
#'
#' This function parsing *_fastqc.zip files from input path.
#' It assumes that input path is a directory or a single *_fastqc.zip
#' file.
#'
#' @description Create html report for alignment
#' @param input mapping_stat.csv files
#' @param output Path to the result file. alignment_stat.html
#'   Don't add the file extension.
#' @param template a character vector specifying the path to an Rmd template.
#'  file.
#' @param preview logical value. If TRUE, shows a preview of the report.
#' @examples
#' \donotrun{
#' # Demo
#'
#' }
#'
#' @import readr
#' @import rmarkdown
#'
#' @return A list of paths of zip files
#' @export
align_report <- function (input, output, template = NULL, preview = TRUE) {
  # list files
  fx <- mapply(function(f){
    if(file.exists(f)) {
      return(f)
    }
  }, input)

  if(length(fx) == 0){
    stop("stat files not detected!")
  }

  if(! dir.exists(output)){
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
  }

  ## output
  output <- normalizePath(output)
  align_html <- file.path(output, "alignment_report.html")

  if (is.null(template)) {
    report_template <- system.file("report_templates",
                                   "alignment_report_template.Rmd",
                                   package = "goldclipReport")
  } else {
    report_template <- template
  }

  oldwd <- getwd()
  setwd(dirname(output))

  rmarkdown::render(input = report_template, output_file = align_html,
                    params = list(stat_files = fx))

  if (preview) {
    utils::browseURL(result.file)
  }
  message("\n--------------------------\nOutput file: ", align_html,
          "\n--------------------------\n")
  setwd(oldwd)
}




