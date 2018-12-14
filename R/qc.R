
#' QCBaseContent
#'
#' Make plots for base-content
#'
#' @param path directory of the base-content
#' @param name name of the output file
#'
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#'
#' @export
QCBaseContent <- function(path, name) {
  f1 <- file.path(path, "base_content.txt")
  f2 <- file.path(path, "length_distribution.txt")

  if (! file.exists(f1) & file.exists(f2)) {
    stop("file not exists: base_content.txt, length_distribution.txt")
  }

  ##----------------------------------------------------------------------------##
  df <- read.table(f1, header = TRUE, sep = " ")
  df[is.na(df)] <- 0 # NA to 0
  df <- tidyr::gather(df, base, count, -position)
  df <- dplyr::filter(df, base %in% c("A", "C", "G", "T")) %>%
    dplyr::mutate(base = factor(base, levels = c("A", "C", "G", "T"))) %>%
    dplyr::group_by(position) %>%
    dplyr::mutate(freq = count / sum(count))

  p1 <- ggplot(df, aes(position, freq, fill = base)) +
    geom_bar(position = "fill", stat = "identity") +
    xlab("Position in read (bp)") +
    ylab("Per base sequence content (%)") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    scale_fill_manual(values = c("green3", "blue3", "grey30", "red3")) +
    theme_classic()

  p2 <- ggplot(df, aes(position, freq, color = base)) +
    geom_line(size = .5) +
    #geom_bar(position = "fill", stat = "identity") +
    xlab("Position in read (bp)") +
    ylab("Per base sequence content (%)") +
    scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
    scale_color_manual(values = c("green3", "blue3", "grey30", "red3")) +
    theme_classic()

  ##----------------------------------------------------------------------------##
  df2 <- read.table(f2, header = FALSE, sep = " ")
  names(df2) <- c("length", "count")
  # if only one record in df2
  if (nrow(df2) == 1) {
    n <- df2[1, "length"]
    df2x <- data.frame(length = c(n - 1, n + 1),
                       count = 0)
    df2 <- rbind(df2, df2x)
  }

  p3 <- ggplot(df2, aes(length, count)) +
    geom_line(size = .5, color = "red3") +
    xlab("Length of reads (nt)") +
    ylab("Number of reads") +
    theme_classic()

  ##----------------------------------------------------------------------------##
  png1 <- file.path(path, paste0(name, ".base_content.bar.png"))
  png2 <- file.path(path, paste0(name, ".base_content.line.png"))
  png3 <- file.path(path, paste0(name, ".length_distribution.line.png"))
  ggsave(png1, p1, width = 5, height = 3, units = "in")
  ggsave(png2, p2, width = 5, height = 3, units = "in")
  ggsave(png3, p3, width = 5, height = 3, units = "in")
  ##----------------------------------------------------------------------------##
  return(list(p1, p2, p3))
}



#' fastqcFiles
#'
#' This function parsing *_fastqc.zip files from input path.
#' It assumes that input path is a directory or a single *_fastqc.zip
#' file.
#'
#' @param qc.path Path to the input directory
#' @return A list of paths of zip files
#' @export
FastqcFiles <- function (qc.path) {
  stopifnot(length(qc.path) == 1)
  if (qc.path == ".")
    qc.path <- getwd()
  if (dir.exists(qc.path)) {
    qc.files <- list.files(qc.path, pattern = "*_fastqc.zip", full.names = TRUE)
    if (length(qc.files) == 0)
      stop("Can't find any file that match: ", theBasename)
  } else if (! dir.exists(qc.path) & file.exists(qc.path)) {
    qc.files <- qc.path
  } else if (! dir.exists(qc.path) & ! file.exists(qc.path)) {
    theDirname  <- dirname(qc.path)
    theBasename <- basename(qc.path)
    qc.files    <- list.files(theDirname,
                              pattern = paste0("^", theBasename, "*_fastqc.zip"),
                              full.names = TRUE)
    if (length(qc.path) == 0)
      stop("Can't find any file that match: ", theBasename)
  } else if (! dir.exists(qc.path) |  ! file.exists(qc.path)) {
    stop("Specified QC path doesn't exist.")
  }
  return(qc.files)
}

#' parse fastQC output
#'
#' @param qc.file the zip file of fastqc ouptut
#'
#' @import fastqcr
#'
#' @export
FastqcPlot <- function (qc.file) {
  stopifnot(length(qc.file) == 1)
  stopifnot(file.exists(qc.file))
  qc <- fastqcr::qc_read(qc.file)
  p1 <- fastqcr::qc_plot(qc, "Per base sequence quality")
  p2 <- fastqcr::qc_plot(qc, "Per base sequence content")
  p3  <- cowplot::plot_grid(p1, p2, labels = c("A", "B"))
  # now add the title
  prefix <- gsub("_fastqc.zip", "", basename(qc.file))
  title  <- cowplot::ggdraw() + cowplot::draw_label(prefix, fontface = 'bold')
  p <- plot_grid(title, p3, ncol = 1, rel_heights = c(0.1, 1))
  # rel_heights values control title margins
  return(p)
}



#' FastqcReport
#'
#' This function parsing *_fastqc.zip files from input path.
#' It assumes that input path is a directory or a single *_fastqc.zip
#' file.
#'
#' @description Create html report for fastqc output
#' @param qc.path Path to the input directory, Default is the
#'   current working directory.
#' @param result.file path to the result file prefix (e.g., path/to/qc-result).
#'   Don't add the file extension.
#' @param experiment text specifying a short description of the experiment. For
#'  example experiment = "goldCLIP seqeuencing for AGO2".
#' @param template a character vector specifying the path to an Rmd template.
#'  file.
#' @param preview logical value. If TRUE, shows a preview of the report.
#' @examples
#' \donotrun{
#' # Demo
#' qc.path <- system.file("fastqc_results", package = "fastqcr")
#'
#' FastqcReport(qc.path, result.file = "demo")
#' }
#'
#' @import fastqcr
#' @import rmarkdown
#'
#' @return A list of paths of zip files
#' @export
FastqcReport <- function (qc.path, result.file, experiment = NULL,
                          template = NULL, preview = TRUE) {
  if (qc.path == ".")
    qc.path <- getwd()
  if (! dir.exists(qc.path) & ! file.exists(qc.path)) {
    theDirname  <- dirname(qc.path)
    theBasename <- basename(qc.path)
    qc.path     <- list.files(theDirname,
                              pattern = paste0("^", theBasename, "*_fastqc.zip"),
                              full.names = TRUE)
    if (length(qc.path) == 0)
      stop("Can't find any file that match: ", theBasename)
    qc.path <- qc.path[1]
  }
  if (! dir.exists(qc.path) | !file.exists(qc.path)) {
    stop("Specified QC path doesn't exist.")
  }
  qc.path <- normalizePath(qc.path)
  if(!dir.exists(dirname(result.file)))
    dir.create(dirname(result.file), showWarnings = FALSE, recursive = TRUE)
  oldwd <- getwd()
  setwd(dirname(result.file))
  result.file <- paste0(basename(result.file), ".html")
  result.file <- file.path(getwd(), result.file)
  if (is.null(template)) {
    report_template <- system.file("report_templates",
                                   "fastqc_report_template.Rmd",
                                   package = "goldclipData")
  } else {
    report_template <- template
  }
  rmarkdown::render(input = report_template, output_file = result.file,
                    params = list(qc.path = qc.path, experiment = experiment))
  if (preview) {
    utils::browseURL(result.file)
  }
  message("\n--------------------------\nOutput file: ", result.file,
          "\n--------------------------\n")
  setwd(oldwd)
}



#
#
# funcA <- function(s) {
#   s   <- gsub("^\n", "", s)
#   s1  <- unlist(strsplit(s, "\n", TRUE))
#   tag <- unlist(strsplit(gsub(">>", "", s1[1]), "\t"))[1] # title
#   header <- unlist(strsplit(gsub("^#", "", s1[2]), "\t")) # header
#   s2 <- s1[-c(1:2)]
#   df <- tidyr::separate(data.frame(name = s1[-c(1:2)]), name, header, sep = "\t")
#   t <- list(tag = df)
#   names(t) <- tag
#   return(t)
# }
#
# f <- "demo_fastqc.zip"
# group <- "Per base sequence quality"
#
#
# # tempdir
# dd <- tempdir()
# #file.path(gsub(".zip$", "", basename(f)), "fastqc_data.txt")
# unzip(f, exdir = dd)
# fd <- list.files(dd, "fastqc_data.txt", recursive = TRUE, full.names = TRUE)
# line <- paste(readLines(fd)[-1], collapse = "\n")
# a <- unlist(strsplit(line, ">>END_MODULE", fixed = TRUE))
# qc <- lapply(a, funcA)
# qc2 <- list()
# for(i in seq(length(qc))){
#   qc2 <- c(qc2, qc[[i]])
# }
# ## title
# if(group %in% names(qc2)) {
#   df <- qc2[[group]]
# }
#
#
# library(ggplot2)
# ## plots
# # 1. base-quality
# df1 <- qc2[["Per base sequence quality"]]
#
# # 2. base-content
# df2 <- qc2[["Per base sequence content"]] %>%
#   tidyr::gather("base", "freq", 2:5) %>%
#   dplyr::rename(position = Base) %>%
#   dplyr::mutate(freq = as.numeric(freq),
#                 position = factor(position, levels = unique(position)),
#                 base = factor(base, levels = c("A", "C", "G", "T")))
# # barplot
# p2a <- ggplot(df2, aes(position, freq, fill = base)) +
#   geom_bar(position = "fill", stat = "identity") +
#   xlab("Position in read (bp)") +
#   ylab("Per base sequence content (%)") +
#   scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
#   scale_fill_manual(values = c("green3", "blue3", "grey30", "red3")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
# # lineplot
# p2b <- ggplot(df2, aes(position, freq, color = base, group = base)) +
#   geom_point(size = .7) +
#   geom_line(size = .7) +
#   xlab("Position in read (bp)") +
#   ylab("Per base sequence content (%)") +
#   scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = seq(0, 100, 20)) +
#   scale_color_manual(values = c("green3", "blue3", "grey30", "red3")) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
#
#







