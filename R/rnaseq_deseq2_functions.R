
# R code for DESeq2 analysis
#
# Author: Ming Wang
# Date: 2018-09-30
#
# DESeq2: v1.18.1



#' run DESeq2 for matrix data
#'
#' @param file Path to the matrix file, featureCount output
#'
#' @param species Name of the species, dm3, dm6, mm9, mm10, hh19, default: dm3
#'
#' @param outdir Path to directory to save results, default: cwd
#'
#' @export
#'
#' @describeIn Suppose multiple BAM files in featureCounts output
#'   require at least four samples, name of file endswith "_rep[1-3]"
#'
DESeq2_for_featureCounts <- function(x, organism = "dm3", outdir = "./",
                                     pvalue_cutoff = 0.05) {
  # read counts
  df_count <- featureCountsReader(x, organism)

  # fetch replicate samples
  col_reps <- grep("_rep[1-9]", names(df_count), ignore.case = TRUE)
  df_reps  <- df_count[, col_reps]
  if (length(df_reps) < 4) {
    stop("less than 4 BAM files of replicate samples")
  }

  # make pairs: Control vs experiment
  g <- unique(gsub("_rep[1-9]$", "", names(df_reps), ignore.case = TRUE))
  g_pairs <- combn(g, 2)

  # DE analysis
  for (i in seq_len(ncol(g_pairs))) {
    smp <- g_pairs[, i]
    rep_names <- names(df_reps)
    g1  <- rep_names[startsWith(rep_names, smp[1])]
    g2  <- rep_names[startsWith(rep_names, smp[2])]
    smp_line <- paste(smp, collapse = "_vs_")
    smp_dir  <- file.path(outdir, smp_line)
    if (! dir.exists(smp_dir)) {
      dir.create(smp_dir, recursive = TRUE, mode = "0750")
    }

    # make countdata
    ma <- as.matrix(df_reps[, c(g1, g2)])

    # make coldata
    coldata <- design_example(colnames(ma), smp)

    # run DESeq2
    deseq2_run(ma, coldata, smp_dir, pvalue_cutoff)

    ##------------------------------------------------------------------------##
    # rename gene id
    fs <- file.path(smp_dir, "transcripts_deseq2.csv")
    fsFix <- file.path(smp_dir, "transcripts_deseq2.fix.xls")
    ## save table
    fName = paste("genelist", organism, "rda", sep = ".")
    f = system.file("extdata", fName, package = "goldclipReport")
    load(f) # genelist

    ## read data
    df <- DESeq2_csv2df(fs)
    df2 <- dplyr::mutate(df, id = as.character(id)) %>%
      dplyr::mutate(id = plyr::mapvalues(id, genelist$gene_id, genelist$gene_name, FALSE)) %>%
      dplyr::filter(! is.na(padj)) %>%
      dplyr::arrange(padj)

    readr::write_delim(df2, fsFix, delim = "\t", col_names = TRUE)
    ##----------------------------------------------------------------------------##
  }
}



#' read featureCount matrix file
#'
#' @param x Path the matrix file
#'
#' @param normalizeTo1M Logical value, whether or not normalize the counts of
#'   each BAM file to 1 million reads. default: FALSE
#'
#' @import dplyr
#' @import readr
#' @import tidyr
#' @import tibble
#'
#'
#' @export
#'
featureCountsReader <- function(x, organism = "dm3", normalizeTo1M = FALSE,
                                fixZero = 0) {
  # parse file
  df_exp <- readr::read_delim(x, "\t",
                              col_types = readr::cols(),
                              comment = "#") %>%
    dplyr::rename(id = Geneid) %>% # "id" column
    dplyr::select(-(2:6)) %>% # remove extra columns
    as.data.frame()

  # convert wide to long table
  df_exp2 <- tidyr::gather(df_exp, key = sample, value = count, -id)

  # fix Zero values
  if (fixZero > 0) {
    df_exp2$count[df_exp2$count == 0] <- fixZero
  }

  # normalize to 1 Million reads (RPM)
  if (isTRUE(normalizeTo1M)) {
    dfNorm <- dplyr::group_by(df_exp2, sample) %>%
      mutate(norm = count / sum(count) * 1e6) %>%
      dplyr::rename(value = norm) %>%
      dplyr::select(-count)
  } else {
    dfNorm <- dplyr::rename(df_exp2, value = count)
  }

  # convert count to int (--fraction option in featureCounts)
  dfNorm$value <- round(dfNorm$value, 0)

  # output
  dfOut <- dplyr::ungroup(dfNorm) %>%
    dplyr::select(id, sample, value) %>%
    tidyr::spread(key = sample, value = value) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("id")

  # reorder columns
  origin_header <- names(df_exp)[-1]
  dfOut <- dfOut[, origin_header]

  # rename colnames
  smps <- basename(colnames(dfOut))
  ext  <- str_common(smps, suffix = TRUE)
  colnames(dfOut) <- gsub(ext, "", smps)

  # rename rownames
  genelist <- get_genelist(organism)
  stopifnot(all(c("gene_id", "gene_name") %in% names(genelist)))
  # rownames(dfOut) <- plyr::mapvalues(rownames(dfOut),
  #                                    from = genelist$gene_id,
  #                                    to = genelist$gene_name,
  #                                    warn_missing = FALSE)
  # # remove temp objects
  rm(df_exp)
  rm(df_exp2)
  rm(dfNorm)

  # report
  return(dfOut)
}



#' run DESeq2 analysis
#'
#' @param ma count data in matrix
#'
#' @param coldata experiment design, in data.frame format
#'
#' @param outdir Directory to save the results
#'
#' @param pvalue_cutoff Cutoff to filt the records, padj for DESeq2 output,
#'   default: 0.05
#'
#'
#' @import DESeq2
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
#'
deseq2_run <- function(ma, coldata, outdir, pvalue_cutoff = 0.05) {
  stopifnot(is.matrix(ma))

  # prepare files
  de_count <- file.path(outdir, "transcripts_deseq2.csv")
  de_plots <- file.path(outdir, c("figure1_MA_plot.png",
                                  "figure2_MA_plot_LFC.png",
                                  "figure3_sample_counts.png",
                                  "figure4_PCA_plot.png",
                                  "figure5_dispersion.png",
                                  "figure6_sample_distance.png",
                                  "figure7_top_genes.png",
                                  "figure8_volcano.png"))

  # prepare experiment design
  countdata <- ma
  smp <- unique(gsub("_rep[1-9]$", "", colnames(countdata), ignore.case = TRUE))
  coldata <- design_example(ids = colnames(countdata), conditions = smp)

  # load matrix
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata, colData = coldata,
                                        design = ~condition)

  # Run the DESeq pipeline
  dds <- DESeq2::DESeq(dds)

  # Get differential expression results
  res <- DESeq2::results(dds)
  resLFC <- DESeq2::lfcShrink(dds, coef = 2, res = res)
  rld <- DESeq2::rlogTransformation(dds, blind = FALSE)
  ntd <- DESeq2::normTransform(dds)

  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  resSig <- subset(as.data.frame(res), padj < pvalue_cutoff)

  # Normalied counts
  ncount <- DESeq2::counts(dds, normalized = TRUE)

  ## Merge with normalized count data
  resdata <- merge(as.data.frame(ncount), as.data.frame(res), by = "row.names",
                   sort = FALSE)
  resdata <- resdata[order(resdata$padj), ]
  names(resdata)[1] <- "Gene"

  # save data to file
  write.csv(resdata, de_count, quote = FALSE, row.names = TRUE)

  # MA
  png(de_plots[1], width = 1200, height = 1200, res = 300)
  DESeq2::plotMA(res, ylim = c(-2, 2))
  dev.off()

  # MA for LFC
  png(de_plots[2], width = 1200, height = 1200, res = 300)
  DESeq2::plotMA(resLFC, ylim = c(-2, 2))
  dev.off()

  # Sample counts
  png(de_plots[3], width = 1200, height = 1200, res = 300)
  DESeq2::plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")
  dev.off()

  # PCA
  png(de_plots[4], width = 2000, height = 2000, res = 300)
  print(DESeq2::plotPCA(rld, intgroup = c("condition")))
  dev.off()

  # Dispersion
  png(de_plots[5],  width = 1500, height = 1500, res = 300)
  DESeq2::plotDispEsts(dds)
  dev.off()

  # Sample distance
  png(de_plots[6], width = 1000, height = 1000, res = 300)
  sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rld$condition
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
                     clustering_distance_rows = sampleDists,
                     clustering_distance_cols = sampleDists,
                     col = colors)
  dev.off()

  # Top genes
  png(de_plots[7], width = 1200, height = 1200, res = 300)
  select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
                  decreasing = TRUE)[1:30]
  ma <- SummarizedExperiment::assay(ntd)[select, ]
  df <- as.data.frame(SummarizedExperiment::colData(dds)[, c("condition")])
  colnames(df) <- "condition"
  rownames(df) <- colnames(ma)
  pheatmap::pheatmap(SummarizedExperiment::assay(ntd)[select,],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     show_rownames = FALSE,
                     annotation_col = df)
  dev.off()
}



#' return the longest common string
#'
#' @param list A list of strings
#'
#' @param suffix A logical value, if TRUE, return the longest suffix of strings
#'   default: FALSE
#'
#' @export
#'
str_common <- function(list, suffix = FALSE) {
  # return the longest common string between list of strings
  # default: prefix
  if (isTRUE(suffix)) {
    list <- stringi::stri_reverse(list)
  }

  # split string
  dl <- lapply(list, function(s) {
    unlist(strsplit(s, ''))
  })

  # min len
  len_min <- min(unlist(lapply(dl, length)))

  # check common
  cs = c()
  for (i in seq_len(len_min)) {
    s <- unique(unlist(lapply(dl, function(ds){ds[i]})))

    if (length(s) > 1){
      break
    }
    cs <- c(cs, s)
  }

  if (isTRUE(suffix)) {
    cs <- rev(cs)
  }

  # output
  ss <- paste0(cs, collapse = "")
  return(ss)
}



#' create experiment design for DESeq2 analysis
#'
#' @param ids A set of variables, the name of samples, default: letters[1:4]
#'
#' @param conditions A set of variables, the levels for design, default:
#'   letters[1:2]
#'
#' @export
design_example <- function(ids = letters[1:4], conditions = letters[1:2]) {
  rep = length(ids) / length(conditions)
  if ( floor(rep) < rep) {
    warning("Number of ids not divided by conditinos")
  }
  df <- data.frame(condition = factor(rep(conditions, each = rep),
                                      levels = conditions))
  rownames(df) <- ids
  return(df)
}



#' import matrix to DESeq2
#'
#' @param df A data.frame of count matrix
#'
#' @param coldata A data.frame of design
#'
#' @export
#'
DESeq2_data_from_matrix <- function(df, coldata) {
  stopifnot(is.data.frame(df))
  if(! "condition" %in% colnames(coldata)) {
    print(colnames(coldata))
    stop("condition - not found in coldata")
  }

  # subset of dataframe
  dfx <- dplyr::select(df, rownames(coldata))
  if (length(dfx) != nrow(coldata)) {
    stop("data.frame and coldata not match")
  }

  # load DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(df, coldata, ~condition)
  dds <- DESeq2::DESeq(dds)
  return(dds)
}



#' fetch gene id and name from Ensembl
#' using biomaRt package
#'
#' @param organism A character of organism, default: dm3
#'
#' @import biomaRt
#'
#' @export
#'
#' @description to-do, support all avilable organism
#'  gene_id, gene_name/external_gene_id
#'
#' load from package data,
#' optional: download from ensembl, update
get_genelist <- function(organism = "dm3") {
  # library(biomaRt)
  print(organism)
  f_name = paste("genelist", organism, "rda", sep = ".")
  f = system.file("extdata", f_name, package = "goldclipReport")
  if (file.exists(f)) {
    print(f)
    load(f) # genelist
    return(genelist)
  } else {
    stop('genelist file not exists')
  }
}



#' cal mean of replicates
#' default, 2 replicates for each condition
#'
#' @param data A data fraem of DESeq2 output, first 5 columns:
#'   <id> <a-rep1> <a-rep2> <b-rep1> <b-rep2> ...
#'
#' @export
DESeq2_add_mean <- function(data){
  # only for DESeq2 output csv
  # insert to: col-2, col-3
  stopifnot(is.data.frame(data))
  ctl_cols <- colnames(data)[2:3]
  tre_cols <- colnames(data)[4:5]
  ctl_name <- str_common(ctl_cols)
  tre_name <- str_common(tre_cols)
  ctl_name <- gsub("_R|_rep", "", ctl_name, ignore.case = TRUE)
  tre_name <- gsub("_R|_rep", "", tre_name, ignore.case = TRUE)
  # calculate mean
  df_mean  <- data.frame(id  = data$id,
                         ctl = rowMeans(data[, 2:3]),
                         tre = rowMeans(data[, 4:5]))
  colnames(df_mean) <- c("id", ctl_name, tre_name)
  # merge
  df <- merge(df_mean, data, by = "id")
  return(df)
}



#' read output of DESeq2_run function
#' csv format
#'
#' @param file Path to csv file
#'
#'
#' @import readr
#' @import dplyr
#'
#'
#'
#' @export
DESeq2_csv2df <- function(x) {
  # read DESeq2 csv file
  # library(readr)
  # library(dplyr)
  df <- readr::read_csv(x, col_types = readr::cols())
  # colnames(df)[1] <- 'id'
  df[, 1] <- NULL # remove the first column
  colnames(df)[1] <- "id"
  df$id <- as.character(df$id)
  # add mean columns
  df2 <- DESeq2_add_mean(df)
  return(df2)
}



#
# de_deseq2_plot <- function(df, withTE = FALSE, add_sigNames = TRUE,
#                            extra_ids = c("CG9754", "nxf2", "piwi")) {
#   s1 <- colnames(df)[2:3]
#   df[, s1[1]] <- log10(df[, s1[1]])
#   df[, s1[2]] <- log10(df[, s1[2]])
#   # gene list
#   dfGene <- filtGeneId(df, gene = TRUE, te = FALSE)
#   # filter > 0
#   # dfGene <- dfGene[is.finite(dfGene$shattp2) & is.finite(dfGene$shNxf2), ]
#   dfGene <- dplyr::filter_(dfGene, s1[1] > 0 & s1[2] > 0)
#   p1 <- de_scatter(dfGene, s1[1], s1[2],
#                    add_sigNames = add_sigNames,
#                    extra_ids = extra_ids)
#   # p2 <- de_ma(dfGene, add_sigNames = add_sigNames, extra_ids = extra_ids)
#   # print(p1)
#   # print(p2)
#   # pout <- list(p1, p2)
#   pout <- list(p1)
#
#   # te list
#   if (isTRUE(withTE)) {
#     dfTE <- filtGeneId(df, gene = FALSE, te = TRUE)
#     p3 <- de_scatter(dfTE, s1[1], s1[2], extra_ids = extra_ids)
#     pout <- c(pout, p3)
#     # print(p3)
#   }
#   return(pout)
# }

#
#
#
# de_scatter <- function(df, xName, yName, add_sigNames = TRUE,
#                        extra_ids = c("piwi", "CG9754"),
#                        tt = "Differential expression", cutoff = 0.05) {
#   library(ggplot2)
#   library(ggrepel)
#   # at most, 10 labels
#   dfSig <- dplyr::filter(df, padj <= cutoff) %>%
#     dplyr::arrange(padj)
#   dfEx <- dplyr::filter(df, id %in% extra_ids)
#   if (isTRUE(add_sigNames)) {
#     dfLabel <- rbind(dfEx, dfSig)
#   } else {
#     dfLabel <- dfEx
#   }
#   if (nrow(dfLabel) > 3) dfLabel <- dfLabel[1:3, ]
#   # tt = paste(xName, yName, sep = "_")
#   tt = glue::glue("DE-analysis: {xName} vs {yName}")
#
#   p <- ggplot(df, aes_string(xName, yName)) +
#     geom_abline(slope = 1, intercept = 0, color = "grey10") +
#     geom_abline(slope = 1, intercept = log10(2),
#                 color = "grey30", linetype = 2) +
#     geom_abline(slope = 1, intercept = -log10(2),
#                 color = "grey30", linetype = 2) +
#     geom_point(size = 1.3,
#                color = ifelse(is.na(df$padj), "grey20",
#                               ifelse(df$padj <= cutoff, "red", "grey20"))) +
#     # ggtitle(tt) +
#     xlab(xName) + ylab(yName) +
#     scale_x_continuous(limits = c(0, 6),
#                        breaks = seq(0, 6, 1),
#                        labels = seq(0, 6, 1),
#                        expand = c(0, 0, 0, 0)) +
#     scale_y_continuous(limits = c(0, 6),
#                        breaks = seq(0, 6, 1),
#                        labels = seq(0, 6, 1),
#                        expand = c(0, 0, 0, 0)) +
#     theme_classic() +
#     theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
#           plot.title = element_text(color = "black", hjust = .5, size = 14),
#           axis.title = element_text(color = "black", size = 12),
#           axis.text = element_text(color = "black", size = 10))
#   if (nrow(dfLabel) > 0) {
#     p <- p + geom_text_repel(
#       data          = dfLabel,
#       label         = dfLabel$id,
#       box.padding   = 1,
#       segment.size  = 0.4,
#       segment.color = "grey50",
#       direction     = "both"
#     )
#   }
#   return(p)
# }



#
# matrix_DESeq2_batch <- function(fn, genelist, outDir) {
#   # only for featureCount output
#   stopifnot(file.exists(fn))
#   stopifnot(is.data.frame(genelist))
#   stopifnot(all(c("gene_id", "gene_name") %in% names(genelist)))
#   df <- fc2df(fn)
#   df <- dplyr::mutate(df,
#                       id = plyr::mapvalues(id,
#                                            genelist$gene_id,
#                                            genelist$gene_name,
#                                            warn_missing = FALSE))
#   # recursive pair
#   smp_ids <- names(df)[-1]
#   smp_ids <- gsub("_rep\\d+", "", smp_ids)
#   smp_ids <- unique(smp_ids)
#   # combinations
#   cc <- combn(smp_ids, 2)
#   for (i in seq(ncol(cc))) {
#     smp <- cc[, i]
#     smp_line <- paste(smp, collapse = "_vs_")
#     g1 <- stringr::str_which(names(df), smp[1])
#     g2 <- stringr::str_which(names(df), smp[2])
#     ## subset of data.frame
#     dfx <- dplyr::select_(df, .dots = c(g1, g2)) %>% as.data.frame()
#     rownames(dfx) <- df$id
#     ## make coldata
#     coldata <- data.frame(row.names = colnames(dfx),
#                           condition = rep(c("wt", "mut"), each = 2))
#     dds <- load_matrix_to_DESeq2(dfx, coldata)
#     ## out_dir
#     dfDir <- file.path(outDir, smp_line)
#     if (! dir.exists(dfDir)) {
#       dir.create(dfDir, recursive = TRUE)
#     }
#     oldDir <- getwd()
#     setwd(dfDir)
#     run_DESeq2(dds)
#     setwd(oldDir)
#     print(smp_line)
#   }
# }
#
#
# kallistoToDESeq2Batch <- function(fn, outDir, tx2gene = NULL,
#                                   gene_level = FALSE) {
#   # only for kallisto output
#   # directory structure
#   # suppose 2 replicates for each sample
#   # path-to-kallisto/prjname/abundance.h5
#   fn <- sort(fn)
#   # combination
#   smp_ids <- gsub("RNAseq_|_rep\\d+", "", basename(dirname(fn)))
#   smp_ids <- unique(smp_ids)
#   # combinations
#   cc <- combn(smp_ids, 2)
#   for (i in seq(ncol(cc))) {
#     smp <- cc[, i]
#     smp_line <- paste(smp, collapse = "_vs_")
#     g1 <- stringr::str_which(fn, smp[1])
#     g2 <- stringr::str_which(fn, smp[2])
#     # sub files
#     fx <- fn[c(g1, g2)]
#     dds <- kallistoToDESeq2(fx, tx2gene, gene_level = TRUE)
#     ## out_dir
#     dfDir <- file.path(outDir, smp_line)
#     if (! dir.exists(dfDir)) {
#       dir.create(dfDir, recursive = TRUE)
#     }
#     oldDir <- getwd()
#     setwd(dfDir)
#     run_DESeq2(dds)
#     setwd(oldDir)
#     print(smp_line)
#   }
# }

#
# tetoolkit2DESeq2 <- function(fn) {
#   #  f <- "../tetranscript/dsCG9754_vs_dsPiwi/dsCG9754_vs_dsPiwi.cntTable"
#   df <- readr::read_delim(fn, "\t", quote = "\"", col_types = cols())
#   colnames(df)[1] <- 'id'
#   dfx <- df[, -1] %>% as.data.frame()
#   n <- gsub(".Aligned.sortedByCoord.out.bam|.C|.T|RNAseq_", "", basename(names(dfx)))
#   rownames(dfx) <- df$id
#   colnames(dfx) <- n
#   coldata <- data.frame(row.names = colnames(dfx),
#                         condition = rep(c("mut", "wt"), each = 2))
#   dds <- load_matrix_to_DESeq2(dfx, coldata)
#   ## out_dir
#   dfDir <- dirname(fn)
#   oldDir <- getwd()
#   setwd(dfDir)
#   run_DESeq2(dds)
#   setwd(oldDir)
#   print(basename(fn))
# }


#
# fc2df <- function(fn, normalize = FALSE,
#                   convertLog10 = FALSE,
#                   convertInt = TRUE) {
#   # parse featureCounts
#   # convert count to INT
#   # convert count to log10
#   # normalize, RPM
#   library(readr)
#   library(dplyr)
#   library(tidyr)
#   df <- read_delim(fn, "\t", col_types = cols(), comment = "#") %>%
#     dplyr::rename(id = Geneid)
#   df1 <- df[, -c(2:6)]
#
#   # normalization
#   # RPM: reads per million
#   if (isTRUE(normalize)) {
#     df2 <- gather(df1, sample, count, -1) %>%
#       dplyr::group_by(sample) %>%
#       mutate(norm = count / sum(count) * 1e6) %>%
#       dplyr::rename(value = norm)
#   } else {
#     df2 <- gather(df1, sample, count, -1) %>%
#       dplyr::rename(value = count)
#   }
#   #
#   if (isTRUE(convertLog10)) {
#     df2 <- dplyr::mutate(df2, value = log10(value))
#   }
#
#   if(isTRUE(convertInt)) {
#     df2 <- dplyr::mutate(df2, value = round(value, 0))
#   }
#
#   df3 <- ungroup(df2) %>%
#     dplyr::select(id, sample, value) %>%
#     mutate(sample = gsub("RNAseq_|.piRNA_clusters|.unique|.merged|.Aligned|.sortedByCoord|.out|.bam", "", basename(sample))) %>%
#     spread(sample, value)
#   return(df3)
# }
#


#
# #' import kallisto to DESeq2
# #'
# #' @param kal_files A set of files
# #'
# #' @param coldata A data.frame of experiment design
# #'
# #'
# #' @import tximport
# #'
# #' @export
# #'
# DESeq2_data_from_kallisto <- function(kal_files, coldata, t2g) {
#
#   library(tximport)
#   kal_files <- kal_files[file.exists(kal_files)]
#   # transcript-level
#   txi <- tximport(files = kal_files, type = "kallisto", txOut = TRUE,
#                   dropInfReps = TRUE)
#   # # gene-level
#   # txi <- tximport(files = kal_files, type = "kallisto", txOut = FALSE,
#   #                 tx2gene = t2g, dropInfReps = TRUE)
#   rownames(txi$counts) <- as.character(rownames(txi$counts))
#   # checkout dataset
#   if(length(kal_files) != nrow(coldata)) {
#     stop('kallisto file not match coldata')
#   }
#   if(! "condition" %in% colnames(coldata)) {
#     print(colnames(coldata))
#     stop("condition - not found in coldata")
#   }
#   dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata,
#                                   design = ~condition)
#   dds <- DESeq(dds)
#   return(dds)
# }



#
# filtGeneId <- function(df, gene = TRUE, te = TRUE) {
#   # filt output of TEtranscripts
#   if (isTRUE(gene) & isTRUE(te)) {
#     return(df) # do not filt
#   } else if (isTRUE(gene)) {
#     # id, do not contain ":"
#     gids <- grepl(":", df$id)
#     return(df[! gids, ])
#   } else if(isTRUE(te)) {
#     # id, only contain ":"
#     teids <- grep(":", df$id)
#     return(df[teids, ])
#   } else {
#     warning("do not return records")
#   }
# }




#' run regular DESeq2 analysis
#'
#' @param dds A variable of dds
#'
#'
#' @import ggplot2
#' @import DESeq2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
#'
# DESeq2_run <- function(dds, path_out, pval_cutoff = 0.05) {
#
#   # prepare files
#   file_count <- file.path(path_out, "transcripts_deseq2.csv")
#   file_plot <- file.path(path_out, c("figure1_MA_plot.png",
#                                      "figure2_MA_plot_LFC.png",
#                                      "figure3_sample_counts.png",
#                                      "figure4_PCA_plot.png",
#                                      "figure5_dispersion.png",
#                                      "figure6_sample_distance.png",
#                                      "figure7_top_genes.png",
#                                      "figure8_volcano.png"))
#   #
#   #   # prepare variables
#   #   library(DESeq2, quietly = TRUE)
#   #   library("ggplot2")
#   #   library(pheatmap)
#   #   library("RColorBrewer")
#
#   # start
#   res <- DESeq2::results(dds)
#   resOrdered <- res[order(res$padj), ]
#   # print(head(resOrdered))
#   # resSig <- dplyr::filter(resOrdered, padj < pval_cutoff)
#   resSig <- subset(as.data.frame(resOrdered), padj < pval_cutoff)
#   resLFC <- DESeq2::lfcShrink(dds, coef = 2, res = res)
#   rld <- DESeq2::rlog(dds, blind = FALSE)
#   #vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
#   #vsd.fast <- DESeq2::vst(dds, blind = FALSE)
#   ntd <- DESeq2::normTransform(dds)
#   #add normalized counts
#   ncount <- DESeq2::counts(dds, normalized = TRUE)
#   exp <- cbind(as.data.frame(ncount), as.data.frame(res))
#   expOrdered <- exp[order(exp$padj), ]
#
#   # save data to file
#   write.csv(expOrdered, file_count, quote = FALSE, row.names = TRUE)
#
#   # make plots
#   # MA plot
#   png(file_plot[1], width = 1200, height = 1200, res = 300)
#   DESeq2::plotMA(res, ylim = c(-2, 2))
#   dev.off()
#
#   png(file_plot[2], width = 1200, height = 1200, res = 300)
#   DESeq2::plotMA(resLFC, ylim = c(-2, 2))
#   dev.off()
#
#   png(file_plot[3], width = 1200, height = 1200, res = 300)
#   DESeq2::plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")
#   dev.off()
#
#   # PCA
#   png(file_plot[4], width = 2000, height = 2000, res = 300)
#   print(DESeq2::plotPCA(rld, intgroup = c("condition")))
#   dev.off()
#
#   png(file_plot[5],  width = 1500, height = 1500, res = 300)
#   DESeq2::plotDispEsts(dds)
#   dev.off()
#
#   # Sample distance
#   png(file_plot[6], width = 1000, height = 1000, res = 300)
#   sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
#   sampleDistMatrix <- as.matrix(sampleDists)
#   rownames(sampleDistMatrix) <- rld$condition
#   colnames(sampleDistMatrix) <- NULL
#   colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
#   pheatmap::pheatmap(sampleDistMatrix,
#                      clustering_distance_rows = sampleDists,
#                      clustering_distance_cols = sampleDists,
#                      col = colors)
#   dev.off()
#
#   # Count matrix of top genes
#   png(file_plot[7], width = 1200, height = 1200, res = 300)
#   select <- order(rowMeans(DESeq2::counts(dds, normalized = TRUE)),
#                   decreasing = TRUE)[1:30]
#   ma <- SummarizedExperiment::assay(ntd)[select, ]
#   df <- as.data.frame(SummarizedExperiment::colData(dds)[, c("condition")])
#   colnames(df) <- "condition"
#   rownames(df) <- colnames(ma)
#   pheatmap::pheatmap(SummarizedExperiment::assay(ntd)[select,],
#                      cluster_rows = FALSE,
#                      cluster_cols = FALSE,
#                      show_rownames = FALSE,
#                      annotation_col = df)
#   dev.off()
# }

