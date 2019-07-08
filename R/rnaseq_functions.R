
# R code for DESeq2 analysis
#
# Author: Ming Wang
# Date: 2018-09-30
#
# DESeq2: v1.18.1



#' run DESeq2 for matrix data
#'
#' @param countA data.frame or list of BAM files for control
#' @param countB data.frame or list of BAM files for treatment
#' @param organism Name of the species, dm3, dm6, mm9, mm10, hh19, default: dm3
#' @param nameA the name of control, default prefix of countA
#' @param nameB the name of treatment, default prefix of countB
#' @param outdir Path to directory to save results, default: cwd
#' @param pvalue_cutoff float value, default: 0.1
#'
#' @description Suppose multiple BAM files in featureCounts output
#'   require at least four samples, name of file endswith "_rep[1-3]"
#' @export
deseqHub <- function(countA, countB, organism = "dm3",
                     nameA = NULL, nameB = NULL,
                     outdir = "./", pvalue_cutoff = 0.1,
                     readable = False) {
  if(inherits(countA, "data.frame") & inherits(countB, "data.frame")) {
    ## at least two replicates
    t1 <- dim(countA) > c(0, 2) & dim(countB) > c(0, 2)
    t2 <- "id" %in% names(countA) & "id" %in% names(countB)
    if(! all(t1, t2)) {
      stop("require: >= 2 replicates, id in header")
    }
    ## nameA
    n1 <- str_common(basename(colnames(countA)))
    n2 <- str_common(basename(colnames(countB)))
    nameA <- ifelse(is.null(nameA), n1, nameA)
    nameB <- ifelse(is.null(nameB), n1, nameB)
    ## combine data.frame
    df <- merge(countA, countB, by = "id")
  } else if(inherits(countA, "character") & inherits(countB, "character")) {
    ## replicates are required for each sample/condition
    df_1 <- featureCountsParser(countA, organism)
    df_2 <- featureCountsParser(countB, organism)
    df   <- merge(df_1, df_2, by = "id")
  }

  ## expression matrix
  ma <- as.matrix(tibble::column_to_rownames(df, "id"))

  ## run deseq2
  deseqNode(ma, outdir, pvalue_cutoff)

  ## update gene names
  cnt_csv <- file.path(outdir, "transcripts_deseq2.csv")
  cnt_xls <- file.path(outdir, "transcripts_deseq2.fix.xls")

  ## add mean values
  df_csv <- deseqCsvParser(cnt_csv)
  df_new <- deseqAddMeanValue(df_csv)

  ## sort by pvalue
  df_filt <- df_new %>%
    dplyr::filter(! is.na(padj)) %>% # remove na rows
    dplyr::arrange(padj)

  ## readable
  if(isTRUE(readable)) {
    df_filt$id <- deseqIdConvertor(df_filt$id,
                                   organism = organism,
                                   from = "gene_id",
                                   to = "gene_name")
  }

  readr::write_delim(df_filt, cnt_xls, delim = "\t", col_names = TRUE)
}


#' run DESeq2 analysis
#'
#' @param ma count data in matrix
#' @param coldata experiment design, in data.frame format
#' @param outdir Directory to save the results
#' @param pvalue_cutoff Cutoff to filt the records, padj for DESeq2 output,
#'   default: 0.05
#'
#' @import DESeq2
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
#'
deseqNode <- function(ma, outdir, pvalue_cutoff = 0.05) {
  stopifnot(is.matrix(ma))
  if(! dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  }

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
  # matrix data
  countdata <- ma

  # prepare experiment design
  smp     <- colnames(countdata)
  cond_1  <- smp[1:(length(smp)/2)]
  cond_2  <- smp[-1:-(length(smp)/2)]
  coldata <- deseqDesign(cond_1, cond_2)

  # load matrix to DESeq2
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = countdata,
                                        colData = coldata,
                                        design = ~condition)

  # Run the DESeq pipeline
  dds <- DESeq2::DESeq(dds)

  # Get differential expression results
  res    <- DESeq2::results(dds)
  resLFC <- DESeq2::lfcShrink(dds, coef = 2, res = res, type = "apeglm")
  rld    <- DESeq2::rlogTransformation(dds, blind = FALSE)
  ntd    <- DESeq2::normTransform(dds)

  # Order by adjusted p-value
  res    <- res[order(res$padj), ]
  resSig <- subset(as.data.frame(res), padj < pvalue_cutoff)

  # Normalied counts
  ncount <- DESeq2::counts(dds, normalized = TRUE)

  ## Merge with normalized count data
  resdata <- merge(as.data.frame(ncount), as.data.frame(res), by = "row.names",
                   sort = FALSE)
  resdata <- resdata[order(resdata$padj), ]
  names(resdata)[1] <- "Gene"

  # save data to file
  write.csv(resdata, de_count, quote = TRUE, row.names = TRUE)

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



#' create coldata for kallisto
#' the dirname of kallisto output
#'
#'
#' @export
kalColdata <- function(file, condition=NULL) {
  fnames <- basename(dirname(file))
  # condition
  if(is.null(condition)) {
    condition <- gsub("_rep\\d+$|_r\\d+$", "", fnames, perl = TRUE)
  }
  stopifnot(length(file) == length(condition))
  # data.frame
  coldata <- data.frame(condition = condition,
                        stringsAsFactors = TRUE)
  rownames(coldata) <- fnames

  return(coldata)
}



#' load from kallisto
#' using tximport
#'
#' @param file h5 file of kallisto output
#' @param coldata  experiment design, data.frame
#' @param format the file format for input, salmon, default: kallisto
#'
#' @imort tximport
#'
#' @export
deseqImportData <- function(file, coldata, format = "kallisto", geneLevel = TRUE,
                            organism = "fruitfly") {
  stopifnot(inherits(coldata, "data.frame"))
  ## t2g
  t2g_file <- system.file("extdata", "t2g.rda", package = "goldclipReport")
  load(t2g_file) # t2g
  organism <- tolower(organism)
  if(organism %in% c("dm3", "fruitfly")) {
    t2g <- t2g_dm3
  } else if(organism %in% c("hg19", "human")) {
    t2g <- t2g_hg19
  } else if(organism %in% c("hg38")) {
    t2g <- t2g_hg38
  } else if(organism %in% c("mm10", "mouse")) {
    t2g <- t2g_mm10
  } else {
    stop(paste0("unknown organism: ", organism))
  }
  t2g <- dplyr::select(t2g, target_id, ext_gene)

  ## import file
  if(isTRUE(geneLevel)) {
    txi <- tximport::tximport(file, type = format, tx2gene = t2g)
  } else {
    txi <- tximport::tximport(file, type = format)
  }
  stopifnot(rownames(coldata) == colnames(txi$counts))

  ## convert to DESeqData
  dds <- DESeqDataSetFromTximport(txi, coldata, ~condition)

  ## remove zeros
  dds <- dds[rowSums(counts(dds)) > 1, ]

  return(dds)
}



#' run DESeq2 for matrix data
#'
#' @param countA data.frame or list of BAM files for control
#' @param countB data.frame or list of BAM files for treatment
#' @param organism Name of the species, dm3, dm6, mm9, mm10, hh19, default: dm3
#' @param nameA the name of control, default prefix of countA
#' @param nameB the name of treatment, default prefix of countB
#' @param outdir Path to directory to save results, default: cwd
#' @param pvalue_cutoff float value, default: 0.1
#'
#' @description Suppose multiple BAM files in featureCounts output
#'   require at least four samples, name of file endswith "_rep[1-3]"
#' @export
deseqHub2 <- function(countA, countB, organism = "dm3",
                      nameA = NULL, nameB = NULL, geneLevel = TRUE,
                      outdir = "./", pvalue_cutoff = 0.1) {
  if(inherits(countA, "data.frame") & inherits(countB, "data.frame")) {
    ## at least two replicates
    t1 <- dim(countA) > c(0, 2) & dim(countB) > c(0, 2)
    t2 <- "id" %in% names(countA) & "id" %in% names(countB)
    if(! all(t1, t2)) {
      stop("require: >= 2 replicates, id in header")
    }
    ## nameA
    n1 <- str_common(basename(colnames(countA)))
    n2 <- str_common(basename(colnames(countB)))
    nameA <- ifelse(is.null(nameA), n1, nameA)
    nameB <- ifelse(is.null(nameB), n1, nameB)
    ## combine data.frame
    df <- merge(countA, countB, by = "id")
    ## expression matrix
    ma <- as.matrix(tibble::column_to_rownames(df, "id"))
    ## run deseq2
    deseqNode(ma, outdir, pvalue_cutoff)
    ## kallisto output
  } else if(all(grepl("abundance.h5$", c(countA, countB)))) {
    kal_files <- c(countA, countB)
    ## coldata
    coldata <- kalColdata(kal_files)
    ## factor levels
    coldata$condition <- factor(coldata$condition, levels = c(nameA, nameB))
    dds <- deseqImportData(kal_files, coldata, "kallisto", geneLevel = geneLevel)
    deseqNode2(dds, outdir, pvalue_cutoff)
  } else if(inherits(countA, "character") & inherits(countB, "character")) {
    ## replicates are required for each sample/condition
    df_1 <- featureCountsParser(countA, organism)
    df_2 <- featureCountsParser(countB, organism)
    df   <- merge(df_1, df_2, by = "id")
    ## expression matrix
    ma <- as.matrix(tibble::column_to_rownames(df, "id"))
    ## run deseq2
    deseqNode(ma, outdir, pvalue_cutoff)
  }

  ## update gene names
  cnt_csv <- file.path(outdir, "transcripts_deseq2.csv")
  cnt_xls <- file.path(outdir, "transcripts_deseq2.fix.xls")

  ## add mean values
  df_csv <- deseqCsvParser(cnt_csv)
  df_new <- deseqAddMeanValue(df_csv)

  ## sort by pvalue
  df_filt <- df_new %>%
    dplyr::filter(! is.na(padj)) %>% # remove na rows
    dplyr::arrange(padj)

  readr::write_delim(df_filt, cnt_xls, delim = "\t", col_names = TRUE)
}


#' run DESeq2 analysis
#' import: dds
#' gene-level
#'
#' @param ma count data in matrix
#' @param coldata experiment design, in data.frame format
#' @param outdir Directory to save the results
#' @param pvalue_cutoff Cutoff to filt the records, padj for DESeq2 output,
#'   default: 0.05
#'
#' @import DESeq2
#' @import ggplot2
#' @import pheatmap
#' @import RColorBrewer
#' @import SummarizedExperiment
#'
#' @export
#'
deseqNode2 <- function(dds, outdir, pvalue_cutoff = 0.05) {
  stopifnot(inherits(dds, "DESeqDataSet"))
  if(! dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  }

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
  # Run the DESeq pipeline
  dds <- DESeq2::DESeq(dds)

  # save table
  res    <- DESeq2::results(dds)

  # model test
  # Get differential expression results
  resLFC <- DESeq2::lfcShrink(dds, coef = 2, res = res, type = "apeglm")
  rld    <- DESeq2::rlogTransformation(dds, blind = FALSE)
  ntd    <- DESeq2::normTransform(dds)

  ## save to file
  res    <- res[order(res$padj), ]   # Order by adjusted p-value
  resSig <- subset(as.data.frame(res), padj < pvalue_cutoff)
  ncount <- DESeq2::counts(dds, normalized = TRUE)   # Normalied counts
  ## Merge with normalized count data
  resdata <- merge(as.data.frame(ncount), as.data.frame(res), by = "row.names",
                   sort = FALSE)
  resdata <- resdata[order(resdata$padj), ]
  names(resdata)[1] <- "Gene"
  # save data to file
  write.csv(resdata, de_count, quote = TRUE, row.names = TRUE)


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



#' import matrix to DESeq2
#'
#' @param df A data.frame of count matrix
#'
#' @param coldata A data.frame of design
#'
#' @export
#'
deseqFromMatrix <- function(x, coldata){
  # DESeq2_data_from_matrix <- function(df, coldata) {
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



#' create experiment design for DESeq2 analysis
#'
#' @param name_control vector, a list of character, names of control
#' @param name_treatment vector, a list of character, names of treatment
#'
#' @export
deseqDesign <- function(name_control, name_treatment) {
  # only support 2 conditions
  conditions <- c("control", "treatment")
  condC <- rep("control", times = length(name_control))
  condT <- rep("treatment", times = length(name_treatment))
  df <- data.frame(condition = factor(c(condC, condT),
                                      levels = conditions))
  ## names
  rownames(df) <- c(name_control, name_treatment)
  ## return
  return(df)
}



#' read output of DESeq2_run function
#' csv format
#'
#' @param x Path to csv file
#' @import readr
#'
#' @export
deseqCsvParser <- function(x) {
  # read DESeq2 csv file
  df <- readr::read_csv(x, col_types = readr::cols())
  # colnames(df)[1] <- 'id'
  df[, 1] <- NULL # remove the first column
  colnames(df)[1] <- "id"
  df$id <- as.character(df$id)
  return(df)
}



#' cal mean of replicates
#'
#' @param x the csv file of DESeq2 output, default
#'
#' @export
deseqAddMeanValue <- function(x) {
  # only for DESeq2 output csv
  if(is.character(x)){
    df <- deseqCsvParser(x) # convert to data.frame
  } else if(is.data.frame(x)){
    df <- x
  } else {
    stop("unknown file format")
  }
  # remove the last 6-columns
  # determine the replicates
  ix <- length(df) - 6
  smp_names <- colnames(df)[2:ix] # split into two groups
  smp_1  <- smp_names[1:(length(smp_names)/2)]
  smp_2  <- smp_names[-1:-(length(smp_names)/2)]
  name_1 <- str_common(smp_1)
  name_2 <- str_common(smp_2)
  # select sub_groups
  df_1   <- dplyr::select(df, .dots = smp_1)
  df_2   <- dplyr::select(df, .dots = smp_2)
  # create data.frame
  df_mean <- data.frame(id    = df$id,
                        nameA = rowMeans(df_1),
                        nameB = rowMeans(df_2))
  colnames(df_mean) <- c("id", name_1, name_2)
  # merge data.frame
  df_new  <- merge(df_mean, df, by = "id")
  return(df_new)
}



#' read featureCount matrix file
#'
#' @param x Path the matrix file
#' @param organism character, the name of organism, default: dm3
#' @param normalizeTo1M Logical value, whether or not normalize the counts of
#'   each BAM file to 1 million reads. default: FALSE
#' @param fixZero float, to replicate the zero values in matrix
#'
#' @import dplyr
#' @import readr
#' @import tidyr
#' @import tibble
#'
#' @export
#'
featureCountsParser <- function(x, organism = "dm3", normalizeTo1M = FALSE,
                                fixZero = 0, readable = FALSE) {
  # count values
  df_exp <- readr::read_delim(x, "\t",
                              col_types = readr::cols(),
                              comment = "#") %>%
    dplyr::rename(id = Geneid) %>% # "id" column
    dplyr::select(-(2:6)) %>% # remove extra columns
    as.data.frame() %>%
    tibble::column_to_rownames("id")

  # fix Zero values
  if (fixZero > 0) {
    df_exp[df_exp == 0] <- fixZero
  }

  # normalize to 1 Million reads (RPM)
  norm_scales <- 1e6 / colSums(df_exp)

  if (isTRUE(normalizeTo1M)) {
    df_norm <- data.frame(mapply(`*`, df_exp, norm_scales, SIMPLIFY=FALSE))
  } else {
    df_norm <- df_exp
  }

  # convert count to int (--fraction option in featureCounts)
  df_norm <- round(df_norm, 0) # int

  # rename colnames
  smps <- basename(colnames(df_norm))
  ext  <- str_common(smps, suffix = TRUE)
  colnames(df_norm) <- gsub(ext, "", smps)

  # rename rownames
  if(isTRUE(readable)) {
    rownames(df_norm) <- deseqIdConvertor(rownames(df_norm),
                                          organism,
                                          from = "gene_id",
                                          to = "gene_name")
  }

  # report
  df_out <- tibble::rownames_to_column(df_norm, "id")
  return(df_out)
}



#' ID convertor
#' from id to name
#' only support goldclipReport
#'
#' @param x vector, a list of gene ids
#' @param organism character, the name of organism, support dm3, hg19, mm10, hg38, dm6
#' @param from the type fo ids, gene_id, refseq_mrna
#'
#' @import plyr
#' @import dplyr
#'
#'
#' @export
deseqIdConvertor <- function(x, organism, from = "gene_id", to = "gene_name") {
  # library(biomaRt)
  # print(organism)
  # retrieve data from goldclipReport package
  f_name = paste0("genelist.", organism, ".rda")
  f = system.file("extdata", f_name, package = "goldclipReport")
  if (file.exists(f)) {
    load(f) # genelist
  } else {
    stop('genelist file not exists')
  }
  # convert
  from_list <- genelist %>% dplyr::pull(from)
  to_list   <- genelist %>% dplyr::pull(to)
  # for gencode gene_id ENSG00000000000.X
  # trim the .x suffix
  x <- gsub("\\.\\d+$", "", x, perl = TRUE)
  y <- plyr::mapvalues(x,
                       from = from_list,
                       to   = to_list,
                       warn_missing = FALSE)

  # return
  return(y)
}



#' read output of DESeq2
#' return data.frame
#'
#' @param x xls or csv
#'
#' @import readr
#' @import dplyr
#'
#' @export
#'
deseqReader <- function(x, sep = ",") {
  df <- readr::read_delim(x, sep, col_types = readr::cols())
  return(df)
}






#' count reads
#' using featureCounts from Rsubread
#' @param bam_files a character vector, BAM files
#' @param organism a character vector, the name of reference genome, dm3, mm9, mm10, hg19, hg38
#' @param outdir a character vector, the directory of output results
#' @param strand an integer vector indicating if strand-specific read count should be performed. 0 (unstranded), 1 (stranded), 2 (reversely stranded). Default: 0
#' @param gtf_feature a character string giving the feature type used to select rows in the GTF annotation
#' @param gtf_attr a character string specifying extra GTF attribute type in the GTF annotation
#' @param threads an integer vector, the number of threads to use, Default: 1
#' @param data_path a character string, the path the genome data files of species, Default: ~/data/genome
#'
#' @import Rsubread
#'
#' @description count reads on features using featureCounts from Rsubread
#'
#' @export
#'
#' -M -O --fraction -g gene_id -t exon -T %s -a %s -s %s
# fc_count <- function(bam_files, organism, outdir, strand = 0,
#                      gtf_feature = "exon", gtf_attr = "gene_id",
#                      threads = 1, data_path = "~/data/genome") {
#   gtf_name <- paste0(organism, ".ensembl.gtf")
#   gtf_file <- file.path(data_path, organism, "annotation_and_repeats", gtf_name)
#   stopifnot(file.exists(gtf_file))
#
#   fx = Rsubread::featureCounts(files                  = bam_files,
#                                annot.ext              = gtf_file, # -a
#                                isGTFAnnotationFile    = TRUE, #
#                                GTF.featureType        = gtf_feature, # -g
#                                GTF.attrType           = gtf_attr, # -t
#                                allowMultiOverlap      = TRUE, # -M
#                                countMultiMappingReads = TRUE, # -O
#                                fraction               = TRUE, # --fraction
#                                strandSpecific         = strand, # -s
#                                nthreads               = threads )# -T
#
#   # save to file
#   ## 1, count, gene counts matrix, gene + bam
#   ## 2, annotation, gene feature data.frame
#   ## 3, targets, bam list
#   ## 4, stat, summary
#   outdir <- "deseq"
#   fc_file_names <- c("count.txt", "annotation.txt", "bam_list.txt", "summary.txt")
#   fc_file <- file.path(outdir, fc_file_names)
#   rda_file <- file.path(outdir, "fc_out.rda")
#
#   # count
#   df1 <- as.data.frame(fx[[1]], row_names = TRUE)
#   fc_count <- tibble::rownames_to_column(df1, "id")
#   readr::write_csv(fc_count, fc_file[1], col_names = TRUE)
#
#   # annotation
#   fc_anno <- fx[[2]]
#   readr::write_csv(fc_anno, fc_file[2], col_names = TRUE)
#
#   # bam_list
#   fc_bam <- fx[[3]]
#   readr::write_lines(fc_bam, fc_file[3], sep = "\t")
#
#   # stat
#   fc_stat <- fx[[4]]
#   readr::write_delim(fc_stat, fc_file[4], delim = "\t")
#
#   # r-object
#   fc_count <- df1x
#   save(fc_count, fc_anno, file = rda_file, compress = "gzip")
#
#   # return the matrix of count
#   return(fc_count) # column "id"
#
#   # sub-functions
#   # get gtf (GoldCLIP project), parameters, or inbuilt
#   # get
# }
