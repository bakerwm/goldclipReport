
#' merge a list of data.frames
#'
#' @param list a list of data.frames
#'
#' @param by column name for function merge()
#'
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
    list <- rlist::list.remove(list, c(1, 2))
    # list_new <- rlist::list.append(list, df)
    list_new <- rlist::list.insert(list, 1, df)
    # list[[1]] <- NULL
    # list[[2]] <- NULL
    # list
    # list_new <- c(list, df)
    # print(length(list_new))
    return(merge_list_of_dataframes(list_new, by = by))
  }
}





#' merge a list of data.frames
#'
#' @param x a list of data.frames
#' @param by column name for function merge()
#' @param all logical
#'
#' @import rlist
#'
#' @export
merge_df <- function(x, by, all = TRUE) {
  stopifnot(inherits(x, "list"))

  # check all data.frame
  if(!all(unlist(lapply(x, is.data.frame)))){
    stop("only data.frames allowed in x, exception detected")
  }

  # check id (by=)
  n <- lapply(x, function(i) {
    by %in% names(i)
  })
  if(!all(unlist(n))) {
    stop(paste0("by=", by, " not present in all data.frame"))
  }

  # merge
  if(length(x) == 1) {
    return(x[[1]])
    #  } else if(length(x) == 2) {
    #    return(merge(x[[1]], x[[2]], by = by, all = all))
  } else {
    df <- merge(x[[1]], x[[2]], by = by, all = all)
    # convert NA to 0
    x  <- rlist::list.remove(x, c(1, 2))
    x2 <- rlist::list.insert(x, 1, df)
    return(merge_df(x2, by = by, all = all))
  }
}




#' return the longest common string
#'
#' @param list A list of strings
#'
#' @param suffix A logical value, if TRUE, return the longest suffix of strings
#'   default: FALSE
#'
#' @description return the longest common string between list of strings
#' default: prefix
#'
#' @import stringi
#'
#' @export
#'
str_common <- function(list, suffix = FALSE) {
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


#' convert json to data.frame
#' only fit for report_json
#'
#' @param x file or URL, report_demx.json file
#'
#' @import jsonlite
#'
#' @return
#'
#' @export
json2dataframe <- function(x){
  xlist <- jsonlite::fromJSON(txt = x, flatten = TRUE)

  tmp <- lapply(names(xlist), function(i){
    data.frame(xlist[[i]], stringsAsFactors = FALSE) %>%
      dplyr::mutate(index = i) %>%
      dplyr::select(index, name, count)
  })

  df <- dplyr::bind_rows(tmp) %>%
    dplyr::arrange(name)
}



#' bacteria report for each lane
#'
#'
#' @param x path or URL, directory of bacteria dir
#'
#'
#' @import RColorBrewer
#' @import pheatmap
#'
#'
#' @export
#'
demxBacteria <- function(x, save2pdf = FALSE,
                         width = 8, height = 8) {
  # pdf file
  pdfFile <- file.path(x, "bacteria_content.pdf")
  csvFile <- file.path(x, "bacteria_content.csv")

  # read file
  stat_files <- list.files(x, "kraken2.stat$",
                           all.files = TRUE,
                           full.names = TRUE)
  dl <- lapply(stat_files, function(i){
    f_name <- gsub(".kraken2.stat", "", basename(i))
    readr::read_delim(i, "\t", col_types = readr::cols()) %>%
      mutate(reads_in_tax = reads_in_tax * 1e6 / reads_total) %>%
      dplyr::select(name, reads_in_tax) %>%
      dplyr::rename(!! f_name := reads_in_tax)
  })

  df <- goldclipReport::merge_list_of_dataframes(dl, by = "name")
  # df <- df %>% dplyr::select(sort(colnames(df)))
  df <- tibble::column_to_rownames(df, "name")
  ma <- log10(df + 1)

  ## make plot
  cc = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7,
                                       name = "RdYlBu")))
  breaks = seq(0, 6, length.out = 100)

  # plot
  p <- pheatmap::pheatmap(ma,
                silent = TRUE,
                cluster_cols = FALSE,
                color = cc(100),
                breaks = breaks,
                border_color = "grey40")

  # return
  if(isTRUE(save2pdf)) {
    pdf(pdfFile, width = width, height = height)
    print(p)
    dev.off()
    ## save to table
    write.csv(df, csvFile, row.names = TRUE, quote = FALSE)
  } else {
    return(p)
  }
}



#' parse reads from demx dir
#'
#' @param x path or URL,
#'
#' @import dplyr
#' @import scales
#'
#' @export
#'
#' the following are the structure of demx directory
#'
#' results
#'   ├── bacteria
#'   ├── barcode_ATCACG
#'   ├── p7_index*
#'   └── qc
#'
#' p7_index: required, saved main data
#' barcode_ATACACTG: optional, 0-N directories
#' qc: fastqc output of all fastq files
#' bacteria: bacteria content
#'
demxParser <- function(x) {
  # the YY number
  # YY35_20190325/results/
  yy <- basename(dirname(x))
  yy <- unlist(strsplit(yy, "_"))[1]
  # search directories
  rdirs <- list.files(x, "^p7_index$|^barcode_[ACGT]+$|^bacteria$",
                      all.files    = TRUE,
                      full.names   = TRUE,
                      include.dirs = TRUE)
  p7_dir  <- rdirs[basename(rdirs) == "p7_index"]
  bc_dirs <- rdirs[grep("^barcode_[ACGT]+$", basename(rdirs))]
  bac_dir <- rdirs[basename(rdirs) == "bacteria"]

  # check p7_dir
  p7_json <- file.path(p7_dir, "report_demx.json")
  p7_df   <- json2dataframe(p7_json) %>%
    dplyr::rename(p7_index = index) %>%
    dplyr::mutate(barcode  = "NULL",
                  yy       = yy) %>%
    dplyr::select(yy, p7_index, barcode, name, count)

  # check barcode dirs
  bc_list <- lapply(bc_dirs, function(i){
    i_json <- file.path(i, "report_demx.json")
    i_df   <- json2dataframe(i_json) %>%
      dplyr::rename(barcode = index) %>%
      dplyr::mutate(p7_index = gsub("barcode_", "", basename(i)),
                    yy = yy) %>%
      dplyr::select(yy, p7_index, barcode, name, count)
  })

  # filt undemx
  bc_df <- dplyr::bind_rows(bc_list)

  # combine p7 and barcode
  p7_df <- p7_df %>%
    dplyr::filter(! grepl("^[ACGTN]{4,10}$", name))

  demx_df <- rbind(p7_df, bc_df) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(count = sum(count)) %>%
    dplyr::rename(id = name) %>%
    dplyr::mutate(yy = yy) %>%
    dplyr::select(yy, id, count)

  # save bacteria report
  tmp <- demxBacteria(bac_dir, save2pdf = TRUE)

  return(demx_df)
}




#' create html report for demx
#'
#'
#' @param x path or URL, reults directory of demx
#'
#'
#' @export
#'
demxReport <- function(x) {
  # template
  report_template <- system.file("report_templates",
                                 "demx_report_template.Rmd",
                                 package = "goldclipReport")
  output_html <- file.path(x, "demx_report.html")
  output_tsv  <- file.path(x, "demx_report.txt")

  # save to table
  df <- demxParser(x)
  write.table(df,
              file  = output_tsv,
              quote = FALSE,
              sep   = "\t",
              row.names = FALSE)

  # render html
  rmarkdown::render(input = report_template,
                    output_file = output_html,
                    params = list(demx_path = x))
}


