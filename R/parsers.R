
#' BedParser
#'
#' read BED file
#'
#' @param file BED file
#' @param type BED3, BED6 or BED12
#' @param extra logical, if TRUE, including extra columns in output
#' @param score.desc logical, if TRUE re-arrange BED records decreasingly
#'
#' @import dplyr
#' @import readr
#'
#' @export
BedParser <- function(file, type = "BED6", extra = TRUE, score.desc = FALSE) {
  # read BED file, default: BED6+N
  df <- readr::read_delim(file, delim = "\t", col_names = FALSE,
                          col_types = readr::cols(.default = "c"))
  if(type == "BED3") {
    if(length(df) < 3) {
      stop("BED3, less than 3 columns detected")
    } else {
      names(df)[1:3] <- c("chr", "start", "end")
    }
  } else if(type == "BED6") {
    if(length(df) < 6) {
      stop("BED6, less than 6 columns detected")
    } else {
      names(df)[1:6] <- c("chr", "start", "end", "name", "score", "strand")
    }
  } else if(type == "BED12") {
    if(length(df) < 12) {
      stop("BED12, less than 12 columns detected")
    }
  } else {
    stop(paste0(type, ": unknown type, BED3, BED6, BED12"))
  }
  # BED3 BED6 BED12 All
  df <- dplyr::mutate(df,
                      chr = as.character(chr),
                      start = as.numeric(start),
                      end = as.numeric(end))
  # sort
  if (type %in% c("BED6", "BED12") & score.desc) {
    df <- dplyr::arrange(df, desc(score))
  }

  # select columns
  if (! extra) {
    df <- dplyr::select(df, - dplyr::starts_with("X"))
  }

  return(df)
}


#' BedLenParser
#'
#' @param file path to the BED file, BED6
#' @param prefix name of the BED file
#'
#' @import dplyr
#'
#' @export
#'
BedLenParser <- function(file, prefix = NULL) {
  # parsing the number, length of BED records
  if( ! is.null(prefix)) {
    stopifnot(length(file) == length(prefix))
  }
  # sub function
  f2d <- function(file, id) {
    bed <- BedParser(file)
    df <- as.data.frame(table(bed$end - bed$start)) #Var1 Freq
    names(df) <- c("length", "count")
    if(is.null(id)) {
      df <- dplyr::mutate(df, id = gsub(".bed$", "", basename(file)))
    } else {
      df <- dplyr::mutate(df, id = id)
    }
    return(df)
  }
  #
  dl <- lapply(seq(length(file)), function(i){
    f2d(file[i], prefix[i])
  })
  df <- do.call("rbind", dl)
}






#' AnnoParser
#'
#' @param file file of the annotation
#' @param name name of the annotation
#'
#' @import dplyr
#'
#' @export
AnnoParser <- function(file, prefix = NULL, dm3_basic = FALSE) {
  # parsing the annotation output
  if(! is.null(prefix)) {
    stopifnot(length(file) == length(prefix))
  }
  # sub function
  f2d <- function(file, id) {
    df <- read.table(file, header = FALSE, sep = " ",
                     col.names = c('id', 'genome', 'type', 'count'))
    # multiple samples in one file
    nsmps <- length(unique(df$id))
    if(nsmps > 1) {
      df <- df
    }else if(is.null(id)) {
      df <- dplyr::mutate(df, id = file)
    } else {
      df <- dplyr::mutate(df, id = id)
    }
    df <- dplyr::select(df, id, type, count) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(pct = round(count / sum(count) * 100, 2))
  }
  dl <- lapply(seq(length(file)), function(i){
    f2d(file[i], prefix[i]) %>%
      ungroup() %>%
      mutate(id = as.character(id))
  })
  tags <- c('tts', 'rRNA', 'pseudo', 'promoters', 'ncRNA', 'introns',
            'intergenic', 'coding', 'utr5', 'utr3')
  labels <- c("TTS", "rRNA", "Pseudo", "Promoters", "ncRNA", "Introns",
              "Intergenic", "CDS", "5' UTR", "3' UTR")
  if(isTRUE(dm3_basic)) {
    tags <- c('TE', 'genicPiRNA', 'nonGenicPiRNA', 'rRNA', 'tRNA', '3u', '5u',
              'exon', 'intron', 'igr')
    labels <- c('TE', 'piRNA', 'nongenic-piRNA', 'rRNA', 'tRNA', "3' UTR", "5' UTR",
                "Exon", "Intron", "Intergenic")
  }
  df1 <- do.call("rbind", dl) %>%
    dplyr::filter(type != "others") %>%
    dplyr::mutate(type = plyr::mapvalues(type, tags, labels)) %>%
    dplyr::mutate(type = factor(type, levels = labels))
  return(df1)
}


#' AnnoParser2
#' parse the output of bed_annotation.py
#' <type> <count> <sample>
#'
#' @param x path,file of the annotation
#' @param prefix name of the annotation
#' @param dm3_basic logical, use dm3 types
#'
#' @import dplyr
#'
#' @export
AnnoParser2 <- function(x, prefix = NULL, dm3_basic = FALSE) {
  # sub function
  f2d <- function(x) {
    df <- read.table(x, header = TRUE, col.names = c("count", "id"))
    df <- tibble::rownames_to_column(df, "type") %>%
      dplyr::select(id, type, count) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(pct = round(count / sum(count) * 100, 2)) %>%
      dplyr::ungroup() %>%
      as.data.frame()
  }
  # parse multiple x files
  da <- lapply(x, f2d)

  # parsing the annotation output
  if(! is.null(prefix)) {
    if (length(da) != length(prefix)) {
      stop("x and prefix are differ in length")
    }
  } else {
    prefix = ""
  }

  # add names
  tags <- c('tts', 'rRNA', 'pseudo', 'promoters', 'ncRNA', 'introns',
            'intergenic', 'coding', 'utr5', 'utr3')
  labels <- c("TTS", "rRNA", "Pseudo", "Promoters", "ncRNA", "Introns",
              "Intergenic", "CDS", "5' UTR", "3' UTR")
  if(isTRUE(dm3_basic)) {
    tags <- c('TE', 'genicPiRNA', 'nonGenicPiRNA', 'rRNA', 'tRNA', '3u', '5u',
              'exon', 'intron', 'igr')
    labels <- c('TE', 'piRNA', 'nongenic-piRNA', 'rRNA', 'tRNA', "3' UTR", "5' UTR",
                "Exon", "Intron", "Intergenic")
  }
  df1 <- do.call("rbind", da) %>%
    dplyr::filter(! type %in% c("other", "others")) %>%
    dplyr::mutate(type = plyr::mapvalues(type, from = tags, to = labels,
                                         warn_missing = FALSE)) %>%
    dplyr::mutate(type = factor(type, levels = labels))
  return(df1)
}



#' AnnoParser3
#' parse the output of bed_annotation.py
#' <type> <count> <sample>
#'
#' @param x path, file of the annotation
#' @param prefix name of the annotation
#' @param dm3_basic logical, use dm3 types
#'
#' @import dplyr
#'
#' @export
AnnoParser3 <- function(x, prefix = NULL, dm3_basic = FALSE) {
  # sub function
  f2d <- function(f) {
    df <- readr::read_delim(f, "\t", col_names = TRUE, col_types = readr::cols())
    # df <- tibble::rownames_to_column(df, "type") %>%
    #   dplyr::select(id, type, count) %>%
    #   dplyr::group_by(id) %>%
    #   dplyr::mutate(pct = round(count / sum(count) * 100, 2)) %>%
    #   dplyr::ungroup() %>%
    #   as.data.frame()
    return(df)
  }
  # parse multiple x files
  da <- lapply(x, f2d)

  # parsing the annotation output
  if(! is.null(prefix)) {
    if (length(da) != length(prefix)) {
      stop("x and prefix are differ in length")
    }
  } else {
    prefix = ""
  }

  # add names
  tags <- c('tts', 'rRNA', 'pseudo', 'promoters', 'ncRNA', 'introns',
            'intergenic', 'coding', 'utr5', 'utr3')
  labels <- c("TTS", "rRNA", "Pseudo", "Promoters", "ncRNA", "Introns",
              "Intergenic", "CDS", "5' UTR", "3' UTR")
  if(isTRUE(dm3_basic)) {
    tags <- c('TE', 'genicPiRNA', 'nonGenicPiRNA', 'rRNA', 'tRNA', '3u', '5u',
              'exon', 'intron', 'igr')
    labels <- c('TE', 'piRNA', 'nongenic-piRNA', 'rRNA', 'tRNA', "3' UTR", "5' UTR",
                "Exon", "Intron", "Intergenic")
  }
  df1 <- do.call("rbind", da) %>%
    dplyr::filter(type %in% tags) %>%
    dplyr::mutate(type = plyr::mapvalues(type, from = tags, to = labels,
                                         warn_missing = FALSE)) %>%
    dplyr::mutate(type = factor(type, levels = labels))
    # do.call("rbind", da) %>%
    # dplyr::filter(! type %in% c("other", "others")) %>%
    # dplyr::mutate(type = plyr::mapvalues(type, from = tags, to = labels,
    #                                      warn_missing = FALSE)) %>%
    # dplyr::mutate(type = factor(type, levels = labels))
  return(df1)
}
















