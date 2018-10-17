# compute peak sensitivity
# 1. peaks overlap with seed-matches
# 2. filt seed-matches by miRNA expression
# 3. output: intervals, peaks, seed-matches numbers
# 4. peaks overlap with miRNAs


##----------------------------------------------------------------------------##
## function

#' isBed6
#'
#' @param data dataframe
#' @param type BED6 or BED12
#'
#' @import GenomicRanges
#' @import IRanges
#'
.isBed6 <- function(data, type = "bed6") {
  if("data.frame" %in% class(data) & length(data) >= 6) {
    if(all(c("chr", "start", "end", "name", "score", "strand") %in% names(data))) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}


#' GenerateInterval
#'
#' @param from integer, left edge of seq
#' @param to integer, right edge of seq
#' @param interval integer,
#'
#' @export
GenerateInterval <- function(from, to, interval) {
  vec <- do.call(what = seq, args = list(from, to, interval))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}


#' DataframeColVal
#'
#' @param data data.frame
#' @param col integer or character, column names in data.frame
#' @param return.values logical, if TRUE return the values of the column
#'
#' @dplyr
#'
#' @export
DataframeColVal <- function(data, col, return.values = FALSE) {
  if(is.data.frame(data) & length(col) == 1) {
    if(is.numeric(col) & length(data) >= col || is.character(col)) {
      if(is.numeric(col)) {
        col <- as.integer(col)
      }
      if(return.values) {
        return(dplyr::pull(data, col))
      } else {
        return(TRUE)
      }
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}


#' Bed2GRange
#'
#' @param data dataframe, in BED6 format
#'
#' @import GenomicRanges
#'
#' @export
Bed2GRange <- function(data) {
  # This function loads a BED-like file and stores it as a GRanges object.
  # The tab-delimited file must be ordered as 'chr', 'start', 'end', 'id', 'score', 'strand'.
  # The minimal BED file must have the 'chr', 'start', 'end' columns.
  # Any columns after the strand column are ignored.
  stopifnot(.isBed6(data)) # type
  data <- as.data.frame(data)
  data$strand <- gsub("[^+-]+", "*", data$strand) # ignore strand !!!!
  gr <- with(data, GenomicRanges::GRanges(seqnames = chr,
                                          ranges = IRanges::IRanges(start, end),
                                          score = score,
                                          strand = strand))
  return(gr)
}


#' DataframeTopN
#'
#' @param data data.frame data
#' @param sortByCol integer or character in data.frame column names
#' @param topN integer, top expressed records
#'
#' @import dplyr
#'
#' @export
DataframeTopN <- function(data, sortByCol, topN = TRUE) {
  stopifnot("data.frame" %in% class(data))
  if(! DataframeColVal(data, sortByCol)) {
    stop("sortByCol not included in data.frame")
  }
  n <- ifelse(is.numeric(sortByCol), names(data)[as.integer(sortByCol)],
              sortByCol)
  # if(! n %in% names(data)) {
  #   stop(paste0(n, ": not in column names"))
  # }
  if(! is.numeric(dplyr::pull(data, n))) {
    stop(paste0(n, ": column are not numeric data"))
  }
  #stopifnot(n %in% names(data))
  #stopifnot(is.numeric(dplyr::pull(data, n)))
  myCol <- paste0("desc(", n, ")") # descending
  df <- dplyr::arrange_(data, .dots = myCol)
  if(isTRUE(topN)){
    return(df2)
  } else if(is.numeric(topN) & topN > 0) {
    topN <- ifelse(nrow(data) < topN, nrow(data), as.integer(topN))
    return(df[seq(topN), ])
  } else {
    return(FALSE)
  }
}


#' MiRNAIdFormater
#'
#' @param x string of miRNA id, eg: gene:miR-1-3p/2-5p
#'
#' @import stringr
#'
#' @export
#'
# MiRNAIdFormater <- function(x) {
#   # convert gene:miR-1-3p/2-5p to miR-1-3p, miR-2-5p
#   GetId <- function(name) {
#     p <- unlist(strsplit(name, split = ":", fixed = TRUE))
#     g <- p[1]
#     n <- p[2]
#     n <- gsub("^miR-", "", n)
#     ns <- unlist(strsplit(n, "/", fixed = TRUE))
#     nids <- paste0("miR-", ns)
#     # trim "miR-124-3p.1" to "miR-124-3p"
#     nids <- gsub("\\.\\d+$", "", nids)
#     return(data.frame(gene = g, name = name, id = nids))
#   }
#   df <- do.call("rbind", lapply(x, GetId))
#   return(df)
# }
MiRNAIdFormater <- function(x) {
  # convert gene:miR-1-3p/2-5p to miR-1-3p, miR-2-5p
  x <- gsub("\\.\\d", "", x) # remove tail: 3p.1 to 3p
  x <- gsub("miR-", "", x) # remove miR-
  # split into: gene + id
  df1 <- tidyr::separate(data.frame(name = x), name, c("gene", "name"), sep = ":")
  # split ids:
  tabs <- paste0("X", seq(max(stringr::str_count(x, "/")) + 1))
  df2 <- tidyr::separate(df1, name, tabs, sep = "/", fill = "right")
  return(df2)
}


#' SeedmatchFilter
#'
#' @param data data.frame of seedmatch, BED6
#' @param miRs string, miRNA ids
#'
#' @import dplyr
#'
#' @export
#'
# SeedmatchFilter <- function(data, miRs) {
#   #df <- BedParser(file, "BED6", extra = TRUE) # seed-match in BED6
#   miRs <- gsub("^\\w+-miR", "miR", miRs) # trim "species" in miRs
#   df   <- MiRNAIdFormater(data$name) # temp data.frame, <gene> <name> <id>
#   df   <- dplyr::filter(df, id %in% miRs) # candidates
#   hit  <- dplyr::filter(data, name %in% df$name) # seedmatch filt
#   rm(df)
#   return(hit)
# }
SeedmatchFilter <- function(data, miRs) {
  miRs <- gsub("^\\w+-miR-", "miR-", miRs) # trim "species" in miRs
  df   <- MiRNAIdFormater(data$name) # temp data.frame, <gene> <name> <id>
  dfHits <- apply(df, 2, function(i) i %in% miRs)
  hits <- as.logical(rowSums(dfHits)) # hit rows
  return(data[hits, ])
}


#'GRintersect
#'
#' compare between two GRange objects
#'
#' @param subject GRange object
#' @param query GRange object
#' @param ignore.strand logical, if TRUE do not consider strand, default TRUE
#' @param gap integer, maximum distance between two intervals
#'
#' @import GenomicRanges
#' @import IRanges
#'
#' @export
BedIntersect <- function(subject, query, ignore.strand = TRUE, gap = 0) {
  #stopifnot("GRanges" %in% class(subject))
  #stopifnot("GRanges" %in% class(query))
  stopifnot("data.frame" %in% class(subject))
  stopifnot("data.frame" %in% class(query))
  stopifnot(is.logical(ignore.strand))
  stopifnot(is.numeric(gap))
  gap <- as.integer(gap) # convert to int
  grQuery <- Bed2GRange(query)
  grSubject <- Bed2GRange(subject)
  hit <- GenomicRanges::countOverlaps(grQuery, grSubject, maxgap = gap,
                                      ignore.strand = ignore.strand)
  return(length(hit[hit>0])) # count overlaps
}


#' GRintersectIntervals
#'
#' @param subject data.frame of subject in BED6 format
#' @param query data.frame of query in BED6 format
#' @param interval integer, number of records for each subgroup
#' @param ignore.strand logical, if TRUE do not consider strand, default TRUE
#' @param gap integer, maximum distance between two intervals
#'
#' @import GenomicRanges
#' @import IRanges
#'
#' @export
BedIntersectIntervals <- function(query, subject, interval = 10, gap = 0,
                                 ignore.strand = FALSE) {
  # count the number of subjects for specific group of querys
  # countOverlaps()
  stopifnot(is.numeric(interval))
  interval <- abs(as.integer(interval)) # int
  interval <- ifelse(interval > length(query), length(query), interval) # range
  queryIvl <- GenerateInterval(interval, length(query), interval)
  queryCnt <- lapply(seq(length(queryIvl)), function(i){
    if(i == 1) {
      s = 1
    } else {
      s <- queryIvl[i - 1] + 1
    }
#    s <- ifelse(s < 0, 1, s) # left-edge
    querySub <- query[s:queryIvl[i], ]
    c <- BedIntersect(subject, querySub, ignore.strand, gap)
    return(data.frame(interval = queryIvl[i], count = c))
  })
  df <- do.call("rbind", queryCnt)
  return(df)
}


##----------------------------------------------------------------------------##
## main
#' PeakSensitivity
#'
#' @param peak peak file in BED6 format
#' @param seedmatch seedmatch file in BED6 format
#' @param miRexp file of miRNA expression, BED6+extra format
#' @param interval integer, interval for the peaks
#'
#'
#' @import dplyr
#' @import GenomicRanges
#'
#' @export
#'
PeakSensitivity <- function(peak, seedmatch, miRexp, miRTopN = 20,
                            interval = 10, gap = 0, ignore.strand = FALSE) {
  dfPeak <- BedParser(peak)
  dfSm   <- BedParser(seedmatch)
  dfmiR  <- BedParser(miRexp)
  # filt topN miRNAs
  dfx   <- apply(dfmiR[, -c(1:6)], 2, as.numeric)
  dfmiR <- dplyr::mutate(dfmiR, exp = rowMeans(dfx))
  dfmiRTopN <- DataframeTopN(dfmiR, "exp", miRTopN)
  dfmiRTopN <- tidyr::separate(dfmiRTopN, name, into = c("name", "id"), sep = ",")
  miRidsTopN <- gsub("^\\w+-miR-|^miR-", "", dfmiRTopN$name)
  # filt seedmatch
  dfSmHits <- SeedmatchFilter(dfSm, miRidsTopN)
  # overlap
  df <- BedIntersectIntervals(peak, seedmatch, interval, gap, ignore.strand)
  return(df)
}


# peak <- "/home/data/projects/goldclip/goldclip_AGO2_20170114_merged/results/peak_bed/clipper/HEK293_NoAGO2_254_CLIP_merged.fixed.bed"
# seedmatch <- "../targetscan/Predicted_Targets.hg19.bed"
# miRexp <- "/home/data/projects/goldclip/goldclip_AGO2_20170114_merged/results/peak_sensitivity/miRNA_exp.txt"


# df <- PeakSensitivity(peak, seedmatch, miRexp, miRTopN = 200, interval = 5000, gap = 0, ignore.strand = FALSE)



















