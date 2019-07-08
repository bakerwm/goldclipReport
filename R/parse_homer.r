
## parse values from motif head line
#' motif1.motif
#' the matrix file of homer output
#'
#' example
#' T:2330.0(21.38%),B:442.9(8.47%),P:1e-374
#'
#' output
#' [1] "21.38%" "8.47%"  "1e-374"
#'
#' homerMotifDescVals(motif)
#'
#' @export
homerMotifDescVals <- function(file, nucl_type = 'RNA') {
  stopifnot(file.exists(file))
  ## parse header and content separately
  dfHeader <- read.table(file, sep = "\t", nrows = 1, stringsAsFactors = FALSE)
  dfScore  <- read.table(file, sep = "\t", skip = 1, stringsAsFactors = FALSE)
  ## parse desc vals
  dseq <- gsub(">", "", dfHeader$V1)
  #dlogPval <- format(dfHeader$V4, scientific = TRUE, digits = 3)
  dlogPval <- formatC(dfHeader$V4, format = "E", digits = 2) # switch to "E"
  dval <- unlist(strsplit(dfHeader$V6, ","))
  dval <- gsub("^[A-Z]\\:", "", dval, perl = TRUE)
  dval <- gsub("^[0-9.]+\\(|\\)", "", dval, perl = TRUE)
  dfMotif <- data.frame(seq = dseq,
                        pval = toupper(dval[3]),
                        logPval = dlogPval,
                        target = dval[1],
                        bg = dval[2],
                        stringsAsFactors = FALSE)
  rm(dfHeader, dseq, dlogPval, dval)
  ## format matrix, row=ACGT, column=position
  proportion <- function(x){
    rsum <- sum(x);
    return(x / rsum);
  }
  dfMatrix <- apply(dfScore, 1, proportion)
  nucl <- c('A', 'C', 'G', 'U')
  if(nucl_type == 'DNA') { nucl = c('A', 'C', 'G', 'T') }
  rownames(dfMatrix) <- nucl
  return(list(desc = dfMotif, matrix = dfMatrix))
}



#' @export
ggseqlogo2 <- function(data, facet = 'wrap', scales = 'free_y',
                       ncol = 1L, nrow = NULL,
                       theTitle = NULL, ...) {
  p <- ggplot() +
    geom_logo(data = data, method = 'prob', col_scheme = 'nucleotide2') +
    ggtitle(theTitle) +
    theme_logo() +
    theme(axis.title   = element_blank(),
          axis.text    = element_blank(),
          #strip.text   = element_blank(),
          strip.text.y = element_text(colour = 'black', face = 'bold',
                                      size = 10, angle = 180),
          plot.title   = element_text(hjust = .5))
  if (! "list" %in% class(data))
    return(p)
  facet_opts = c('grid', 'wrap')
  pind = pmatch(facet, facet_opts)
  facet = facet_opts[pind]
  if (is.na(facet))
    stop("facet option must be set to 'wrap' or 'grid'")
  if (facet == "grid") {
    p = p + facet_grid(~seq_group, scales = scales)
  }
  else if (facet == "wrap") {
    p = p + facet_wrap(~seq_group, strip.position = 'l',
                       ncol = ncol, nrow = nrow, scales = scales)
  }
  return(p)
}




#' @export
homerMotifDescPlot <- function(motifs, motif.only = FALSE, title = NULL, ...) {
  if(! "character" %in% class(motifs))
    #return(0)
    stop("motifs: is not character")
  if(length(motifs) < 1)
    #return(0)
    stop("motifs: characters not found")
  ## parsing motif files
  da <- lapply(motifs, function(x){
    d <- homerMotifDescVals(x)
  })
  ## matrix
  dm <- lapply(1:length(da), function(x){
    ma <- da[[x]][['matrix']]
  })
  ## description
  ddes <- lapply(1:length(da), function(x){
    de <- da[[x]][['desc']]
  })
  ## format description
  ddf <- reshape::merge_all(ddes) %>%
    select(- seq) %>%
    tibble::rownames_to_column(var = 'rank') %>%
    arrange(rank) %>%
    mutate(x = 1, y = length(bg):1) %>%
    gather(key = 'type', value = 'value', 2:5) %>%
    mutate(type = factor(type,
                         levels = c('pval', 'logPval', 'target', 'bg'),
                         labels = c('P-value', 'log P-value', 'Targets%', 'Background%')))
  ## plots
  p1 <- ggseqlogo2(dm, ncol = 1) + theme(plot.margin = margin(10, 2, 5, 2))
  p2 <- ggplot(ddf, aes(x = x, y = y, label = value)) +
    geom_blank() +
    geom_text(size = 3) +
    facet_wrap(~type, nrow = 1) +
    #theme_nothing()
    theme_minimal() +
    theme(text         = element_text(color = "black", size = 8),
          axis.title   = element_blank(),
          axis.text    = element_blank(),
          axis.line    = element_blank(),
          panel.grid   = element_blank(),
          panel.border = element_blank(),
          plot.title   = element_text(hjust = .5, face = 'bold'),
          plot.margin  = margin(0, 5, 15, 5))
  ## combine plots
  suppressPackageStartupMessages(library(cowplot))
  mainTitle <- ggdraw() + draw_label(title, fontface = 'bold')
  if(motif.only)
    if(is.null(title)) {
      return(p1)
    } else {
      p1x <- plot_grid(mainTitle, p1, ncol = 1, rel_heights = c(0.1, 1)) # rel_heights values control title margins
      return(p1t)
    }
  p <- plot_grid(p1, p2, align = 'h', axis = 't', nrow = 1)
  if(is.null(title)) {
    return(p)
  } else {
    px <- plot_grid(mainTitle, p, ncol = 1, rel_heights = c(0.1, 1))
    return(px)
  }
}


#' homerid2table
#'
#' @param x, the input data.frame
#'
#' @export
homerid2table <- function(x) {
  stopifnot(any(class(x) == "data.frame"))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  names(x) <- c('seq', 'target', 'n1', 'logPval', 'n2', 'n3')
  df <- x %>%
    mutate(target = gsub("^.*BestGuess:", "", target, perl = TRUE),
           logPval = format(logPval, scientific = TRUE, digits = 2)) %>%
    separate(col = n3, into = c('tp', 'bp', 'pval'), sep = ",") %>%
    select(pval, logPval, tp, bp) %>%
    mutate(pval = gsub("\\w\\:", "", pval, perl = TRUE),
           tp = gsub("\\w\\:", "", tp, perl = TRUE),
           bp = gsub("\\w\\:", "", bp, perl = TRUE)) %>%
    mutate(tp = gsub("^\\d+\\.?\\d|\\(|\\)", "", tp),
           bp = gsub("^\\d+\\.?\\d|\\(|\\)", "", bp))
  return(df)
}

#' fig2url
#'
#' convert figure to url in markdown
#'
#' @export
fig2url <- function(path) {
  f <- paste0("![](", path, ")")
  return(f)
}


#' parseHomer
#'
#' @param path path to the HOMER (findMotifs.pl) output
#' @param topN top N motifs
#'
#' @export
parseHomer <- function(path, topN = 5) {
  suppressPackageStartupMessages(library(dplyr))
  library(readr)
  library(digest)
  path <- normalizePath(path)
  stopifnot(dir.exists(path))
  subpath <- file.path(path, 'homerResults')
  if(! dir.exists(subpath)) return(0)
  # motif files, file2table
  m <- paste0('motif', 1:topN, ".motif")
  mfile <- file.path(subpath, m)
  dl <- lapply(mfile, function(x){
    if(! file.exists(x)) return(NA) #next
    id <- read.table(x, sep = "\t", nrows = 1)
    #id <- read_delim(x, delim = "\t", n_max = 1, col_names = FALSE)
    return(id)
  })
  df1 <- do.call('rbind', dl)
  ## weblogo (svg2png)
  pg <- mapply(function(x){
    logo  <- paste0('motif', x, ".logo.svg")
    mlogo <- file.path(subpath, logo)
    fig   <- gsub('.svg$', '.png', logo, perl = TRUE)
    ## create dir
    n1 <- digest(path, algo = "md5", serialize = FALSE)
    mdir  <- file.path("./", "images", n1)
    dir.create(mdir, recursive = TRUE, showWarnings = FALSE, mode = "0755")
    mfig  <- file.path(mdir, fig)
    if(! file.exists(mlogo)) return(NA)
    rsvg::rsvg_png(mlogo, mfig)
    return(mfig)
  }, 1:topN)
  ## output table
  df2 <- homerid2table(df1) %>%
    mutate(Rank = 1:topN,
           Motif = fig2url(pg)) %>%
    select(Rank, Motif, pval, logPval, tp, bp)
  colnames(df2) <- c('Rank', 'Motif', 'P-value', 'log P-value', '% of Target', '% of Background')
  df2 <- na.omit(df2)
  return(df2)
}


## render HOMER results
#' @export
homerDirToHtml <- function(dirlist,
                           outputRmd = "demo.Rmd",
                           RBP    = "PTB",
                           topN   = "All",
                           render = FALSE,
                           template = "~/work/yu_2017/bin/homer_report.template/homer_report.template.v2.Rmd") {
  stopifnot(file.exists(template)) # template
  stopifnot(file.exists(dirlist)) # file contain HOMER dir
  dirlist <- normalizePath(dirlist) # absolute path
  ## checkout dir existence
  dirs <- readr::read_delim(dirlist,
                            delim = " ",
                            col_names = 'id',
                            col_types = cols()) %>% pull(id)
  dirs <- dirs[dir.exists(dirs)]
  stopifnot(length(dirs) > 0) ##
  ## change temple
  lines     <- readLines(template)
  lines[2]  <- paste0("title: \"", RBP, ": HOMER analysis report top ", topN, "\"")
  lines[14] <- paste0("homerDirList <- \"", dirlist, "\" ")
  ## save Rmd file
  outdir <- dirname(outputRmd)
  if(! dir.exists(outdir)) dir.create(outdir, mode = "0755")
  status <- writeLines(lines, outputRmd)
  if(render) rmarkdown::render(outputRmd, quiet = TRUE) ## render Rmd to html file
}


## function to create summary table
#' @export
homerListToTable <- function(
  inputDirs,
  project = "HEK293_PTB",
  RBP     = 'PTB',
  tool    = 'pyicoclip',
  topN    = 'All',
  subdir  = 'data',
  render  = FALSE) {
  stopifnot(length(inputDirs) > 1) ##
  ## save sub-html files to data/
  subRmd <- file.path(subdir, paste("HOMER", project, topN, "Rmd", sep = "."))
  subList <- file.path(subdir, paste("HOMER", project, topN, "txt", sep = "."))
  ## parse and save HOMER dirs
  fileConn <- file(subList)
  writeLines(inputDirs, fileConn)
  close(fileConn)
  ## create Rmd file
  homerDirToHtml(dirlist   = subList,
                 outputRmd = subRmd,
                 topN      = topN,
                 render    = render,
                 RBP       = RBP)
  ## save data.frame
  dn <- data.frame(RBP = RBP, topN = topN, rmd = subRmd, peakCaller = tool) %>%
    mutate(html = gsub(".Rmd$", ".html", rmd)) %>%
    mutate(link = url2mdlink(html, title = rep('link', times = length(html)))) %>%
    select(RBP, topN, peakCaller, link)
  return(dn)
}



