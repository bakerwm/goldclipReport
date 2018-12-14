

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


