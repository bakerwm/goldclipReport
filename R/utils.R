
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
    list_new <- rlist::list.append(list, df)
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



