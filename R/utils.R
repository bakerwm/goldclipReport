
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
    list <- list.remove(list, c(1, 2))
    list_new <- list.append(list, df)
    # list[[1]] <- NULL
    # list[[2]] <- NULL
    # list
    # list_new <- c(list, df)
    # print(length(list_new))
    return(merge_list_of_dataframes(list_new, by = by))
  }
}
