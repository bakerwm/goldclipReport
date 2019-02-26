
#' demx_stat
#'
#' @param x path to the directory of demx
#'
#' @export
demx_stat <- function(x) {

}


library(rjson)
f <- "/nas/yulab/seq_data/Yu_2019/YY31_20190121/demo/p7_index/report_demx.json"
d <- rjson::fromJSON(file = f)
dn <- lapply(d, function(i) {as.data.frame(i, stringsAsFactors = FALSE)})
df <- dplyr::bind_rows(dn)
