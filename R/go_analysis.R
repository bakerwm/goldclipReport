

##----------------------------------------------------------------------------##
## GO analysis



#' go_group
#'
#' This function perform GO analysis using clusterProfiler package
#' It creates table and plots
#'
#' @description Perform GO analysis
#' @param gene_input list of genes
#' @param out_dir path to the results
#' @param orgdb the org.Hs.eg.db
#' @param show_num the number of categories to show in plot
#' @param level the level of GO
#' @param text_width the max width of labels in plot
#' @param overwrite whether overwrite exist file
#'
#'  @examples
#'  \donotrun{
#'    # Demo
#'  }
#'
#'  @import clusterProfiler
#'  @import stringr
#'  @import cowplot
#'
#'  @export
go_group <- function(gene_input, out_dir, orgdb,
                     show_num = 12, level = 3,
                     text_width = 30,
                     overwrite = FALSE) {
  ## group of GO
  library(clusterProfiler)
  ##
  if(! dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ## check exists file
  go_group_table <- file.path(out_dir, "go_group.xlsx")
  if(overwrite) {
    if(file.exists(go_group_table)){
      old_table <- file.path(out_dir, paste0(format(Sys.time( ), "%Y%m%d%H%M%S"), ".", "go_group.xlsx"))
      file.rename(go_group_table, old_table)
    }
  } else {
    if(file.exists(go_group_table)){
      stop("file exists, return go_group() with overwrite=TRUE")
    }
  }
  ## all-in-one
  plist <- lapply(c("BP", "CC", "MF"), function(g){
    print(paste0("GO group: ", g))
    ggo <- groupGO(gene     = gene_input,
                   keyType  = "ENSEMBL",
                   OrgDb    = orgdb,
                   ont      = g,
                   level    = level,
                   readable = TRUE)
    #head(ggo)

    ## save table
    if(nrow(ggo) > 0) {
      xlsx::write.xlsx(ggo, go_group_table, sheetName=g,
                       col.names=TRUE, row.names=TRUE, append=TRUE)
    }
    ## create barplot
    library(stringr)
    p1 <- barplot(ggo, drop = TRUE, showCategory = 12, order = TRUE) +
      scale_x_discrete(labels=function(x) str_wrap(x, width=text_width)) +
      ylab("Number of genes") +
      ggtitle(g) +
      theme(plot.title = element_text(hjust = 0.5))

    ## return plot
    return(p1)
  })

  ## save plot
  library(cowplot)
  go_group_plot  <- file.path(out_dir, "go_group.pdf")
  p <- cowplot::plot_grid(plotlist = plist, ncol = 1)
  cowplot::ggsave(go_group_plot, p, width = 7, height = 10)
}


#' go_enrich analysis
#'
#' This function perform GO analysis using clusterProfiler package
#' It creates table and plots
#'
#' @description Perform GO analysis
#' @param gene_hit a list of gene
#' @param gene_all a list of genes
#' @param out_dir path to the results
#' @param orgdb the org.Hs.eg.db
#' @param show_num the number of categories to show in plot
#' @param level the level of GO
#' @param text_width the max width of labels in plot
#' @param pvalueCutoff the cutoff to filt GO enrichment analysis
#' @param qvalueCutoff the cutoff to filt GO enrichment enrichment analysis
#' @param overwrite whether overwrite exist file#'
#'  @examples
#'  \donotrun{
#'    # Demo
#'  }
#'
#'  @import clusterProfiler
#'  @import stringr
#'  @import cowplot
#'
#'  @export
go_enrich <- function(gene_hit, gene_all, out_dir, orgdb,
                      show_num = 12, level = 3,
                      text_width = 30,
                      overwrite = FALSE,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05) {
  ## group of GO
  library(clusterProfiler)
  library(stringr)
  ##
  if(! dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ## check exists file
  go_enrich_table <- file.path(out_dir, "go_enrich.xlsx")
  if(overwrite) {
    if(file.exists(go_enrich_table)){
      old_table <- file.path(out_dir, paste0(format(Sys.time( ), "%Y%m%d%H%M%S"), ".", "go_enrich.xlsx"))
      file.rename(go_enrich_table, old_table)
    }
  } else {
    if(file.exists(go_enrich_table)){
      stop("file exists, return go_enrich() with overwrite=TRUE")
    }
  }

  ##
  plist <- lapply(c("BP", "CC", "MF"), function(g){
    print(paste0("GO enrich: ", g))
    ## enrich
    ego <- enrichGO(gene          = gene_hit,
                    universe      = gene_all,
                    keyType       = "ENSEMBL",
                    OrgDb         = orgdb,
                    ont           = g,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pvalueCutoff,
                    qvalueCutoff  = qvalueCutoff,
                    readable      = TRUE)
    # head(ego)

    ## save table
    if(nrow(ego) > 0) {
      xlsx::write.xlsx(ego, go_enrich_table, sheetName=g,
                       col.names=TRUE, row.names=TRUE, append=TRUE)
    }

    ## create barplot
    p1 <- dotplot(ego, showCategory = 12) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_width)) +
      xlab("Gene Ratio") +
      ggtitle(g) +
      theme(plot.title = element_text(hjust = 0.5))

    ## return
    return(p1)
  })

  ## save plot
  library(cowplot)
  go_enrich_plot  <- file.path(out_dir, "go_enrich.pdf")
  p <- cowplot::plot_grid(plotlist = plist, ncol = 1)
  cowplot::ggsave(go_enrich_plot, p, width = 7, height = 12)
}

# out_dir <- "go_group"
# go_group(gene_input, out_dir, orgdb,
#          show_num = 12,
#          level = 3,
#          text_width = 30,
#          overwrite = FALSE)

# out_dir <- "go_enrich"
# go_enrich(gene_hit, gene_all,
#           out_dir,
#           orgdb = org.Hs.eg.db,
#           show_num = 12,
#           level = 3,
#           text_width = 30,
#           overwrite = TRUE)
#
# setwd("~/work/wmlib/wmlib/go_analysis/")
# ## group go
# ##
# ## barplot
# data(geneList, package="DOSE")
# gene <- names(geneList)[abs(geneList) > 2]
# gene.df <- bitr(gene,
#                 fromType = "ENTREZID",
#                 toType   = c("ENSEMBL", "SYMBOL"),
#                 OrgDb    = org.Hs.eg.db)
#
# gene_all     <- names(geneList)
# gene_hit     <- gene
# show_num     <- 12
# text_width   <- 30
# level        <- 3
# pvalueCutoff <- 0.01
# qvalueCutoff <- 0.05
# organism     <- "human"
# orgdb <- org.Hs.eg.db


