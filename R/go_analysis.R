

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
                     keytype = "ENSEMBL",
                     show_num = 12, level = 3,
                     text_width = 30,
                     overwrite = FALSE,
                     readable = TRUE) {
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
                   keyType  = keytype,
                   OrgDb    = orgdb,
                   ont      = g,
                   level    = level,
                   readable = readable)
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
                      keytype = "ENSEMBL",
                      show_num = 12, level = 3,
                      text_width = 30,
                      overwrite = FALSE,
                      readable = TRUE,
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05) {
  ## group of GO
  # library(clusterProfiler)
  # library(stringr)
  ##
  if(! dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ## check exists file
  go_enrich_table <- file.path(out_dir, "go_enrich.xlsx")
  if(overwrite) {
    if(file.exists(go_enrich_table)){
      old_table <- file.path(out_dir,
                             paste0("go_enrich.",
                                    format(Sys.time( ), "%Y%m%d%H%M%S"), ".xlsx"))
      file.rename(go_enrich_table, old_table)
    }
  } else {
    if(file.exists(go_enrich_table)){
      stop("file exists, return go_enrich() with overwrite=TRUE")
    }
  }

  # save plot
  # library(cowplot)
  go_enrich_plot  <- file.path(out_dir, "go_enrich.pdf")
  pdf(go_enrich_plot, width = 10, height = 16)

  plist <- lapply(c("BP", "CC", "MF"), function(g){
    print(paste0("GO enrich: ", g))
    ## GO enrich
    ego <- enrichGO(gene          = gene_hit,
                    universe      = gene_all,
                    keyType       = keytype,
                    OrgDb         = orgdb,
                    ont           = g,
                    pAdjustMethod = "BH",
                    pvalueCutoff  = pvalueCutoff,
                    qvalueCutoff  = qvalueCutoff,
                    readable      = readable)
    # head(ego)
#
#     ## GO GSEA
#     ego3 <- gseGO(geneList     = gene_hit,
#                   keyType      = keytype,
#                   OrgDb        = orgdb,
#                   ont          = g,
#                   nPerm        = 1000,
#                   minGSSize    = 100,
#                   maxGSSize    = 500,
#                   pvalueCutoff = 0.05,
#                   verbose      = FALSE)
#
    ## save table
    pa <- NULL
    if(is.null(ego)) {
      print(paste0("No enrich genes detected: ", g))
    } else if(nrow(ego) > 0) {
      ## simplify GO terms
      ego2  <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

      ## write to table
      xlsx::write.xlsx(ego2, go_enrich_table, sheetName=g,
                       col.names=TRUE, row.names=TRUE, append=TRUE)

      ## create barplot
      p1 <- barplot(ego2, showCategory = 12) +
        scale_x_discrete(labels=function(x) str_wrap(x, width=text_width)) +
        ylab("Gene Ratio") +
        ggtitle(g) +
        theme(plot.title = element_text(hjust = 0.5))

      ## dotplot
      p2 <- dotplot(ego, showCategory = 12) +
        scale_x_discrete(labels=function(x) str_wrap(x, width=text_width)) +
        ylab("Gene Ratio") +
        ggtitle(g) +
        theme(plot.title = element_text(hjust = 0.5))

      ## create netplot
      p3 <- cnetplot(ego2)

      ## combine plot
      pa <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO",
                               align = "v", ncol = 1,
                               rel_heights = c(1, 1, 3))
    } else {
      print(paste0("No enrich genes detected: ", g))
    }

    ## return
    # return(pa)
    print(pa)
  })

  dev.off()

  # p <- cowplot::plot_grid(plotlist = plist, ncol = 1)
  # cowplot::ggsave(go_enrich_plot, p, width = 7, height = 12)
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


