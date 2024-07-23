
# ------- plot count num. of editing types ---------------
# barplot
#' Title
#'
#' @param df rna_type_count
#' @param save whether to save
#' @param file FilePath to save
#'
#' @return barplot of rna type count
#' @export
#'
#' @examples
#' plot_edit_type_count(rna_type_count)
#'
#' plot_edit_type_count(rna_type_count, save=TRUE, file="plot_rna_type_count.pdf")

plot_edit_type_count = function(df, save=FALSE, file=NA) {

  # TODO
  # format check of df (ie. rna_type_count)


  rna_type_count = df
  rna_type_count_long = rna_type_count %>%
    tidyr::pivot_longer(c(2:ncol(rna_type_count)), names_to = "sampleID", values_to = "counts")

  y_list=c("counts")

  for(item in y_list) {
    print(paste0(item, " is plotting ..."))
    library(ggpubr)
    p1 = ggbarplot(rna_type_count_long %>% na.omit, x = "type", y = "counts",
                   fill = "skyblue",
                   add = "mean_se", error.plot = "upper_errorbar") +
      scale_y_continuous(n.breaks = 10)

    # library(gg.gap)
    # gg.gap(plot = p1,
    #        segments = c(200, 10000),
    #        tick_width = c(40,20000),
    #        rel_heights = c(0.1, 0, 0.3),
    #        ylim = c(0, 130000)
    # )

    if(save) {
      print(paste0(item, " plot saves to ..."))
      ggsave(filename=file, width = 6, height = 5)
    }
  }
  return(p1)
}


