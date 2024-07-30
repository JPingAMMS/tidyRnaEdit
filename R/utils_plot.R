
# ------- plot count num. of editing types ---------------
# barplot
#' Title
#'
#' @param data dataframe, rna_type_count
#' @param save whether to save
#' @param file FilePath to save
#'
#' @return barplot of rna editing type count
#' @export
#'
#' @examples
#' plot_edit_type_count(rna_type_count)
#'
#' plot_edit_type_count(rna_type_count, save=TRUE, file="plot_rna_type_count.pdf")

plot_edit_type_count = function(data=rna_type_count, save=FALSE, file=NA) {

  # TODO
  # format check of data (ie. rna_type_count)


  # rna_type_count = data
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


#' Title
#'
#' @param data
#' @param save
#' @param file
#'
#' @return
#' @export
#'
#' @examples
plot_edit_rate_mean = function(data=mat_mean, save=FALSE, file=NA){
  # TODO
  # format check of data (ie. mat_mean)

  # test
  # mat_mean=read.csv(gzfile(paste0(path_save,"/mat_mean.csv.gz")),header=TRUE)
  # plot_mat_mean(data=mat_mean,file=paste0(path_figure,"/mean_edit_group.pdf"))

  data_long <- data %>%
    pivot_longer(cols = -symbol, names_to = "group", values_to = "mean_edit_rate")

  #print(head(data_long))
  pairs=combn(unique(mat_mean_long$group),2,simplify=FALSE)
  p <- ggboxplot(data_long, x="group", y="mean_edit_rate", color = "group",
                 palette = c("#00AFBB", "#E7B800", "#FC4E07"))+
    #add = "jitter", shape="group")+
    stat_compare_means(comparisons = pairs,method = 't.test')

  if(save) {
    ggplot2::ggsave(p,filename=file, width = 5, height = 5)
  }

  return(p)
}



#' Title
#'
#' @param data
#' @param FDR_thres
#' @param delta_thres
#' @param show_up_down
#' @param save
#' @param file
#'
#' @return
#' @export
#'
#' @examples
plot_volcano = function(data,
                        FDR_thres = 0.05,
                        delta_thres = 0,
                        show_up_down = FALSE,
                        save=FALSE, file=NA){
  # TODO
  # 1. format check of data (ie. mat_mean); check FDR_thres/delta_thres...
  # 2. if exists `FDR` `delta`, then skip `rename cols`
  # 3. if show_up_down, then up=red, down=green. if not show_up_down, sig=green.

  #rename cols
  data_rename=data
  names(data_rename) <- sub(".*\\.", "", names(data_rename))
  data_rename=data_rename[,-c(1,2)]
  print(head(data_rename))

  # mark sig sites
  data_rename$sig=as.factor(ifelse(data_rename$FDR < FDR_thres,'with_sig','nonsense'))
  print(paste0('number of sites significant:',sum(data_rename$FDR < FDR_thres)))
  print(paste0('number of sites nonsense:   ',sum(data_rename$FDR >= FDR_thres)))

  # calculate log
  data_rename$log10FDR=-log10(data_rename$FDR)

  # volcano & save
  print('plot volcano ...')
  p1 <- ggplot(data_rename,
               aes(x =delta, y=lgFDR,colour=sig)
               ) +
    geom_point(alpha=0.65, size=0.1) +
    scale_color_manual(values=c("#d2dae2", "green"))  +
    geom_vline(xintercept=c(-delta_thres, delta_thres),
               lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(FDR_thres),
               lty=4,col="black",lwd=0.8) +
    labs(title = "Volcano plot",
         x="delta editing rate",
         y="-log10(FDR)") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="right",
          legend.title = element_blank()
          )

  if(save) {
    print('saving vocanol...')
    ggplot2::ggsave(p1,filename=file, width = 7, height = 5)
  }

  return(p1)
}



