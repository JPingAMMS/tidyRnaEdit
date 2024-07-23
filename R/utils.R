#######################
##     Read data     ##
#######################
# ------ read-in data ------

#' read_redit
#'
#' @param sample_info df of 3 cols(filename, sampleid, group)
#' @param soft c("reditools_known") currently supporting reditools results.
#'
#' @return matrix of RNA editing sites
#' @importFrom readr read_tsv
#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr full_join
#' @importFrom tidyr separate
#' @export
#'
#' @examples
#' read_redit(sample_info, soft="reditools_known")

read_redit = function(sample_info, soft="reditools_known") {

  # TODO:
  # format check on sample_info

  # soft = "reditools_known"
  if (soft == "reditools_known") {
    importer <- function(x) suppressMessages(readr::read_tsv(x, progress=FALSE))
  }

  # Editing data Read-In & Merge
  merged_mat = NULL
  for (i in seq(nrow(sample_info))) {
    # i=1
    id=sample_info[i,"sampleid"]

    time=format(Sys.time(), "%Y-%m-%d__%H:%M:%S")
    message(time, "   Read-in Sample ", i,"/", nrow(sample_info), " --> ",id, " \n",appendLF=FALSE)

    # import
    tmp = as.data.frame(importer(as.character(sample_info[i,"filename"])))

    # # check for columns
    stopifnot(all(c("Region", "Position", "AllSubs", "Frequency") %in% names(tmp)))

    tmp = tmp %>% dplyr::filter(Frequency >= 0.001) %>%
      dplyr::mutate(symbol = paste0(Region,"_",Position,":",AllSubs)) %>%
      dplyr::select(c(symbol, Frequency))
    colnames(tmp) = c("symbol", id)

    # out matrix
    if(i == 1) {
      merged_mat = tmp
    } else {
      merged_mat = merged_mat %>% dplyr::full_join(tmp, by = "symbol")
    }
  }

  # look at results
  # head(merged_mat)[,1:5]

  merged_data = merged_mat %>% tidyr::separate(col="symbol",into=c("chr","pos","type"), sep="[_:]")
  # look at results
  head(merged_data)[,1:5]

  return(merged_data)
}

# --------- save csv && compressed as gzfile --------------
save_csv_gz = function(df, file){
  write.csv(df, file=gzfile(file), row.names = FALSE)
}


# ------- count num. of editing types ---------------
edit_type_count = function(df) {

  # TODO:
  # format check on df (i.e. merged_data)

  merged_data = df

  # lib
  library(dplyr)
  # init
  rna_type_count = NULL
  for(i in 4:ncol(merged_data)){
    # i=3
    id=colnames(merged_data)[i]
    tmp = merged_data[,c(3,i)] %>% na.omit()
    tmp_count = tmp %>%
      group_by(type) %>%
      summarise(cnt = n())
    colnames(tmp_count) = c("type", id)

    # out matrix
    if(i == 4) {
      rna_type_count = tmp_count
    } else {
      rna_type_count = rna_type_count %>% dplyr::full_join(tmp_count, by = "type")
    }
  }

  # look at results
  head(rna_type_count)

  return(rna_type_count)
}

#######################
##        QC         ##
#######################
# ------- Filter variant ------
filter_variant = function(merged_data,
                          sample_info,
                          edit_rate_indi = 0.05,
                          edit_rate_mean = 0.05,
                          missing_rate = NA,
                          left_n_total = NA,
                          left_n_in_each_group = 3
){

  # TODO
  # format check of merged_data and sample_info

  # ------ get args ------
  n_group = sample_info$group %>% unique %>% length
  name_group = sample_info$group %>% unique

  thres_edit_rate_indi = edit_rate_indi
  thres_edit_rate_mean = edit_rate_mean
  thres_miss_rate = missing_rate
  thres_left_n_total = left_n_total

  # thres_miss_n_in_each_group
  if( (left_n_in_each_group %>% length) == n_group) {
    thres_left_n_in_each_group = left_n_in_each_group
  }
  if( (left_n_in_each_group %>% length) != n_group) {
    thres_left_n_in_each_group = (rep(left_n_in_each_group[1], n_group) %>% as.vector())
  }

  # ------------ df 2 mat -----------
  merged_mat = merged_data %>% tidyr::unite(symbol, chr, pos, type, sep = '_')

  # ------------ QC --------------
  cat(paste0("######    Quanlity Control  ######",
             "\n", "----------------------------------", "\n",
             "Sites when QC beginning: ", nrow(merged_mat),"\n",
             "----------------------------------", "\n"))

  cat(paste0("\n",
             "Params input:", "\n",
             "(1) edit_rate_indi: ", thres_edit_rate_indi, "\n",
             "(2) left_n_total:", thres_left_n_total, "\n",
             "(3) left_n_in_each_group: c(",
             paste(thres_left_n_in_each_group, collapse = ','),") ","\n",
             "(4) edit_rate_mean: ", thres_edit_rate_mean, "\n",
             "(5) miss_rate: ", thres_miss_rate, "\n",
             "----------------------------------", "\n"
  ))

  # QC: edit_rate_indi
  if( !is.na(thres_edit_rate_indi)){
    mat_nosymbol = merged_mat %>% tibble::column_to_rownames("symbol") %>% as.matrix

    # num of NA
    sites_raw_na = sum(is.na(mat_nosymbol))
    # num of NA newly added (ie. not passing thres)
    sites_removed = sum((mat_nosymbol < thres_edit_rate_indi) %>% as.vector %>% na.omit)

    # set sites, which is not passing thres, to NA
    mat_nosymbol[mat_nosymbol < thres_edit_rate_indi] = NA

    # num of NA now.
    sites_left = sum(is.na(mat_nosymbol))

    merged_mat = mat_nosymbol %>% as.data.frame() %>% tibble::rownames_to_column("symbol")
  }

  # # check dups
  # merged_mat[duplicated(merged_mat$symbol) == T,]

  cat(paste0(  "\n",
               "running QC by edit_rate_indi >= ",thres_edit_rate_indi,
               ";\n  Sites set to NA: ",sites_removed,"\n"))

  # QC: left_n_total
  sites_raw=sites_left=NULL
  sites_raw=nrow(merged_mat)
  if( !is.na(thres_left_n_total)){
    # calc miss_n_total
    left_n_total = apply(merged_mat, 1, function(x) {sum(!is.na(x))})
    # FALSE if not passing
    qc_left_n_total = (left_n_total >= thres_left_n_total) %>% as.logical()
    merged_mat = merged_mat[qc_left_n_total, ]
  }
  sites_left=nrow(merged_mat)

  # # check dups
  # merged_mat[duplicated(merged_mat$symbol) == T,]

  cat(paste0( "\n",
              "running QC by left_n_total >= ",thres_left_n_total,
              ";\n  Sites filtered: ",sites_raw, " --> ", sites_left,"\n"))

  # QC: left_n_in_each_group
  sites_raw=sites_left=NULL
  sites_raw=nrow(merged_mat)
  if( !is.na(thres_left_n_in_each_group[1])){
    # sampleID in each group

    left_n_in_each_group = NULL
    qc_left_n_in_each_group = NULL
    qc_left_n_in_group = NULL

    # calc miss_n_in_each_group
    for(i in seq(1:n_group)) {
      # i=1

      s_list = NULL
      s_list = sample_info %>% dplyr::filter(group == name_group[i]) %>% pull(sampleid)
      left_n_in_each_group = apply(merged_mat %>% dplyr::select(symbol, all_of(s_list)),
                                   1, function(x) {sum(!is.na(x))})
      qc_left_n_in_each_group = (left_n_in_each_group >= thres_left_n_in_each_group[i])
      if(i == 1){
        qc_left_n_in_group = qc_left_n_in_each_group
      } else {
        qc_left_n_in_group = qc_left_n_in_each_group * qc_left_n_in_group
      }
    }

    # FALSE if not passing
    qc_left_n_in_group = qc_left_n_in_group %>% as.logical()
    # print(qc_left_n_in_group[1:50])
    merged_mat = merged_mat[qc_left_n_in_group, ]
  }
  sites_left=nrow(merged_mat)

  # # check dups
  # merged_mat[duplicated(merged_mat$symbol) == T,]

  cat(paste0("\n",
             "running QC by qc_left_n_in_group >= c(",
             paste(thres_left_n_in_each_group, collapse = ','),") in each group",
             ";\n  Sites filtered: ",sites_raw, " --> ",sites_left,"\n"))

  # QC: edit_rate_mean
  sites_raw=sites_left=NULL
  sites_raw=nrow(merged_mat)
  if( !is.na(thres_edit_rate_mean)){
    rownames(merged_mat) = NULL
    # calc mean_edit_rate
    mean_edit_rate = rowMeans(merged_mat %>% tibble::column_to_rownames("symbol") %>% as.matrix, na.rm=TRUE)
    # FALSE if not passing
    thres_edit_rate_mean = ifelse(thres_edit_rate_mean <= 0.5, thres_edit_rate_mean, 1-thres_edit_rate_mean)

    # sites edit rate >= thres_edit_rate_mean
    qc_edit_rate_dn = (mean_edit_rate >= thres_edit_rate_mean)
    # merged_mat = merged_mat[qc_edit_rate_dn, ]

    # sites edit rate <= (1-thres_edit_rate_mean)
    qc_edit_rate_up = (mean_edit_rate <= (1-thres_edit_rate_mean))
    qc_edit_rate = (qc_edit_rate_dn * qc_edit_rate_up) %>% as.logical()

    merged_mat = merged_mat[qc_edit_rate, ]

    # num of qc
    sites_qc_edit_rate_dn = sum(!qc_edit_rate_dn)
    sites_qc_edit_rate_up = sum(!qc_edit_rate_up)

    # qc_edit_rate = (qc_edit_rate_dn * qc_edit_rate_up) %>% as.logical()
    # merged_mat = merged_mat[qc_edit_rate, ]
  }
  sites_left=nrow(merged_mat)

  # # check dups
  # merged_mat[duplicated(merged_mat$symbol) == T,]

  cat(paste0("\n",
             "running QC by edit_rate_mean >= ",thres_edit_rate_mean,
             "\n  Sites filtered: ",sites_raw, " --> ",sites_left,
             "; including,",
             "\n    sites removed by edit_rate_mean <= ",thres_edit_rate_mean,": ",sites_qc_edit_rate_dn,
             "\n    sites removed by edit_rate_mean >= ",(1-thres_edit_rate_mean),": ",sites_qc_edit_rate_up,"\n"))

  # QC: miss_rate
  sites_raw=sites_left=NULL
  sites_raw=nrow(merged_mat)
  if( !is.na(thres_miss_rate)){
    # calc miss_rate
    miss_rate = apply(merged_mat, 1, function(x) {sum(is.na(x))/length(x)})
    # FALSE if not passing
    qc_miss_rate = (miss_rate <= thres_miss_rate) %>% as.logical()
    merged_mat = merged_mat[qc_miss_rate, ]
  }
  sites_left=nrow(merged_mat)

  cat(paste0("\n",
             "running QC by miss_rate <= ",thres_miss_rate,
             "\n  Sites filtered: ",sites_raw, " --> ", sites_left,"\n"))

  cat(paste0("\n", "----------------------------------", "\n",
             "Sites when QC finished: ", nrow(merged_mat),"\n",
             "----------------------------------", "\n"))

  return(merged_mat)
}


#######################
##        diff       ##
#######################
calc_mean = function(df, sample_info){

  # TODO
  # format check of df
  # if (merged_data) then stop(need qc)
  # if (merged_mat) then ok

  library(dplyr)
  merged_mat = df
  rownames(merged_mat) = NULL
  # ------------ df 2 mat -----------
  # merged_mat = merged_data %>% tidyr::unite(symbol, chr, pos, type, sep = '_')
  symbol = merged_mat$symbol
  merged_df_nosymbol = merged_mat %>% tibble::column_to_rownames("symbol")

  # get args
  n_group = sample_info$group %>% unique %>% length
  name_group = sample_info$group %>% unique

  # calc means
  for(i in seq(1:n_group)) {
    # i=1

    s_list = NULL
    s_list = sample_info %>% dplyr::filter(group == name_group[i]) %>% pull(sampleid)

    symbol_mean = merged_df_nosymbol %>% dplyr::select(all_of(s_list)) %>%
      rowMeans(na.rm = T) %>% cbind() %>% as.data.frame() %>%
      tibble::rownames_to_column("symbol")
    colnames(symbol_mean) =  c("symbol", name_group[i])

    if(i == 1){
      df_mean = symbol_mean
    } else {
      df_mean = full_join(df_mean, symbol_mean, by = c("symbol"))
    }
  }

  return(df_mean)
}


diff_edit_site = function(df, sample_info, comparison) {

  # TODO
  # format check of df
  # if (merged_data) then stop(need qc)
  # if (merged_mat) then ok


  # merged_mat = df
  # df = mat_qc
  # comparison = c("Han","Tibetan")  # NOTE: cont = comparison[1]

  rownames(df) = NULL


  library(dplyr)
  library(ezlimma)

  # 重排样本顺序
  symbol = df$symbol
  mat = df %>% dplyr::select(-symbol)
  samples = colnames(mat)

  ID_A = sample_info %>% dplyr::filter(group==comparison[1]) %>% pull(sampleid)
  ID_B = sample_info %>% dplyr::filter(group==comparison[2]) %>% pull(sampleid)

  mat_comparison = mat %>% dplyr::select(all_of(ID_A),
                                         all_of(ID_B))
  rownames(mat_comparison) = symbol

  # sort by c(ID_A, ID_B)
  sample_order = data.frame(sampleid=c(ID_A,ID_B))
  sample_info = sample_order %>% dplyr::left_join(sample_info, by= c("sampleid"))

  # define group
  grp = sample_info$group
  # define matrix
  M = mat_comparison %>% as.matrix

  # define comparison
  contr = paste0(comparison[2],"-",comparison[1])
  contr.v <- c(comparison[1], comparison[2], contr)
  names(contr.v) <- c(comparison[1], comparison[2], contr)

  # run limma & tidy results
  diff_res <- limma_contrasts(M, grp=grp, contrast.v = contr.v)

  # add delta for percentage data, such as methylation, editing rate.
  diff_res[,paste0(comparison[2],"-",comparison[1],".delta")] = diff_res[2] - diff_res[1]
  diff_res = diff_res[,c(1,2,15,14,11,12)]

  return(diff_res)
}

#######################
##       Anno        ##
#######################
anno_edit_site = function(df, annoFile) {

  # TODO
  # format check of df
  # if (symbol split to 3 piece) then check(need chr pos type)
  # if (chr pos type) then ok
  # other; then stop (need chr pos type)

  #
  library(dplyr)
  diff_sig = df

  # Read-in anno
  anno_symbol = anno %>% dplyr::mutate(symbol = paste0(Region,"_",Position,"_",Ref,Ed))

  # check
  # anno %>% dplyr::filter(Region=="chr11") %>% dplyr::filter(Position==65441384)

  df_anno = diff_sig %>%
    tibble::rownames_to_column("symbol") %>%
    # tidyr::separate(col="symbol",into=c("chr","pos","type"), sep="[_:]") %>%
    left_join(anno_symbol, by = c("symbol"))

  # merged_data = merged_mat %>% tidyr::separate(col="symbol",into=c("chr","pos","type"), sep="[_:]")

  return(df_anno)
}
