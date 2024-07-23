# Example run

```r
######################################################################
# packages
######################################################################
library(dplyr)
library(readr)
library(readxl)
source("E:/pj/gy习服.RNA编辑/R/utils.R")
source("E:/pj/gy习服.RNA编辑/R/utils_plot.R")


######################################################################
# Path
######################################################################
#--------- work space ------------
setwd("E:/pj/gy.RNA编辑/a2i_known_gz_T15H10")

#--------- path_data ------------
# sample_info (3 cols: filename, sampleid, group)
path_sample = "E:/pj/gy.RNA编辑/reditools_samples.xlsx"

# anno (download from REDIportal)
path_anno = "F:/AMS2018_RNAedit/TABLE1_hg38.txt.gz"

#--------- path_save ------------
# root
path_save = "E:/pj/gy.RNA编辑/T15H10"
# editing_type
path_save_type        = paste0(path_save, "")
# OverallEditing
path_save_overalledit = paste0(path_save, "")


######################################################################
# Read-in data
######################################################################
# sample info
sample_info = readxl::read_excel(path_sample)

# anno
anno = read.table(gzfile(path_anno), sep ='\t', fill=TRUE, header = TRUE)

# editing data
merged_data = read_redit(sample_info, soft = "reditools_known")

# save
write.csv(merged_data, file=gzfile(paste0(path_save,"/mat_raw.csv.gz")), row.names = FALSE)

merged_data = read.csv(gzfile(paste0(path_save,"/mat_raw.csv.gz")),header=TRUE)
save(merged_data, sample_info, file="T15H10.RData")

######################################################################
# Count num. of editing types
######################################################################
# count
rna_type_count = edit_type_count(merged_data)
# save
write.csv(rna_type_count, file=paste0(path_save_type,"/mat_raw_count.csv"), row.names = FALSE)

# plot
plot_edit_type_count(rna_type_count, save = T,
                     file=paste0(path_save_type,"/plot_rna_type_count.pdf"))

######################################################################
# Filtering variants
######################################################################
mat_qc = filter_variant(merged_data, sample_info,
                        edit_rate_indi = 0.01,
                        edit_rate_mean = 0.05,
                        missing_rate = NA,
                        left_n_total = NA,
                        left_n_in_each_group = 5)

write.csv(mat_qc, file=gzfile(paste0(path_save,"/mat_qc.csv.gz")), row.names = FALSE)

mat_qc[duplicated(mat_qc$symbol) == T,]

######################################################################
# diff_analy_limma for genes
######################################################################
mat_mean = calc_mean(mat_qc, sample_info = sample_info)

diff_res = diff_edit_site(mat_qc,
                          sample_info = sample_info,
                          comparison = c("Han", "Tibetan"))  # control_group = "Han"

# threshold
FDR_thres = 0.05/nrow(diff_res)
deltaEdit_thres = 0.30

# filter sig symbol
diff_sig = diff_res %>%
  dplyr::filter(abs(`Tibetan-Han.delta`)> deltaEdit_thres) %>%
  dplyr::filter(`Tibetan-Han.FDR` < FDR_thres)

# # diff_summary
# TODO
# omic6_diff_limma_num = meth_diff_limma_num(omic6_diff, deltaB_thres = deltaEdit_thres)
# omic6_diff_limma_num_total = meth_diff_limma_num_total(omic6_diff, deltaB_thres = deltaEdit_thres)

######################################################################
# annotate
######################################################################
diff_sig_anno = anno_edit_site(diff_sig, anno)

write.csv(diff_sig_anno, file=paste0(path_save,"/diff_sig_anno.csv"), row.names = FALSE)
```