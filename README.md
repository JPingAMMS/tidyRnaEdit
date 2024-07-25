# Example run

```r
######################################################################
# packages
######################################################################
library(tidyRnaEdit)

######################################################################
# Read-in example data
######################################################################
# Read-in the example dataset.
data(T15H10,package="tidyRnaEdit")

# This example dataset contains two files. 
# 1. editing_rate matrix: merged_data
# 2. sample_info: sample_info

# get a look at this dataset
merged_data %>% head
sample_info %>% head

# Read-in annotation file.
data(anno_hg38_REDIportal,package="tidyRnaEdit")

# get a look at annotation infomation.
anno %>% head

######################################################################
# Count num. of editing types
######################################################################
# count
edit_type_count = edit_type_count(merged_data)
# save
write.csv(edit_type_count, file=paste0(path_save_type,"/mat_raw_count.csv"), row.names = FALSE)

# plot without save
plot_edit_type_count(edit_type_count)

# plot and save
path_save_type="."  # change to your save_directory_path
plot_edit_type_count(edit_type_count, save = TRUE,
                     file=paste0(path_save_type,"/plot_edit_type_count.pdf"))

######################################################################
# Filtering variants
######################################################################
mat_qc = filter_variant(merged_data, sample_info,
                        edit_rate_indi = 0.01,
                        edit_rate_mean = 0.05,
                        missing_rate = NA,
                        left_n_total = NA,
                        left_n_in_each_group = 5)

path_save="."  # change to your save_directory_path
write.csv(mat_qc, file=gzfile(paste0(path_save,"/mat_qc.csv.gz")), row.names = FALSE)

mat_qc[duplicated(mat_qc$symbol) == TRUE,]

######################################################################
# diff_analy_limma for genes
######################################################################
mat_mean = calc_mean(mat_qc, sample_info = sample_info)

diff_res = diff_edit_site(mat_qc,
                          sample_info = sample_info,
                          comparison = c("Han", "Tibetan"))  # control_group = "Han"

# threshold
FDR_thres = 0.05
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

path_save="."  # change to your save_directory_path
write.csv(diff_sig_anno, file=paste0(path_save,"/diff_sig_anno.csv"), row.names = FALSE)

# recoding site
diff_sig_recoding_site = diff_sig_anno %>% dplyr::filter(Func.wgEncodeGencodeBasicV34=="exonic")

diff_sig_recoding_site %>% head
```

# Run with your data

## Download the example data of REDItools outTables.

We provided an example data at [github](https://github.com/JPingAMMS/tidyRnaEdit_exampledata), which you can download at your PC.
Then, unzip the `tidyRnaEdit_exampledata.zip`.

The files named `*outTable.gz` are REDItools outTables. Only chr17 is provided for test.
The file `reditools_sample.xlsx` is a dataframe of 3 columns, including filepath, sampleid, group.

NOTE:
1. `*outTable.gz` is the standard output by REDItools.
2. When you generate your `reditools_sample.xlsx` file, You should KEEP (not modify) the column names.

## Download the annotation file.

A curated annotated file is available at [REDIportal](http://srv00.recas.ba.infn.it/atlas/download.htmcl).

The annotation file with hg38 positions is already to use in this package, namely, the `anno`.

## Tutorials with the example data.

```r
######################################################################
# packages
######################################################################
library(tidyRnaEdit)


######################################################################
# Path
######################################################################
#--------- work space ------------
setwd("/data/pingj/soft/tidyRnaEdit_exampledata/")  # change to your_path_where_saved_the_outTables

#--------- path_data ------------
# sample_info (3 cols: filename, sampleid, group)
path_sample = "reditools_samples.xlsx"  # change to your_path_where_saved_the_outTables

# anno (using `anno` in our package, or you can download from REDIportal)
path_anno = "TABLE1_hg38.txt.gz"  # change to your_path_where_saved_the_outTables

#--------- path_save ------------
# root
path_save = "."  # change to your save_directory_path
# editing_type
path_save_type        = paste0(path_save, "")   # change to your save_directory_path
# OverallEditing
path_save_overalledit = paste0(path_save, "")   # change to your save_directory_path

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

# Then, you can do down-stream analysis, follow the `example run` section. Enjoy!
```