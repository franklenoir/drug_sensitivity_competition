rm(list = ls())
if (is.integer(dev.list())){dev.off()}
cat("\014")
set.seed(1)
setwd("/export/wflenoir/dream_competition/")

head2 <- function(x) head(x)[,1:5]
`%!in%` <- Negate(`%in%`)
library(tidyverse)
library(factoextra)
library(cowplot)
library(edgeR)
library(fgsea)
library(preprocessCore)

### this is dataset provided by merve of genes that are never expressed, used to filter out genes that we are not interested in
never_exp <- read.table("/export/wflenoir/dream_competition/5044_never_exp_GMM_dream_comp.txt", quote="\"", comment.char="")
#sometimes_exp <- read.table("/export/wflenoir/dream_competition/6625_sometimes_exp_GMM_dream_comp.txt", quote="\"", comment.char="")

### commented out chunk calculated average dmso and untreated expression of original 11 lines
# dat1 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/ASPC1-RNAseq-Perturbations.csv")
# dat2 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/DU145-RNAseq-Perturbations.csv")
# dat3 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/EFO21-RNAseq-Perturbations.csv")
# dat4 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/H1793-RNAseq-Perturbations.csv")
# dat5 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/HCC1143-RNAseq-Perturbations.csv")
# dat6 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/HF2597-RNAseq-Perturbations.csv")
# dat7 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/HSTS-RNAseq-Perturbations.csv")
# dat8 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/KRJ1-RNAseq-Perturbations.csv")
# dat9 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/LNCAP-RNAseq-Perturbations.csv")
# dat10 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/PANC1-RNAseq-Perturbations.csv")
# dat11 <- read.csv("/export/wflenoir/dream_competition/rnaseq_concealed/U87-RNAseq-Perturbations.csv")
# 
# df.list <- list(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11)
# rm(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11)
# 
# data_process <- function(dat,never_exp) {
#   rownames(dat) <- dat$hgnc_symbol
#   dat$hgnc_symbol = NULL
#   cell_name <-
#     unique(str_extract(
#       colnames(dat),
#       "(?<=[\\d]_[\\d]{1,2}_)[[:alnum:]]+(?=_.+)"
#     ))
#   print(cell_name)
#   sample_info <-
#     tibble(sample_id = colnames(dat), sample_sum = colSums(dat))
#   n_samples <- nrow(sample_info)
#   sample_info$cmpd <- str_extract(colnames(dat), "cmpd_[[:alpha:]]+")
#   sample_info$type <- "drug"
#   sample_info$type[sample_info$cmpd == "cmpd_untreated"] <-
#     "untreated"
#   sample_info$type[sample_info$cmpd == "cmpd_dmso"] <- "dmso"
#   dat$gene <- rownames(dat)
#   dat %>% filter(gene %!in% never_exp$V1)
#   dat <- dat %>% gather(sample_id, exp, -gene)
#   dat <- merge(dat, sample_info)
#   dat <- dat %>% filter(type != "drug")
#   dat_agg <- aggregate(
#     dat$exp,
#     by = list(dat$gene),
#     FUN = mean,
#     na.rm = TRUE
#   )
#   colnames(dat_agg) <- c("gene", "mean_exp")
#   dat_agg$cell = cell_name
#   return(dat_agg)
# }
# 
# res <- lapply(df.list, function(x) data_process(x,never_exp))
# browser()
# 
# res <- bind_rows(res, .id = "gene")
# colnames(res)[1] <- "num"
# res <- res %>% dplyr::select(-num)
# 
# try <- res %>% spread(cell,mean_exp)
# 
# write_delim(try,"./mean_exp_known.txt",delim = "\t")

mean_exp_known <- read.delim("/export/wflenoir/dream_competition/mean_exp_known.txt") ### original cell lines WT expression
ccle_dat <- read.csv("/export/wflenoir/dream_competition/data_ccle_RNAseq_DREAMv2_FIXED.csv") ### CCLE expression

###filtering for overlapping genes
ccle_dat <- ccle_dat %>% filter(X %!in% never_exp$V1) %>% filter(X %in% mean_exp_known$gene)
#ccle_dat <- ccle_dat %>% filter(X %in% sometimes_exp$V1) %>% filter(X %in% mean_exp_known$gene)
mean_exp_known <- mean_exp_known %>% filter(gene %in% ccle_dat$X)
#browser()
###adding string to mean exp to avoid name overlap
colnames(mean_exp_known) <- paste("orig", colnames(mean_exp_known), sep = "_")

###this is sensitivity date of 11 cell lines against 30 drugs. 
sensitivity_pred <- read.csv("/export/wflenoir/dream_competition/sens_v2.txt")

colnames(ccle_dat)[1] <- "gene"
colnames(mean_exp_known)[1] <- "gene"
ccle_cells <- colnames(ccle_dat)[-1]

rna_seq_merge <- merge(ccle_dat,mean_exp_known)
rows <- rna_seq_merge$gene
rna_seq_merge$gene <- NULL
cols <- colnames(rna_seq_merge)

rna_seq_merge <- normalize.quantiles(as.matrix(rna_seq_merge)) %>% as.data.frame()

colnames(rna_seq_merge) <- cols
rownames(rna_seq_merge) <- rows

head(rna_seq_merge)

dist_matrix <- as.matrix(dist(as.matrix(t(rna_seq_merge))))

cell_orig <- c()
cell_closest <- c()
cell_2_closest <- c()

temp_results <- sensitivity_pred %>% filter(X == "") 

sensitivity_values <- sensitivity_pred  %>% column_to_rownames(var = "X")
for (cell in ccle_cells){
  cell_orig <- c(cell_orig,cell)
  temp <- as.data.frame(dist_matrix[cell,colnames(mean_exp_known)[-1]]) %>% rownames_to_column("Cell_Line")#%>% dplyr::select(colnames(mean_exp_known)[-1])
  temp <- temp %>% arrange(`dist_matrix[cell, colnames(mean_exp_known)[-1]]`)
  
  close_cells <- temp$Cell_Line[1:2]
  distances <- temp$`dist_matrix[cell, colnames(mean_exp_known)[-1]]`[1:2]
  distances <- 1/distances 
  distances <- distances/sum(distances)

  close_cells <- str_remove_all(close_cells, "orig_")
  
  cell_closest <- c(cell_closest,close_cells[1])
  cell_2_closest <- c(cell_2_closest,close_cells[2])

  sens_dat <- sensitivity_values[close_cells,]
  sens_dat <- sens_dat * distances
  res_sen <- colSums(sens_dat) %>% t()  %>% as.data.frame() 
  res_sen$X <- cell
  temp_results <- rbind(temp_results,res_sen)
}

result_close_distance <- tibble(cell_orig,cell_closest,cell_2_closest)

browser()

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
results <- temp_results
rownames(results) <- results$X
results$X <- NULL
res <- range01(results)

browser()

###ordering to the template

template_final <- read.csv("/export/wflenoir/dream_competition/template_leaderboard.csv")
res$cell_line <- rownames(res)
rownames(res) <- NULL

try_res <- res[colnames(template_final)]
try_res <- try_res[match(template_final$cell_line, try_res$cell_line),]

write.csv(try_res,"./submission_v5.csv",row.names = FALSE)

