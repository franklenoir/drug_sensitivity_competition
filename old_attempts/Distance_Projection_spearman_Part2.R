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
library(uwot)

### this is dataset provided by merve of genes that are never expressed, used to filter out genes that we are not interested in
#never_exp <- read.table("/export/wflenoir/dream_competition/5044_never_exp_GMM_dream_comp.txt", quote="\"", comment.char="")
sometimes_exp <- read.table("/export/wflenoir/dream_competition/6625_sometimes_exp_GMM_dream_comp.txt", quote="\"", comment.char="")

###creating a matching cell_id name
cell_id <- c("ASPC1","DU145","EFO21","NCIH1793","HCC1143","LNCAPCLONEFGC","U87MG")

### CCLE expression
ccle_dat <- read.csv("/export/wflenoir/dream_competition/data_ccle_RNAseq_DREAMv2_FIXED.csv")

###filtering for overlapping genes
#ccle_dat <- ccle_dat %>% filter(X %!in% never_exp$V1) #%>% filter(X %in% mean_exp_known$gene)
ccle_dat <- ccle_dat %>% filter(X %in% sometimes_exp$V1) #%>% filter(X %in% mean_exp_known$gene)

colnames(ccle_dat)[1] <- "gene"
ccle_cells <- colnames(ccle_dat)[-1]


### At this point ccle_dat and mean_exp_known have both been put in a form where both can be used, pushed through umap
spearman_pipe <- function(rna_seq_dat){
  
  rows <- rna_seq_dat$gene
  rna_seq_dat$gene <- NULL
  cols <- colnames(rna_seq_dat)
  rna_seq_dat <- normalize.quantiles(as.matrix(rna_seq_dat)) %>% as.data.frame()
  
  colnames(rna_seq_dat) <- cols
  rownames(rna_seq_dat) <- rows
  
  # dat_umap <- uwot::umap(as.data.frame(t(rna_seq_dat)),n_neighbors = 30)
  # dat_umap_df <- as.data.frame(dat_umap)
  # rownames(dat_umap_df) <- cols
  # 
  dist_df <- cor(rna_seq_dat,method = "spearman")
  return(dist_df)
}

#perturb_dist <- umap_pipe(mean_exp_known)
#browser()
ccle_dist <- spearman_pipe(ccle_dat)

###this is sensitivity date of 11 cell lines against 30 drugs. 
sensitivity_pred <- read.csv("/export/wflenoir/dream_competition/sens_v2.txt") %>% filter(X %!in% c("PANC1","HSTS","KRJ1","HF2597")) 
#sensitivity_pred$X <- c("ASPC1","DU145","EFO21","NCIH1793","HCC1143","HF2597","HSTS","KRJ1","LNCAPCLONEFGC","PANC1","U87MG")
sensitivity_pred$X <- c("ASPC1","DU145","EFO21","NCIH1793","HCC1143","LNCAPCLONEFGC","U87MG")
sensitivity_pred <- sensitivity_pred %>% column_to_rownames("X")

drugs <- colnames(sensitivity_pred)

ccle_dist_dat <- ccle_dist %>% as.matrix() %>% as.data.frame()
ccle_dist_dat <- ccle_dist_dat %>% dplyr::select(ASPC1,DU145,EFO21,NCIH1793,HCC1143,LNCAPCLONEFGC,U87MG)
ccle_dist_dat$cell = rownames(ccle_dist_dat)
ccle_dist_dat <- ccle_dist_dat %>% filter(cell %in% colnames(ccle_dist_dat))

sensitivity_pred$cell = rownames(sensitivity_pred)

lm_perturb_dat <- merge(ccle_dist_dat,sensitivity_pred)

results <- ccle_dist %>% as.matrix() %>% as.data.frame() %>% dplyr::select()

for (drug in drugs){
  sensitivity_pred_models <- lm(eval(parse(text = drug))~ASPC1+DU145+EFO21+NCIH1793+HCC1143+LNCAPCLONEFGC+U87MG,lm_perturb_dat)
  temp <- ccle_dist %>% as.matrix() %>% as.data.frame()
  temp <- temp %>% dplyr::select(ASPC1,DU145,EFO21,NCIH1793,HCC1143,LNCAPCLONEFGC,U87MG)
  temp$cell = rownames(temp)
  lm_estimates <- predict(sensitivity_pred_models,temp)
  results[drug] <- lm_estimates
}

#browser()

head(results)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#results <- temp_results
#rownames(results) <- results$
#results$X <- NULL
res <- range01(results)
head(res)
#browser()

###ordering to the template

template_final <- read.csv("/export/wflenoir/dream_competition/template_leaderboard.csv")
res$cell_line <- rownames(res)
rownames(res) <- NULL

try_res <- res[colnames(template_final)]
try_res <- try_res[match(template_final$cell_line, try_res$cell_line),]

write.csv(try_res,"./submission_v8.csv",row.names = FALSE)

