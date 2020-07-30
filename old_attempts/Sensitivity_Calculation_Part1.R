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

calc_cmpd_path <- function(cell_df) {
  
  # getting cmpds to iterate
  cpmds <- cell_df %>% filter(compound_id %!in% c("cmpd_dmso", "cmpd_untreated")) %>% 
    distinct(compound_id) %>% 
    pull()
  
  diffexp_list <- lapply(cpmds,
                         function(i) {
                           
                           
                           wide <- cell_df %>%
                             filter(compound_id %in% c("cmpd_dmso",
                                                       i))  %>%
                             select(ids, count, gene) %>%
                             pivot_wider(names_from = ids,
                                         values_from = count,
                                         id_cols = gene) %>%
                             column_to_rownames("gene")
                           
                           groups <- str_extract(colnames(wide), "cmpd_[[:alpha:]]+")
                           
                           fit <- lmFit(wide, design = as.numeric(as.factor(groups)))
                           fit <- eBayes(fit)
                           tt <- topTable(fit, number = Inf)
                           
                           return(tt)
                         })
  
  names(diffexp_list) <- cpmds
  return(diffexp_list)
  
}

run_gsea <- function(diffexp, path =  gmtPathways("./c2.cp.kegg.v7.0.symbols.gmt")) {
  lapply(diffexp, function(x) {
    
    res <- x %>%
      rownames_to_column("gene") %>%
      select(gene, t) %>%
      deframe()
    
    fgseaRes <- fgsea(pathways = path,
                      stats = res,
                      nperm = 500)
    
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
      arrange(desc(NES), padj)
    
    return(fgseaResTidy)
  })
}

plot_gsea <- function(gsea_list,cell_name = NULL) {
  apop <- lapply(gsea_list, function(x) {
    x %>%
      select(pathway, NES, padj) %>%
      filter(pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450"))
      #filter(pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450","KEGG_APOPTOSIS"))
  })
  
  apop_df <- bind_rows(apop, .id = "compound") %>%
    mutate(rank = rank(NES))
  #apop_df$NES <- abs(apop_df$NES)
  apop_df <- apop_df %>% group_by(pathway) %>% mutate(normed_NES = scale(NES))
  # apop_df %>% spread(pathway,normed_NES) -> try
  # try1 <- try %>% dplyr::select(compound,KEGG_APOPTOSIS)
  # try2 <- try %>% dplyr::select(compound,KEGG_DRUG_METABOLISM_CYTOCHROME_P450)
  # try <- merge(try1,try2)
  # try %>% ggplot(aes(KEGG_APOPTOSIS,KEGG_DRUG_METABOLISM_CYTOCHROME_P450)) + geom_point() + theme_cowplot()
  # 
  
  apop_df <- aggregate(
    apop_df$NES,
      by = list(apop_df$compound),
      FUN = sum,
      na.rm = TRUE
    )
  colnames(apop_df) <- c("compound","NES")
  return(apop_df)
  
  # ggplot(apop_df) +
  #   geom_point(aes(fct_reorder(compound, NES), NES,
  #                  size = padj)) +
  #   coord_flip() +
  #   theme_classic() +
  #   ylab(paste("KEGG_apoptosis enrichment scores for cell line", cell_name)) +
  #   xlab("")
  
}

data_process <- function(dat){
  dat_tidy <- dat %>% 
    rownames_to_column("gene") %>% 
    gather(key = "ids",
           value = "count",
           -gene)
  # adding labels
  dat_tidy <- dat_tidy %>% 
    mutate(cell_line = str_extract(ids, "(?<=[\\d]_[\\d]{1,2}_)[[:alnum:]]+(?=_.+)"),
           compound_id = str_extract(ids, "cmpd_[[:alpha:]]+"),
           concentration = str_extract(ids, "(?<=[[:alpha:]]_)[\\d\\.]+(?=_\\d)"),
           treatment_time = str_extract(ids, "(?<=[\\d\\.]_)[\\d]+(?=_.+)"),
           replicate = str_extract(ids, "(?<=_)\\d+$"))
  cell_type <- dat_tidy %>% pull(cell_line) %>% unique()
  print(cell_type)
  dat_pathways <- calc_cmpd_path(dat_tidy)
  dat_gsea <- run_gsea(dat_pathways)
  result <- plot_gsea(dat_gsea) %>% dplyr::select(compound,NES)
  result$cell <- cell_type
  return(result)
  }

dat1 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/ASPC1-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat2 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/DU145-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat3 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/EFO21-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat4 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/HCC1143-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat5 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/HF2597-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat6 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/HSTS-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat7 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/KRJ1-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat8 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/H1793-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat9 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/LNCAP-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat10 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/PANC1-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)
dat11 <- read.table("/export/wflenoir/dream_competition/rnaseq_concealed/U87-RNAseq-Perturbations.csv",sep = ",",header = T,row.names = 1)

df.list <- list(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11)
rm(dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9,dat10,dat11)
res <- lapply(df.list, function(x) data_process(x))
head(res)
res <- bind_rows(res, .id = "num")
temp <- res %>% dplyr::select(-num) %>% spread(cell,NES)
cols <- temp$compound
temp$compound <- NULL
temp <- t(temp)
colnames(temp) <- cols


write.csv(temp,"./sens_v2.txt")
