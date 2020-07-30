rm(list = ls())
if (is.integer(dev.list())) {
  dev.off()
}
cat("\014")
set.seed(1)

head2 <- function(x)
  head(x)[, 1:5]
`%!in%` <- Negate(`%in%`)

library(readr)
library(dplyr)
library(tidyr)

library(stringr)
library(edgeR)
library(fgsea)

print("step_2")

###Function that identifies gene expression differences between compound vs "normal" conditions
calc_cmpd_path <- function(cell_df) {
  cpmds <-
    cell_df %>% filter(compound_id %!in% c("cmpd_dmso", "cmpd_untreated")) %>%
    distinct(compound_id) %>%
    pull()
  
  diffexp_list <- lapply(cpmds,
                         function(i) {
                           wide <- cell_df %>%
                             filter(compound_id %in% c("cmpd_dmso",
                                                       i))  %>%
                             dplyr::select(ids, count, gene) %>%
                             pivot_wider(names_from = ids,
                                         values_from = count,
                                         id_cols = gene) %>% as.data.frame()
                             #column_to_rownames("gene")
                          
                             row.names(wide) <- wide$gene
                             wide$gene <- NULL
                           
                           groups <-
                             str_extract(colnames(wide), "cmpd_[[:alpha:]]+")
                           
                           fit <-
                             lmFit(wide, design = as.numeric(as.factor(groups)))
                           fit <- eBayes(fit)
                           tt <- topTable(fit, number = Inf)
                           
                           return(tt)
                         })
  
  names(diffexp_list) <- cpmds
  return(diffexp_list)
  
}

###Function that identifies pathway differences based on gene expression differences
run_gsea <-
  function(diffexp,
           path =  gmtPathways("/usr/local/bin/input/c2.cp.kegg.v7.0.symbols.gmt")) {
    lapply(diffexp, function(x) {
	x$gene <- row.names(x)
	res <- x %>% dplyr::select(gene, t) #%>% deframe()
	#browser()
  t_vec <- res$t %>% as.vector()
  names(t_vec) <- res$gene
  res <- t_vec
      fgseaRes <- fgsea(pathways = path,
                        stats = res,
                        nperm = 500)
      
      fgseaResTidy <- fgseaRes %>%
        as_data_frame() %>%
        dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
        arrange(desc(NES), padj)
      
      return(fgseaResTidy)
    })
  }

###function originally used to create plots, but now just used to score drugs based on Kegg p450 pathway
plot_gsea <- function(gsea_list, cell_name = NULL) {
  apop <- lapply(gsea_list, function(x) {
    x %>%
      dplyr::select(pathway, NES, padj) %>%
      filter(pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450"))
  })
  
  apop_df <- bind_rows(apop, .id = "compound") %>%
    mutate(rank = rank(NES))
  apop_df <-
    apop_df %>% group_by(pathway) %>% mutate(normed_NES = scale(NES))
  
  apop_df <- aggregate(
    apop_df$NES,
    by = list(apop_df$compound),
    FUN = sum,
    na.rm = TRUE
  )
  colnames(apop_df) <- c("compound", "NES")
  return(apop_df)
  
}

###function putting it all together
data_process <- function(dat) {
	
  dat$gene <- row.names(dat)

  dat_tidy <- dat %>%
   # rownames_to_column("gene") %>%
   gather(key = "ids",
           value = "count", -gene)
  # adding labels
  dat_tidy <- dat_tidy %>%
    mutate(
      cell_line = str_extract(ids, "(?<=[\\d]_[\\d]{1,2}_)[[:alnum:]]+(?=_.+)"),
      compound_id = str_extract(ids, "cmpd_[[:alpha:]]+"),
      concentration = str_extract(ids, "(?<=[[:alpha:]]_)[\\d\\.]+(?=_\\d)"),
      treatment_time = str_extract(ids, "(?<=[\\d\\.]_)[\\d]+(?=_.+)"),
      replicate = str_extract(ids, "(?<=_)\\d+$")
    )
  cell_type <- dat_tidy %>% pull(cell_line) %>% unique()
  print(cell_type)
  dat_pathways <- calc_cmpd_path(dat_tidy)
  dat_gsea <- run_gsea(dat_pathways)
  result <- plot_gsea(dat_gsea) %>% dplyr::select(compound, NES)
  result$cell <- cell_type
  return(result)
}

###loading in data

dat1 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/ASPC1-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat2 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/DU145-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat3 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/EFO21-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat4 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/HCC1143-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat5 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/HF2597-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat6 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/HSTS-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat7 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/KRJ1-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat8 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/H1793-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat9 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/LNCAP-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat10 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/PANC1-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )
dat11 <-
  read.table(
    "/usr/local/bin/input/RNAseq_Perturb/U87-RNAseq-Perturbations.csv",
    sep = ",",
    header = T,
    row.names = 1
  )

df.list <-
  list(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11)
rm(dat1, dat2, dat3, dat4, dat5, dat6, dat7, dat8, dat9, dat10, dat11)

###processing and scoring each dataset/drug

res <- lapply(df.list, function(x)
  data_process(x))

head(res)

res <- bind_rows(res, .id = "num")
temp <- res %>% dplyr::select(-num) %>% spread(cell, NES)
cols <- temp$compound
temp$compound <- NULL
temp <- t(temp)
colnames(temp) <- cols
cells <- row.names(temp)
temp <- temp %>% as_data_frame()
temp$X <- cells

###outputting sensitivity scores.

write_delim(temp, "/usr/local/bin/input/sens_v2.txt")
