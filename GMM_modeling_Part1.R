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

library(ClusterR)

print("step_1")

exp <- read.csv("/usr/local/bin/input/data_ccle_RNAseq_DREAMv2_FIXED.csv") 

exp$mean <- rowMeans(exp[-1])
exp$min <- apply(exp[-1], 1, FUN=min)
exp$max <- apply(exp[-1], 1, FUN=max)
exp$STD <- apply(exp[-1], 1, FUN=sd)

sub_exp <- exp %>% dplyr::select(mean,STD) %>% droplevels()

gmm = GMM(sub_exp, 3)
pr = predict_GMM(sub_exp, gmm$centroids, gmm$covariance_matrices, gmm$weights)
sub_exp$cluster_labels <- pr$cluster_labels
exp$cluster_labels <- pr$cluster_labels

agg = aggregate(exp$mean,
                by = list(exp$cluster_labels),
                FUN = mean)

group_fil <- agg %>% filter(x == min(agg$x)) %>% droplevels() %>% pull(Group.1) %>% as.numeric()

low_expressed_genes <- exp %>% filter(cluster_labels == group_fil) %>% droplevels() %>% dplyr::select(X)

write_delim(low_expressed_genes,"/usr/local/bin/input/never_exp_GMM_dream_comp.txt")
