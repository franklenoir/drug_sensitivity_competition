# Docker inheritance
FROM bioconductor/bioconductor_docker:devel@sha256:0a7f63c1e462c6f6372f4b33efc6c61b8fdc8ef34633d2169b7a8b334884bb72

# Update apt-get
RUN apt-get update \
	&& apt-get install -y --no-install-recommends build-essential r-base \
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Install required Bioconductor package
RUN R -e 'BiocManager::install("edgeR")'
RUN R -e 'BiocManager::install("fgsea")'
RUN R -e 'BiocManager::install("preprocessCore")'

#install reps
RUN R -e "install.packages(c('dplyr','tidyR','factoextra','data.table','readr','ClusterR'),  repos='http://cran.rstudio.com/')"

RUN mkdir /input
RUN mkdir /output

COPY GMM_modeling_Part1.R /usr/local/bin/GMM_modeling_Part1.R
COPY Sensitivity_Calculation_Part2.R /usr/local/bin/Sensitivity_Calculation_Part2.R
COPY Distance_Projection_spearman_Part3.R /usr/local/bin/Sensitivity_Calculation_Part3.R
COPY RNAseq_Perturb/ /usr/local/bin/input/RNAseq_Perturb/
COPY data_ccle_RNAseq_DREAMv2_FIXED.csv /usr/local/bin/input/data_ccle_RNAseq_DREAMv2_FIXED.csv 
COPY template_final.csv /usr/local/bin/input/template_final.csv
COPY run.sh /run.sh
COPY c2.cp.kegg.v7.0.symbols.gmt /usr/local/bin/input/c2.cp.kegg.v7.0.symbols.gmt

RUN chmod 775 /usr/local/bin/GMM_modeling_Part1.R
RUN chmod 775 /run.sh
RUN chmod 775 /usr/local/bin/Sensitivity_Calculation_Part2.R
RUN chmod 775 /usr/local/bin/Sensitivity_Calculation_Part3.R

ENTRYPOINT ["/bin/bash", "/run.sh"]
