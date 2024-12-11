FROM mambaorg/micromamba:latest
LABEL authors="Eike Wacker" \
      description="Docker image containing R with the packages data.table and tidyverse for gwas-assoc pipeline plotting processes"


#RUN conda install -c conda-forge mamba && conda clean -a

COPY plotting_environment.yml /
RUN micromamba env create -f /plotting_environment.yml && micromamba clean --all --yes
ENV PATH /opt/conda/envs/eikematthias-rplotting/bin:$PATH
