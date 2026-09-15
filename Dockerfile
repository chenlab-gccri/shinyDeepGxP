# Use the correct base image (match R version in renv.lock)
FROM rocker/shiny:4.4.2

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    bzip2 \
    ca-certificates \
    build-essential \
    software-properties-common \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libxt-dev \
    libhdf5-dev \
    patch \
    libglpk40 \
    libglpk-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

## 3. Install Miniconda
#RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
#    bash /tmp/miniconda.sh -b -p /opt/conda && \
#    rm /tmp/miniconda.sh
#
#ENV PATH="/opt/conda/bin:$PATH"
#
## 4. Create clean Conda env and install Python packages
#RUN conda create -n r-DL_env python=3.10 -y && \
#    conda install -n r-DL_env -c conda-forge \
#        tensorflow-cpu=2.10 \
#        numpy=1.21 \
#        pandas \
#        scikit-learn -y && \
#    conda run -n r-DL_env pip install magic-impute
	
# 3. Install Miniforge
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh \
    -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh

ENV PATH="/opt/conda/bin:$PATH"	
	
# 4. Create clean Conda env and install Python packages
RUN conda create -n r-DL_env python=3.10 -y && \
    conda install -n r-DL_env \
        tensorflow-cpu=2.10 \
        numpy=1.21 \
        pandas \
        scikit-learn -y && \
    conda run -n r-DL_env pip install magic-impute

# 5. R packages
RUN Rscript -e "install.packages('reticulate')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Rmagic/Rmagic_2.0.3.tar.gz', repos = NULL, type = 'source')"
RUN Rscript -e "install.packages('shinyjs')"
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "install.packages('shinyWidgets')"
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('Seurat')"
RUN Rscript -e "install.packages('hdf5r')"
RUN Rscript -e "install.packages('remotes')"
#RUN Rscript -e "remotes::install_github('mojaveazure/seurat-disk')"
RUN Rscript -e "remotes::install_github('mojaveazure/seurat-disk#198')"
RUN Rscript -e "install.packages('keras')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('DOSE')"
RUN Rscript -e "install.packages('rintrojs')"
RUN Rscript -e "install.packages('shinyBS')"
RUN Rscript -e "install.packages('shinydashboardPlus')"
RUN Rscript -e "install.packages('shinycssloaders')"
RUN Rscript -e "install.packages('shinydashboard')"
RUN Rscript -e "install.packages('DT')"
RUN Rscript -e "install.packages('tibble')"
RUN Rscript -e "install.packages('methods')"
RUN Rscript -e "install.packages('zip')"
RUN Rscript -e "install.packages('stats')"
RUN Rscript -e "install.packages('VennDiagram')"
RUN Rscript -e "install.packages('arrow')"
RUN Rscript -e "install.packages('data.table')"


RUN conda clean -afy

#6. Copy all app files, including renv infrastructure
COPY . /srv/shiny-server/
WORKDIR /srv/shiny-server/

# Set ownership for shiny
RUN chown -R shiny:shiny /srv/shiny-server

## Expose Shiny port
#EXPOSE 3838

## Launch Shiny server
#CMD ["/usr/bin/shiny-server"]

USER shiny
CMD ["shiny-server"]