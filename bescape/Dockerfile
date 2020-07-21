## Start from this Docker image
FROM ubuntu:19.10

# Set locale && Set environment variable as noninteractive
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y locales && \
    sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8 

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app

COPY docker_files/requirements.txt ./

# ubuntu shipped with python --> update and install pip3
# Install any needed packages specified in requirements.txt
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y python3-pip && \
    pip3 install --trusted-host pypi.python.org -r requirements.txt


# Instal R and its dependiencies
RUN apt-get update && apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev libssl-dev r-base=3.6.1-4 pandoc

# Run R scripts to install R packages
RUN R -e "install.packages(c('xml2', 'roxygen2', 'rversions', 'devtools', 'tidyverse', 'BiocManager', 'nnls', 'hdf5r'))" && \
    R -e "BiocManager::install('Biobase')" && \
    R -e "devtools::install_github('xuranw/MuSiC')" && \
    R -e "devtools::install_github('meichendong/SCDC')" && \
    R -e "devtools::install_github('asrinivasan-oa/ganalyse')" && \
    R -e "devtools::install_github('GfellerLab/EPIC', build_vignettes=FALSE)" && \
    R -e "devtools::install_github('meichendong/SCDC')" && \
    R -e "install.packages('BisqueRNA')"

COPY docker_files/bescape.py docker_files/_helper_deconv.py docker_files/musicpy.py docker_files/script_music.R docker_files/epicpy.py docker_files/epicpy_ref.py docker_files/SCDC_tkt_simplified.Rmd docker_files/scdcpy.py docker_files/scdc_input.R docker_files/epic.R ./

# Make script executable
RUN chmod a+x bescape.py && \
    chmod a+x musicpy.py && \
    chmod a+x script_music.R && \
    chmod a+x epicpy.py && \
    chmod a+x epicpy_ref.py && \
    chmod a+x SCDC_tkt_simplified.Rmd && \
    chmod a+x scdcpy.py && \
    chmod a+x scdc_input.R

# Make Docker container executable
CMD ["python", "/app/musicpy.py"]
