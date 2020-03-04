## Start from this Docker image
FROM ubuntu:19.10

# Set locale && Set environment variable as noninteractive
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y locales
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8 

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app

COPY docker_files/requirements.txt ./

# ubuntu shipped with python --> update and install pip3
RUN apt-get update &&  apt-get -y upgrade && apt-get install -y python3-pip

# Install any needed packages specified in requirements.txt
RUN pip3 install --trusted-host pypi.python.org -r requirements.txt

# Instal R and r-based deconvolution methods
RUN apt-get update && apt-get install -y libxml2 libxml2-dev libcurl4-gnutls-dev libssl-dev r-base=3.6.1-4

# Run R scripts to install R packages
RUN R -e "install.packages(c('xml2', 'roxygen2', 'rversions', 'devtools', 'tidyverse', 'BiocManager'))"
RUN R -e "BiocManager::install('Biobase')"
RUN R -e "devtools::install_github('xuranw/MuSiC')"
RUN R -e "devtools::install_github('meichendong/SCDC')"
RUN R -e "devtools::install_github('asrinivasan-oa/ganalyse')"
RUN R -e "install.packages('BisqueRNA')"

COPY docker_files/bescape.py ./
COPY docker_files/_helper_deconv.py ./
COPY docker_files/musicpy.py ./
COPY docker_files/script_music.R ./
COPY docker_files/epic.R ./


# Make script executable
RUN chmod a+x bescape.py
RUN chmod a+x musicpy.py
RUN chmod a+x script_music.R

# Make Docker container executable
CMD ["python", "/app/bescape.py"]
