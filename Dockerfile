############################################################
# Dockerfile to build PrimeDesign CLI and WebApp
############################################################

# Set the base image to anaconda
FROM continuumio/miniconda3

# File Author / Maintainer
MAINTAINER Jonathan Y. Hsu

ENV SHELL bash

#RUN conda install r-base
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

# Update packages of the docker system
RUN apt-get update && apt-get install gsl-bin libgsl0-dev -y && apt-get install libgomp1 -y && apt-get clean

# Install crispritz package (change the dafault version of python to 3.8)
RUN conda update -n base -c defaults conda
RUN conda install python=3.8 -y
RUN conda install crispritz -y && conda clean --all -y

# Add website dependencies
RUN pip install dash==1.9.1  # Dash core
RUN pip install dash-bio==0.4.8 # Dash bio
RUN pip install dash_daq
RUN pip install dash-bootstrap-components
RUN pip install seqfold
RUN pip install gunicorn
RUN pip install biopython

# Create environment
COPY PrimeDesign /PrimeDesign
#RUN mkdir /tmp/UPLOADS_FOLDER
#RUN mkdir /tmp/RESULTS_FOLDER

# Reroute to enable the PrimeDesign CLI and WebApp
WORKDIR /PrimeDesign
EXPOSE 9994
RUN chmod +x /PrimeDesign/web_app/start_server_docker.sh
ENTRYPOINT ["python", "/PrimeDesign/primedesign_router.py"]
