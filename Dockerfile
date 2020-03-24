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

# Add website dependencies
RUN pip install dash==1.9.1  # Dash core
RUN pip install dash-bio==0.4.8 # Dash bio
RUN pip install gunicorn

# Create environment
COPY primeDesign /primeDesign
#RUN mkdir /tmp/UPLOADS_FOLDER
#RUN mkdir /tmp/RESULTS_FOLDER

# Reroute to enable the PrimeDesign CLI and WebApp
WORKDIR /primeDesign
EXPOSE 9993
RUN chmod +x /primeDesign/web_app/start_server_docker.sh
ENTRYPOINT ["python", "/primeDesign/primedesign_router.py"]