FROM rocker/binder:3.6.0
LABEL maintainer='Yanxian Li'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}


RUN apt-get update -y\
&& apt-get install -y dpkg-dev zlib1g-dev libssl-dev libffi-dev\
&& apt-get install -y curl libcurl4-openssl-dev\
&& apt-get install -y git\
&& R -e "install.packages('BiocManager')"\
&& R -e "BiocManager::install('ComplexHeatmap')" 

RUN wget https://github.com/KrogdahlLab/AquaFly-SeawaterGutHealth-Aquaculture-2018/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-08-15/')); devtools::install_deps()"
