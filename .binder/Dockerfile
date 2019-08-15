FROM rocker/binder:3.6.0
LABEL maintainer='Yanxian Li'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}
RUN wget https://github.com/KrogdahlLab/AquaFly-SeawaterGutHealth-Aquaculture-2018/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-08-15/')); devtools::install_deps()"
