FROM rocker/rstudio-stable:3.4.0 
# This dockerfile is based on:
#  rocker/rstudio-stable:3.4.0 (based on debian:jessie) 
#  combinelab (available at https://github.com/COMBINE-lab/salmon/blob/master/docker/Dockerfile)

MAINTAINER Jakub Niwa "jakub.niwa@roche.com"

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

### <R installation including required packages> ###
# Set the CRAN mirror
RUN echo "local({\n  r <- getOption(\"repos\")\n\
  r[\"CRAN\"] <- \
  \"https://cran.r-project.org\"\n\
  options(repos = r)\n\
  })\n" >> /etc/R/Rprofile.site

RUN apt-get update
RUN apt-get install -y \
  libxml2-dev libpng-dev openssh-client

# Set MRAN repo to date
RUN sed -i /usr/local/lib/R/etc/Rprofile.site -e 's/2017-04-26/2017-05-23/'
COPY install_libraries.R /code/install_libraries.R

RUN R -f /code/install_libraries.R
### </R installation including required packages> ###

### <chipVis installation including required packages> ###
ENV PACKAGES git gcc make g++ cmake libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf

WORKDIR /home

RUN apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean

CMD ["/init"]


