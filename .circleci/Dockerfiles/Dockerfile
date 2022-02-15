# Taken from (and slightly modified) https://gist.github.com/pangyuteng/f5b00fe63ac31a27be00c56996197597
#
#   ref https://github.com/tebeka/pythonwise/blob/master/docker-miniconda/Dockerfile
#
#   miniconda vers: http://repo.continuum.io/miniconda
#   sample variations:
#     Miniconda3-latest-Linux-armv7l.sh
#     Miniconda3-latest-Linux-x86_64.sh
#     Miniconda3-py38_4.10.3-Linux-x86_64.sh
#     Miniconda3-py37_4.10.3-Linux-x86_64.sh
#
#

ARG UBUNTU_VER=18.04
ARG CONDA_VER=latest
ARG OS_TYPE=x86_64

FROM ubuntu:${UBUNTU_VER}

# System packages 
RUN apt-get update && apt-get install -yq curl wget jq vim gcc g++ make cmake git ssh bash tar gzip ca-certificates

# Use the above args during building https://docs.docker.com/engine/reference/builder/#understand-how-arg-and-from-interact
ARG CONDA_VER
ARG OS_TYPE
# Install miniconda to /miniconda
RUN curl -LO "http://repo.continuum.io/miniconda/Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh"
RUN bash Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh -p /miniconda -b
RUN rm Miniconda3-${CONDA_VER}-Linux-${OS_TYPE}.sh


# some devs like to let it hang https://stackoverflow.com/a/42873832/868736
# ENTRYPOINT ["tail", "-f", "/dev/null"]
