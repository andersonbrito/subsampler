FROM continuumio/miniconda3

LABEL maintainer "Daniel Park <dpark@broadinstitute.org>"

# Set default locale to en_US.UTF-8
ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="C"

# apt packages
RUN apt-get update
RUN apt-get install -y -qq --no-install-recommends \
        lsb-release ca-certificates wget rsync curl python-is-python3 python3-pip \
        less nano vim git locales make \
        dirmngr parallel gnupg file \
        liblz4-tool pigz bzip2 lbzip2 zip unzip zstd

# bring in our files
COPY . /opt/subsampler

# install conda dependencies
WORKDIR /app
RUN conda env create -f /opt/subsampler/config/subsampler.yaml

# Make RUN commands use the new environment:
RUN echo "conda activate subsampler" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# set up entrypoint
WORKDIR /opt/subsampler
CMD ["/bin/bash"]

