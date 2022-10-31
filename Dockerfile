FROM ubuntu:jammy-20221020

LABEL maintainer "Daniel Park <dpark@broadinstitute.org>"

COPY . /opt/subsampler

# Set default locale to en_US.UTF-8
ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8" MINICONDA_PATH="/opt/miniconda"
ENV PATH="/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# apt packages
RUN apt-get update
RUN apt-get install -y -qq --no-install-recommends \
        lsb-release ca-certificates wget rsync curl python-is-python3 python3-pip \
        less nano vim git locales make \
        dirmngr parallel gnupg file \
        liblz4-tool pigz bzip2 lbzip2 zip unzip zstd

# install conda
RUN /opt/subsampler/scripts/install-miniconda.sh

# install conda dependencies
WORKDIR /opt/subsampler/config
RUN conda env create -f subsampler.yaml
RUN conda activate subsampler

# set up entrypoint
CMD ["/bin/bash"]

