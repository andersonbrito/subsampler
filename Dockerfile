FROM ubuntu:jammy-20221020

LABEL maintainer "Daniel Park <dpark@broadinstitute.org>"

COPY . /opt/subsampler

RUN apt-get update

# Set default locale to en_US.UTF-8
ENV LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8" MINICONDA_PATH="/opt/miniconda"

# install conda
RUN /opt/subsampler/scripts/install-miniconda.sh

# install conda dependencies
WORKDIR /opt/subsampler/config
RUN conda env create -f subsampler.yaml
RUN conda activate subsampler

# set up entrypoint
ENV PATH="$MINICONDA_PATH/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
CMD ["/bin/bash"]

