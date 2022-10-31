#!/bin/bash
set -e -o pipefail

MINICONDA_VERSION="4.12.0"
MINICONDA_URL="https://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh"

# download and run miniconda installer script
curl -sSL $MINICONDA_URL > "/tmp/Miniconda3-${MINICONDA_VERSION}-x86_64.sh"
chmod a+x "/tmp/Miniconda3-${MINICONDA_VERSION}-x86_64.sh"
/tmp/Miniconda3-${MINICONDA_VERSION}-x86_64.sh -b -f -p "$MINICONDA_PATH"
rm /tmp/Miniconda3-${MINICONDA_VERSION}-x86_64.sh

PATH="$MINICONDA_PATH/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda config --add channels r
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set auto_update_conda false
conda install -y mamba -c conda-forge # compatible CLI with faster solver: https://github.com/mamba-org/mamba
conda clean -y --all
