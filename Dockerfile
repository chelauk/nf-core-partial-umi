FROM nfcore/base:1.9
LABEL authors="Chela James" \
      description="Docker image containing all software requirements for the umi_preprocessing"

# Install the conda environment
RUN conda install -c conda-forge mamba
COPY environment.yml /
RUN mamba env create -f /environment.yml && mamba clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-umipreprocessing-1.0.0/bin:$PATH

# These are needed for R packages
RUN apt-get update && apt-get install build-essential -y

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-umipreprocessing-1.0.0 > nf-core-umipreprocessing-1.0.0.yml
