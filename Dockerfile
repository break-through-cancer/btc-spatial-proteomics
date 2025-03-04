FROM debian:bullseye-slim

ARG DEBIAN_FRONTEND=noninteractive

LABEL maintainer="Break Through Cancer"
LABEL description="Docker image for SPACEc pipeline"
LABEL version="cpu"

# Update the package list and install dependencies
RUN apt-get update && \
    apt-get install -y \
    fontconfig \
    fonts-dejavu-core \
    fonts-droid-fallback \
    fonts-liberation \
    fonts-noto-mono \
    ttf-bitstream-vera \
    git \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Install Miniforge
RUN apt-get update && \
    apt-get install -y wget && \
    wget "https://github.com/conda-forge/miniforge/releases/download/23.11.0-0/Miniforge3-23.11.0-0-Linux-x86_64.sh" && \
    bash Miniforge3-*.sh -b -p /miniforge && \
    rm Miniforge3-*.sh

# Add Miniforge to PATH
ENV PATH="/miniforge/bin:${PATH}"

# Create a new conda environment with Python 3.10
RUN conda install -c conda-forge python=3.10 -y

# Install mamba
RUN conda install -c conda-forge mamba -y

# Install libxml2
RUN mamba install -c conda-forge libxml2=2.13.5 -y

# Install other dependencies
RUN mamba install -c conda-forge graphviz libvips pyvips openslide-python -y

# Install gcc
RUN apt-get update && \
    apt-get install -y build-essential gcc && \
    rm -rf /var/lib/apt/lists/*

# Install SPACEc
RUN mamba run pip install spacec

# Install additional Python packages
RUN mamba run pip install networkx==3.2.* protobuf==3.20.0 numpy==1.24.*

# Clean up
RUN mamba clean --all -f -y && \
    rm -rf /root/.cache/pip
