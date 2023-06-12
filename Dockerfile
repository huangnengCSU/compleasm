FROM continuumio/miniconda3
USER root

WORKDIR /app

SHELL [ "/bin/bash", "--login", "-c" ]

RUN conda create -n compleasm -c conda-forge -c bioconda compleasm --yes --all

RUN conda init bash

RUN echo "conda activate compleasm" >> ~/.bashrc

ENV PATH=/opt/conda/envs/compleasm/bin:$PATH

