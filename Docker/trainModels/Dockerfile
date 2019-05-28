# Basics
FROM gcc:5.4
FROM continuumio/anaconda

MAINTAINER Andrew Bailey, andbaile@ucsc.edu

# apt-get installs
RUN apt-get update -qq && apt-get install -y git make zlib1g-dev g++ autotools-dev autoconf cmake libhdf5-serial-dev
WORKDIR /home/

# signalAlign
RUN git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign
WORKDIR /home/signalAlign/
RUN conda env create -f requirements.yml python=3.6
RUN echo "source activate signalalign" > ~/.bashrc
ENV PATH /opt/conda/envs/signalalign/bin:$PATH
ENV LIBRARY_PATH /opt/conda/envs/signalalign/lib:$LIBRARY_PATH
ENV LD_LIBRARY_PATH /opt/conda/envs/signalalign/lib:$LD_LIBRARY_PATH
ENV C_INCLUDE_PATH /opt/conda/envs/signalalign/include:$C_INCLUDE_PATH

RUN make
WORKDIR /opt/conda/envs/signalalign/lib
RUN ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
#RUN samtools
WORKDIR /home/signalAlign/
# set PATH variables
ENV ROOT /
ENV PATH "$ROOT/bwa:$PATH"
ENV PATH "$ROOT/bin:$PATH"

RUN make test

# set PATH variables
ENV ROOT /
ENV PATH "$ROOT/bwa:$PATH"
ENV PATH "$ROOT/bin:$PATH"

# set signalAlign bin as workDir
WORKDIR /home/signalAlign/bin/

COPY train_wrapper.sh /opt/signalAlign/

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/opt/signalAlign/train_wrapper.sh"]
