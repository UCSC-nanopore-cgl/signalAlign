# syntax = docker/dockerfile:1.0-experimental
FROM baileybase/embed_dependencies:v0.0.1 AS build
MAINTAINER Andrew Bailey, andbaile@ucsc.edu

ENV VIRTUAL_ENV=/opt/venv
RUN python3.7 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# cryptography needed for moto[s3] but new versions need rust installed
RUN python3.7 -m pip install --upgrade pip
RUN python3.7 -m pip -q install wheel
#RUN python3.7 -m pip -q install wheel cryptography==2.3.0 moto[s3]==1.3.14

ENV BUILD_SHARED_LIBS="OFF"
RUN git clone --recursive --single-branch https://github.com/adbailey4/embed_fast5.git && \
    cd /root/embed_fast5 && \
    python3.7 -m pip install .

# install bwa
WORKDIR /root/
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && \
    cp bwa /opt/venv/bin

COPY . /home/signalAlign
RUN mkdir /home/signalAlign/build && \
    cd /home/signalAlign/sonLib && \
    make clean && \
    cd /home/signalAlign/externalTools && \
    make clean && \
    cd /home/signalAlign/build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=. -DCMAKE_VERBOSE_MAKEFILE=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=RELEASE && \
    make -j 4 && \
    cd /home/signalAlign && \
    python3.7 -m pip install . && \
    mkdir /root/signalAlign && mkdir /root/signalAlign/build && \
    cp -r /home/signalAlign/build /root/signalAlign && \
    mkdir /root/signalAlign/bin && \
    cp -r /home/signalAlign/bin /root/signalAlign && \
    cp -r /home/signalAlign/Docker/run_wrapper.sh /root/signalAlign/run_wrapper.sh


FROM ubuntu:18.04 as runtime
RUN apt-get update && \
    apt-get install -y --no-install-recommends python3.7 python3.7-venv libgomp1 libcurl4-openssl-dev samtools hdf5-tools && \
    apt-get clean && \
    apt-get purge


COPY --from=build /root/signalAlign/ /root/signalAlign/
COPY --from=build /opt/venv /opt/venv

ENV VIRTUAL_ENV=/opt/venv
ENV PATH="$VIRTUAL_ENV/bin:/root/signalAlign/bin:$PATH"

RUN mkdir /data
WORKDIR /data

ENTRYPOINT ["sh", "/root/signalAlign/run_wrapper.sh"]

FROM runtime as test
RUN apt-get -qq install -y --no-install-recommends python3-distutils
RUN python3.7 -m pip -q install pytest
COPY ./tests /root/signalAlign/tests
COPY ./models /root/signalAlign/models
COPY ./pytest.ini /root/signalAlign/pytest.ini

COPY ./src/signalalign/tests /root/signalAlign/src/signalalign/tests

RUN cd /root/signalAlign && python3.7 -m pytest
WORKDIR /root/signalAlign
