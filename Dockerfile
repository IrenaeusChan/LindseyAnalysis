FROM debian:bullseye-slim
MAINTAINER "Irenaeus Chan <chani@wustl.edu>"

# Volumes
VOLUME /build

ARG TAG=v0.0.0
ARG PYTHON_VERSION=3.12.6
ARG HTSLIB_VERSION=1.21
ARG BEDTOOLS_VERSION=2.31.1

ARG BOWTIE2_VERSION=2.5.4

# Install build dependencies
RUN apt-get update -qq \
    && apt-get -y install apt-transport-https ca-certificates \
    && apt-get update -qq \
    && apt-get -y install \
    build-essential \
    git \
    less \
    bc \
    vim \
    libnss-sss \
    libcurl4-openssl-dev \
    curl \
    wget \
    unzip \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libffi-dev \
    libncurses5-dev \
    libncursesw5-dev \
    liblzma-dev \
    libtool \
    libtool-bin \
    libltdl7 \
    libltdl-dev \
    autoconf \
    m4 \
    automake \
    pkg-config \
    ca-certificates \
    --no-install-recommends

RUN apt-get clean all

# Install Bowtie2
ENV BOWTIE2_INSTALL_DIR=/opt/bowtie2

WORKDIR /tmp
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    mv bowtie2-${BOWTIE2_VERSION}-linux-x86_64 ${BOWTIE2_INSTALL_DIR} && \
    rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    ln -s ${BOWTIE2_INSTALL_DIR}/bowtie2 /usr/bin/bowtie2 && \
    ln -s ${BOWTIE2_INSTALL_DIR}/bowtie2-build /usr/bin/bowtie2-build && \
    ln -s ${BOWTIE2_INSTALL_DIR}/bowtie2-inspect /usr/bin/bowtie2-inspect

# Download and install HTSlib
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2 && \
    tar --bzip2 -xvf htslib-$HTSLIB_VERSION.tar.bz2 && \
    cd /tmp/htslib-$HTSLIB_VERSION && \
    ./configure --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix && \
    ln -s $HTSLIB_INSTALL_DIR/bin/bgzip /usr/bin/bgzip

# Bcftools #
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools

WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/$HTSLIB_VERSION/bcftools-$HTSLIB_VERSION.tar.bz2 && \
    tar --bzip2 -xvf bcftools-$HTSLIB_VERSION.tar.bz2 && \
    cd /tmp/bcftools-$HTSLIB_VERSION && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$BCFTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/bcftools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/bcftools-$HTSLIB_VERSION

# Samtools #
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/$HTSLIB_VERSION/samtools-$HTSLIB_VERSION.tar.bz2 && \
    tar --bzip2 -xvf samtools-$HTSLIB_VERSION.tar.bz2 && \
    cd /tmp/samtools-$HTSLIB_VERSION && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-$HTSLIB_VERSION

# Download and install Python
RUN wget https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz \
    && tar -xzf Python-$PYTHON_VERSION.tgz \
    && cd Python-$PYTHON_VERSION \
    && ./configure --prefix=/opt/python \
    && make \
    && make install

# Update PATH
ENV PATH="/opt/python:$PATH"
ENV PYTHONPATH="/opt/python/lib/python$PYTHON_VERSION/site-packages:$PYTHONPATH"
ENV PYTHON_EXECUTABLE="/opt/python/bin/python3"

# Create symlink for python3
RUN ln -s /opt/python/bin/python3 /usr/local/bin/python3
RUN ln -s /opt/python/bin/pip3 /usr/local/bin/pip3

# Install pip and the package
RUN /opt/python/bin/pip3 install --upgrade pip
RUN /opt/python/bin/pip3 install pysam numpy biopython pandas intervaltree

# Install cutadapt
RUN /opt/python/bin/pip3 install cutadapt && \
    ln -s /opt/python/bin/cutadapt /usr/local/bin/cutadapt

# Bedtools #
ENV BEDTOOLS_INSTALL_DIR=/opt/bedtools

WORKDIR /tmp
RUN wget https://github.com/arq5x/bedtools2/releases/download/v$BEDTOOLS_VERSION/bedtools-$BEDTOOLS_VERSION.tar.gz && \
    tar -zxvf bedtools-$BEDTOOLS_VERSION.tar.gz && \
    cd bedtools2 && \
    make && \
    cp ./bin/* /usr/local/bin/

# Pandaseq
ENV PANDASEQ_INSTALL_DIR=/opt/pandaseq
RUN git clone http://github.com/neufeld/pandaseq.git/ && \
    cd pandaseq && \
    ./autogen.sh && ./configure && make && make install && \
    ldconfig && \
    cd / && \
    rm -rf pandaseq 
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
