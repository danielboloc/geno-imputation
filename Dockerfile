FROM nfcore/base:2.1

LABEL authors="Daniel Boloc" \
       description="Docker image containing all software requirements, without Michigan Imputation Server"

# Install build tools
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    git \
    pkg-config \
    wget \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev p7zip-full && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -xvf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && rm -rf bcftools-1.9*

# Install Eagle
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    gunzip Eagle_v2.4.1.tar.gz && \
    tar xvf Eagle_v2.4.1.tar && \
    mv Eagle_v2.4.1/eagle /usr/local/bin/ && \
    rm -rf Eagle_v2.4.1

RUN useradd --create-home --shell /bin/bash ubuntu && \
    chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]