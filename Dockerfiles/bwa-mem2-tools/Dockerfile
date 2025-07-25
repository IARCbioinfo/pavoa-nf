FROM ubuntu:24.04

# Éviter les prompts interactifs lors de l'installation
ENV DEBIAN_FRONTEND=noninteractive

# Mise à jour du système et installation des dépendances
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    python3 \
    python3-pip \
    unzip \
    autoconf \
    automake \
    libtool \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Créer un répertoire de travail
WORKDIR /opt

# Installation de bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.3/bwa-mem2-2.3_x64-linux.tar.bz2 \
  | tar jxf - && \
  mv bwa-mem2-2.3_x64-linux/* /usr/local/bin/ && \
  rm -rf bwa-mem2-2.3_x64-linux

#RUN git clone --recursive https://github.com/bwa-mem2/bwa-mem2.git && \
#    cd bwa-mem2 && \
#    make && \
#    cp bwa-mem2* /usr/local/bin/ && \
#    cd .. && rm -rf bwa-mem2

# Installation de sambamba
RUN wget https://github.com/biod/sambamba/releases/download/v1.0.1/sambamba-1.0.1-linux-amd64-static.gz && \
    gunzip sambamba-1.0.1-linux-amd64-static.gz && \
    mv sambamba-1.0.1-linux-amd64-static /usr/local/bin/sambamba && \
    chmod +x /usr/local/bin/sambamba

# Installation de samblaster
RUN git clone https://github.com/GregoryFaust/samblaster.git && \
    cd samblaster && \
    make && \
    cp samblaster /usr/local/bin/ && \
    cd .. && rm -rf samblaster

# Installation de k8
RUN wget https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 && \
    tar -jxf k8-1.2.tar.bz2 k8-1.2/k8-x86_64-Linux && \
    cp k8-1.2/k8-x86_64-Linux /usr/local/bin/k8 && rm -rf k8-1.2

RUN wget https://github.com/lh3/bwa/raw/master/bwakit/bwa-postalt.js

# Créer un répertoire de travail pour les données
WORKDIR /data

# Point d'entrée par défaut
CMD ["/bin/bash"]
