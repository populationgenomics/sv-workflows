# use the ubuntu base image
FROM ubuntu:18.04 AS builder

LABEL maintainer="Tobias Rausch rausch@embl.de"

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install specific delly branch
ARG DELLY_VERSION="v0.8.7"
RUN cd /opt \
    && git clone \
      --recurse-submodules \
      --branch $DELLY_VERSION \
      https://github.com/dellytools/delly \
    && cd /opt/delly/ \
    && make STATIC=1 all \
    && make install


# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/delly/bin
WORKDIR /opt/delly/bin
COPY --from=builder /opt/delly/bin/delly .

# Workdir
WORKDIR /root/

# Add Delly to PATH
ENV PATH="/opt/delly/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/sh"]
