FROM ubuntu:17.04

RUN apt-get -qq update && apt-get install -yq \
  autoconf \
  automake \
  build-essential \
  cmake \
  libboost-all-dev \
  libfreetype6-dev \
  liblzma-dev \
  libpng-dev \
  libtool \
  m4 \
  python3-pip \
  python-pkgconfig \
  software-properties-common \
  zlib1g-dev \
  libbz2-dev \
  valgrind \
  cppcheck \
  clang \
  clang-tidy \
  samtools \
  wget curl zip \
  &&  \
  apt-get clean -y && \
  rm -rf /var/lib/apt/lists/*

RUN pip3 install nose mock intervaltree pysam jsonschema pylint pep8
RUN mkdir /opt/paragraph-dependencies
WORKDIR /opt/paragraph-dependencies
RUN wget http://github.com/samtools/htslib/archive/1.5.tar.gz -O htslib.tar.gz
RUN tar xf htslib.tar.gz && cd htslib-1.5 && \
    make DESTDIR=/opt/paragraph-dependencies/htslib-install prefix='' install && \
    rm -rf /opt/paragraph-dependencies/htslib.tar.gz /opt/paragraph-dependencies/htslib-1.5
RUN wget http://github.com/google/protobuf/archive/v3.3.0.tar.gz -O protobuf.tar.gz
RUN tar xf protobuf.tar.gz && cd protobuf-3.3.0 && ./autogen.sh && \
    PKG_CONFIG_PATH=/opt/paragraph-dependencies/protobuf-install/lib/pkgconfig ./configure --prefix=/opt/paragraph-dependencies/protobuf-install && \
    make && \
    make install && \
    echo /opt/paragraph-dependencies/protobuf-install/lib > /etc/ld.so.conf.d/protobuf.conf && \
    rm -rf /opt/paragraph-dependencies/protobuf.tar.gz /opt/paragraph-dependencies/protobuf-3.3.0

RUN mkdir /opt/paragraph-dependencies/data
WORKDIR /opt/paragraph-dependencies/data
ADD make_references.sh /opt/paragraph-dependencies/data/
RUN sh make_references.sh

WORKDIR /

ENV PROTOBUF_INSTALL_PATH /opt/paragraph-dependencies/protobuf-install
ENV HTSLIB_INSTALL_PATH /opt/paragraph-dependencies/htslib-install
ENV HG19 /opt/paragraph-dependencies/data/hg19.fa
ENV HG38 /opt/paragraph-dependencies/data/hg38.fa
RUN samtools faidx ${HG19}
RUN samtools faidx ${HG38}
