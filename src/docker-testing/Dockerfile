FROM ubuntu:17.10

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

RUN mkdir -p /opt/paragraph-dependencies/data
WORKDIR /opt/paragraph-dependencies/data
ADD make_references.sh /opt/paragraph-dependencies/data/
RUN sh make_references.sh

WORKDIR /

ENV HG19 /opt/paragraph-dependencies/data/hg19.fa
ENV HG38 /opt/paragraph-dependencies/data/hg38.fa
RUN samtools faidx ${HG19}
RUN samtools faidx ${HG38}
