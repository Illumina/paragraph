FROM ubuntu:18.04

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
  libbz2-dev &&  \
  apt-get clean -y && \
  rm -rf /var/lib/apt/lists/*

ADD . /opt/paragraph-source
RUN pip3 install -r /opt/paragraph-source/requirements.txt

RUN mkdir /opt/paragraph-build
WORKDIR /opt/paragraph-build
RUN cmake /opt/paragraph-source -DCMAKE_INSTALL_PREFIX=/opt/paragraph && make && make install
RUN rm -rf /opt/paragraph-source

ENTRYPOINT ["/usr/bin/python3"]
CMD ["/opt/paragraph/bin/multigrmpy.py", "-h"]
