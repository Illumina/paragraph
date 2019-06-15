# Installation of Paragraph

* [System Requirements](#SystemRequirements)
	* [Hardware](#Hardware)
	* [Operating systems](#Operatingsystems)
	* [Third-party libraries](#ThirdPartyLibraries)
* [Static Build](#StaticBuild)
* [Installation](#Installation)
	* [Native build](#NativeBuild)
	* [From Docker image](#FromDockerImage)

## <a name='SystemRequirements'></a>System Requirements

### <a name='Hardware'></a>Hardware

A standard workstation with at least 8GB of RAM should be sufficient for compilation and testing of the program.

### <a name='Operatingsystems'></a>Operating systems

Paragrpah is supported on the following systems:

- Ubuntu 16.04 and CentOS 5-7,
- macOS 10.11+,

Python 3.6+ is required.

We recommend using g++ (6.0+), or a recent version of Clang.

We use the C++11 standard, any Posix compliant compiler supporting this standard
should be usable.

### <a name='ThirdPartyLibraries'></a>Third-party libraries

Please check [requirements](../requirements.txt) for required python modules.

We have included copies of other dependent libraries in external/. They are:
- Google Test and Google Mock (v1.8.0)
- Htslib (v1.9)
- Spdlog

## <a name='Static Build'></a>Static Build

We provide a static build that works for GCC 5.2+ under linux environment. No installation is required for the static build.

Download the static build under "release" tag of the github repo.

## <a name='Installation'></a>Installation

### <a name='NativeBuild'></a>Native buid

[Boost libraries](http://www.boost.org) version >= 1.5 is required.
- We prefer to statically link Boost libraries to Paragraph executables:

  ```bash
  cd ~
  wget http://downloads.sourceforge.net/project/boost/boost/1.65.0/boost_1_65_0.tar.bz2
  tar xf boost_1_65_0.tar.bz2
  cd boost_1_65_0
  ./bootstrap.sh
  ./b2 --prefix=$HOME/boost_1_65_0_install link=static install
  ```

- To point Cmake to your version of Boost use the `BOOST_ROOT` environment variable:

  ```bash
  export BOOST_ROOT=$HOME/boost_1_65_0_install
  ```

Once you have boost installed, checkout the repository like so:

  ```bash
  git clone https://github.com/Illumina/paragraph.git
  cd paragraph-tools
  ```

  Then create a new directory for the program and compile it there:

  ```bash
  # Create a separate build folder.
  cd ..
  mkdir paragraph-tools-build
  cd paragraph-tools-build

  # Configure
  # optional:
  # export BOOST_ROOT=<path-to-boost-installation>
  cmake ../paragraph-tools
  # if this doesn't work, run this instead:
  # cmake ../paragraph-tools -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc` -DBOOST_ROOT=$BOOST_ROOT

  # Make, use -j <n> to use n parallel jobs to build, e.g. make -j4
  make
  ```

### <a name='FromDockerImage'></a>From Docker Image
We also provide a [Dockerfile](Dockerfile). To build a Docker image, run the following command inside the source
  checkout folder:

  ```bash
  docker build .
  ```

  Once the image is built you can find out its ID like this:

  ```bash
  docker images
  ```
  ```
  REPOSITORY                             TAG                 IMAGE ID            CREATED             VIRTUAL SIZE
  <none>                                 <none>              54c7d4015330        16 seconds ago      1.76 GB
  ```
 
  Check the below section for how to run Paragraph, and execute this before running:

  ```bash
  sudo docker run -v `pwd`:/data 54c7d4015330
  ```

  The current directory can be accessed as `/data` inside the Docker container.
  
  The default entry point is `multigrmpy.py`.
  
  To override the default entrypoint and get an interactive shell, run:

  ```bash
  sudo docker run --entrypoint /bin/bash -it 54c7d4015330
  ```