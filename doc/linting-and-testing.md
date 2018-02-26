# Linting and testing for developers

Paragraph contains scripts to help assessing code quality and
to run our integration tests. All testing and linting can be 
carried out inside a [Docker image](Dockerfile). A build of this image is 
available at [Docker Hub]( 
https://cloud.docker.com/swarm/ilmncgrpmi/repository/docker/ilmncgrpmi/paragraph-testing-dev/general) as 
`ilmncgrpmi/paragraph-testing-dev:master` and
can be obtained as follows:

```bash
# (run with sudo if necessary)
docker pull ilmncgrpmi/paragraph-testing-dev:master
```

Using this image we automate the following test and 
code analysis tasks:

<!-- vscode-markdown-toc -->
* [Python Linting](#PythonLinting)
* [C++ Static analysis using Cppcheck](#CStaticanalysisusingCppcheck)
* [C++ Builds](#CBuilds)
* [Unit, Blackbox, Python tests](#UnitBlackboxPythontests)
* [Valgrind Dynamic Code Analysis](#ValgrindDynamicCodeAnalysis)
* [Clang Formatting and Static Analysis](#ClangFormattingandStaticAnalysis)

<!-- vscode-markdown-toc-config
	numbering=false
	autoSave=true
	/vscode-markdown-toc-config -->
<!-- /vscode-markdown-toc -->

## <a name='PythonLinting'></a>Python Linting

We use [pylint](https://www.pylint.org/) and [pep8](https://pypi.python.org/pypi/pep8).

Assuming we are in a paragraph checkout directory, 
we can run these checks as follows:

```bash
./src/docker-testing/python-linting.sh `pwd`
```
```
master: Pulling from ilmncgrpmi/paragraph-testing-dev
Digest: ...
Status: Image is up to date for ilmncgrpmi/paragraph-testing-dev:master

------------------------------------
Your code has been rated at 10.00/10

...
```

The return code is zero when no errors or warnings were found.

## <a name='CStaticanalysisusingCppcheck'></a>C++ Static analysis using Cppcheck

We can run static code analysis with [cppcheck](http://cppcheck.sourceforge.net/)
using the [cppcheck.sh](../src/docker-testing/cppcheck.sh) script.

**Warning**: this command will overwrite / erase the file `cppcheck.xml` and the
directory `cppcheck-result` in the current working directory.

```bash
src/docker-testing/cppcheck.sh `pwd`
master: Pulling from ilmncgrpmi/paragraph-testing-dev
Digest: ...
Status: Image is up to date for ilmncgrpmi/paragraph-testing-dev:master
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp ...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: CYGWIN...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: CYGWIN;_DEBUG...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: HTS_ENDIAN_NEUTRAL...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: _DEBUG...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: _DEBUG;__GNUC__...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: _MIPSEB;__AAARCHEB__;__ARMEB__;__BIG_ENDIAN__;__BYTE_ORDER__;__MIPSEB;__MIPSEB__;__THU
MBEB__...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: _MIPSEL;__AARCH64EL__;__ARMEL__;__BYTE_ORDER__;__LITTLE_ENDIAN__;__MIPSEL;__MIPSEL__;_
_THUMBEL__...
Checking /opt/paragraph-source/src/c++/lib/common/Alignment.cpp: __GNUC__...
...
```

The results can be viewed by opening index.html in the cppcheck-result directory. If no
issues were found, no such file will be created.

## <a name='CBuilds'></a>C++ Builds

We can create a build of paragraph which can be run inside our Docker image
(or compatible systems running Ubuntu) using the [make_gcc_build.sh](../src/docker-testing/make_gcc_build.sh) and [make_clang_build.sh](../src/docker-testing/make_clang_build.sh) scripts. These scripts produce
build folders and release tarballs which enable further code analysis
and integration testing.

The two scripts can be used in the same way, here is an example using gcc:

```bash
src/docker-testing/make_gcc_build.sh `pwd` test
```
This should produce two files: `test-gcc-buildfolder.tar.gz` and `test-gcc-release.tar.gz`.

Note that our builds are run with `-Werror`, so code changes that introduce
warnings will make the build fail.

## <a name='UnitBlackboxPythontests'></a>Unit, Blackbox, Python tests

When we have a paragraph build, we can run tests. We have 
the following sets of tests available:

* C++ Unit tests (test_grm executable) - these test the basic functionality
  of our C++ APIs. All test data is included in with the Paragraph installation, 
  and no large / external reference files are required.
* C++ Blackbox tests (test_blackbox executable) - these are C++-level integration
  tests. Some test data is included in with the Paragraph installation, 
  but some reference FASTA files are required which are discovered using 
  the `HG19` and `HG38` environment variables inside the Docker container.
* Python integration tests -- these are workflow-level integration tests. 
  Some reference FASTA files are required which are discovered using 
  the `HG19` and `HG38` environment variables inside the Docker container.

All tests are run using the [testing.sh](../src/docker-testing/testing.sh) script.

```bash
src/docker-testing/testing.sh `pwd` test-gcc-release.tar.gz
```
```
master: Pulling from ilmncgrpmi/paragraph-testing-dev
Digest: sha256:8853fdcb6b6d4dc5bd215bc60df1a2c9b528a1ebf5ab496c58e4a2334d6d8b97
Status: Image is up to date for ilmncgrpmi/paragraph-testing-dev:master
[==========] Running 121 tests from 35 test cases.
[----------] Global test environment set-up.
[----------] 1 test from Alignment
[ RUN      ] Alignment.KlibBasic
Testing klib alignment.
...
```

Any test that fails will cause a non-zero return code.

## <a name='ValgrindDynamicCodeAnalysis'></a>Valgrind Dynamic Code Analysis

We run a subset of our unit and black box tests through valgrind. We exclude some
of the more time-consuming tests for older code which has not been updated for a 
while.

All tests are run using the [valgrind-testing.sh](../src/docker-testing/valgrind-testing.sh) 
script, which has a non-zero return code when one of the tests fails. To work out which
errors occurred, you may inspect the files `valgrind*.xml` which will be written into the 
current working directory (**Warning** these files will be overwritten if they already exist).

```bash
src/docker-testing/valgrind-testing.sh `pwd` test-gcc-release.tar.gz
```
```
master: Pulling from ilmncgrpmi/paragraph-testing-dev
Digest: sha256:8853fdcb6b6d4dc5bd215bc60df1a2c9b528a1ebf5ab496c58e4a2334d6d8b97
Status: Image is up to date for ilmncgrpmi/paragraph-testing-dev:master
...
```

## <a name='ClangFormattingandStaticAnalysis'></a>Clang Formatting and Static Analysis

We can use clang-format to achieve consistent formatting of our code, 
and clang-tidy for warnings.

To check for formatting issues, you may run the [clang-formatting.sh](../src/docker-testing/clang-formatting.sh) shell script:

```bash
src/docker-testing/clang-formatting.sh `pwd`
clang-master: Pulling from ilmncgrpmi/paragraph-testing-dev
Digest: sha256:3d9416a7d8f535effda0dfa3a4887039e2b7a9cbef0df11c31aec1aa594b853b
Status: Image is up to date for ilmncgrpmi/paragraph-testing-dev:clang-master
No formatting issues found.
```

Passing all clang-tidy warnings is not yet a requirement for our code, but 
will be in the near future. In order to use clang-tidy we require a clang
build folder that matches the current source tree.

```bash
src/docker-testing/make_clang_build.sh `pwd` build
src/docker-testing/clang-linting.sh `pwd` build-clang-buildfolder.tar.gz
```
```
clang-master: Pulling from ilmncgrpmi/paragraph-testing-dev
Digest: sha256:3d9416a7d8f535effda0dfa3a4887039e2b7a9cbef0df11c31aec1aa594b853b
Status: Image is up to date for ilmncgrpmi/paragraph-testing-dev:clang-master
...
```