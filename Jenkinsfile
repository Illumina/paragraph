node ('ussd-prd-rdln01') {
    script {
        env.PROTOBUF_INSTALL_PATH = "/illumina/development/grmpy/protobuf/protobuf-3.3.0-install"
    }

    stage ('Checkout') {
        cleanWs()
        checkout scm
    }

    stage ('Pylint') {
        sh 'cd ${WORKSPACE} && . /bioinfoSD/users/pkrusche/py3/bin/activate && PYTHONPATH=${WORKSPACE}/src/python/lib:${WORKSPACE}/src/python/bin python3 `which pylint` --rcfile=.pylintrc --output-format=parseable --reports=no ${WORKSPACE}/src/python/lib || exit 0'
        sh 'cd ${WORKSPACE} && . /bioinfoSD/users/pkrusche/py3/bin/activate && PYTHONPATH=${WORKSPACE}/src/python/lib:${WORKSPACE}/src/python/bin python3 `which pylint` --rcfile=.pylintrc --output-format=parseable --reports=no ${WORKSPACE}/src/python/test || exit 0'
        sh 'cd ${WORKSPACE} && . /bioinfoSD/users/pkrusche/py3/bin/activate && PYTHONPATH=${WORKSPACE}/src/python/lib:${WORKSPACE}/src/python/bin python3 `which pylint` --rcfile=.pylintrc --output-format=parseable --reports=no ${WORKSPACE}/src/python/bin || exit 0'
    }

    stage ('PEP8') {
        sh 'cd ${WORKSPACE} && . /bioinfoSD/users/pkrusche/py3/bin/activate && PYTHONPATH=${WORKSPACE}/src/python/lib  python3 `which pep8` --ignore=E126 --max-line-length=160 ${WORKSPACE}/src/python || exit 0'
    }

    stage ('Pylint publish') {
        step([$class: 'WarningsPublisher', canComputeNew: false, canResolveRelativePaths: false, consoleParsers: [[parserName: 'PyLint']], defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', unHealthy: '0'])
    }

    stage ('pep8 publish') {
        step([$class: 'WarningsPublisher', canComputeNew: false, canResolveRelativePaths: false, consoleParsers: [[parserName: 'Pep8']], defaultEncoding: '', excludePattern: '', healthy: '', includePattern: '', messagesPattern: '', unHealthy: '0'])
    }

    stage ('Python Tests') {
        sh 'cd ${WORKSPACE} && sh ${WORKSPACE}/src/sh/python-tests.sh'
    }

    stage ('Configure') {
        sh 'cd ${WORKSPACE} && rm -rf build install && mkdir build && cd build && . ${WORKSPACE}/src/sh/illumina-setup.sh && cmake .. -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc` -DBOOST_ROOT=/illumina/development/grmpy/boost/boost_1_65_0_install_gcc6 -DCMAKE_INSTALL_PREFIX=${WORKSPACE}/install'
    }

    stage ('Build') {
        sh 'cd ${WORKSPACE}/build && . ${WORKSPACE}/src/sh/illumina-setup.sh && make -j4 &> ${WORKSPACE}/compile.log'

        step([$class: 'WarningsPublisher',
              canComputeNew: false,
              canResolveRelativePaths: false, defaultEncoding: '', excludePattern: '',
              failedTotalAll: '0', healthy: '0',
              includePattern: '', messagesPattern: '',
              parserConfigurations: [[parserName: 'GNU Make + GNU C Compiler (gcc)', pattern: '**/compile.log']],
              unHealthy: '0', unstableTotalAll: '0'])
    }

    stage ('Install') {
        sh 'cd ${WORKSPACE}/build && . ${WORKSPACE}/src/sh/illumina-setup.sh && make -j4 install'
    }

    stage ('Test') {
        sh 'cd ${WORKSPACE}/install && . ${WORKSPACE}/src/sh/illumina-setup.sh && bin/test_grm && bin/test_blackbox'
    }

    stage ('Python Tests (post-build)') {
        sh 'export GRMPY_INSTALL=${WORKSPACE}/install && cd ${WORKSPACE} && sh ${WORKSPACE}/src/sh/python-tests.sh'
    }

    stage ('Valgrind Check') {
        sh 'cd ${WORKSPACE} && . ${WORKSPACE}/src/sh/valgrind-test.sh'
    }
}
