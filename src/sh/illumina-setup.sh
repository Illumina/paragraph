#!/bin/bash
#
# Load modules to support compiling
#
# Author: Peter Krusche <pkrusche@illumina.com>

if [[ -z "$(module --version 2>&1 | grep Lua)" ]]; then
    . /etc/profile.d/modules.sh
    module purge
    unset MODULEPATH
    module use /illumina/sync/software/thirdparty/HPCBIOS/modules/all
    module use /illumina/development/grmpy/easybuild_centos65/modules/all
    module load CMake/3.2.1-GCC-4.9.2 &> /dev/null
    module load XZ
    module unload GCC
    module load GCC/5.2.0
    module load libtool &> /dev/null
    module load Automake &> /dev/null
    module load  /illumina/sync/software/thirdparty/HPCBIOS.20150401/modules/all/Valgrind/3.8.1-goolf-1.4.10 &> /dev/null

    echo "Note: older versions of Centos aren't supported for building anymore, please switch to Centos 6.8 or higher."
else
    module purge
    module unload use.own
    module unload use.own.eb
    module load HPCBIOS/2017q2-el7 &> /dev/null
    module load Valgrind
    module load CMake/3.7.2-GCCcore-6.3.0
    module load Python/3.6.1-intel-2017a
    module load libtool
    module load Automake

    # add virtualenv modules
    VIRTUAL_ENV="/illumina/sync/software/groups/graphtools/py3"
    export VIRTUAL_ENV
    PATH="$VIRTUAL_ENV/bin:$PATH"
    export PATH

    # unset PYTHONHOME if set
    if ! [ -z "${PYTHONHOME+_}" ] ; then
        _OLD_VIRTUAL_PYTHONHOME="$PYTHONHOME"
        unset PYTHONHOME
    fi

    if [ -z "${VIRTUAL_ENV_DISABLE_PROMPT-}" ] ; then
        _OLD_VIRTUAL_PS1="$PS1"
        if [ "x" != x ] ; then
            PS1="$PS1"
        else
            PS1="(paragraph-dev) $PS1"
        fi
        export PS1
    fi

    # Make sure to unalias pydoc if it's already there
    alias pydoc 2>/dev/null >/dev/null && unalias pydoc

    pydoc () {
        python -m pydoc "$@"
    }
fi

