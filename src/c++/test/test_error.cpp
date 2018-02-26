// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2017 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 *  \brief Test error handling function
 *
 * \file test_error.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/Error.hh"
#include "gtest/gtest.h"

#include <sstream>

struct cerr_redirect
{
    explicit cerr_redirect(std::streambuf* new_buffer)
        : old(std::cerr.rdbuf(new_buffer))
    {
        tmp[0] = 0;
        strcpy(tmp, "/tmp/cerrXXXXXXX");
        fd = mkstemp(tmp);

        fflush(stderr);
        fgetpos(stderr, &pos);
        fd = dup(fileno(stderr));
        assert(freopen(tmp, "w", stderr) != nullptr);
    }

    ~cerr_redirect()
    {
        fflush(stderr);
        dup2(fd, fileno(stderr));
        close(fd);

        FILE* new_stderr = fopen(tmp, "r");

        // capture stderr and redirect
        char buf[1025];
        size_t len = 0;
        while ((len = fread(buf, 1, 1024, new_stderr)) != 0)
        {
            buf[len] = 0;
            std::cerr << buf;
        }
        fclose(new_stderr);
        unlink(tmp);
        std::cerr.flush();
        std::cerr.rdbuf(old);
    }

private:
    std::streambuf* old;
    int fd;
    fpos_t pos;
    char tmp[L_tmpnam];
};

TEST(Error, Message)
{
    std::ostringstream oss;
    {
        volatile cerr_redirect r(oss.rdbuf());
        ASSERT_THROW({ error("Test error (ignore me and my traceback)!"); }, std::runtime_error);
    }
#ifdef _DEBUG
    // the error macro writes this to stderr in debug mode only
    ASSERT_TRUE(strstr(oss.str().c_str(), "Test error") != NULL);
#endif
}
