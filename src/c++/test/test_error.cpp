// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Paragraph
// Copyright (c) 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// You may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//		http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
// See the License for the specific language governing permissions and limitations
//
//

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
