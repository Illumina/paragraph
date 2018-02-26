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

#include "genotyping/GenotypingParameters.hh"
#include "gtest/gtest.h"

using namespace genotyping;

TEST(GenotypingParameters, setFromJson)
{
    const std::vector<std::string> alleles = { "REF", "ALT1", "ALT2" };
    GenotypingParameters param(alleles);

    Json::Value js;
    js["allele_names"] = Json::arrayValue;
    js["allele_names"].append("ALT1");
    js["allele_names"].append("REF");
    js["allele_names"].append("ALT2");
    js["allele_error_rates"] = Json::arrayValue;
    js["allele_error_rates"].append(0.1);
    js["allele_error_rates"].append(0.04);
    js["allele_error_rates"].append(0.1);
    js["het_haplotype_fractions"] = Json::arrayValue;
    js["het_haplotype_fractions"].append(0.33);
    js["het_haplotype_fractions"].append(0.33);
    js["het_haplotype_fractions"].append(0.33);
    param.setFromJson(js);

    EXPECT_EQ((int)param.possibleGenotypes().size(), 6);

    std::vector<double> new_allele_error_rates = { 0.04, 0.1, 0.1 };
    EXPECT_EQ(param.alleleErrorRates(), new_allele_error_rates);
}