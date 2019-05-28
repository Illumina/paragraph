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