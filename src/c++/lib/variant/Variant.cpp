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

#include "variant/Variant.hh"

#include <map>

namespace variant
{

struct VariantCandidateList::VariantCandidateListImpl
{
    std::string reference;

    std::map<std::string, std::unique_ptr<Variant>> variants;
    std::map<std::string, PileupData> variant_pileups;
    std::vector<PileupData> reference_pileups;
    std::vector<PileupData> nonreference_pileups;
};

VariantCandidateList::VariantCandidateList()
    : _impl(new VariantCandidateList::VariantCandidateListImpl())
{
}

VariantCandidateList::VariantCandidateList(std::string const& reference)
    : _impl(new VariantCandidateList::VariantCandidateListImpl())
{
    _impl->reference = reference;
    PileupData zero_pileup;
    _impl->reference_pileups.resize(reference.size(), zero_pileup);
    _impl->nonreference_pileups.resize(reference.size(), zero_pileup);
}

VariantCandidateList::VariantCandidateList(VariantCandidateList const& rhs) {}

VariantCandidateList::VariantCandidateList(VariantCandidateList&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

VariantCandidateList::~VariantCandidateList() = default;

VariantCandidateList& VariantCandidateList::operator=(VariantCandidateList const& rhs)
{
    if (this != &rhs)
    {
        _impl->reference = rhs._impl->reference;
        _impl->nonreference_pileups = rhs._impl->nonreference_pileups;
        _impl->reference_pileups = rhs._impl->reference_pileups;
        _impl->variant_pileups = rhs._impl->variant_pileups;
        _impl->variants.clear();
        for (auto const& v : rhs._impl->variants)
        {
            _impl->variants.emplace(v.first, std::unique_ptr<Variant>(new Variant(*v.second)));
        }
    }
    return *this;
}

VariantCandidateList& VariantCandidateList::operator=(VariantCandidateList&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

/**
 * @return the reference sequence passed at construction
 */
std::string const& VariantCandidateList::getReference() const { return _impl->reference; }

/**
 * Add observation for RefVar allele from a read.
 *
 * @param rv the allele
 * @param is_rev true if observation is on reverse strand
 * @param left_boundary position of preceding variant to the left (if any)
 * @param pqual phred-scaled base quality for rv
 * @return the rightmost position this RefVar could be placed
 */
int64_t VariantCandidateList::addRefVarObservation(RefVar rv, bool is_rev, int64_t left_boundary, int pqual)
{
    // update reference information
    if (rv.end >= rv.start && rv.alt == ".")
    {
        for (auto pos = (size_t)rv.start; pos < std::min(_impl->reference.size(), (const unsigned long)rv.end + 1);
             ++pos)
        {
            _impl->reference_pileups[pos].addObs(is_rev, pqual);
        }
        return rv.end;
    }

    int64_t rightmost = std::max(rv.start, rv.end);
    if (rv.alt != ".")
    {
        variant::rightShift(_impl->reference, rv);
        rightmost = std::max(rv.start, rv.end);
        variant::leftShift(_impl->reference, rv, left_boundary);
        variant::trimLeft(_impl->reference, rv, false);

        // add as non-ref obs
        for (auto pos = (size_t)rv.start; pos < std::min(_impl->reference.size(), (size_t)rightmost + 1); ++pos)
        {
            _impl->nonreference_pileups[pos].addObs(is_rev, pqual);
        }

        // add to canonical variant pileup
        const std::string key = rv.repr();
        auto var_it = _impl->variants.find(key);
        if (var_it == _impl->variants.end())
        {
            var_it = _impl->variants.emplace(key, std::unique_ptr<Variant>(new Variant())).first;
            var_it->second->set_start((int)rv.start);
            var_it->second->set_end((int)rv.end);
            var_it->second->set_alt(rv.alt);
            var_it->second->set_leftmost((int)rv.start);
            var_it->second->set_rightmost((int)rightmost);
        }
        auto pileup_it = _impl->variant_pileups.find(key);
        if (pileup_it == _impl->variant_pileups.end())
        {
            pileup_it = _impl->variant_pileups.emplace(key, PileupData()).first;
        }
        pileup_it->second.addObs(is_rev, pqual);
    }
    return rightmost;
}

/**
 * Get piled up depth summary by position.
 */
PileupData const& VariantCandidateList::getRefPileup(int pos) const { return _impl->reference_pileups[pos]; }

/**
 * Get piled up depth summary by position.
 */
PileupData const& VariantCandidateList::getNonrefPileup(int pos) const { return _impl->nonreference_pileups[pos]; }

/**
 * Get variants
 *
 * @return consolidated list of variants + depths
 */
std::list<Variant*> VariantCandidateList::getVariants() const
{
    std::list<Variant*> result;
    for (auto const& v : _impl->variants)
    {
        auto pileup_it = _impl->variant_pileups.find(v.first);
        if (pileup_it != _impl->variant_pileups.end())
        {
            PileupData ref_pile;
            PileupData other_pile;
            int start_pos = v.second->leftmost();
            int end_pos = v.second->rightmost();
            if (end_pos < start_pos)
            {
                // fully trimmed insertion: use reference bases before and after
                std::swap(start_pos, end_pos);
            }
            start_pos = std::max(0, start_pos);
            for (int pos = start_pos; pos <= end_pos; ++pos)
            {
                if (pos >= (int)_impl->reference.size())
                {
                    break;
                }
                ref_pile += _impl->reference_pileups[pos];
                other_pile += _impl->nonreference_pileups[pos];
            }
            const int reflen = end_pos - start_pos + 1;

            if (reflen > 1)
            {
                ref_pile /= reflen;
                other_pile /= reflen;
            }
            other_pile -= pileup_it->second;

            v.second->set_adr_forward(ref_pile.stranded_DP[0]);
            v.second->set_adr_backward(ref_pile.stranded_DP[1]);
            v.second->set_wadr_forward(ref_pile.qual_weighted_DP[0]);
            v.second->set_wadr_backward(ref_pile.qual_weighted_DP[1]);

            v.second->set_ado_forward(other_pile.stranded_DP[0]);
            v.second->set_ado_backward(other_pile.stranded_DP[1]);
            v.second->set_wado_forward(other_pile.qual_weighted_DP[0]);
            v.second->set_wado_backward(other_pile.qual_weighted_DP[1]);

            v.second->set_ada_forward(pileup_it->second.stranded_DP[0]);
            v.second->set_ada_backward(pileup_it->second.stranded_DP[1]);
            v.second->set_wada_forward(pileup_it->second.qual_weighted_DP[0]);
            v.second->set_wada_backward(pileup_it->second.qual_weighted_DP[1]);
        }
        result.push_back(v.second.get());
    }
    return result;
}

/**
 * Append coverage values to JSON value
 *
 * This will append to arrays like this:
 *
 * {
 *    "ref": <ref-matching depth>
 *    "ref:FWD": <...>
 *    "other": non-ref-matching
 *    ...
 * }
 *
 * @param coords coordinate system on the graph to add position information
 * @param node_name which node name to write + use for the Graph coordinates
 * @param coverage coverage JSON. If entries are present already, this function will append
 */
void VariantCandidateList::appendCoverage(
    graphs::GraphCoordinates const& coords, std::string const& node_name, Json::Value& coverage) const
{
    static const std::list<std::string> fields = {
        "cpos",      "node",      "offset", "base",     "ref",      "ref:FWD", "ref:REV",    "other",
        "other:FWD", "other:REV", "wref",   "wref:FWD", "wref:REV", "wother",  "wother:FWD", "wother:REV",
    };
    // make sure sizes are the same + arrays are present
    auto s = static_cast<size_t>(-1);
    for (auto const& f : fields)
    {
        if (coverage.isMember(f))
        {
            if (s == static_cast<size_t>(-1))
            {
                s = coverage.size();
            }
            else
            {
                assert(s == coverage.size());
            }
        }
        else
        {
            coverage[f] = Json::Value(Json::ValueType::arrayValue);
            if (s == static_cast<size_t>(-1))
            {
                s = 0;
            }
            else
            {
                assert(s == 0);
            }
        }
    }

    auto node_start = coords.canonicalPos(node_name);
    for (auto pos = static_cast<size_t>(0); pos < _impl->reference.size(); ++pos)
    {
        auto const& ref_pileup = getRefPileup(static_cast<int>(pos));
        auto const& other_pileup = getNonrefPileup(static_cast<int>(pos));

        coverage["cpos"].append((Json::Value::UInt64)node_start + pos);
        coverage["node"].append(node_name);
        coverage["offset"].append((Json::UInt64)pos);
        coverage["base"].append(std::string(1, (char)_impl->reference[pos]));
        coverage["ref"].append(ref_pileup.stranded_DP[0] + ref_pileup.stranded_DP[1]);
        coverage["ref:FWD"].append(ref_pileup.stranded_DP[0]);
        coverage["ref:REV"].append(ref_pileup.stranded_DP[1]);
        coverage["other"].append(other_pileup.stranded_DP[0] + other_pileup.stranded_DP[1]);
        coverage["other:FWD"].append(other_pileup.stranded_DP[0]);
        coverage["other:REV"].append(other_pileup.stranded_DP[1]);

        coverage["wref"].append(ref_pileup.qual_weighted_DP[0] + ref_pileup.qual_weighted_DP[1]);
        coverage["wref:FWD"].append(ref_pileup.qual_weighted_DP[0]);
        coverage["wref:REV"].append(ref_pileup.qual_weighted_DP[1]);
        coverage["wother"].append(other_pileup.qual_weighted_DP[0] + other_pileup.qual_weighted_DP[1]);
        coverage["wother:FWD"].append(other_pileup.qual_weighted_DP[0]);
        coverage["wother:REV"].append(other_pileup.qual_weighted_DP[1]);
    }
}
}
