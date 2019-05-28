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
 *  \brief Implementation of BCF file format helpers
 *
 * \file BCFHelpers.cpp
 * \author Peter Krusche
 * \email pkrusche@illumina.com
 *
 */

#include "common/BCFHelpers.hh"

#include "common/Error.hh"

#include <cstdio>
#include <limits>
#include <set>
#include <sstream>

#include <htslib/vcf.h>
#include <memory>

/**
 * @brief Helper to get out GT fields
 */

namespace bcfhelpers
{
const int MAX_GT = 5;

namespace _impl
{
    /** C++-ified version of bcf_get_format_values */
    template <typename target_type_t> struct bcf_get_numeric_format
    {
        bcf_get_numeric_format() {}

        /* return true when successful
         */
        void operator()(
            const bcf_hdr_t* hdr, bcf1_t* line, const char* tag, int isample, std::vector<target_type_t>& dest) const
        {
            dest.clear();
            int i;
            int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);

            if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id))
            {
                return;
            }

            int nsmpl = bcf_hdr_nsamples(hdr);
            if (isample >= nsmpl)
            {
                return;
            }

            bcf_unpack(line, BCF_UN_FMT);

            for (i = 0; i < line->n_fmt; i++)
            {
                if (line->d.fmt[i].id == tag_id)
                {
                    break;
                }
            }

            if (i == line->n_fmt)
            {
                return;
            }

            bcf_fmt_t* fmt = &line->d.fmt[i];
            int type = fmt->type;
            dest.resize((unsigned long)fmt->n);

            for (i = 0; i < fmt->n; ++i)
            {
                if (type == BCF_BT_FLOAT)
                {
                    static const auto make_missing_float = []() -> float {
                        float f;
                        int32_t _mfloat = 0x7F800001;
                        memcpy(&f, &_mfloat, 4);
                        return f;
                    };
                    static const float bcf_missing_float = make_missing_float();
                    // cppcheck-suppress invalidPointerCast
                    float res = ((float*)(fmt->p + isample * fmt->size))[i];
                    if (res != bcf_missing_float)
                    {
                        dest[i] = target_type_t(res);
                    }
                }
                else if (type == BCF_BT_INT8)
                {
                    int8_t r = ((int8_t*)(fmt->p + isample * fmt->size))[i];
                    if (r == bcf_int8_vector_end)
                    {
                        break;
                    }
                    if (r != bcf_int8_missing)
                    {
                        dest[i] = target_type_t(r);
                    }
                }
                else if (type == BCF_BT_INT16)
                {
                    int16_t r = ((int16_t*)(fmt->p + isample * fmt->size))[i];
                    if (r == bcf_int16_vector_end)
                    {
                        break;
                    }
                    if (r != bcf_int16_missing)
                    {
                        dest[i] = target_type_t(r);
                    }
                }
                else if (type == BCF_BT_INT32)
                {
                    int32_t r = ((int32_t*)(fmt->p + isample * fmt->size))[i];
                    if (r == bcf_int32_vector_end)
                    {
                        break;
                    }
                    if (r != bcf_int32_missing)
                    {
                        dest[i] = target_type_t(r);
                    }
                }
                else
                {
                    auto logger = LOG();
                    logger->warn("String format field ignored when looking for numeric formats!");
                    dest[i] = -1;
                }
            }
        }
    };

    template <typename type_t> struct bcf_get_gts
    {

        bcf_get_gts() {}

        /**
         * @brief Extract GT from Format field
         *
         */
        void operator()(bcf_fmt_t* gt, int i, int* igt, bool& phased) const
        {
            phased = false;
            type_t* p = (type_t*)(gt->p + i * gt->size);
            int ial;
            for (ial = 0; ial < MAX_GT; ial++)
            {
                igt[ial] = -1;
            }
            for (ial = 0; ial < gt->n; ial++)
            {
                if (p[ial] == vector_end)
                {
                    break;
                } /* smaller ploidy */
                if (!(p[ial] >> 1) || p[ial] == missing)
                    continue; /* missing allele */
                int al = (p[ial] >> 1) - 1;
                igt[ial] = al;
                phased = phased || ((p[ial] & 1) != 0);
            }
        }

        static const type_t missing;
        static const type_t vector_end;
    };

    template <> const int8_t bcf_get_gts<int8_t>::missing(bcf_int8_missing);

    template <> const int8_t bcf_get_gts<int8_t>::vector_end(bcf_int8_vector_end);

    template <> const int16_t bcf_get_gts<int16_t>::missing(bcf_int16_missing);

    template <> const int16_t bcf_get_gts<int16_t>::vector_end(bcf_int16_vector_end);

    template <> const int32_t bcf_get_gts<int32_t>::missing(bcf_int32_missing);

    template <> const int32_t bcf_get_gts<int32_t>::vector_end(bcf_int32_vector_end);

    template <typename type_t> struct bcf_get_info
    {
        bcf_get_info() {}

        type_t operator()(bcf_info_t* field, int which = 0) const;
    };

    template <> int bcf_get_info<int>::operator()(bcf_info_t* field, int which) const
    {
        if (field->type != BCF_BT_CHAR && field->len <= which)
        {
            error("Cannot extract int from non-scalar INFO field (len = %i, requested: %i).", field->len, which);
        }
        void* dataptr = nullptr;
        if (which == 0)
        {
            dataptr = &field->v1;
        }
        else
        {
            dataptr = field->vptr;
        }
        assert(dataptr);
        switch (field->type)
        {
        case BCF_BT_NULL:
            return 0;
        case BCF_BT_INT8:
            return (int)*(((int8_t*)dataptr) + which);
        case BCF_BT_INT16:
            return (int)*(((int16_t*)dataptr) + which);
        case BCF_BT_INT32:
            return (int)*(((int32_t*)dataptr) + which);
        case BCF_BT_FLOAT:
            return (int)*(((float*)dataptr) + which);
        case BCF_BT_CHAR:
            if (which > 0)
            {
                error("Cannot extract int %i from string INFO field", which);
            }
            return atoi((const char*)field->vptr);
        default:
            break;
        }
        return -1;
    }

    template <> float bcf_get_info<float>::operator()(bcf_info_t* field, int which) const
    {
        switch (field->type)
        {
        case BCF_BT_NULL:
            return std::numeric_limits<float>::quiet_NaN();
        case BCF_BT_INT8:
        case BCF_BT_INT16:
        case BCF_BT_INT32:
        {
            const bcf_get_info<int> gi;
            return (float)gi(field, which);
        }
        case BCF_BT_FLOAT:
            if (field->len <= which)
            {
                error("Cannot extract int from non-scalar INFO field (len = %i, requested: %i).", field->len, which);
            }
            return *(((float*)field->vptr) + which);
        case BCF_BT_CHAR:
            if (which > 0)
            {
                error("Cannot extract int %i from string INFO field", which);
            }
            return (float)atof((const char*)field->vptr);
        default:
            break;
        }
        return std::numeric_limits<float>::quiet_NaN();
    }

    template <> std::string bcf_get_info<std::string>::operator()(bcf_info_t* field, int) const
    {
        char num[256];
        switch (field->type)
        {
        case BCF_BT_NULL:
            return std::string("NULL");
        case BCF_BT_INT8:
        case BCF_BT_INT16:
        case BCF_BT_INT32:
            snprintf(num, 256, "%i", field->v1.i);
            return std::string(num);
        case BCF_BT_FLOAT:
            snprintf(num, 256, "%g", field->v1.f);
            return std::string(num);
        case BCF_BT_CHAR:
            if (field->vptr && field->len > 0)
            {
                return std::string((const char*)field->vptr, (unsigned long)field->len);
            }
        default:
            break;
        }
        return std::string();
    }
}

void bcfHeaderHG19(bcf_hdr_t* header)
{
    bcf_hdr_append(header, "##reference=hg19");
    bcf_hdr_append(header, "##contig=<ID=chr1,length=249250621>");
    bcf_hdr_append(header, "##contig=<ID=chr2,length=243199373>");
    bcf_hdr_append(header, "##contig=<ID=chr3,length=198022430>");
    bcf_hdr_append(header, "##contig=<ID=chr4,length=191154276>");
    bcf_hdr_append(header, "##contig=<ID=chr5,length=180915260>");
    bcf_hdr_append(header, "##contig=<ID=chr6,length=171115067>");
    bcf_hdr_append(header, "##contig=<ID=chr7,length=159138663>");
    bcf_hdr_append(header, "##contig=<ID=chr8,length=146364022>");
    bcf_hdr_append(header, "##contig=<ID=chr9,length=141213431>");
    bcf_hdr_append(header, "##contig=<ID=chr10,length=135534747>");
    bcf_hdr_append(header, "##contig=<ID=chr11,length=135006516>");
    bcf_hdr_append(header, "##contig=<ID=chr12,length=133851895>");
    bcf_hdr_append(header, "##contig=<ID=chr13,length=115169878>");
    bcf_hdr_append(header, "##contig=<ID=chr14,length=107349540>");
    bcf_hdr_append(header, "##contig=<ID=chr15,length=102531392>");
    bcf_hdr_append(header, "##contig=<ID=chr16,length=90354753>");
    bcf_hdr_append(header, "##contig=<ID=chr17,length=81195210>");
    bcf_hdr_append(header, "##contig=<ID=chr18,length=78077248>");
    bcf_hdr_append(header, "##contig=<ID=chr19,length=59128983>");
    bcf_hdr_append(header, "##contig=<ID=chr20,length=63025520>");
    bcf_hdr_append(header, "##contig=<ID=chr21,length=48129895>");
    bcf_hdr_append(header, "##contig=<ID=chr22,length=51304566>");
    bcf_hdr_append(header, "##contig=<ID=chrX,length=155270560>");
    bcf_hdr_append(header, "##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">");
    bcf_hdr_append(
        header,
        "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that "
        "could not be imported.\">");
    bcf_hdr_append(header, "##FORMAT=<ID=AGT,Number=1,Type=String,Description=\"Genotypes at ambiguous locations\">");
    bcf_hdr_append(header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(header, "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">");
    bcf_hdr_append(header, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
    bcf_hdr_append(header, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">");
    bcf_hdr_append(
        header, "##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">");
}

/** extract chrom / pos / length */
void getLocation(bcf_hdr_t* hdr, bcf1_t* rec, int64_t& refstart, int64_t& refend)
{
    bcf_unpack(rec, BCF_UN_STR);
    refstart = rec->pos;
    refend = refstart;

    int endfield = getInfoInt(hdr, rec, "END", -1);

    if (endfield > 0)
    {
        // if there is an end field, don't validate the ref allele
        refend = endfield - 1;
    }
    else
    {
        refend = refstart + strlen(rec->d.allele[0]) - 1;
        if (strchr(rec->d.allele[0], '.') || strchr(rec->d.allele[0], '-'))
        {
            // length might be inaccurate now
            refend = refstart;
            throw bcfhelpers::importexception(
                std::string("[W] Unsupported REF allele with undefined length: ") + std::string(rec->d.allele[0]));
        }
    }
}

/**
 * @brief Retrieve an info field as an integer
 *
 * @param result the default to return if the field is not present
 */
std::string getInfoString(bcf_hdr_t* header, bcf1_t* line, const char* field, const char* def_result)
{
    std::string result = def_result;
    bcf_info_t* info_ptr = bcf_get_info(header, line, field);
    if (info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<std::string> i;
        result = i(info_ptr);
    }
    return result;
}

/**
 * @brief Retrieve an info field as an integer
 *
 * @param result the default to return if the field is not present
 */
int getInfoInt(bcf_hdr_t* header, bcf1_t* line, const char* field, int result)
{
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t* info_ptr = bcf_get_info(header, line, field);
    if (info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<int> i;
        result = i(info_ptr);
    }
    return result;
}

std::vector<int> getInfoInts(bcf_hdr_t* header, bcf1_t* line, const char* field)
{
    std::vector<int> result;
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t* info_ptr = bcf_get_info(header, line, field);
    if (info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<int> i;
        for (int j = 0; j < info_ptr->len; ++j)
        {
            result.push_back(i(info_ptr, j));
        }
    }
    return result;
}

/**
 * @brief Retrieve an info field as a double
 *
 * @return the value or NaN
 */
float getInfoFloat(bcf_hdr_t* header, bcf1_t* line, const char* field)
{
    float result = std::numeric_limits<float>::quiet_NaN();
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t* info_ptr = bcf_get_info(header, line, field);
    if (info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<float> i;
        result = i(info_ptr);
    }
    return result;
}

std::vector<float> getInfoFloats(bcf_hdr_t* header, bcf1_t* line, const char* field)
{
    std::vector<float> result;
    bcf_unpack(line, BCF_UN_INFO);
    bcf_info_t* info_ptr = bcf_get_info(header, line, field);
    if (info_ptr)
    {
        static const bcfhelpers::_impl::bcf_get_info<float> i;
        for (int j = 0; j < info_ptr->len; ++j)
        {
            result.push_back(i(info_ptr, j));
        }
    }
    return result;
}

/**
 * @brief Retrieve an info flag
 *
 * @return true of the flag is set
 */
bool getInfoFlag(bcf_hdr_t* hdr, bcf1_t* line, const char* field)
{
    bcf_unpack(line, BCF_UN_INFO);
    return bcf_get_info_flag(hdr, line, field, nullptr, 0) == 1;
}

/**
 * @brief Read the GT field
 */
void getGT(bcf_hdr_t* header, bcf1_t* line, int isample, int* gt, int& ngt, bool& phased)
{
    bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t* fmt_ptr = bcf_get_fmt(header, line, "GT");
    if (fmt_ptr)
    {
        ngt = fmt_ptr->n;
        if (ngt > MAX_GT)
        {
            auto logger = LOG();
            logger->warn("Found a variant with {} > {} (max) alt alleles. These become no-calls.", ngt, MAX_GT);

            for (int i = 0; i < MAX_GT; ++i)
            {
                gt[i] = -1;
            }
            phased = false;
            return;
        }

        // check if all alleles are populated
        switch (fmt_ptr->type)
        {
        case BCF_BT_INT8:
        {
            static const bcfhelpers::_impl::bcf_get_gts<int8_t> b;
            b(fmt_ptr, isample, gt, phased);
        }
        break;
        case BCF_BT_INT16:
        {
            static const bcfhelpers::_impl::bcf_get_gts<int16_t> b;
            b(fmt_ptr, isample, gt, phased);
        }
        break;
        case BCF_BT_INT32:
        {
            static const bcfhelpers::_impl::bcf_get_gts<int32_t> b;
            b(fmt_ptr, isample, gt, phased);
        }
        break;
        default:
            error(
                "Unsupported GT type: %d at %s:%d\n", fmt_ptr->type, header->id[BCF_DT_CTG][line->rid].key,
                line->pos + 1);
            break;
        }
    }
    else
    {
        ngt = 0;
        phased = false;
    }
}

/** read GQ(X) -- will use in this order: GQX, GQ, -1
 *
 * This is somewhat Illumina-specific, TODO: refactor to avoid using this function in a general setting
 * */
void getGQ(const bcf_hdr_t* header, bcf1_t* line, int isample, float& gq)
{
    using namespace _impl;
    static const bcf_get_numeric_format<float> gf;

    std::vector<float> values;
    gf(header, line, "GQX", isample, values);
    if (values.empty())
    {
        gf(header, line, "GQ", isample, values);
    }
    if (values.size() > 1)
    {
        auto logger = LOG();
        logger->warn("Too many GQ fields at {}:{}", header->id[BCF_DT_CTG][line->rid].key, line->pos);
    }
    if (values.empty())
    {
        gq = -1;
    }
    else
    {
        gq = values[0];
    }
}

/** read AD */
void getAD(const bcf_hdr_t* header, bcf1_t* line, int isample, int* ad, int max_ad)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;

    std::vector<int> values;
    gf(header, line, "AD", isample, values);
    if (max_ad < (int)values.size())
    {
        auto logger = LOG();
        logger->warn(
            "Too many AD fields at {}:{} max_ad = {} retrieved: {}", header->id[BCF_DT_CTG][line->rid].key, line->pos,
            max_ad, values.size());
    }
    for (size_t q = 0; q < values.size(); ++q)
    {
        ad[q] = values[q];
    }
}

/** read DP(I) -- will use in this order: DP, DPI, -1 */
void getDP(const bcf_hdr_t* header, bcf1_t* line, int isample, int& dp)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;

    std::vector<int> values;
    gf(header, line, "DP", isample, values);
    if (values.empty())
    {
        gf(header, line, "DPI", isample, values);
    }
    if (values.size() > 1)
    {
        auto logger = LOG();
        logger->warn("Too many DP fields at {}:{}", header->id[BCF_DT_CTG][line->rid].key, line->pos);
    }
    if (values.empty())
    {
        dp = 0;
    }
    else
    {
        dp = values[0];
    }
}

/** read a format field as a single int. */
int getFormatInt(const bcf_hdr_t* header, bcf1_t* line, const char* field, int isample, int defaultresult)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;

    std::vector<int> values;
    gf(header, line, field, isample, values);
    if (values.size() > 1)
    {
        std::ostringstream os;
        os << "[W] too many " << field << " fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos;
        throw importexception(os.str());
    }
    if (values.size() == 1)
    {
        return values[0];
    }
    return defaultresult;
}

std::vector<int> getFormatInts(const bcf_hdr_t* header, bcf1_t* line, const char* field, int isample)
{
    using namespace _impl;
    static const bcf_get_numeric_format<int> gf;
    std::vector<int> result;
    gf(header, line, field, isample, result);
    return result;
}

/** read a format field as a single float. default return value is NaN */
float getFormatFloat(const bcf_hdr_t* header, bcf1_t* line, const char* field, int isample)
{
    using namespace _impl;
    float result = std::numeric_limits<float>::quiet_NaN();
    static const bcf_get_numeric_format<float> gf;

    std::vector<float> values;
    gf(header, line, field, isample, values);
    if (values.size() > 1)
    {
        std::ostringstream os;
        os << "[W] too many " << field << " fields at " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos;
        throw importexception(os.str());
    }
    if (values.size() == 1)
    {
        result = values[0];
    }
    return result;
}

std::vector<float> getFormatFloats(const bcf_hdr_t* header, bcf1_t* line, const char* field, int isample)
{
    using namespace _impl;
    static const bcf_get_numeric_format<float> gf;
    std::vector<float> result;
    gf(header, line, field, isample, result);
    return result;
}

/** read a format field as a single double. result will not be overwritten on failure */
std::string getFormatString(const bcf_hdr_t* hdr, bcf1_t* line, const char* field, int isample, const char* result)
{
    int nsmpl = bcf_hdr_nsamples(hdr);
    if (isample >= nsmpl)
    {
        return result;
    }

    bcf_unpack(line, BCF_UN_FMT);
    int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, field);

    if (!bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, tag_id))
    {
        return result;
    }

    if (!(line->unpacked & BCF_UN_FMT))
    {
        bcf_unpack(line, BCF_UN_FMT);
    }

    // index in format fields
    int i = 0;
    for (i = 0; i < line->n_fmt; i++)
    {
        if (line->d.fmt[i].id == tag_id)
        {
            break;
        }
    }

    if (i == line->n_fmt)
    {
        return result;
    }

    bcf_fmt_t* fmt = &line->d.fmt[i];

    if (fmt == NULL)
    {
        return result;
    }

    int type = fmt->type;

    if (fmt->n < 1)
    {
        return result;
    }

    std::string str_result = result;
    if (type == BCF_BT_FLOAT)
    {
        static const float bcf_missing_float = missing_float();
        // cppcheck-suppress invalidPointerCast
        float res = *((float*)(fmt->p + isample * fmt->size));
        if (res != bcf_missing_float)
        {
            str_result = std::to_string(res);
        }
    }
    else if (type == BCF_BT_INT8)
    {
        int8_t r = *((int8_t*)(fmt->p + isample * fmt->size));
        if (r != bcf_int8_missing && r != bcf_int8_vector_end)
        {
            str_result = std::to_string(r);
        }
    }
    else if (type == BCF_BT_INT16)
    {
        int16_t r = *((int16_t*)(fmt->p + isample * fmt->size));
        if (r != bcf_int16_missing && r != bcf_int16_vector_end)
        {
            str_result = std::to_string(r);
        }
    }
    else if (type == BCF_BT_INT32)
    {
        int32_t r = *((int32_t*)(fmt->p + isample * fmt->size));
        if (r != bcf_int32_missing && r == bcf_int32_vector_end)
        {
            str_result = std::to_string(r);
        }
    }
    else
    {
        const char* src = (const char*)fmt->p + isample * fmt->size;
        if (src && fmt->size)
        {
            str_result = std::string(src, (unsigned long)fmt->size);
            // deal with 0 padding
            str_result.resize(strlen(str_result.c_str()));
        }
    }

    return str_result;
}

/** update format string for a single sample.  */
void setFormatStrings(const bcf_hdr_t* hdr, bcf1_t* line, const char* field, const std::vector<std::string>& formats)
{
    // TODO this can probably be done faster / better
    std::unique_ptr<const char* []> p_fmts = std::unique_ptr<const char* []>(new const char*[line->n_sample]);
    bcf_unpack(line, BCF_UN_FMT);
    bool any_nonempty = false;
    if (formats.size() != line->n_sample)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " " << hdr->id[BCF_DT_ID][line->rid].key << ":" << line->pos
           << " -- we have " << formats.size() << " for " << line->n_sample << " samples";
        throw importexception(os.str());
    }
    for (int si = 0; si < line->n_sample; ++si)
    {
        p_fmts.get()[si] = formats[si].c_str();
        if (!formats[si].empty())
        {
            any_nonempty = true;
        }
    }
    int res;
    if (any_nonempty)
    {
        res = bcf_update_format_string(hdr, line, field, p_fmts.get(), line->n_sample);
    }
    else
    {
        res = bcf_update_format_string(hdr, line, field, NULL, 0);
    }

    if (res != 0)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " " << hdr->id[BCF_DT_ID][line->rid].key << ":" << line->pos
           << " -- we have " << formats.size() << " for " << line->n_sample << " samples";
        throw importexception(os.str());
    }
}

/** update format with single float values.  */
void setFormatFloats(const bcf_hdr_t* header, bcf1_t* line, const char* field, const std::vector<float>& value)
{
    if (value.empty())
    {
        int res = bcf_update_format(header, line, field, NULL, 0, 0);
        if (res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":"
               << line->pos << " -- we have " << value.size() << " for " << line->n_sample << " samples";

            throw importexception(os.str());
        }
        return;
    }
    std::unique_ptr<float[]> p_dbl = std::unique_ptr<float[]>(new float[line->n_sample]);

    for (size_t i = 0; i < line->n_sample; ++i)
    {
        if (i < value.size())
        {
            p_dbl.get()[i] = value[i];
        }
        else
        {
            p_dbl.get()[i] = std::numeric_limits<float>::quiet_NaN();
        }
    }

    int res = bcf_update_format_float(header, line, field, p_dbl.get(), line->n_sample);
    if (res != 0)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
           << " -- we have " << value.size() << " for " << line->n_sample << " samples";
        throw importexception(os.str());
    }
}

void setFormatInts(const bcf_hdr_t* header, bcf1_t* line, const char* field, const std::vector<int>& value)
{
    if (value.empty())
    {
        int res = bcf_update_format(header, line, field, NULL, 0, 0);
        if (res != 0)
        {
            std::ostringstream os;
            os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":"
               << line->pos << " -- we have " << value.size() << " for " << line->n_sample << " samples";

            throw importexception(os.str());
        }
        return;
    }
    std::unique_ptr<int[]> p_dbl = std::unique_ptr<int[]>(new int[line->n_sample]);

    for (size_t i = 0; i < line->n_sample; ++i)
    {
        if (i < value.size())
        {
            p_dbl.get()[i] = value[i];
        }
        else
        {
            p_dbl.get()[i] = bcf_int32_missing;
        }
    }

    int res = bcf_update_format_int32(header, line, field, p_dbl.get(), line->n_sample);
    if (res != 0)
    {
        std::ostringstream os;
        os << "[W] cannot update format " << field << " " << header->id[BCF_DT_CTG][line->rid].key << ":" << line->pos
           << " -- we have " << value.size() << " for " << line->n_sample << " samples";
        throw importexception(os.str());
    }
}

/** return sample names from header */
std::list<std::string> getSampleNames(const bcf_hdr_t* hdr)
{
    std::list<std::string> l;
    for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i)
    {
        std::string samplename = hdr->samples[i];
        if (samplename == "*")
        {
            auto logger = LOG();
            logger->warn("Skipping sample named '*'");
            continue;
        }
        l.push_back(samplename);
    }
    return l;
}

/** return number of reference padding bases */
int addRefPad(bcf_hdr_t* hdr, bcf1_t* rec, common::FastaFile const& ref, int npad)
{
    if (npad <= 0)
    {
        error_stop("npad<=0", 0, "npad<=0");
    }
    int64_t start, end;
    getLocation(hdr, rec, start, end);
    std::string pad = ref.query(bcf_hdr_int2id(hdr, BCF_DT_CTG, rec->rid), start - npad, start - 1);
    rec->pos -= npad;
    std::string new_alleles = "";
    for (int i = 0; i < rec->n_allele; i++)
    {
        if (rec->d.allele[i][0] != '<')
        {
            if (i == 0)
            {
                new_alleles += pad;
                new_alleles += (std::string)rec->d.allele[i];
            }
            else
            {
                new_alleles += ",";
                new_alleles += pad;
                new_alleles += (std::string)rec->d.allele[i];
            }
        }
    }
    bcf_update_alleles_str(hdr, rec, new_alleles.c_str());
    return (npad);
}

/** return number of reference padding bases */
int isRefPadded(bcf1_t* line)
{
    bcf_unpack(line, BCF_UN_SHR);
    if (line->n_allele == 1)
    {
        return 0;
    }

    const char* ref = line->d.allele[0];
    const int reflen = (int)strlen(ref);

    int max_match = reflen;
    for (int al = 1; al < line->n_allele; ++al)
    {
        const char* alt = line->d.allele[al];
        // symbolic or missing ALT
        if (strcmp(alt, ".") == 0 || *alt == '<')
        {
            max_match = 0;
            break;
        }
        int rpos = 0;
        for (rpos = 0; rpos < reflen; ++rpos, ++alt)
        {
            const char rb = *(ref + rpos);
            if (*alt == 0 || *alt != rb)
            {
                break;
            }
        }
        max_match = std::min(rpos, max_match);
    }
    return max_match;
}

// everything after here is modified from bcftools 1.3.1. vcfnorm.c

static void split_info_numeric(
    const bcf_hdr_t* hdr, bcf1_t* src, bcf_info_t* info, int ialt, bcf1_t* dst, uint8_t*& tmp_arr1, int& ntmp_arr1)
{
#define BRANCH_NUMERIC(type, type_t)                                                                                   \
    {                                                                                                                  \
        const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);                                                   \
        int ntmp = ntmp_arr1 / sizeof(type_t);                                                                         \
        int ret = bcf_get_info_##type(hdr, src, tag, &tmp_arr1, &ntmp);                                                \
        ntmp_arr1 = ntmp * sizeof(type_t);                                                                             \
        assert(ret > 0);                                                                                               \
        type_t* vals = (type_t*)tmp_arr1;                                                                              \
        int len = bcf_hdr_id2length(hdr, BCF_HL_INFO, info->key);                                                      \
        if (len == BCF_VL_A)                                                                                           \
        {                                                                                                              \
            assert(ret == src->n_allele - 1);                                                                          \
            bcf_update_info_##type(hdr, dst, tag, vals + ialt, 1);                                                     \
        }                                                                                                              \
        else if (len == BCF_VL_R)                                                                                      \
        {                                                                                                              \
            assert(ret == src->n_allele);                                                                              \
            if (ialt != 0)                                                                                             \
                vals[1] = vals[ialt + 1];                                                                              \
            bcf_update_info_##type(hdr, dst, tag, vals, 2);                                                            \
        }                                                                                                              \
        else if (len == BCF_VL_G)                                                                                      \
        {                                                                                                              \
            assert(ret == src->n_allele * (src->n_allele + 1) / 2);                                                    \
            if (ialt != 0)                                                                                             \
            {                                                                                                          \
                vals[1] = vals[bcf_alleles2gt(0, ialt + 1)];                                                           \
                vals[2] = vals[bcf_alleles2gt(ialt + 1, ialt + 1)];                                                    \
            }                                                                                                          \
            bcf_update_info_##type(hdr, dst, tag, vals, 3);                                                            \
        }                                                                                                              \
        else                                                                                                           \
            bcf_update_info_##type(hdr, dst, tag, vals, ret);                                                          \
    }
    switch (bcf_hdr_id2type(hdr, BCF_HL_INFO, info->key))
    {
    case BCF_HT_INT:
        BRANCH_NUMERIC(int32, int32_t);
        break;
    case BCF_HT_REAL:
    {
        // cppcheck-suppress invalidPointerCast
        BRANCH_NUMERIC(float, float);
        break;
    }
    }
#undef BRANCH_NUMERIC
}
// Find n-th field in a comma-separated list and move it to dst.
// The memory areas may overlap.
#define STR_MOVE_NTH(dst, src, end, nth, len)                                                                          \
    {                                                                                                                  \
        char *ss = src, *se = src;                                                                                     \
        int j = 0;                                                                                                     \
        while (*se && se < (end))                                                                                      \
        {                                                                                                              \
            if (*se == ',')                                                                                            \
            {                                                                                                          \
                if (j == nth)                                                                                          \
                    break;                                                                                             \
                j++;                                                                                                   \
                ss = se + 1;                                                                                           \
            }                                                                                                          \
            se++;                                                                                                      \
        }                                                                                                              \
        if (j == nth)                                                                                                  \
        {                                                                                                              \
            int n = se - ss;                                                                                           \
            memmove((dst), ss, n);                                                                                     \
            src = se;                                                                                                  \
            len += n;                                                                                                  \
        }                                                                                                              \
        else                                                                                                           \
            len = -1;                                                                                                  \
    }

static void split_info_string(
    const bcf_hdr_t* hdr, bcf1_t* src, bcf_info_t* info, int ialt, bcf1_t* dst, uint8_t*& tmp_arr1, int& ntmp_arr1)
{
    const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);
    int ret = bcf_get_info_string(hdr, src, tag, &tmp_arr1, &ntmp_arr1);
    assert(ret > 0);

    kstring_t str;
    str.m = (size_t)ntmp_arr1;
    str.l = (size_t)ret;
    str.s = (char*)tmp_arr1;

    int len = bcf_hdr_id2length(hdr, BCF_HL_INFO, info->key);
    if (len == BCF_VL_A)
    {
        char* tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s, tmp, str.s + str.l, ialt, len);
        if (len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(hdr, dst, tag, str.s);
    }
    else if (len == BCF_VL_R)
    {
        char* tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s, tmp, str.s + str.l, 0, len);
        str.s[len] = ',';
        tmp++;
        len++;
        STR_MOVE_NTH(&str.s[len], tmp, str.s + str.l, ialt, len);
        if (len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(hdr, dst, tag, str.s);
    }
    else if (len == BCF_VL_G)
    {
        // cppcheck-suppress duplicateExpression
        int i0a = bcf_alleles2gt(0, ialt + 1), iaa = bcf_alleles2gt(ialt + 1, ialt + 1);
        char* tmp = str.s;
        int len = 0;
        STR_MOVE_NTH(str.s, tmp, str.s + str.l, 0, len);
        str.s[len] = ',';
        tmp++;
        len++;
        STR_MOVE_NTH(&str.s[len], tmp, str.s + str.l, i0a - 1, len);
        if (len < 0)
        {
            return;
        } // wrong number of fields: skip
        str.s[len] = ',';
        tmp++;
        len++;
        STR_MOVE_NTH(&str.s[len], tmp, str.s + str.l, iaa - i0a - 1, len);
        if (len < 0)
            return; // wrong number of fields: skip
        str.s[len] = 0;
        bcf_update_info_string(hdr, dst, tag, str.s);
    }
    else
        bcf_update_info_string(hdr, dst, tag, str.s);
}

static void split_info_flag(
    const bcf_hdr_t* hdr, bcf1_t* src, bcf_info_t* info, int ialt, bcf1_t* dst, uint8_t*& tmp_arr1, int& ntmp_arr1)
{
    const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, info->key);
    int ret = bcf_get_info_flag(hdr, src, tag, &tmp_arr1, &ntmp_arr1);
    bcf_update_info_flag(hdr, dst, tag, NULL, ret);
}

static void split_format_genotype(
    const bcf_hdr_t* hdr, bcf1_t* src, bcf_fmt_t* fmt, int ialt, bcf1_t* dst, uint8_t*& tmp_arr1, int& ntmp_arr1)
{
    int ntmp = ntmp_arr1 / 4;
    int ngts = bcf_get_genotypes(hdr, src, &tmp_arr1, &ntmp);
    ntmp_arr1 = ntmp * 4;
    assert(ngts > 0);

    int32_t* gt = (int32_t*)tmp_arr1;
    int i, j, nsmpl = bcf_hdr_nsamples(hdr);
    ngts /= nsmpl;
    for (i = 0; i < nsmpl; i++)
    {
        for (j = 0; j < ngts; j++)
        {
            if (gt[j] == bcf_int32_vector_end)
            {
                break;
            }
            // cppcheck-suppress clarifyCalculation
            if (bcf_gt_is_missing(gt[j]) || bcf_gt_allele(gt[j]) == 0)
            {
                continue; // missing allele or ref: leave as is
            }
            if (bcf_gt_allele(gt[j]) == ialt + 1)
            {
                gt[j] = bcf_gt_unphased(1) | bcf_gt_is_phased(gt[j]); // set to first ALT
            }
            else
            {
                gt[j] = bcf_gt_unphased(0) | bcf_gt_is_phased(gt[j]); // set to REF
            }
        }
        gt += ngts;
    }
    bcf_update_genotypes(hdr, dst, tmp_arr1, ngts * nsmpl);
}

static void split_format_numeric(
    const bcf_hdr_t* hdr, bcf1_t* src, bcf_fmt_t* fmt, int ialt, bcf1_t* dst, uint8_t*& tmp_arr1, int& ntmp_arr1)
{
#define BRANCH_NUMERIC(type, type_t, is_vector_end, set_vector_end)                                                    \
    {                                                                                                                  \
        const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);                                                     \
        int ntmp = ntmp_arr1 / sizeof(type_t);                                                                         \
        int nvals = bcf_get_format_##type(hdr, src, tag, &tmp_arr1, &ntmp);                                            \
        ntmp_arr1 = ntmp * sizeof(type_t);                                                                             \
        assert(nvals > 0);                                                                                             \
        type_t* vals = (type_t*)tmp_arr1;                                                                              \
        int len = bcf_hdr_id2length(hdr, BCF_HL_FMT, fmt->id);                                                         \
        int i, nsmpl = bcf_hdr_nsamples(hdr);                                                                          \
        if (nvals == nsmpl) /* all values are missing */                                                               \
        {                                                                                                              \
            bcf_update_format_##type(hdr, dst, tag, vals, nsmpl);                                                      \
            return;                                                                                                    \
        }                                                                                                              \
        if (len == BCF_VL_A)                                                                                           \
        {                                                                                                              \
            assert(nvals == (src->n_allele - 1) * nsmpl);                                                              \
            nvals /= nsmpl;                                                                                            \
            type_t *src_vals = vals, *dst_vals = vals;                                                                 \
            for (i = 0; i < nsmpl; i++)                                                                                \
            {                                                                                                          \
                dst_vals[0] = src_vals[ialt];                                                                          \
                dst_vals += 1;                                                                                         \
                src_vals += nvals;                                                                                     \
            }                                                                                                          \
            bcf_update_format_##type(hdr, dst, tag, vals, nsmpl);                                                      \
        }                                                                                                              \
        else if (len == BCF_VL_R)                                                                                      \
        {                                                                                                              \
            assert(nvals == src->n_allele * nsmpl);                                                                    \
            nvals /= nsmpl;                                                                                            \
            type_t *src_vals = vals, *dst_vals = vals;                                                                 \
            for (i = 0; i < nsmpl; i++)                                                                                \
            {                                                                                                          \
                dst_vals[0] = src_vals[0];                                                                             \
                dst_vals[1] = src_vals[ialt + 1];                                                                      \
                dst_vals += 2;                                                                                         \
                src_vals += nvals;                                                                                     \
            }                                                                                                          \
            bcf_update_format_##type(hdr, dst, tag, vals, nsmpl * 2);                                                  \
        }                                                                                                              \
        else if (len == BCF_VL_G)                                                                                      \
        {                                                                                                              \
            if (nvals != src->n_allele * (src->n_allele + 1) / 2 * nsmpl && nvals != src->n_allele * nsmpl)            \
                error(                                                                                                 \
                    "Error at %s:%d, the tag %s has wrong number of fields\n", bcf_seqname(hdr, src), src->pos + 1,    \
                    bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id));                                                          \
            nvals /= nsmpl;                                                                                            \
            int all_haploid = nvals == src->n_allele ? 1 : 0;                                                          \
            type_t *src_vals = vals, *dst_vals = vals;                                                                 \
            for (i = 0; i < nsmpl; i++)                                                                                \
            {                                                                                                          \
                int haploid = all_haploid;                                                                             \
                if (!haploid)                                                                                          \
                {                                                                                                      \
                    int j;                                                                                             \
                    for (j = 0; j < nvals; j++)                                                                        \
                        if (is_vector_end)                                                                             \
                            break;                                                                                     \
                    if (j != nvals)                                                                                    \
                        haploid = 1;                                                                                   \
                }                                                                                                      \
                dst_vals[0] = src_vals[0];                                                                             \
                if (haploid)                                                                                           \
                {                                                                                                      \
                    dst_vals[1] = src_vals[ialt + 1];                                                                  \
                    if (!all_haploid)                                                                                  \
                        set_vector_end;                                                                                \
                }                                                                                                      \
                else                                                                                                   \
                {                                                                                                      \
                    dst_vals[1] = src_vals[bcf_alleles2gt(0, ialt + 1)];                                               \
                    dst_vals[2] = src_vals[bcf_alleles2gt(ialt + 1, ialt + 1)];                                        \
                }                                                                                                      \
                dst_vals += all_haploid ? 2 : 3;                                                                       \
                src_vals += nvals;                                                                                     \
            }                                                                                                          \
            bcf_update_format_##type(hdr, dst, tag, vals, all_haploid ? nsmpl * 2 : nsmpl * 3);                        \
        }                                                                                                              \
        else                                                                                                           \
            bcf_update_format_##type(hdr, dst, tag, vals, nvals);                                                      \
    }
    switch (bcf_hdr_id2type(hdr, BCF_HL_FMT, fmt->id))
    {
    case BCF_HT_INT:
        BRANCH_NUMERIC(int32, int32_t, src_vals[j] == bcf_int32_vector_end, dst_vals[2] = bcf_int32_vector_end);
        break;
    case BCF_HT_REAL:
    {
        // cppcheck-suppress invalidPointerCast
        BRANCH_NUMERIC(float, float, bcf_float_is_vector_end(src_vals[j]), bcf_float_set_vector_end(dst_vals[2]));
        break;
    }
    }
#undef BRANCH_NUMERIC
}

static void squeeze_format_char(char* str, int src_blen, int dst_blen, int n)
{
    int i, isrc = 0, idst = 0;
    for (i = 0; i < n; i++)
    {
        memmove(str + idst, str + isrc, dst_blen);
        idst += dst_blen;
        isrc += src_blen;
    }
}

static void split_format_string(
    const bcf_hdr_t* hdr, bcf1_t* src, bcf_fmt_t* fmt, int ialt, bcf1_t* dst, uint8_t*& tmp_arr1, int& ntmp_arr1)
{
    const char* tag = bcf_hdr_int2id(hdr, BCF_DT_ID, fmt->id);
    int ret = bcf_get_format_char(hdr, src, tag, &tmp_arr1, &ntmp_arr1);
    assert(ret > 0);

    kstring_t str;
    str.m = ntmp_arr1;
    str.l = ret;
    str.s = (char*)tmp_arr1;

    int nsmpl = bcf_hdr_nsamples(hdr);
    int len = bcf_hdr_id2length(hdr, BCF_HL_FMT, fmt->id);
    if (len == BCF_VL_A)
    {
        int i, blen = ret / nsmpl, maxlen = 0;
        char* ptr = str.s;
        for (i = 0; i < nsmpl; i++)
        {
            char* tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(tmp, tmp, ptr + blen, ialt, len);
            if (len < 0)
            {
                return;
            } // wrong number of fields: skip
            if (maxlen < len)
                maxlen = len;
            ptr += blen;
        }
        if (maxlen < blen)
        {
            squeeze_format_char(str.s, blen, maxlen, nsmpl);
        }
        bcf_update_format_char(hdr, dst, tag, str.s, nsmpl * maxlen);
    }
    else if (len == BCF_VL_R)
    {
        int i, blen = ret / nsmpl, maxlen = 0;
        char* ptr = str.s;
        for (i = 0; i < nsmpl; i++)
        {
            char* tmp = ptr;
            int len = 0;
            STR_MOVE_NTH(ptr, tmp, ptr + blen, 0, len);
            ptr[len] = ',';
            tmp++;
            len++;
            STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, ialt, len);
            if (len < 0)
            {
                return;
            } // wrong number of fields: skip
            if (maxlen < len)
                maxlen = len;
            ptr += blen;
        }
        if (maxlen < blen)
        {
            squeeze_format_char(str.s, blen, maxlen, nsmpl);
        }
        bcf_update_format_char(hdr, dst, tag, str.s, nsmpl * maxlen);
    }
    else if (len == BCF_VL_G)
    {
        int i, blen = ret / nsmpl, maxlen = 0, i0a = bcf_alleles2gt(0, ialt + 1);
        // cppcheck-suppress duplicateExpression
        int iaa = bcf_alleles2gt(ialt + 1, ialt + 1);
        char* ptr = str.s;
        for (i = 0; i < nsmpl; i++)
        {
            char *se = ptr, *sx = ptr + blen;
            int nfields = 1;
            while (*se && se < sx)
            {
                if (*se == ',')
                {
                    nfields++;
                }
                se++;
            }
            assert(nfields == src->n_allele * (src->n_allele + 1) / 2 || nfields == src->n_allele);
            int len = 0;
            if (nfields == src->n_allele) // haploid
            {
                char* tmp = ptr;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, 0, len);
                ptr[len] = ',';
                tmp++;
                len++;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, ialt, len);
                if (len < 0)
                {
                    return;
                } // wrong number of fields: skip
            }
            else // diploid
            {
                char* tmp = ptr;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, 0, len);
                ptr[len] = ',';
                tmp++;
                len++;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, i0a - 1, len);
                if (len < 0)
                {
                    return;
                } // wrong number of fields: skip
                ptr[len] = ',';
                tmp++;
                len++;
                STR_MOVE_NTH(&ptr[len], tmp, ptr + blen, iaa - i0a - 1, len);
                if (len < 0)
                    return; // wrong number of fields: skip
            }
            if (maxlen < len)
            {
                maxlen = len;
            }
            ptr += blen;
        }
        if (maxlen < blen)
        {
            squeeze_format_char(str.s, blen, maxlen, nsmpl);
        }
        bcf_update_format_char(hdr, dst, tag, str.s, nsmpl * maxlen);
    }
    else
        bcf_update_format_char(hdr, dst, tag, str.s, str.l);
}

void splitMultiAllelics(bcf_hdr_t* hdr, bcf1_t* line, std::vector<p_bcf1>& out)
{
    bcf_unpack(line, BCF_UN_ALL);

    // Init the target biallelic lines
    int ntmp_lines = line->n_allele - 1;
    kstring_t tmp = { 0, 0, 0 };
    kputs(line->d.allele[0], &tmp);
    kputc(',', &tmp);
    int rlen = tmp.l;
    int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");

    // work arrays - this should probably be moved outside the function.
    uint8_t* tmp_arr1 = NULL;
    int ntmp_arr1 = 0;

    out.clear();
    for (int i = 0; i < ntmp_lines; i++) // for each ALT allele
    {
        p_bcf1 p_dst = p(bcf_init1());
        bcf1_t* dst = p_dst.get();
        out.push_back(p_dst);

        dst->rid = line->rid;
        dst->pos = line->pos;
        dst->qual = line->qual;

        // Not quite sure how to handle IDs, they can be assigned to a specific
        // ALT.  For now we leave the ID unchanged for all.
        bcf_update_id(hdr, dst, line->d.id ? line->d.id : ".");

        tmp.l = rlen;
        kputs(line->d.allele[i + 1], &tmp);
        bcf_update_alleles_str(hdr, dst, tmp.s);

        if (line->d.n_flt)
        {
            bcf_update_filter(hdr, dst, line->d.flt, line->d.n_flt);
        }

        for (int j = 0; j < line->n_info; j++)
        {
            bcf_info_t* info = &line->d.info[j];
            int type = bcf_hdr_id2type(hdr, BCF_HL_INFO, info->key);
            if (type == BCF_HT_INT || type == BCF_HT_REAL)
            {
                split_info_numeric(hdr, line, info, i, dst, tmp_arr1, ntmp_arr1);
            }
            else if (type == BCF_HT_FLAG)
            {
                split_info_flag(hdr, line, info, i, dst, tmp_arr1, ntmp_arr1);
            }
            else
            {
                split_info_string(hdr, line, info, i, dst, tmp_arr1, ntmp_arr1);
            }
        }

        dst->n_sample = line->n_sample;
        for (int j = 0; j < line->n_fmt; j++)
        {
            bcf_fmt_t* fmt = &line->d.fmt[j];
            int type = bcf_hdr_id2type(hdr, BCF_HL_FMT, fmt->id);
            if (fmt->id == gt_id)
            {
                split_format_genotype(hdr, line, fmt, i, dst, tmp_arr1, ntmp_arr1);
            }
            else if (type == BCF_HT_INT || type == BCF_HT_REAL)
            {
                split_format_numeric(hdr, line, fmt, i, dst, tmp_arr1, ntmp_arr1);
            }
            else
            {
                split_format_string(hdr, line, fmt, i, dst, tmp_arr1, ntmp_arr1);
            }
        }
    }
    free(tmp.s);
    free(tmp_arr1);
}

bcf1_t* extractRefFromMNP(bcf_hdr_t* hdr, bcf1_t* rec, int i)
{
    assert(rec->n_allele > 1);

    int nvali = 1, nvalf = 1;
    int32_t* worki = (int*)malloc(sizeof(int32_t) * nvali);
    float* workf = (float*)malloc(sizeof(float) * nvalf);

    bcf1_t* new_rec = bcf_init();
    new_rec->pos = rec->pos + i;
    new_rec->rid = rec->rid;
    bcf_float_set_missing(new_rec->qual);

    char alleles[4] = "X,.";
    alleles[0] = rec->d.allele[0][i];
    bcf_update_alleles_str(hdr, new_rec, alleles);

    bcf_update_filter(hdr, new_rec, rec->d.flt, rec->d.n_flt);

    int homref_gt[2];
    homref_gt[0] = homref_gt[1] = bcf_gt_unphased(0);
    bcf_update_genotypes(hdr, new_rec, homref_gt, 2);

    if (bcf_get_format_float(hdr, rec, "GQ", &workf, &nvalf) == 1)
    {
        worki[0] = (int32_t)workf[0];
        assert(bcf_update_format_int32(hdr, new_rec, "GQX", worki, 1) == 0);
    }

    if (bcf_get_format_int32(hdr, rec, "DP", &worki, &nvali) != 1)
    {
        // TODO should this actually be an error?
        bcf_destroy1(new_rec);
        return NULL;
    }
    assert(bcf_update_format_int32(hdr, new_rec, "DP", worki, 1) == 0);
    if (bcf_get_format_int32(hdr, rec, "DPF", &worki, &nvali) != 1)
    {
        bcf_destroy1(new_rec);
        return NULL;
    }
    bcf_update_format_int32(hdr, new_rec, "DPF", worki, 1);
    free(workf);
    free(worki);

    return new_rec;
}
} // namespace bcfhelpers
