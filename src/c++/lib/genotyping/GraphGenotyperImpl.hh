
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
 * Internal data structure for graph genotyper
 *
 * \author Sai Chen & Egor Dolzhenko & Peter Krusche
 * \email schen6@illumina.com & pkrusche@illumina.com & edolzhenko@illumina.com
 *
 */

#include "genotyping/BreakpointFinder.hh"
#include "graphcore/Graph.hh"
#include "json/json.h"

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#pragma once

namespace genotyping
{

using graphtools::Graph;

struct GraphGenotyper::GraphGenotyperImpl
{
    /**
     * Our graph
     */
    Graph const* graph;

    /**
     *  depths per sample
     */
    std::vector<std::pair<double, int>> depths;

    /**
     * sex information per sample
     */
    std::vector<SampleInfo::Sex> sexes;

    /**
     *  Read counts for each sample
     */
    std::vector<BreakpointMap> breakpoint_maps;

    /**
     * Name of each sample
     */
    std::vector<std::string> samplenames;

    /**
     * sample name -> index
     */
    std::unordered_map<std::string, size_t> samplenameindex;

    /**
     *  whole-variant and breakpoint genotypes
     */
    std::map<std::pair<std::string, std::string>, Genotype> graph_genotypes;

    /**
     * Extract basic event info from alignment JSONs
     */
    Json::Value basic_info;

    /**
     * Names of breakpoints, alleles and edges
     */
    std::list<std::string> breakpointnames;
    std::vector<std::string> allelenames;
};
};