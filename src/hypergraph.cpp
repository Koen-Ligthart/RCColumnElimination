#include <algorithm>
#include <cassert>
#include <ranges>
#include "hypergraph.h"
#include "profiler.h"

ExplicitIndependentSet::ExplicitIndependentSet(std::vector<VertexList> const & newE) : E{ newE } {}

ExplicitIndependentSet::~ExplicitIndependentSet() {}

void ExplicitIndependentSet::add(int const i) {
    PROFILE(PROFILE_IndependenceOracle);
    points.push_back(i);
}

void ExplicitIndependentSet::removeLastAdded() {
    PROFILE(PROFILE_IndependenceOracle);
    points.pop_back();
}

bool ExplicitIndependentSet::isIndependent() const {
    PROFILE(PROFILE_IndependenceOracle);
    VertexList sortedPoints{points};
    std::ranges::sort(sortedPoints);
    for (VertexList const & e : E)
        if (std::ranges::includes(sortedPoints, e))
            return false;
    return true;
}

VertexList ExplicitIndependentSet::findContainedCircuit() const {
    PROFILE(PROFILE_IndependenceOracle);
    
    std::vector<std::pair<int, int>> indexedSortedPoints; // contains entries of the form (v, i) so that points[i] = v
    indexedSortedPoints.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        indexedSortedPoints.push_back({points[i], i});
    // sort on ascending v
    std::ranges::sort(indexedSortedPoints, [](std::pair<int, int> const & a, std::pair<int, int> const & b){ return a.first < b.first; });

    std::vector<bool> inConflict(points.size()); // inConflict[i] indicates whether points[i] is part of a conflict
    std::optional<VertexList> sortedConflict; // candidate return value
    int minI = -1;
    for (VertexList const & conflict : E) {
        // check whether this set is a subset of E and update inConflict
        std::fill(inConflict.begin(), inConflict.end(), false);
        bool isConflict = true;
        size_t i = 0;
        for (int const v : conflict) {
            // scan for next element v of the candidate conflicting edge in this set
            while (i < points.size() && indexedSortedPoints[i].first != v)
                ++i;
            if (i == points.size()) { // v is not in this set
                isConflict = false;
                break;
            }
            inConflict[indexedSortedPoints[i].second] = true;
        }

        if (isConflict) {
            // extract the conflict in order of insertion
            VertexList newSortedConflict;
            int newMinI = std::numeric_limits<int>::max();
            for (size_t j = 0; j < points.size(); ++j) {
                if (inConflict[j]) {
                    newSortedConflict.push_back(points[j]);
                    if (newMinI == std::numeric_limits<int>::max())
                        newMinI = j;
                }
            }
            // only consider the conflict with maximum minimum index
            if (newMinI > minI) {
                sortedConflict = newSortedConflict;
                minI = newMinI;
            }
        }
    }

    if (sortedConflict)
        return sortedConflict.value();
    throw std::runtime_error("");
}

VertexList const & ExplicitIndependentSet::getElements() const { return points; }