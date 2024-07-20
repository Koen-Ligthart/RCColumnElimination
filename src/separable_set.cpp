#include <algorithm>
#include <ios>
#include <iostream>
#include <sstream>
#include "profiler.h"
#include "scip_exception_util.hpp"
#include "separable_set.h"

// creates empty LP with empty constraints
SCIP_LPI * buildLPWithConstraints(int const d) {
    SCIP_LPI * scipLpi;
    // create LP
    SCIP_CALL_EXC( SCIPlpiCreate(&scipLpi, NULL, "separable_set", SCIP_OBJSEN_MAXIMIZE) );

    // add constraints
    double * const rhs = new double[d + 1](); // = 0
    char * * const rownames = new char *[d + 1];
    std::ostringstream namebuf{std::ios_base::ate};
    for (int i = 0; i <= d; ++i) {
        namebuf.str("c");
        namebuf << i;
        rownames[i] = strdup(namebuf.str().c_str());
    }
    SCIP_CALL_EXC( SCIPlpiAddRows(scipLpi, d + 1, rhs, rhs, rownames, 0, NULL, NULL, NULL) );
    delete[] rhs;
    for (int i = 0; i <= d; ++i)
        free(rownames[i]);
    delete[] rownames;

    return scipLpi;
}

enum LPColType { LPColType_x, LPColType_y };

// add columns to the LP for the specified set of points
template<LPColType colType>
void LPAddCols(SCIP_LPI * scipLpi, int const d, std::vector<Point const *> const & points, AlgorithmSettings const & settings) {
    std::ostringstream namebuf{std::ios_base::ate};
    int const n = points.size();
    double * const obj = new double[n];
    std::fill(obj, obj + n, colType == LPColType_x ? 0.0 : 1.0);
    double * const lb = new double[n](); // = 0
    double * const ub = new double[n];
    std::fill(ub, ub + n, SCIPlpiInfinity(scipLpi));
    char ** const colnames = new char *[n];
    for (int i = 0; i < n; ++i) {
        namebuf.str(colType == LPColType_x ? "x" : "y");
        namebuf << i;
        colnames[i] = strdup(namebuf.str().c_str());
    }
    int nnonz = 0;
    std::vector<int> beg(n);
    std::vector<int> ind;
    std::vector<double> val;
    constexpr double consCoeffMult = colType == LPColType_x ? 1.0 : -1.0;
    for (int i = 0; i < n; ++i) {
        beg[i] = nnonz;
        for (int j = 0; j <= d; ++j) {
            double coordinate = j == d ? consCoeffMult : (consCoeffMult * (*points[i])[j]);
            if (!settings.isZero(coordinate)) {
                ind.push_back(j);
                val.push_back(coordinate);
                ++nnonz;
            }
        }
    }
    SCIP_CALL_EXC( SCIPlpiAddCols(scipLpi, n, obj, lb, ub, colnames, nnonz, &beg[0], &ind[0], &val[0]) );
    delete[] obj;
    delete[] lb;
    delete[] ub;
    for (int i = 0; i < n; ++i)
        free(colnames[i]);
    delete[] colnames;
}

// add columns to the LP for the specified set of points
template<LPColType colType>
void LPAddCols(SCIP_LPI * scipLpi, int const d, std::vector<Point> const & points, AlgorithmSettings const & settings) {
    std::vector<Point const *> pointers(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        pointers[i] = &points[i];
    LPAddCols<LPColType_x>(scipLpi, d, pointers, settings);
}

SeparableSet::SeparableSet(RCInstance const & newInstance, AlgorithmSettings const & newSettings) : instance{ newInstance }, settings{ newSettings } {
    PROFILE(PROFILE_IndependenceOracle);

    scipLpi = buildLPWithConstraints(instance.d);

    // add initial columns corresponding to X
    LPAddCols<LPColType_x>(scipLpi, instance.d, instance.X, settings);
}

SeparableSet::~SeparableSet() {
    SCIP_CALL_NOEXC( SCIPlpiFree(&scipLpi) );
}

void SeparableSet::add(const int i) {
    PROFILE(PROFILE_IndependenceOracle);

    points.push_back(i);
    LPAddCols<LPColType_y>(scipLpi, instance.d, {&instance.Y[i]}, settings);
    SCIP_CALL_EXC( SCIPlpiSolvePrimal(scipLpi) );
}

void SeparableSet::removeLastAdded() {
    PROFILE(PROFILE_IndependenceOracle);

    assert( !points.empty() );
    const int c = instance.X.size() + points.size() - 1;
    SCIP_CALL_EXC( SCIPlpiDelCols(scipLpi, c, c) );
    points.pop_back();
    SCIP_CALL_EXC( SCIPlpiSolvePrimal(scipLpi) );
}

bool SeparableSet::isIndependent() const {
    PROFILE(PROFILE_IndependenceOracle);

    if (points.empty() || SCIPlpiIsOptimal(scipLpi)) {
        return true;
    } else {
        assert( SCIPlpiHasPrimalRay(scipLpi) );
        return false;
    }
}

VertexList SeparableSet::findContainedCircuit() const {
    PROFILE(PROFILE_IndependenceOracle);

    assert( SCIPlpiHasPrimalRay(scipLpi) );
    // detect which points contribute to the unbounded ray
    double * const ray = new double[instance.X.size() + points.size()];
    SCIP_CALL_EXC( SCIPlpiGetPrimalRay(scipLpi, ray) );
    VertexList conflictCandidates;
    std::vector<Point const *> YInRay;
    for (size_t l = 0; l < points.size(); ++l) {
        if (!settings.isZero(ray[instance.X.size() + l])) {
            conflictCandidates.push_back(points[l]);
            YInRay.push_back(&instance.Y[points[l]]);
        }
    }
    delete[] ray;

    // create new LP with only the points contributing to the ray
    SCIP_LPI * secondaryScipLpi = buildLPWithConstraints(instance.d);
    LPAddCols<LPColType_x>(secondaryScipLpi, instance.d, instance.X, settings);
    LPAddCols<LPColType_y>(secondaryScipLpi, instance.d, YInRay, settings);
    
    // try to remove each point if possible
    // observe that is guaranteed that our last point added remains in points as before its addition there were no conflicts
    std::vector<bool> essential(conflictCandidates.size());
    essential[conflictCandidates.size() - 1] = true;
    for (size_t p = 0; p < conflictCandidates.size() - 1; ++p) {
        SCIP_CALL_EXC( SCIPlpiDelCols(secondaryScipLpi, instance.X.size(), instance.X.size()) );
        SCIP_CALL_EXC( SCIPlpiSolvePrimal(secondaryScipLpi) );
        if (SCIPlpiIsOptimal(secondaryScipLpi)) { // this point was essential for the conflict, add it back
            essential[p] = true;
            LPAddCols<LPColType_y>(secondaryScipLpi, instance.d, {&instance.Y[conflictCandidates[p]]}, settings);
        } else { // this point was not essential for the conflict, we can leave it out
            assert( SCIPlpiHasPrimalRay(secondaryScipLpi) );
        }
    }
    SCIP_CALL_NOEXC( SCIPlpiFree(&secondaryScipLpi) );
    VertexList conflict;
    for (size_t p = 0; p < conflictCandidates.size(); ++p)
        if (essential[p])
            conflict.push_back(conflictCandidates[p]);
    return conflict;
}

VertexList const & SeparableSet::getElements() const {
    return points;
}