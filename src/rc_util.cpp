#include "hiding_sets.h"
#include "scip_exception.hpp"
#include "profiler.h"
#include "rc_util.h"
#include "separable_set.h"

extern "C" {
    #include "problem_rc.c"
}

/** allocate a Datapoints struct with the given points using SCIP */
Datapoints * allocateDatapoints(SCIP * scip, int const d, std::vector<Point> const & points) {
    int const n = points.size();
    Datapoints * datapoints;
    SCIP_CALL_EXC( SCIPallocBlockMemory(scip, &datapoints) );
    datapoints->ndatapoints = n;
    datapoints->dimension = d;
    SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip, &(datapoints->points), n) );
    for (int i = 0; i < n; ++i)
        SCIP_CALL_EXC( SCIPallocBlockMemoryArray(scip, &(datapoints->points[i]), d) );
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < d; ++j)
            datapoints->points[i][j] = points[i][j];
    return datapoints;
}

OrdinaryGraph computeHidingGraph(RCInstance const & instance, AlgorithmSettings const & settings) {
    PROFILE(PROFILE_ComputeHidingGraph);
    int const n = instance.Y.size();
    Hypergraph H(n, {});
    if (instance.d <= settings.maxConvexHullDimension) {
        SCIP * scip;
        SCIP_CALL_EXC( SCIPcreate(&scip) );
        int * hidingsetidx;
        int nhidingsetidx;
        int maxnhidingsetidx;
        Datapoints * const X = allocateDatapoints(scip, instance.d, instance.X);
        Datapoints * const Y = allocateDatapoints(scip, instance.d, instance.Y);
        SCIP_CALL_EXC( computeHidingSetCuts(scip, X, Y, &hidingsetidx, &nhidingsetidx, &maxnhidingsetidx) );
        
        assert( nhidingsetidx % 2 == 0 );
        for (int i = 0; 2 * i < nhidingsetidx; ++i) {
            int const y = hidingsetidx[2 * i];
            int const z = hidingsetidx[2 * i + 1];
            H.E.push_back({y, z}); // y < z
        }

        if (maxnhidingsetidx != 0)
            SCIPfreeBlockMemoryArray(scip, &hidingsetidx, maxnhidingsetidx);
        SCIP_CALL_EXC( SCIPfreeModel(scip, X, Y ) ); // free the Datapoint clones
        SCIP_CALL_EXC( SCIPfree(&scip) );
    } else {
        // use LPs to determine cardinality 2 hiding sets in polynomial time
        SeparableSet E(instance, settings);
        for (int v1 = 0; v1 < n; ++v1) {
            E.add(v1);
            for (int v2 = v1 + 1; v2 < n; ++v2) {
                E.add(v2);
                if (!E.isIndependent())
                    H.E.emplace_back(std::vector<int>{v1, v2});
                E.removeLastAdded();
            }
            E.removeLastAdded();
        }
    }
    return H;
}

RCInstance filterNonobservers(RCInstance const & instance, AlgorithmSettings const & settings) {
    std::optional<RCInstance> result;
    // attempt to filter non-observers using convex hull computation
    if (instance.d <= settings.maxConvexHullDimension) {
        SCIP * scip;
        SCIP_CALL_EXC( SCIPcreate(&scip) );
        Datapoints * const X = allocateDatapoints(scip, instance.d, instance.X);
        Datapoints * const Y = allocateDatapoints(scip, instance.d, instance.Y);
        SCIP_RETCODE code = SCIPremoveNonobservers(scip, X, Y);
        if (code == SCIP_OKAY) {
            // extract points that were not filtered
            std::vector<Point> vY;
            vY.reserve(Y->ndatapoints);
            for (int i = 0; i < Y->ndatapoints; ++i) {
                std::vector<double> point(instance.d);
                for (int j = 0; j < instance.d; ++j)
                    point[j] = Y->points[i][j];
                vY.emplace_back(std::move(point));
            }
            result = {instance.d, instance.X, vY};
        } else if (code != SCIP_ERROR) {
            SCIP_CALL_EXC( code ); // attempt alternative method if failure because of a generic error
        }
        SCIP_CALL_EXC( SCIPfreeModel(scip, X, Y ) ); // free the Datapoint clones
        SCIP_CALL_EXC( SCIPfree(&scip) );
    }
    if (result)
        return result.value();

    // use LPs to determine observers in polynomial time
    result = instance;
    std::vector<bool> filtered(instance.Y.size());
    for (size_t y = 0; y < instance.Y.size(); ++y) {
        result.value().X.push_back(result.value().Y[y]);
        SeparableSet I(result.value(), settings);
        // check whether conv(X U {y}) contains any elements of Y\{y}
        for (size_t o = 0; o < instance.Y.size(); ++o) {
            if (o != y && !filtered[o]) {
                I.add(o);
                if (I.isIndependent()) {
                    I.removeLastAdded();
                } else { // o is in conv(X U {y})
                    filtered[y] = true;
                    break;
                }
            }
        }
        result.value().X.pop_back();
    }
    // update point set Y
    result.value().Y.clear();
    for (size_t y = 0; y < instance.Y.size(); ++y)
        if (!filtered[y])
            result.value().Y.push_back(instance.Y[y]);
    return result.value();
}