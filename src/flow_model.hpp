#pragma once

#include <algorithm>
#include <ios>
#include <ranges>
#include <sstream>
#include <unordered_set>
#include "scip/scipdefplugins.h"
#include "scip_exception_util.hpp"
#include "bdd_base.hpp"
#include "primitives.h"
#include "settings.h"

/** result of a flow decomposition algorithm invocation */
struct FlowDecompositionResult {

    std::optional<std::vector<VertexList>>      partition;                  /** optional optimal color partitioning if it was found */
    std::vector<ConflictSolution>               conflicts;                  /** list of identified conflicts */

};

/** container for the generalized flow model over a BDD */
template<IsBDD BDD>
class FlowModel {

private:
    
    SCIP *                                      scip;                       /** pointer to SCIP model */
    BDD const &                                 bdd;                        /** BDD */
    bool const                                  relaxed;                    /** whether the integrality constraints of the SCIP model are relaxed */
    std::unordered_map<Arc, std::pair<SCIP_Var *, double>, ArcHash>
                                                arcVars;                    /** map from arcs to arc flow SCIP variables
                                                                             *  and value in optimal solution if it exists during decompose()
                                                                             */
    std::vector<SCIP_Cons *>                    coverConss;                 /** SCIP constraints ensuring that every variable is included on some path */
    SCIP_Cons *                                 objCut;                     /** objective bound constraint */
    std::vector<SCIP_Cons *>                    cliqueCuts;                 /** constraints restricting the solution based on maximal cliques for each vertex */
    std::unordered_map<NodePtr, SCIP_Cons *, NodePtrHash>
                                                flowConss;                  /** SCIP constraints modelling flow conservation over internal nodes */
    IndependenceSystem const &                  independenceSystem;         /** interface to the underlying hypergraph */
    AlgorithmSettings const &                   settings;                   /** settings */

    FlowModel(FlowModel const &) = delete;
    FlowModel(FlowModel &&) = delete;
    FlowModel & operator=(FlowModel const &) = delete;
    FlowModel & operator=(FlowModel &&) = delete;

    /** copies the SCIP solution into arcVars */
    void copyFlowSolution() {
        // we cannot decompose infeasible models (i.e. also not when r == t0)
        assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
        // load flow solution into variables
        SCIP_Sol * const sol = SCIPgetBestSol(scip);
        for (auto &[_, arcVar] : arcVars)
            arcVar.second = relaxed ? SCIPgetSolVal(scip, sol, arcVar.first) : std::round(SCIPgetSolVal(scip, sol, arcVar.first));
    }

    /** follows the flow using the chooser to decompose a single path from the solution
     *  the decomposed path with its corresponding flow is returned
     *  if the flow in the network is zero, this does not return a path
     */
    std::optional<std::pair<Path, double>> decomposeFlowPath(
        std::function<bool(double, double)> const &
                                                flowChooser                 /** given f_0, f_1 (flow to corresponding children over x-arcs) outputs x
                                                                             *  to declare that the path should move to the child at the end of the x-arc
                                                                             */
    ) {
        if (!settings.isPositive(arcVars[{bdd.getRoot(), 0}].second))
            return {};

        Path path{{bdd.getRoot()}, {}};
            
        // decompose flow along solution
        double flow = SCIP_DEFAULT_INFINITY;
        NodePtr p = bdd.getRoot();
        int lP;
        while ((lP = bdd.getLayer(p)) < bdd.getN()) {
            NodePtrPair C = bdd.getChildren(p);
            double F[2];
            for (bool const d : BOOLS) // consider d-arc (p, c)
                F[d] = C[d] == bdd.getTerminal(0) ? 0.0 : arcVars[{p, d}].second;
            bool const d = flowChooser(F[0], F[1]);
            double const f = F[d];
            assert( relaxed || settings.isPositive(f) );
            if (!settings.isPositive(f)) {
                // this only happens due to numerical inaccuracies where flow conservation is not preserved
                // as this only happens in a relaxed flow, missing out on a potential conflict is not problematic
                // we do need to reduce the flow along the path found so far to prevent finding the same path again
                for (size_t k = 0; k < path.length(); ++k)
                    arcVars[{path.nodes[k], path.arcLabels[k]}].second -= flow;
                return decomposeFlowPath(flowChooser); // try to decompose another path
            }
            NodePtr const c = C[d];
            path.nodes.push_back(c);
            path.arcLabels.push_back(d);
            flow = std::min(flow, f);
            p = c;
        }

        // reduce flow along path
        for (size_t k = 0; k < path.length(); ++k)
            arcVars[{path.nodes[k], path.arcLabels[k]}].second -= flow;

        return {{path, flow}};
    }

    /** decomposes the flow solution */
    std::optional<std::vector<VertexList>> decompose(
        std::function<bool(std::unordered_set<int> const &, std::vector<std::unique_ptr<IndependentSet>> const &)> const &
                                                continueCondition,          /** whether to continue decomposing the next path (assuming it exists)
                                                                             *  the first argument is the collection of vertices that have already been treated
                                                                             *  the second argument is the current partition
                                                                             */          
        std::function<bool(double, double)> const &
                                                flowChooser,                /** given f_0, f_1 (flow to corresponding children over x-arcs) outputs x
                                                                             *  to declare that the path should move to the child at the end of the x-arc
                                                                             */
        std::function<std::unique_ptr<IndependentSet>(std::unordered_set<int> &, Path const &, double)> const &
                                                pathFinisher,               /** performs operations given a path that was just obtained from the flow
                                                                             *  the first argument is the collection of vertices that have already been treated
                                                                             *  and the second argument is the set of vertices associated with this path
                                                                             *  the third argument is the minimum flow that was present on the path
                                                                             *  this is expected to add the vertices that are able to be chosen on this path to
                                                                             *  be added to an independent set that is returned
                                                                             */         
        std::function<bool(std::unordered_set<int> &, std::vector<std::unique_ptr<IndependentSet>> &)> const &
                                                finisher                    /** performs work after decomposing all paths and if this returns false
                                                                             *  decompose will return the empty optional to avoid constructing a final partition
                                                                             *  the first argument is the collection of vertices that have already been treated
                                                                             *  the second argument is the current partition
                                                                             */
    ) {
        copyFlowSolution();

        std::unordered_set<int> partitioned;
        std::vector<std::unique_ptr<IndependentSet>> partition;

        std::optional<std::pair<Path, double>> decomposedPath;
        while (continueCondition(partitioned, partition) && (decomposedPath = decomposeFlowPath(flowChooser))) {
            std::unique_ptr<IndependentSet> independentSet = pathFinisher(partitioned, decomposedPath.value().first, decomposedPath.value().second);
            if (!independentSet->getElements().empty())
                partition.emplace_back(std::move(independentSet));
        }

        if (!finisher(partitioned, partition))
            return {};
        
        std::vector<VertexList> result;
        result.reserve(partition.size());
        for (std::unique_ptr<IndependentSet> const & independentSet : partition) {
            result.emplace_back(independentSet->getElements());
            std::ranges::sort(result.back());
        }
        return result;
    }

public:

    /** constructs flow model in SCIP for given DD */
    FlowModel(
        BDD const &                             newBdd,                     /** BDD */
        bool const                              newRelaxed,                 /** whether this should be the LP relaxation of the flow model or the full IP */
        IndependenceSystem const &              newIndependenceSystem,      /** interface to the underlying hypergraph */
        size_t const                            lb,                         /** lower bound on the chromatic number */
        std::optional<std::vector<size_t>> const &
                                                cliqueSizes,                /** size of a maximal clique for each vertex */
        AlgorithmSettings const &               newSettings,                /** settings */
        std::optional<std::vector<VertexList>> const &
                                                sol                         /** valid heuristic coloring to be used to extract a heuristic feasible flow model solution */
    ) :
        bdd{ newBdd },
        relaxed { newRelaxed },
        independenceSystem{ newIndependenceSystem },
        settings{ newSettings }
    {
        PROFILE(PROFILE_FlowModelConstruction);

        // create SCIP model
        SCIP_CALL_EXC( SCIPcreate(&scip) );
        SCIP_CALL_EXC( SCIPincludeDefaultPlugins(scip) );
        if (settings.scipSet)
            SCIP_CALL_EXC( SCIPreadParams(scip, settings.scipSet.value().c_str()) );
        if (!settings.verboseSCIP)
            SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);
        SCIPmessagehdlrSetLogfile(SCIPgetMessagehdlr(scip), "./out/scip.log");
        SCIP_CALL_EXC( SCIPcreateProb(scip, "flow_model", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

        // create constraints for each vertex
        coverConss.reserve(bdd.getN());
        if (settings.cliqueCuts)
            cliqueCuts.reserve(bdd.getN());
        std::ostringstream namebuf(std::ios_base::ate);
        for (int v = 0; v < bdd.getN(); ++v) {
            // create cover constraint
            namebuf.str("c");
            namebuf << v;
            SCIP_Cons * coverCons;
            SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &coverCons, namebuf.str().c_str(), 0, NULL, NULL, 1.0, SCIP_DEFAULT_INFINITY, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            coverConss.push_back(coverCons);
            SCIP_CALL_EXC( SCIPaddCons(scip, coverCons) );

            // create clique cut
            if (settings.cliqueCuts) {
                namebuf.str("cli");
                namebuf << v;
                SCIP_Cons * cliqueCut;
                SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &cliqueCut, namebuf.str().c_str(), 0, NULL, NULL, cliqueSizes.value()[v] - 1, SCIP_DEFAULT_INFINITY, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                cliqueCuts.push_back(cliqueCut);
                SCIP_CALL_EXC( SCIPaddCons(scip, cliqueCut) );
            }
        }

        // create objective lower bound
        if (settings.flowObjCut) {
            SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &objCut, "obj", 0, NULL, NULL, lb, SCIP_DEFAULT_INFINITY, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL_EXC( SCIPaddCons(scip, objCut) );
        }

        // BFS through the BDD to add node and arc related constraints in a suitable order
        std::queue<NodePtr> Q;
        if (bdd.getRoot() != bdd.getTerminal(0))
            Q.push(bdd.getRoot());
        int arcVarNumber = 0;
        while (!Q.empty()) {
            NodePtr const p = Q.front();
            Q.pop();
            NodePtrPair const C = bdd.getChildren(p);
            for (bool const d : BOOLS) { // consider d-arc (p, c)
                NodePtr const c = C[d];
                if (c == bdd.getTerminal(0)) // flow to the 0-terminal is always 0
                    continue;
                
                int const lS = bdd.getLayer(p);
                int const lT = bdd.getLayer(c);

                // create variable for arc
                namebuf.str("f");
                namebuf << arcVarNumber;
                ++arcVarNumber;
                SCIP_Var * arcVar;
                SCIP_CALL_EXC( SCIPcreateVar(scip, &arcVar, namebuf.str().c_str(), 0.0, SCIP_DEFAULT_INFINITY, p == bdd.getRoot() ? 1.0 : 0.0, relaxed ? SCIP_VARTYPE_CONTINUOUS : SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
                arcVars.insert({{p, d}, {arcVar, 0.0}});
                SCIP_CALL_EXC( SCIPaddVar(scip, arcVar) );

                if (p == bdd.getRoot()) {
                    if (settings.flowObjCut) // add to objective bound constraint
                        SCIP_CALL_EXC( SCIPaddCoefLinear(scip, objCut, arcVar, 1.0) );
                } else { // add to flow conservation for p
                    SCIP_CALL_EXC( SCIPaddCoefLinear(scip, flowConss[p], arcVar, -1.0) );
                }

                // add to cover constraint for all covered layers
                for (int l = std::max(lS, 0); l < lT; ++l) {
                    int const v = bdd.getVariableByLayer(l);
                    if (l > lS || d == 1)
                        SCIP_CALL_EXC( SCIPaddCoefLinear(scip, coverConss[v], arcVar, 1.0) );
                    if (settings.cliqueCuts && (l > lS || d == 0))
                        SCIP_CALL_EXC( SCIPaddCoefLinear(scip, cliqueCuts[v], arcVar, 1.0) );
                }

                if (lT < bdd.getN()) { // add to flow conservation for c as long as c is no terminal
                    SCIP_Cons * & flowConsC = flowConss[c];
                    if (flowConsC == NULL) { // create constraint as it does not yet exist
                        namebuf.str("f");
                        namebuf << c;
                        SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &flowConsC, namebuf.str().c_str(), 0, NULL, NULL, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
                        SCIP_CALL_EXC( SCIPaddCons(scip, flowConsC) );
                        Q.push(c);
                    }

                    SCIP_CALL_EXC( SCIPaddCoefLinear(scip, flowConsC, arcVar, 1.0) ); 
                }
            }
        }

        // compute heuristic solution to flow model based on heuristic coloring
        if (settings.warmStart && sol) {
            SCIP_Sol * scipSol;
            SCIP_CALL_EXC( SCIPcreateSol(scip, &scipSol, NULL) );
            for (VertexList const & c : sol.value()) {
                // increase flow over path associated with class c by 1
                VertexList sortedC = c;
                std::ranges::sort(sortedC, [this](int const a, int const b){ return bdd.getLayerByVariable(a) < bdd.getLayerByVariable(b); });
                NodePtr u = bdd.getRoot();
                size_t i = 0;
                // invariant i is the least i s.t. layer(u) <= layer(sortedC[i])
                while (u != bdd.getTerminal(1)) {
                    bool b = false;
                    if (i < sortedC.size() && bdd.getVariable(u) == sortedC[i]) {
                        ++i;
                        b = true;
                    }
                    SCIP_Var * const var = arcVars[{u, b}].first;
                    SCIP_CALL_EXC( SCIPsetSolVal(scip, scipSol, var, SCIPgetSolVal(scip, scipSol, var) + 1.0) );
                    u = bdd.getChild(u, b);
                    assert(u != bdd.getTerminal(0));
                }
            }
            int unsigned stored;
            SCIP_CALL_EXC( SCIPaddSolFree(scip, &scipSol, &stored) );
        }
    }

    ~FlowModel() {
        PROFILE(PROFILE_FlowModelDestruction);
        for (auto pair : flowConss)
            SCIP_CALL_NOEXC( SCIPreleaseCons(scip, &pair.second) );
        if (settings.flowObjCut)
            SCIP_CALL_NOEXC( SCIPreleaseCons(scip, &objCut) );
        if (settings.cliqueCuts)
            for (SCIP_Cons * cons : cliqueCuts)
                SCIP_CALL_NOEXC( SCIPreleaseCons(scip, &cons) );
        for (SCIP_Cons * cons : coverConss)
            SCIP_CALL_NOEXC( SCIPreleaseCons(scip, &cons) );
        for (auto pair : arcVars)
            SCIP_CALL_NOEXC( SCIPreleaseVar(scip, &pair.second.first) );
        SCIP_CALL_NOEXC( SCIPfree(&scip) );
    }

    /** returns the number of active nodes contributing to the size of the flow model excluding the terminals and root */
    size_t activeNodeCount() const { return flowConss.size() + 2 + (bdd.getRoot() != bdd.getTerminal(0)); }

    /** writes IP / LP model to a .lp file */
    void writeToLP(std::string const & file) {
        SCIP_CALL_EXC( SCIPwriteOrigProblem(scip, file.c_str(), NULL, FALSE) );
    }

    /** solves the SCIP flow model and returns the value of an optimal solution if it exists */
    std::optional<double> solve() {
        PROFILE(PROFILE_FlowModelSolving);
        SCIP_CALL_EXC( SCIPsolve(scip) );
        SCIP_STATUS const status = SCIPgetStatus(scip);
        if (status == SCIP_STATUS_INFEASIBLE)
            return {};
        assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
        return relaxed ? SCIPgetPrimalbound(scip) : std::round(SCIPgetPrimalbound(scip));
    }

    /** obtain a mapping from arcs to associated flow values over it in the found flow model solution */
    std::unordered_map<Arc, double, ArcHash> getFlow() const {
        assert( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL );
        std::unordered_map<Arc, double, ArcHash> flow;
        SCIP_Sol * const sol = SCIPgetBestSol(scip);
        for (auto const &[arc, arcVar] : arcVars)
            flow.insert({arc, SCIPgetSolVal(scip, sol, arcVar.first)});
        return flow;
    }

    /** decomposes the flow solution and finds conflicts */
    FlowDecompositionResult findConflicts() {
        PROFILE(PROFILE_FindConflicts);
        assert( settings.maxConflictCount > 0 );

        std::vector<ConflictSolution> conflicts;
        auto const continueCondition = [&conflicts, this](std::unordered_set<int> const & partitioned, std::vector<std::unique_ptr<IndependentSet>> const & parition){
            return conflicts.size() < settings.maxConflictCount;
        };
        auto const pathFinisher = [&conflicts, this](std::unordered_set<int> & partitioned, Path const & path, double const f){
            if (relaxed && settings.isPositive(settings.minConflictFlow - f)) // only detect conflicts for flows at least a threshold, to prevent LP relaxations from generating too many useless conflicts and blowing up the diagram
                return independenceSystem();

            std::unique_ptr<IndependentSet> independentSet = independenceSystem();
            std::vector<bool> partialSolution;

            // trace path
            for (size_t i = 0; i < path.length(); ++i) {
                int const lP = bdd.getLayer(path.nodes[i]);
                int const lC = bdd.getLayer(path.nodes[i + 1]);
                for (int l = lP; l < lC; ++l) {
                    // extend partial solution
                    if (l != -1) {
                        assert(l == (int) partialSolution.size());
                        partialSolution.push_back(path.arcLabels[i] == 1 || l > lP);
                    }

                    if (path.arcLabels[i] == 1 || l > lP) {
                        // attempt to add the vertex associated to layer l to this independent set
                        independentSet->add(bdd.getVariableByLayer(l));
                        if (!independentSet->isIndependent()) {
                            VertexList const & conflictVertices = independentSet->findContainedCircuit();
                            std::vector<int> conflict;
                            conflict.reserve(conflictVertices.size());
                            for (int const vConflict : conflictVertices)
                                conflict.push_back(bdd.getLayerByVariable(vConflict));
                            assert( std::ranges::is_sorted(conflict) );
                            conflicts.emplace_back(partialSolution, conflict);
                            return independentSet;
                        }
                        // vertex is not added to partitioned because findConflicts() does not use this
                    }
                }
            }
            return independentSet;
        };
        auto const finisher = [&conflicts](std::unordered_set<int> & partitioned, std::vector<std::unique_ptr<IndependentSet>> & partition) {
            return conflicts.empty(); // return partition only if we encountered no conflicts
        };
        
        auto const partition = decompose(continueCondition, settings.getFlowChooserByType(settings.flowChooserConflictFinder), pathFinisher, finisher);
        if (partition) {
            assert( conflicts.empty() );
            return {partition, {}};
        } else {
            assert( !conflicts.empty() );
            return {{}, conflicts};
        }
    }

    /** attempts to find a heuristic feasible coloring with
     *  strictly less than ub colors using the flow solution
     *  assumes that the flow problem is not infeasible
     */
    std::optional<std::vector<VertexList>> primalHeuristic(size_t const ub) {
        PROFILE(PROFILE_PrimalHeuristic);
        if (!getProfiler().currentNodeWithinTimeBudget(settings.heuristicTimeBudget))
            return {};

        auto const continueCondition = [this, ub](std::unordered_set<int> const & partitioned, std::vector<std::unique_ptr<IndependentSet>> const & partition){
            // stop if we have already partitioned all vertices or partition still needs to have at least ub classes
            return partitioned.size() < (size_t) bdd.getN() && (settings.heuristicClassGrowth ? (partition.size() < ub - 1) : (partition.size() < ub));
        };
        auto const pathFinisher = [this, ub](std::unordered_set<int> & partitioned, Path const & path, double const f){
            std::unique_ptr<IndependentSet> independentSet = independenceSystem();

            for (size_t i = 0; i < path.length(); ++i) {
                int const lP = bdd.getLayer(path.nodes[i]);
                int const lC = bdd.getLayer(path.nodes[i + 1]);
                for (int l = lP; l < lC; ++l) {
                    if (path.arcLabels[i] == 1 || l > lP) {
                        int const v = bdd.getVariableByLayer(l);
                        // attempt to add v to this class
                        if (!partitioned.contains(v)) {
                            independentSet->add(v);
                            if (independentSet->isIndependent()) {
                                partitioned.insert(v);
                            } else {
                                // undo addition
                                independentSet->removeLastAdded();
                            }
                        }
                    }
                }
            }

            if (settings.heuristicClassGrowth && !independentSet->getElements().empty()) {
                // greedily add other vertices to this color class if possible
                for (int v = 0; v < bdd.getN(); ++v) {
                    if (!partitioned.contains(v)) {
                        independentSet->add(v);
                        if (independentSet->isIndependent()) {
                            partitioned.insert(v);
                        } else {
                            independentSet->removeLastAdded();
                        }
                    }
                }
            }

            return independentSet;
        };
        auto const finisher = [this, ub](std::unordered_set<int> & partitioned, std::vector<std::unique_ptr<IndependentSet>> & partition) {
            // if we have done heuristic class growing after separating every path, we need not consider these classes again here
            size_t const unsaturatedClassStartIndex = settings.heuristicClassGrowth ? partition.size() : 0;
            // greedily add remaining vertices
            for (int v = 0; v < bdd.getN(); ++v) {
                if (!partitioned.contains(v)) {
                    if (partition.size() >= ub - 1)
                        return false;
                    // attempt to place v in one of the already existing independent sets
                    for (size_t c = unsaturatedClassStartIndex; c < partition.size(); ++c) {
                        std::unique_ptr<IndependentSet> const & independentSet = partition[c];
                        independentSet->add(v);
                        if (independentSet->isIndependent()) {
                            partitioned.insert(v);
                            break;
                        } else {
                            independentSet->removeLastAdded();
                        }
                    }
                    // if that is not possible, create a new independent set with v in it
                    if (!partitioned.contains(v)) {
                        partition.emplace_back(independenceSystem());
                        partition.back()->add(v);
                        if (!partition.back()->isIndependent())
                            return false; // infeasible because this singleton is dependent
                        // partitioned.insert(v); not needed as partitioned is no longer checked for inclusion of v
                    }
                }
            }
            return true;
        };
        
        return decompose(continueCondition, settings.getFlowChooserByType(settings.flowChooserPrimalHeuristic), pathFinisher, finisher);
    }

};



namespace {

/** recursively computes the length of the longest path from u to the 1-terminal and updates the dist map */
template<IsBDD BDD>
size_t maxDist(
    BDD const &                                 bdd,                        /** the BDD in which to search */
    std::unordered_map<NodePtr, std::pair<size_t, bool>, NodePtrHash> &
                                                dist,                       /** map of node to the length of the longest path from it to the 1-terminal
                                                                             *  and the label of the first arc to follow to find such longest path
                                                                             */
    NodePtr const                               u                           /** the node of which the maxDist is to be evaluated */
) {
    auto const lookup = dist.find(u);
    if (lookup == dist.end()) {
        // length via 0-arc
        std::pair<size_t, bool> result = {maxDist(bdd, dist, bdd.getChild(u, 0)) + bdd.getLayer(bdd.getChild(u, 0)) - bdd.getLayer(u) - 1, 0};
        // length via 1-arc
        size_t const distC1 = bdd.getChild(u, 1) == bdd.getTerminal(0) ? 0 : (maxDist(bdd, dist, bdd.getChild(u, 1)) + bdd.getLayer(bdd.getChild(u, 1)) - bdd.getLayer(u));

        if (distC1 > result.first)
            result = {distC1, 1};
        dist.emplace(u, result);
        return result.first;
    } else {
        return (*lookup).second.first;
    }
}

}

/** find the longest r-1-terminal path in the BDD */
template<IsBDD BDD>
Path findLongestPath(BDD const & bdd) {
    PROFILE(PROFILE_LongestPath);

    // initialize dist map and compute maxDist for the root
    std::unordered_map<NodePtr, std::pair<size_t, bool>, NodePtrHash> dist;
    dist.emplace(bdd.getTerminal(1), std::pair<size_t, bool>{0, 0}); // maxDist of 1-terminal is 0
    maxDist(bdd, dist, bdd.getRoot());

    // recover path
    NodePtr u = bdd.getRoot();
    Path path({u}, {});
    while (u != bdd.getTerminal(1)) {
        bool const d = dist[u].second;
        u = bdd.getChild(u, d);
        path.nodes.push_back(u);
        path.arcLabels.push_back(d);
    }
    return path;
}

template<IsBDD BDD>
std::optional<ConflictSolution> findConflictOnPath(BDD const & bdd, Path const & path, IndependenceSystem const & independenceSystem) {
    PROFILE(PROFILE_FindConflicts);
    std::unique_ptr<IndependentSet> const independentSet = independenceSystem();
    std::vector<bool> partialSolution;
    // trace path
    for (size_t i = 0; i < path.length(); ++i) {
        int const lP = bdd.getLayer(path.nodes[i]);
        int const lC = bdd.getLayer(path.nodes[i + 1]);
        for (int l = lP; l < lC; ++l) {
            // attempt to add the vertex associated to layer l to this independent set
            if (l != -1) {
                assert(l == (int) partialSolution.size());
                partialSolution.push_back(path.arcLabels[i] == 1 || l > lP);
            }

            if (path.arcLabels[i] == 1 || l > lP) {
                // attempt to add the vertex associated to layer l to this independent set
                independentSet->add(bdd.getVariableByLayer(l));
                if (!independentSet->isIndependent()) {
                    VertexList const & conflictVertices = independentSet->findContainedCircuit();
                    std::vector<int> conflict;
                    conflict.reserve(conflictVertices.size());
                    for (int const vConflict : conflictVertices)
                        conflict.push_back(bdd.getLayerByVariable(vConflict));
                    assert( std::ranges::is_sorted(conflict) );
                    return ConflictSolution{partialSolution, conflict};
                }
            }
        }
    }
    return {};
}