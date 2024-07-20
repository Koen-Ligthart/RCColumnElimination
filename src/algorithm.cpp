#include "algorithm.h"
#include "hypergraph_node_state.hpp"
#include "iterative_refinement.hpp"
#include "separable_set.h"

/** color a hypergraph with the given number of vertices in a state oblivious manner */
std::optional<std::vector<VertexList>> color(int const n, IndependenceSystem const & independenceSystem, std::vector<int> const & variableOrder, std::optional<std::vector<size_t>> const & cliqueSizes, AlgorithmSettings const & settings) {
    assert(settings.separator != CONFLICT_SEPARATOR_LocalStateful);
    if (settings.separator == CONFLICT_SEPARATOR_Local) {
        UniqueTableBDD<SimpleParentListedNode> bdd(n, variableOrder);
        return std::get<0>(iterativelyRefine(bdd, independenceSystem, {}, 0, cliqueSizes, settings));
    } else { // settings.separator == CONFLICT_SEPARATOR_Global
        UniqueTableBDD<SimpleReferenceCountedNode> bdd(n, variableOrder);
        return std::get<0>(iterativelyRefine(bdd, independenceSystem, {}, 0, cliqueSizes, settings));
    }
}

std::optional<std::vector<VertexList>> color(Hypergraph const & H, std::vector<int> const & variableOrder, std::optional<std::vector<size_t>> const & cliqueSizes, AlgorithmSettings const & settings) {
    std::vector<VertexList> sortedE;
    sortedE.reserve(H.E.size());
    // filter out nonminimal edges (additionally filters out duplicate edges)
    for (VertexList const & e : H.E) {
        VertexList sortede{e};
        std::ranges::sort(sortede);
        bool superfluous = false;
        // check whether this is contained in another edge
        for (VertexList const & other : sortedE) {
            if (std::ranges::includes(sortede, other)) {
                superfluous = true;
                break;
            }
        }
        if (!superfluous)
            sortedE.emplace_back(std::move(sortede));
    }
    auto const independenceSystem = [&sortedE]{
        return std::make_unique<ExplicitIndependentSet>(sortedE);
    };
    return color(H.n, independenceSystem, variableOrder, cliqueSizes, settings);
}

std::optional<std::vector<VertexList>> color(OrdinaryGraph const & G, std::vector<int> const & variableOrder, std::optional<std::vector<size_t>> const & cliqueSizes, AlgorithmSettings const & settings) {
    assert( settings.separator == CONFLICT_SEPARATOR_LocalStateful );
    IndependenceSystem independenceSystem = [&G]{
        return std::make_unique<OrdinaryGraphIndependentSet>(G);
    };
    StatefulBDD<OrdinaryGraphStatefulNode, OrdinaryGraphNodeStateManager> bdd(G.n, variableOrder, G);
    return std::get<0>(iterativelyRefine(bdd, independenceSystem, {}, 0, cliqueSizes, settings));
}

std::optional<std::vector<VertexList>> computeRC(RCInstance const & instance, std::vector<int> const & variableOrder, std::optional<std::vector<size_t>> const & cliqueSizes, AlgorithmSettings const & settings) {
    IndependenceSystem independenceSystem = [&instance, &settings]{
        return std::make_unique<SeparableSet>(instance, settings);
    };
    return color(instance.Y.size(), independenceSystem, variableOrder, cliqueSizes, settings);
}

std::optional<std::vector<VertexList>> computeRCTwoStep(RCInstance const & instance, std::vector<int> const & variableOrder, std::optional<std::vector<size_t>> const & cliqueSizes, OrdinaryGraph const & hidingGraph, AlgorithmSettings const & settings) {
    assert( settings.separator == CONFLICT_SEPARATOR_LocalStateful );
    // color hiding graph
    IndependenceSystem hidingGraphIndependenceSystem = [&hidingGraph]{
        return std::make_unique<OrdinaryGraphIndependentSet>(hidingGraph);
    };
    StatefulBDD<OrdinaryGraphStatefulNode, OrdinaryGraphNodeStateManager> hidingGraphBdd(hidingGraph.n, variableOrder, hidingGraph);
    auto const result = iterativelyRefine(hidingGraphBdd, hidingGraphIndependenceSystem, {}, 0, cliqueSizes, settings);
    if (!std::get<0>(result))
        return {}; // return if infeasible
    
    // modify settings
    AlgorithmSettings newSettings(settings);
    newSettings.longestPathIterations = 0;

    // copy and transform BDD to use hypergraph states
    StatefulBDD<HypergraphStatefulNode, HypergraphNodeStateManager> hypergraphBdd(hidingGraphBdd, instance.Y.size(), variableOrder);
    // learn all hiding graph edges
    for (VertexList const & E : hidingGraph.E)
        hypergraphBdd.registerEdge(E);
    // color hypergraph
    IndependenceSystem independenceSystem = [&instance, &settings]{
        return std::make_unique<SeparableSet>(instance, settings);
    };
    std::cout << "proceeding to find hypergraph coloring" << std::endl;
    return std::get<0>(iterativelyRefine(hypergraphBdd, independenceSystem, std::get<2>(result), std::get<1>(result), cliqueSizes, newSettings));
}

std::optional<std::vector<VertexList>> computeRCHypergraph(RCInstance const & instance, std::vector<int> const & variableOrder, std::optional<std::vector<size_t>> const & cliqueSizes, OrdinaryGraph const & hidingGraph, AlgorithmSettings const & settings) {
    assert( settings.separator == CONFLICT_SEPARATOR_LocalStateful );
    IndependenceSystem independenceSystem = [&instance, &settings]{
        return std::make_unique<SeparableSet>(instance, settings);
    };
    StatefulBDD<HypergraphStatefulNode, HypergraphNodeStateManager> bdd(instance.Y.size(), variableOrder, instance.Y.size(), variableOrder);
    // learn all hiding graph edges
    for (const VertexList &E : hidingGraph.E)
        bdd.registerEdge(E);
    auto const result = iterativelyRefine(bdd, independenceSystem, {}, 0, cliqueSizes, settings);
    return std::get<0>(result);
}
