#pragma once

#include "ordinary_graph.h"
#include "rc_instance.h"
#include "settings.h"

/** color a hypergraph using color by explicitly enumerating all edges beforehand */
std::optional<std::vector<VertexList>> color(
    Hypergraph const &                          H,                  /** container of all the edges of the hypergraph */
    std::vector<int> const &                    variableOrder,      /** initial variable order */
    std::optional<std::vector<size_t>> const &  cliqueSizes,        /** size of a maximal clique for each vertex */
    AlgorithmSettings const &                   settings            /** algorithm settings */
);

/** color an ordinary graph using W.-J. van Hoeves algorithm by explicitly enumerating all edges beforehand */
std::optional<std::vector<VertexList>> color(
    OrdinaryGraph const &                       G,                  /** container of all the edges of the graph */
    std::vector<int> const &                    variableOrder,      /** initial variable order */
    std::optional<std::vector<size_t>> const &  cliqueSizes,        /** size of a maximal clique for each vertex */
    AlgorithmSettings const &                   settings            /** algorithm settings */
);

/** compute rc(X, Y) using iterative refinement in a state oblivious manner */
std::optional<std::vector<VertexList>> computeRC(
    RCInstance const &                          instance,           /** the problem instance */
    std::vector<int> const &                    variableOrder,      /** initial variable order */
    std::optional<std::vector<size_t>> const &  cliqueSizes,        /** size of a maximal clique for each vertex */
    AlgorithmSettings const &                   settings            /** algorithm settings */
);

/** compute rc(X, Y) using the two step approach, that is,
 *  first color the hiding graph as ordinary graph using ordinary graph states
 *  and then translate that BDD to the hypergraph setting
 *  and color it with a stateful approach using hypergraph states
 */
std::optional<std::vector<VertexList>> computeRCTwoStep(
    RCInstance const &                          instance,           /** the problem instance */
    std::vector<int> const &                    variableOrder,      /** initial variable order */
    std::optional<std::vector<size_t>> const &  cliqueSizes,        /** size of a maximal clique for each vertex */
    OrdinaryGraph const &                       hidingGraph,        /** hiding graph corresponding to the instance (X, Y) */
    AlgorithmSettings const &                   settings            /** algorithm settings */
);

/** compute rc(X, Y) using iterative refinement and hypergraph node states */
std::optional<std::vector<VertexList>> computeRCHypergraph(
    RCInstance const &                          instance,           /** the problem instance */
    std::vector<int> const &                    variableOrder,      /** initial variable order */
    std::optional<std::vector<size_t>> const &  cliqueSizes,        /** size of a maximal clique for each vertex */
    OrdinaryGraph const &                       hidingGraph,        /** hiding graph corresponding to the instance (X, Y) */
    AlgorithmSettings const &                   settings            /** algorithm settings */
);