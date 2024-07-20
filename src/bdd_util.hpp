#pragma once

#include <fstream>
#include <sstream>
#include "bdd_stateful.hpp"
#include "bdd_unique_table.hpp"
#include "primitives.h"



/** garbage collect to exhaustively remove all non super root non terminal nodes with indegree 0
 *  the start Q should be the set of all nodes that currently have indegree 0
 */
template<IsBDD BDD>
void garbageCollect(BDD & bdd, std::vector<NodePtr> & Q) requires IsReferenceCountedNode<class BDD::NodeType> {
    PROFILE(PROFILE_BDDGarbageCollection);
    // inspired by Kahn's algorithm
    while (!Q.empty()) {
        NodePtr const u = Q.back();
        Q.pop_back();
        if (bdd.isTerminal(u))
            continue;
        NodePtrPair const cs = bdd.getChildren(u);
        if constexpr (IsUniqueTableBDD<BDD>)
            bdd.removeFromUniqueTable(u);
        bdd.removeNode(u);
        for (NodePtr const c : cs) {
            bdd.removeArc(u, c);
            if (bdd.hasZeroIndegree(c))
                Q.push_back(c);
        }
    }
}

/** garbage collect to exhaustively remove all non super root non terminal nodes with indegree 0
 *  the start Q should be the set of all nodes that currently have indegree 0
 */
template<IsBDD BDD>
void garbageCollect(BDD &bdd, std::vector<NodePtr> &Q) requires IsParentListedNode<class BDD::NodeType> {
    PROFILE(PROFILE_BDDGarbageCollection);
    // clear unreachable nodes layer-by-layer
    auto const minLayer = [&bdd](NodePtr const a, NodePtr const b){ return bdd.getLayer(a) > bdd.getLayer(b); };
    std::priority_queue<NodePtr, std::vector<NodePtr>, decltype(minLayer)> pQ(minLayer); // this may be possible to do with topological sort
    std::unordered_set<NodePtr, NodePtrHash> inQ;
    for (NodePtr const u : Q) {
        inQ.insert(u);
        pQ.push(u);
    }
    while (!pQ.empty()) {
        NodePtr const u = pQ.top();
        pQ.pop();
        bdd.clearDeadInParentList(u);
        if (bdd.hasZeroIndegree(u) && !bdd.isTerminal(u)) {
            for (bool const d : BOOLS) {
                NodePtr const c = bdd.getChild(u, d);
                if (inQ.insert(c).second)
                    pQ.push(c);
            }
            if constexpr (IsUniqueTableBDD<BDD>)
                bdd.removeFromUniqueTable(u);
            bdd.removeNode(u);
        }
    }
}

/** performs the sifting heuristic by R. Rudell to heuristically improve the variable order of the diagram */
template<IsUniqueTableBDD BDD>
void sift(BDD & bdd) requires IsReferenceCountedNode<class BDD::NodeType> {
    for (int v = 0; v < bdd.getN(); ++v) {
        int const initialL = bdd.getLayerByVariable(v);
        int bestL = initialL;
        size_t bestSize = bdd.nodeCount();
        if (bestL <= bdd.getN() / 2) {
            // go up first then down
            for (int l = initialL - 1; l >= 0; --l) {
                bdd.swapLayers(l);
                size_t const newSize = bdd.nodeCount();
                if (newSize < bestSize) {
                    bestL = l;
                    bestSize = newSize;
                }
            }
            for (int l = 0; l < bdd.getN() - 1; ++l) {
                bdd.swapLayers(l);
                size_t const newSize = bdd.nodeCount();
                if (newSize < bestSize || (newSize == bestSize && abs(l + 1 - initialL) < abs(bestL - initialL))) {
                    bestL = l + 1;
                    bestSize = newSize;
                }
            }
            // move to up to locally optimal layer
            for (int l = bdd.getN() - 2; l >= bestL; --l)
                bdd.swapLayers(l);
        } else {
            // go down first then up
            for (int l = initialL; l < bdd.getN() - 1; ++l) {
                bdd.swapLayers(l);
                size_t const newSize = bdd.nodeCount();
                if (newSize < bestSize) {
                    bestL = l + 1;
                    bestSize = newSize;
                }
            }
            for (int l = bdd.getN() - 2; l >= 0; --l) {
                bdd.swapLayers(l);
                size_t const newSize = bdd.nodeCount();
                if (newSize < bestSize || (newSize == bestSize && abs(l - initialL) < abs(bestL - initialL))) {
                    bestL = l;
                    bestSize = newSize;
                }
            }
            // move to down to locally optimal layer
            for (int l = 0; l < bestL; ++l)
                bdd.swapLayers(l);
        }
    }
}

/** performs the sifting heuristic by R. Rudell to heuristically improve the variable order of the diagram */
template<IsUniqueTableBDD BDD>
void sift(BDD & bdd) requires IsParentListedNode<class BDD::NodeType> {
    // copy BDD to one using reference counting to reduce overhead
    UniqueTableBDD<SimpleReferenceCountedNode> referenceCountedBDD(bdd);
    sift(referenceCountedBDD);
    // restore the parent lists
    bdd = referenceCountedBDD;
}

/** count the number of subsets of {v(u), v(u) + 1, ..., n} represented by u */
template<IsBDD BDD>
double long countSat(BDD const & bdd, std::unordered_map<NodePtr, double long, NodePtrHash> & computedTable, NodePtr const u) {
    auto const cachedResult = computedTable.find(u);
    if (cachedResult != computedTable.end())
        return cachedResult->second;

    if (u == bdd.getTerminal(0))
        return 0;
    if (u == bdd.getTerminal(1))
        return 1;
    
    double long result = 0;
    for (bool const d : BOOLS)
        result += pow(2.0l, bdd.getLayer(bdd.getChild(u, d)) - bdd.getLayer(u) - 1) * countSat(bdd, computedTable, bdd.getChild(u, d));
    computedTable.insert({u, result});
    return result;
}

/** compute the number of subsets represented by the root node of the BDD */
template<IsBDD BDD>
double long countSat(BDD const & bdd) {
    std::unordered_map<NodePtr, double long, NodePtrHash> computedTable;
    return countSat(bdd, computedTable, bdd.getRoot());
}

/** redirect arcs in a global bottom-up fashion by iterating over all layers and attempting all possible redirection steps */
template<IsStatefulBDD BDD>
void redirectArcs(BDD & bdd, AlgorithmSettings const & settings) {
    assert(settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_BottomUpDirectToExactState || settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_BottomUpDirectToRelaxedState);
    if (!getProfiler().currentNodeWithinTimeBudget(settings.arcRedirectionTimeBudget))
        return;
    std::unordered_map<NodePtr, double long, NodePtrHash> countSatCache;
    for (int l = bdd.getN() - 2; l >= 0; --l) {
        int const v = bdd.getVariableByLayer(l);
        int const nextV = bdd.getVariableByLayer(l + 1);
        std::vector<NodePtr> nextLayer = bdd.getNodesOnLayerInStateTable(nextV);
        for (NodePtr const s : bdd.getNodesOnLayerInStateTable(v)) {
            for (bool const b : BOOLS) {
                auto const exactState = bdd.decide(bdd.getState(s), v, b);
                // it is pointless to redirect arcs going into the infeasible state
                if (bdd.isInfeasible(exactState))
                    continue;
                auto const newDivergentNode = bdd.findInStateTable(nextV, exactState);
                NodePtr newC = bdd.getChild(s, b);
                if (newDivergentNode) { // if it is possible to make this arc exact, do so
                    newC = newDivergentNode.value();
                } else if (settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_BottomUpDirectToRelaxedState) {
                    // find "best" relaxing candidate based on sat count
                    double long minCountSat = countSat(bdd, countSatCache, newC);
                    for (NodePtr const u : nextLayer) {
                        // find a relaxing state that has less satisfying solutions than the off-path node
                        double long const newCountSat = countSat(bdd, countSatCache, u);
                        if (newCountSat < minCountSat && bdd.relaxes(bdd.getState(u), exactState)) {
                            std::cout << "found relaxed state with sat count " << newCountSat << " down from " << minCountSat << " new arc = (" << s << ", " << u << ", " << b << ")" << std::endl;
                            newC = u;
                            minCountSat = newCountSat;
                        }
                    }
                }
                // redirect the arc
                if (newC != bdd.getChild(s, b)) {
                    bdd.removeArc(s, bdd.getChild(s, b));
                    bdd.setChild(s, newC, b);
                    bdd.addArc(s, newC);
                }
            }
        }
    }
}



namespace {

/** computes the conjunction of functions represented by a and b using the given memoization computed table */
template<IsUniqueTableBDD BDD>
NodePtr conjunction(BDD & bdd, std::unordered_map<NodePtrPair, NodePtr, NodePtrPairHash> & computedTable, NodePtr const a, NodePtr const b) {
    NodePtrPair const computedTableKey = canonicallyOrderedPair(a, b); // exploit symmetry of conjunction
    auto const cachedResult = computedTable.find(computedTableKey);
    if (cachedResult != computedTable.end())
        return cachedResult->second;
    NodePtr u;

    // exploit trivial conjunction identities
    if (a == bdd.getTerminal(0) || b == bdd.getTerminal(0))
        return bdd.getTerminal(0);
    if (a == bdd.getTerminal(1))
        return b;
    if (b == bdd.getTerminal(1))
        return a;

    if (bdd.getLayer(a) == bdd.getLayer(b)) {
        int const v = bdd.getVariable(a);
        NodePtrPair const cA = bdd.getChildren(a);
        NodePtrPair const cB = bdd.getChildren(b);
        NodePtr const c0 = conjunction(bdd, computedTable, cA[0], cB[0]);
        NodePtr const c1 = conjunction(bdd, computedTable, cA[1], cB[1]);
        u = bdd.uniqueE(v, {c0, c1});
    } else {
        NodePtr high = a;
        NodePtr low = b;
        if (bdd.getLayer(a) > bdd.getLayer(b)) {
            high = b;
            low = a;
        }
        // layer of getLayer(high) < getLayer(low)
        int const v = bdd.getVariable(high);
        NodePtrPair const cH = bdd.getChildren(high);
        NodePtr const c0 = conjunction(bdd, computedTable, cH[0], low);
        NodePtr const c1 = conjunction(bdd, computedTable, cH[1], low);
        u = bdd.uniqueE(v, {c0, c1});
    }
    computedTable.emplace(computedTableKey, u);
    return u;
}

}



/** compute the conjunction of a and b using the apply algorithm
 *  this creates a new empty computed table dedicated to only this single operation
 */
template<IsUniqueTableBDD BDD>
NodePtr conjunction(BDD & bdd, NodePtr const a, NodePtr const b) {
    PROFILE(PROFILE_BDDConjunction);
    std::unordered_map<NodePtrPair, NodePtr, NodePtrPairHash> computedTable;
    return conjunction<BDD>(bdd, computedTable, a, b);
}

/** writes a DAG representation of a diagram to a .dot file */
template<IsBDD BDD>
void writeToDot(
    BDD const &             bdd,                            /** the BDD to represent */
    std::string const &     file,                           /** name of the file */
    bool const              showAllLayers = false,          /** whether to show layers that have no nodes in them */
    bool const              showT0 = false,                 /** whether to show the 0-terminal and its incident arcs */
    bool const              showState = true,               /** whether to show the associated state to the nodes of the diagram */
    std::unordered_map<Arc, double, ArcHash> const &
                            flow = {}                       /** map of arcs to flow in an associated flow model solution */
) {
    // perform a BFS to discover all nodes
    std::vector<std::vector<NodePtr>> layers(bdd.getN()); // per layer, the list of nodes that are on that layer
    std::vector<std::tuple<NodePtr, NodePtr, bool, double>> arcs; // all s-t arcs with label d and flow f represented as tuples (s, t, d, f)
    std::unordered_set<NodePtr, NodePtrHash> visited; // whether a node is visited in the BFS
    std::queue<NodePtr> Q;
    for (NodePtr s : bdd.getAliveNodes()) {
        if (visited.insert(s).second) {
            // start a BFS in every unexplored node
            Q.push(s);
            while (!Q.empty()) {
                NodePtr const u = Q.front();
                Q.pop();
                int const l = bdd.getLayer(u);
                if (l != bdd.getN()) {
                    if (l != -1)
                        layers[l].push_back(u);
                    for (bool const d : BOOLS) {
                        NodePtr const c = bdd.getChild(u, d);
                        arcs.push_back({u, c, d, flow.contains({u, d}) ? flow.at({u, d}) : 0.0});
                        if (visited.insert(c).second)
                            Q.push(c);
                    }
                }
            }
        }
    }
    std::ofstream out(file);
    out << "digraph \"BDD\" {\n";
    out << "    // layer structure\n";
    out << "    { node [shape = plaintext]; edge [style = invis]; ";
    for (int l = -1; l <= bdd.getN(); ++l) {
        if (showAllLayers || l == -1 || l == bdd.getN() || !layers[l].empty()) {
            if (l != -1)
                out << " -> ";
            out << "\"l" << l << "\"";
        }
    }
    out << "; }\n";
    out << "    // layer labels\n";
    out << "    ";
    for (int l = -1; l <= bdd.getN(); ++l) {
        if (showAllLayers || l == -1 || l == bdd.getN() || !layers[l].empty()) {
            if (l != -1)
                out << " ";
            out << "\"l" << l << "\"" << " [label = \"" << bdd.getVariableByLayer(l) << " (" << l << ")\"];";
        }
    }
    out << "\n";
    out << "    // layer -1\n";
    out << "    { rank = same; \"l-1\"; \"" << bdd.getRoot() << "\"; }\n";
    for (int l = 0; l < bdd.getN(); ++l) {
        if (showAllLayers || !layers[l].empty()) {
            out << "    // layer " << l << "\n";
            out << "    { rank = same; \"l" << l << "\";";
            for (NodePtr const u : layers[l])
                out << " \"" << u << "\";";
            out << " }\n";
        }
    }
    out << "    // layer " << bdd.getN() << "\n";
    out << "    { rank = same; node [shape = box]; \"l" << bdd.getN() << "\";";
    if (showT0)
        out << " \"0\";";
    out << " \"1\"; }\n";
    if constexpr (IsStatefulBDD<BDD>) {
        bool firstDisplay = true;
        for (NodePtr const u : visited) {
            if ((u == bdd.getTerminal(0) && !showT0))
                continue;
            std::string const display = bdd.getState(u).display();
            if (display != "") {
                if (firstDisplay) {
                    out << "    // node labels\n";
                    out << "   ";
                    firstDisplay = false;
                }
                out << " \"" << u << "\"" << " [label = \"" << u << " " << display << "\"];";
            }
        }
        if (!firstDisplay)
            out << "\n";
    }
    out << "    // arcs\n";
    std::ostringstream namebuf(std::ios_base::ate);
    for (std::tuple<NodePtr, NodePtr, bool, double> const & arc : arcs) {
        if (!showT0 && std::get<1>(arc) == bdd.getTerminal(0))
            continue;
        out << "    \"" << std::get<0>(arc) << "\" -> \"" << std::get<1>(arc) << "\"";
        std::vector<std::pair<std::string, std::string>> properties;
        if (!std::get<2>(arc))
            properties.push_back({"style", "dashed"});
        if (std::get<3>(arc)) {
            properties.push_back({"color", "green"});
            namebuf.str("\"");
            namebuf << std::get<3>(arc) << "\"";
            properties.push_back({"label", namebuf.str()});
        }
        if (!properties.empty()) {
            out << " [";
            for (size_t i = 0; i < properties.size(); ++i) {
                if (i != 0)
                    out << "; ";
                out << properties[i].first << " = " << properties[i].second;
            }
            out << "]";
        }
        out << ";\n";
    }
    out << "}";
    out.flush();
}

/** extract the conflicting subpath associated to the provided conflict partial solution if it is still present in the BDD */
template<IsBDD BDD>
std::optional<ConflictPath> extractConflictPath(BDD const & bdd, ConflictSolution const & conflict) {
    ConflictPath conflictPath{{{}, {}}, conflict.edgeLayers};
    NodePtr u = bdd.getRoot();
    int l;
    while ((l = bdd.getLayer(u)) <= conflict.edgeLayers.back()) {
        bool const d = l == -1 ? 0 : conflict.partialSolution[l];
        NodePtr c = bdd.getChild(u, d);
        if (bdd.getLayer(c) > conflict.edgeLayers[0]) {
            conflictPath.path.nodes.push_back(u);
            conflictPath.path.arcLabels.push_back(d);
        }
        u = c;
    }

    if (u == bdd.getTerminal(0))
        return {}; // conflict no longer exists in bdd because the partial solution is not accepted
    conflictPath.path.nodes.push_back(u); // add the last node on the path

    assert( std::ranges::is_sorted(conflictPath.edgeLayers) );
    assert( std::ranges::is_sorted(conflictPath.path.nodes, std::less(), [&bdd](const NodePtr w){ return bdd.getLayer(w); }) );
    assert( bdd.getLayer(conflictPath.path.nodes[0]) <= conflictPath.edgeLayers[0] && conflict.edgeLayers[0] < bdd.getLayer(conflictPath.path.nodes[1]) );
    assert( bdd.getLayer(conflictPath.path.nodes[conflictPath.path.length() - 1]) <= conflictPath.edgeLayers.back() && conflictPath.edgeLayers.back() < bdd.getLayer(conflictPath.path.nodes.back()) );

    return conflictPath;
}