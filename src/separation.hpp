#pragma once

#include <algorithm>
#include <ranges>
#include "bdd_util.hpp"

/** locally remove the conflict from the diagram using the local separation procedure */
template<IsUniqueTableBDD BDD>
requires IsParentListedNode<class BDD::NodeType>
void separateLocal(BDD & bdd, ConflictPath const & conflict) {
    std::cout << "separating";
    for (int const l : conflict.edgeLayers)
        std::cout << " " << l;
    std::cout << std::endl;

    NodePtr const u0 = conflict.path.nodes[0];

    // stage 1: separation
    NodePtr u = bdd.getTerminal(0);
    int iP = conflict.path.length() - 1;
    int iE = conflict.edgeLayers.size() - 1;
    while (iP >= 0) {
        NodePtr c0;
        NodePtr c1;
        const int lP = bdd.getLayer(conflict.path.nodes[iP]);
        int l;
        if (iE < 0 || conflict.edgeLayers[iE] <= lP) {
            // node corresponds to a choice on the conflicting path
            l = lP;
            if (conflict.path.arcLabels[iP]) {
                c0 = bdd.getChild(conflict.path.nodes[iP], 0);
                c1 = u;
            } else {
                c0 = u;
                c1 = bdd.getChild(conflict.path.nodes[iP], 1);
            }
            if (iE >= 0 && conflict.edgeLayers[iE] == lP)
                --iE;
            --iP;
        } else {
            // node corresponds to a choice that is not made on the conflicting path but corresponds to a vertex of the conflicting edge
            l = conflict.edgeLayers[iE];
            c0 = conflict.path.nodes[iP + 1];
            c1 = u;
            --iE;
        }
        u = bdd.uniqueE(bdd.getVariableByLayer(l), {c0, c1});
    }

    // stage 2: merging and garbage collection
    auto const maxLayer = [&bdd](NodePtr const a, NodePtr const b){ return bdd.getLayer(a) < bdd.getLayer(b); };
    std::priority_queue<NodePtr, std::vector<NodePtr>, decltype(maxLayer)> Q(maxLayer); // might be possible to do with topological sort
    std::unordered_set<NodePtr, NodePtrHash> inQ;
    std::vector<NodePtr> gcQ;
    assert((bdd.getRoot() == u0) <= (bdd.getLayer(u0) == bdd.getLayer(u))); // if u0 is the root, then we cannot have that u is on a lower layer
    bdd.removeFromUniqueTable(u0);
    inQ.insert(u0); // guaranteed to return true
    Q.push(u0);
    // in-place modify u0 to match u
    for (bool const b : BOOLS) {
        NodePtr const oldChild = bdd.getChild(u0, b);
        bdd.removeArc(u0, oldChild);
        bdd.setChild(u0, bdd.getChild(u, b), b);
        bdd.addArc(u0, bdd.getChild(u, b));
        if (bdd.hasZeroIndegree(oldChild))
            gcQ.push_back(oldChild);
    }
    bdd.setVariable(u0, bdd.getVariable(u));
    // garbage collect
    if (bdd.hasZeroIndegree(u) && std::ranges::find(gcQ, u) == gcQ.end()) // required check as u may be one of the children of u0
        gcQ.push_back(u);
    garbageCollect(bdd, gcQ);

    // stage 3: cascading unification
    while (!Q.empty()) {
        NodePtr const c = Q.top();
        Q.pop();
        std::optional<NodePtr> newC;
        if (bdd.getChild(c, 0) == bdd.getChild(c, 1)) {
            newC = bdd.getChild(c, 0);
        } else {
            newC = bdd.findInUniqueTable(bdd.getVariable(c), bdd.getChildren(c));
        }
        if (newC) {
            for (ParentListIterator it = bdd.getParents(c).begin(); it != bdd.getParents(c).end(); ++it) {
                NodePtr const p = *it;
                // redirect arcs from c to newC
                for (bool const d : BOOLS) {
                    if (bdd.getChild(p, d) == c) {
                        if (inQ.insert(p).second) {
                            bdd.removeFromUniqueTable(p);
                            Q.push(p);
                        }
                        it.remove();
                        bdd.setChild(p, newC.value(), d);
                        bdd.addArc(p, newC.value());
                    }
                }
            }
            // remove c
            bdd.removeNode(c);
            bdd.clearDeadInParentList(bdd.getChild(c, 0));
            bdd.clearDeadInParentList(bdd.getChild(c, 1));
            assert( !bdd.hasZeroIndegree(bdd.getChild(c, 0)) );
            assert( !bdd.hasZeroIndegree(bdd.getChild(c, 1)) );
        } else {
            bdd.insertIntoUniqueTable(c);
        }
    }
}

/** globally remove the conflict from the diagram by
 *  taking the conjunction with the constraint induced by the edge and the current root
 *  to form the new root of the diagram
 */
template<IsUniqueTableBDD BDD>
requires IsDisconnectedCheckableNode<class BDD::NodeType>
void separateGlobal(BDD & bdd, ConflictPath const & conflict) {
    std::cout << "separating";
    for (int const l : conflict.edgeLayers)
        std::cout << " " << l;
    std::cout << std::endl;

    // construct BDD structure that forbids the conflicting edge
    NodePtr cons = bdd.getTerminal(0);
    for (int const l : std::ranges::views::reverse(conflict.edgeLayers))
        cons = bdd.uniqueE(bdd.getVariableByLayer(l), {bdd.getTerminal(1), cons});
    
    // compute conjunction and update BDD root
    NodePtr const newRoot = conjunction(bdd, bdd.getRoot(), cons);
    assert(newRoot != bdd.getRoot());
    std::vector<NodePtr> gcQ = {bdd.getRoot()};
    if (bdd.hasZeroIndegree(cons))
        gcQ.push_back(cons);
    garbageCollect(bdd, gcQ);
    bdd.setRoot(newRoot);
}

template<IsStatefulBDD BDD>
void separateLocalStateful(BDD & bdd, ConflictPath const & conflict, AlgorithmSettings const & settings) {
    std::cout << "separating";
    for (int const l : conflict.edgeLayers)
        std::cout << " " << l;
    std::cout << std::endl;

    // compute exact states top-down
    std::vector<class BDD::NodeState> newStates;
    newStates.emplace_back(bdd.getState(conflict.path.nodes[0]));
    for (size_t i = 0; i < conflict.path.length() - 1; ++i) // last state need not be generated as it is the trivial infeasible state associated to the 0-terminal
        newStates.emplace_back(bdd.decide(newStates.back(), bdd.getVariable(conflict.path.nodes[i]), conflict.path.arcLabels[i]));
    NodePtr nextOnPath = bdd.getTerminal(0);

    std::unordered_map<NodePtr, double long, NodePtrHash> countSatCache;
    for (long i = conflict.path.length() - 1; i >= 0; --i) {
        // redirect arc to become exact
        int const v = bdd.getVariable(conflict.path.nodes[i]);
        bool const l = conflict.path.arcLabels[i];
        std::optional<NodePtr> newU = bdd.findInStateTable(v, newStates[i]);
        if (newU) { // node with target state already exists
            NodePtr const u = newU.value();
            bdd.removeArc(u, bdd.getChild(u, l));
            bdd.setChild(u, nextOnPath, l);
            bdd.addArc(u, nextOnPath);
            nextOnPath = u;
        } else { // create a new node with target state
            NodePtrPair c;
            c[!l] = bdd.getChild(conflict.path.nodes[i], !l);
            auto const stateDivergent = bdd.decide(newStates[i], v, !l);
            if (bdd.isInfeasible(stateDivergent)) {
                // diverging arc may be redirected to 0-terminal
                c[!l] = bdd.getTerminal(0);
            } else {
                if (settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_DuringSplitDirectToExactState || settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_DuringSplitDirectToRelaxedState) {
                    PROFILE(PROFILE_RedirectArcs);
                    const int nextV = bdd.getVariable(conflict.path.nodes[i + 1]);
                    if (nextV < bdd.getN()) {
                        auto const newDivergentNode = bdd.findInStateTable(nextV, stateDivergent);
                        if (newDivergentNode) { // if it is possible to make this arc exact, do so
                            c[!l] = newDivergentNode.value();
                        } else if (settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_BottomUpDirectToRelaxedState) {
                            // find "best" relaxing candidate based on sat count
                            double long minCountSat = countSat(bdd, countSatCache, c[!l]);
                            for (NodePtr const u : bdd.getNodesOnLayerInStateTable(nextV)) {
                                // find a relaxing state that has less satisfying solutions than the off-path node
                                if (bdd.relaxes(bdd.getState(u), stateDivergent)) {
                                    double long const newCountSat = countSat(bdd, countSatCache, u);
                                    if (newCountSat < minCountSat) {
                                        std::cout << "found relaxed state with sat count " << newCountSat << " down from " << minCountSat << std::endl;
                                        c[!l] = u;
                                        minCountSat = newCountSat;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            c[l] = nextOnPath;
            // construct the new node
            NodePtr const u = bdd.emplaceNode(v, c);
            bdd.setState(u, std::move(newStates[i]));
            bdd.insertIntoStateTable(u);
            bdd.addArc(u, c[0]);
            bdd.addArc(u, c[1]);
            nextOnPath = u;
        }
        newStates.pop_back(); // free this state
    }
    assert(newStates.empty());
}