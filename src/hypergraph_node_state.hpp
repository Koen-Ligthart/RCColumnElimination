#pragma once

#include <algorithm>
#include <ranges>
#include "cuddInt.h"
#include "ordinary_graph_node_state.hpp"
#include "zdd_util.h"

/** state related to a partial independent set in a hypergraph
 *  states are represented by a ZDD node that represents a set of residual edges
 */
struct HypergraphNodeState {

    DdManager *             zdd;                /** ZDD manager of w */
    DdNode *                w;                  /** associated ZDD node */

    HypergraphNodeState(DdManager * const newZdd, DdNode * const newW) : zdd{ newZdd }, w{ newW } {
        if (w != NULL)
            Cudd_Ref(w);
    };

    /** create a trivial state and minimal memory usage */
    HypergraphNodeState() : HypergraphNodeState(NULL, NULL) {};

    HypergraphNodeState(HypergraphNodeState const & state) : HypergraphNodeState(state.zdd, state.w) {}

    HypergraphNodeState & operator=(HypergraphNodeState const & state) {
        if (state.w != NULL)
            Cudd_Ref(state.w);
        if (w != NULL)
            Cudd_RecursiveDerefZdd(zdd, w);
        zdd = state.zdd;
        w = state.w;
        return *this;
    }

    ~HypergraphNodeState() {
        if (w != NULL)
            Cudd_RecursiveDerefZdd(zdd, w);
    }

    /** creates a textual representation of this state */
    std::string display() const {
        return ""; // representation is visualized through dot files in manager writeState()
    }

    bool operator==(HypergraphNodeState const & state) const { return w == state.w; }

    size_t hash() const { return std::hash<DdNode *>()(w); };

};

class HypergraphNodeStateManager {

private:

    int const               n;                  /** number of variables in the ZDD */
    DdManager *             zdd;                /** ZDD manager */
    std::vector<DdNode *>   edgeCollections;    /** edgeCollections[v] is the ZDD node representing
                                                 *  the known edges anchored in v, i.e. {hyperedge E : min E = v}
                                                 */

    HypergraphNodeStateManager(HypergraphNodeStateManager const &) = delete;
    HypergraphNodeStateManager(HypergraphNodeStateManager &&) = delete;
    HypergraphNodeStateManager & operator=(HypergraphNodeStateManager const&) = delete;
    HypergraphNodeStateManager & operator=(HypergraphNodeStateManager const &&) = delete;

public:

    /** create a manager for hypergraph node states and initialize a ZDD with the specified variable order */
    HypergraphNodeStateManager(int const newN, std::vector<int> const & variableOrder) : n{ newN } {
        zdd = Cudd_Init(0, n, CUDD_UNIQUE_SLOTS, CUDD_CACHE_SLOTS, 0);
        assert( Cudd_GarbageCollectionEnabled(zdd) ); // have garbage collecting enabled
        #ifndef NDEBUG
        Cudd_ReorderingType method;
        #endif
        assert( !Cudd_ReorderingStatus(zdd, &method) ); // no automatic dynamic reordering
        assert( !Cudd_ReorderingStatusZdd(zdd, &method) ); // no automatic dynamic reordering
        assert((size_t) n == variableOrder.size());
        int * const order = new int[variableOrder.size()];
        edgeCollections.reserve(n);
        for (int i = 0; i < n; ++i) {
            order[i] = variableOrder[i];
            edgeCollections.push_back(Cudd_ReadZero(zdd));
            Cudd_Ref(Cudd_ReadZero(zdd));
        }
        Cudd_zddShuffleHeap(zdd, order);
        delete[] order;
    }

    ~HypergraphNodeStateManager() {
        for (DdNode * const edgeCollection : edgeCollections)
            Cudd_RecursiveDerefZdd(zdd, edgeCollection);
        Cudd_Quit(zdd);
    }

    HypergraphNodeState create0TerminalState() { return {zdd, DD_ONE(zdd)}; }
    HypergraphNodeState create1TerminalState() { return {zdd, Cudd_ReadZero(zdd)}; }
    HypergraphNodeState createRootState() { return {zdd, Cudd_ReadZero(zdd)}; }

    /** adds a residual edge to the residual edge collection of the given state */
    HypergraphNodeState introduce(HypergraphNodeState const & state, DdNode * residual) {
        DdNode * const newState = Cudd_zddUnion(zdd, state.w, residual);
        zddVerifyNonNull(zdd, newState);
        Cudd_Ref(newState);
        DdNode * const trimmed = zddMinMal(zdd, newState);
        HypergraphNodeState const result = {zdd, trimmed};
        Cudd_RecursiveDerefZdd(zdd, newState);
        return result;
    }

    /** register the existence an edge by updating the respective edgeCollection
     *  additionally this tries to introduce the corresponding residual to many states if possible
     */
    template<IsStatefulBDD BDD>
    requires std::is_same_v<class BDD::NodeStateManagerType, HypergraphNodeStateManager> && IsParentListedNode<class BDD::NodeType>
    void registerEdge(BDD & bdd, VertexList const & E) {
        // remove min E from E to obtain the initial residual
        VertexList sortedE = E;
        std::ranges::sort(sortedE, [this](int const a, int const b){ return Cudd_ReadPermZdd(zdd, a) < Cudd_ReadPermZdd(zdd, b); });
        int const minV = sortedE[0];
        VertexList element;
        for (int const v : sortedE)
            if (v != minV)
                element.push_back(v);
        DdNode* newResidual = zddSingleton(zdd, element);
        Cudd_Ref(newResidual);

        // register the edge by adding it to the respective edgeCollection
        DdNode* newEdgeCollection = Cudd_zddUnion(zdd, edgeCollections[minV], newResidual);
        zddVerifyNonNull(zdd, newEdgeCollection);
        Cudd_Ref(newEdgeCollection);
        if (newEdgeCollection == edgeCollections[minV]) {
            // this did not have an actual effect as the edge was already registered
            Cudd_RecursiveDerefZdd(zdd, newEdgeCollection);
        } else {
            // we added a new edge
            Cudd_RecursiveDerefZdd(zdd, edgeCollections[minV]);
            edgeCollections[minV] = newEdgeCollection;

            // attempt to introduce it to as many states as possible
            
            std::unordered_set<NodePtr, NodePtrHash> candidateNodes;
            for (NodePtr const u : bdd.getNodesOnLayerInStateTable(minV))
                if (bdd.getChild(u, 1) != bdd.getTerminal(0))
                    candidateNodes.insert(bdd.getChild(u, 1));
            
            size_t i = 1; // least index i s.t. l <= pi(sortedE[i])
            DdNode * residual = newResidual;
            for (int l = bdd.getLayerByVariable(minV) + 1; i < sortedE.size(); ++l) {
                assert((int) Cudd_NodeReadIndex(residual) == sortedE[i]);

                int const v = bdd.getVariableByLayer(l);
                bool require1Arc = false;

                if (v == sortedE[i]) {
                    ++i;
                    require1Arc = true; // only states completely including this edge should be considered for updating
                    residual = Cudd_T(residual); // update residual edge to become smaller
                }

                std::vector<NodePtr> candidateNodesSorted(candidateNodes.begin(), candidateNodes.end());
                std::sort(candidateNodesSorted.begin(), candidateNodesSorted.end());
                std::unordered_set<NodePtr, NodePtrHash> newCandidateNodes;
                for (NodePtr const u : candidateNodesSorted) {
                    HypergraphNodeState const newState = introduce(bdd.getState(u), residual);
                    // ensure that we actually update the state and do not create two nodes with the same state in the same layer
                    if (newState == bdd.getState(u) || bdd.findInStateTable(v, newState))
                        continue;
                    bool updateState = true;
                    for (bool const b : BOOLS) {
                        for (NodePtr const p : bdd.getParents(u)) {
                            // this order of iteration takes care of the case where p is a double parent
                            if (bdd.getChild(p, b) == u) {
                                if (!bdd.relaxes(newState, bdd.decide(bdd.getState(p), bdd.getVariable(p), b))) {
                                    updateState = false;
                                    break;
                                }
                            }
                        }
                        if (!updateState)
                            break;
                    }
                    if (updateState) {
                        // add child nodes that could potentially have the new edge be introduced to the list
                        for (bool const b : BOOLS)
                            if ((b == 1 || !require1Arc) && bdd.getChild(u, b) != bdd.getTerminal(0))
                                newCandidateNodes.insert(bdd.getChild(u, b));

                        // modify state of this node
                        bdd.removeFromStateTable(u);
                        if (isInfeasible(newState)) {
                            // if this node is infeasible, delete it and redirect all incoming arcs to the 0 terminal
                            for (bool const b : BOOLS) {
                                for (ParentListIterator it = bdd.getParents(u).begin(); it != bdd.getParents(u).end(); ++it) {
                                    NodePtr const p = *it;
                                    if (bdd.getChild(p, b) == u) {
                                        it.remove();
                                        bdd.setChild(p, bdd.getTerminal(0), b);
                                        bdd.addArc(p, bdd.getTerminal(0));
                                    }
                                }
                                bdd.removeArc(u, bdd.getChild(u, b));
                            }
                            bdd.removeNode(u);
                        } else {
                            // update the state of this node
                            bdd.setState(u, std::move(newState));
                            bdd.insertIntoStateTable(u);
                        }
                    }
                }
                candidateNodes = newCandidateNodes;
            }
        }
        Cudd_RecursiveDerefZdd(zdd, newResidual);
    }

    /** whether this state corresponds to a dependent set */
    bool isInfeasible(HypergraphNodeState const & state) const { return state.w == DD_ONE(zdd); }

    /** computes the logically following state after deciding that v = d starting at the given state */
    HypergraphNodeState decide(HypergraphNodeState const & state, const int v, const bool d) {
        if (v == -1) {
            assert(d == 0);
            return state;
        }
        if (d) {
            // descend to union of 0- and 1-child if this variable is essential
            DdNode * const existingResiduals = Cudd_NodeReadIndex(state.w) == (int unsigned) v ? Cudd_zddUnion(zdd, Cudd_E(state.w), Cudd_T(state.w)) : state.w;
            zddVerifyNonNull(zdd, existingResiduals);
            Cudd_Ref(existingResiduals);
            // add all residuals anchored in v
            DdNode * const newState = Cudd_zddUnion(zdd, existingResiduals, edgeCollections[v]);
            zddVerifyNonNull(zdd, newState);
            Cudd_Ref(newState);
            Cudd_RecursiveDerefZdd(zdd, existingResiduals);
            // compute collection of minimal residuals
            DdNode * const trimmed = zddMinMal(zdd, newState);
            HypergraphNodeState result = {zdd, trimmed};
            Cudd_RecursiveDerefZdd(zdd, newState);
            return result;
        } else {
            // descend to 0-child if this variable is essential
            if (Cudd_NodeReadIndex(state.w) == (int unsigned) v)
                return {zdd, Cudd_E(state.w)};
            assert( Cudd_ReadPermZdd(zdd, Cudd_NodeReadIndex(state.w)) > Cudd_ReadPermZdd(zdd, v) );
            return {zdd, state.w};
        }
    }

    /** returns whether relaxed relaxes strengthened */
    bool relaxes(HypergraphNodeState const & relaxed, HypergraphNodeState const & strengthened) {
        return zddWeakIncl(zdd, relaxed.w, strengthened.w);
    }

    /** convert an set of eligible vertices to set of forbidden vertices
     *  in the form of a collection of singleton residual edges
     */
    HypergraphNodeState convert(int const v, OrdinaryGraphNodeState const & state) {
        VertexList singletons;
        for (int l = v == -1 ? 0 : Cudd_ReadPermZdd(zdd, v); l < n; ++l) {
            int const i = Cudd_ReadInvPermZdd(zdd, l);
            if (!state.eligibleVertices[i])
                singletons.push_back(i);
        }
        return {zdd, zddSingletons(zdd, singletons)};
    }

    /** to write a representation of the state of all nodes to a file */
    void writeState(std::string const & file, std::vector<std::pair<NodePtr, HypergraphNodeState *>> const & states) {
        std::vector<std::pair<std::string, DdNode *>> labelledNodes;
        std::ostringstream namebuf{std::ios_base::ate};
        // label all pointers from BDD nodes to ZDD nodes
        for (const auto &[u, state] : states) {
            namebuf.str("u");
            namebuf << u;
            labelledNodes.emplace_back(namebuf.str(), state->w);
        }
        // label all edgeCollections
        for (int v = 0; v < n; ++v) {
            namebuf.str("c");
            namebuf << v;
            labelledNodes.emplace_back(namebuf.str(), edgeCollections[v]);
        }
        writeZDDToDot(n, zdd, file, labelledNodes);
    }

};

/** BDD node that additionally has an associated state in terms of a pointer to a ZDD node
 *  which represents a set of residual edges
 *  additionally, it contains a field for parent list tracking
 */
struct HypergraphStatefulNode : public NodeBase {

    HypergraphNodeState     state;
    size_t                  parentListIndex = NULL_INDEX;

    HypergraphStatefulNode(int const newV, NodePtrPair const newC) : NodeBase(newV, newC) {}

    void setDead() {
        NodeBase::setDead();
        state = {};
    }

};