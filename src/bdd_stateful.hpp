#pragma once

#include "bdd_base.hpp"
#include <unordered_set>

template<IsStatefulNode Node, class NodeStateManager>
requires IsNodeStateManager<NodeStateManager, Node>
class StatefulBDD;

/** equality comparator for node states managed by a state manager owned by a BDD */
template<IsStatefulNode Node, class NodeStateManager>
requires IsNodeStateManager<NodeStateManager, Node>
struct NodeStateEqual {

    using is_transparent = void;

    StatefulBDD<Node, NodeStateManager> const *
                                bdd;                            /** reference to BDD that owns the state manager owning the states */

    NodeStateEqual(StatefulBDD<Node, NodeStateManager> const * newBdd) : bdd{ newBdd } {}

    bool operator()(NodePtr const a, NodePtr const b) const { return bdd->getState(a) == bdd->getState(b); }
    bool operator()(NodePtr const a, StatefulBDD<Node, NodeStateManager>::NodeState const &b) const { return bdd->getState(a) == b; }
    bool operator()(StatefulBDD<Node, NodeStateManager>::NodeState const & a, NodePtr const b) const { return a == bdd->getState(b); }
    bool operator()(StatefulBDD<Node, NodeStateManager>::NodeState const & a, StatefulBDD<Node, NodeStateManager>::NodeState const & b) const { return a == b; }
    
};

/** hash for node states managed by a state manager owned by a BDD */
template<IsStatefulNode Node, class NodeStateManager>
requires IsNodeStateManager<NodeStateManager, Node>
struct NodeStateHash {

    using is_transparent = void;

    StatefulBDD<Node, NodeStateManager> const *
                                bdd;                            /** reference to BDD that owns the state manager owning the states */

    NodeStateHash(StatefulBDD<Node, NodeStateManager> const * newBdd) : bdd{ newBdd } {}

    size_t operator()(NodePtr const u) const { return bdd->getState(u).hash(); }
    size_t operator()(StatefulBDD<Node, NodeStateManager>::NodeState const & state) const { return state.hash(); };
    
};

template<IsStatefulNode Node, class NodeStateManager>
requires IsNodeStateManager<NodeStateManager, Node>
class StatefulBDD;

template<class T>
concept IsStatefulBDD = std::is_same_v<StatefulBDD<class T::NodeType, class T::NodeStateManagerType>, T>;

template<class T, class Node, class OtherBDD>
concept CanCopyStatesFrom = IsNodeStateManager<T, Node> && IsStatefulBDD<OtherBDD> && requires(
    T t,
    int const v,
    OtherBDD::NodeState const state
) {
    { t.convert(v, state) } -> std::same_as<decltype(Node::state)>;
};

/** BDD that associates every node with a unique state */
template<IsStatefulNode Node, class NodeStateManager>
requires IsNodeStateManager<NodeStateManager, Node>
class StatefulBDD : public BDDBase<Node> {

private:

    NodeStateManager            stateManager;                   /** manages the node states by providing operations on node states */

    std::vector<std::unordered_set<NodePtr, NodeStateHash<Node, NodeStateManager>, NodeStateEqual<Node, NodeStateManager>>>
                                stateTable;                     /** table by variable that maps a state to the unique node associated to that variable and the given state
                                                                 *  all alive nodes except the super root and terminals are represented in the table
                                                                 */

    /** construct a new BDD with a unique table for the given number of variables
     *  together with the two terminal nodes with state determined by the state manager,
     *  which is constructed using the specified arguments
     *  the root is left undefined
     */
    template<class ... NodeStateManagerConstructorArgs>
    StatefulBDD(int const newN, std::vector<int> const & newPi, std::vector<int> const & newPiInv, NodeStateManagerConstructorArgs && ... stateManagerConstructorArgs) :
        BDDBase<Node>(newN, newPi, newPiInv),
        stateManager(stateManagerConstructorArgs...),
        stateTable(this->n, std::unordered_set<NodePtr, NodeStateHash<Node, NodeStateManager>, NodeStateEqual<Node, NodeStateManager>>(0, NodeStateHash<Node, NodeStateManager>(this), NodeStateEqual<Node, NodeStateManager>(this)))
    {
        this->nodes[0].state = stateManager.create0TerminalState();
        this->nodes[1].state = stateManager.create1TerminalState();
    }

    /** copies entire subtree of bdd rooted in o to this bdd and returns the corresponding root
     *  node states are converted by the state manager of this BDD
     */
    template<IsStatefulBDD BDD>
    requires CanCopyStatesFrom<NodeStateManager, Node, BDD>
    NodePtr copy(std::unordered_map<NodePtr, NodePtr, NodePtrHash> & copyMap, BDD const & bdd, NodePtr const o) {
        if (o == bdd.getTerminal(0))
            return this->getTerminal(0);
        if (o == bdd.getTerminal(1))
            return this->getTerminal(1);
        auto result = copyMap.find(o);
        if (result == copyMap.end()) {
            NodePtr const c0 = copy(copyMap, bdd, bdd.getChild(o, 0));
            NodePtr const c1 = copy(copyMap, bdd, bdd.getChild(o, 1));
            NodePtr const u = this->emplaceNode(bdd.getVariable(o), {c0, c1});
            this->addArc(u, c0);
            this->addArc(u, c1);
            copyMap.emplace(o, u);
            setState(u, stateManager.convert(bdd.getVariable(o), bdd.getState(o)));
            insertIntoStateTable(u);
            return u;
        } else {
            return (*result).second;
        }
    }

public:

    typedef decltype(Node::state) NodeState;
    typedef NodeStateManager NodeStateManagerType;



    /** construct a new BDD with a unique table for the given number of variables
     *  together with the two terminal nodes and root with state determined by the state manager,
     *  which is constructed using the specified arguments
     *  the 0-child and 1-child of the root are respectively 1- and 0-terminal
     */
    template<class ... NodeStateManagerConstructorArgs>
    StatefulBDD(int const newN, std::vector<int> const & newPi, NodeStateManagerConstructorArgs && ... stateManagerConstructorArgs) :
        BDDBase<Node>(newN, newPi),
        stateManager(stateManagerConstructorArgs...),
        stateTable(this->n, std::unordered_set<NodePtr, NodeStateHash<Node, NodeStateManager>, NodeStateEqual<Node, NodeStateManager>>(0, NodeStateHash<Node, NodeStateManager>(this), NodeStateEqual<Node, NodeStateManager>(this)))
    {
        this->nodes[0].state = stateManager.create0TerminalState();
        this->nodes[1].state = stateManager.create1TerminalState();
        this->removeArc(this->getRoot(), this->getTerminal(1));
        this->nodes[2].state = stateManager.createRootState();
        for (int l = 0; l < this->n; ++l) {
            // create relaxed node at layer l
            NodePtr const origin = l + 2;
            NodePtr const newNode = l + 3;
            this->emplaceNode(this->getVariableByLayer(l), {this->getTerminal(1), this->getTerminal(1)});
            this->nodes[newNode.index].state = stateManager.decide(this->getState(origin), this->getVariableByLayer(l - 1), 0);
            if (l > 0) {
                this->setChild(origin, newNode, 1);
                this->addArc(origin, newNode);
            }
            this->setChild(origin, newNode, 0);
            this->addArc(origin, newNode);
            insertIntoStateTable(newNode);
        }
        this->addArc(this->n + 2, this->getTerminal(1));
        this->addArc(this->n + 2, this->getTerminal(1));
    }

    /** copy construct a BDD from another BDD
     *  node states are converted by the state manager of this BDD
     *  which is constructed with the specified arguments
     */
    template<IsStatefulBDD BDD, class ... NodeStateManagerConstructorArgs>
    requires CanCopyStatesFrom<NodeStateManager, Node, BDD>
    StatefulBDD(BDD const & bdd, NodeStateManagerConstructorArgs && ... stateManagerConstructorArgs) :
        StatefulBDD(bdd.n, bdd.pi, bdd.piInv, stateManagerConstructorArgs...)
    {
        std::unordered_map<NodePtr, NodePtr, NodePtrHash> copyMap;
        this->setRoot(copy(copyMap, bdd, bdd.getRoot()));
    }

    ~StatefulBDD() {
        for (NodePtr u = 0; u.index < this->nodes.size(); u = u.index + 1)
            if (!this->nodes[u.index].isDead())
                this->nodes[u.index].setDead();
    }

    /** insert a node into the state table if it is not at layer -1 */
    void insertIntoStateTable(NodePtr const u) {
        int const v = this->getVariable(u);
        if (v == -1) return;
        assert(v < this->n);
        stateTable[v].emplace(u);
    }

    /** remove an entry of the state table corresponding to this node if it is not at layer -1 */
    void removeFromStateTable(NodePtr const u) {
        int const v = this->getVariable(u);
        if (v == -1) return;
        assert(v < this->n);
        stateTable[v].erase(u);
    }

    /** returns the node stored with the given state and variable if it exists */
    std::optional<NodePtr> findInStateTable(int const v, NodeState const & state) const {
        assert(v != -1);
        if (v == this->n) {
            if (isInfeasible(state))
                return this->getTerminal(0);
            return this->getTerminal(1);
        }
        auto const u = stateTable[v].find(state);
        if (u == stateTable[v].end()) return {};
        return *u;
    }

    /** get a list of all alive nodes that are associated with v */
    std::vector<NodePtr> getNodesOnLayerInStateTable(int const v) const {
        assert(v >= 0 && v < this->n);
        std::vector<NodePtr> result;
        for (NodePtr const u : stateTable[v])
            result.push_back(u);
        std::sort(result.begin(), result.end()); // required for determinism
        return result;
    }

    /** set the state of a node */
    void setState(NodePtr const u, NodeState const && state) { this->nodes[u.index].state = state; }

    /** get the state of a node */
    NodeState const & getState(NodePtr const u) const { return this->nodes[u.index].state; }

    /** invokes the state manager to register a new hyperedge */
    void registerEdge(VertexList const & E) { return stateManager.registerEdge(*this, E); }

    /** returns whether the state is infeasible */
    bool isInfeasible(NodeState const & state) const { return stateManager.isInfeasible(state); }

    /** computes the logically following state after deciding that v = d starting at the given state */
    NodeState decide(NodeState const & state, int const v, bool const d) { return stateManager.decide(state, v, d); }

    /** returns whether relaxed relaxes strengthened */
    bool relaxes(NodeState const & relaxed, NodeState const & strengthened) { return stateManager.relaxes(relaxed, strengthened); }

    /** invokes the state manager to write a representation of its state to a file */
    void writeStateManagerState(std::string const & file) {
        std::vector<std::pair<NodePtr, NodeState *>> states;
        for (NodePtr u = 0; u.index < this->nodes.size(); u = u.index + 1)
            if (!this->nodes[u.index].isDead())
                states.emplace_back(u, &(this->nodes[u.index].state));
        stateManager.writeState(file, states);
    }

};