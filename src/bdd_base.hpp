#pragma once

#include <cassert>
#include <functional>
#include <optional>
#include <queue>
#include <unordered_map>
#include <variant>
#include "profiler.h"

#define NULL_INDEX (size_t) -1

constexpr bool BOOLS[2] = {0, 1};

class NodeBase;

class NodePtr;

/** pair of NodePtrs */
typedef std::array<NodePtr, 2> NodePtrPair;

template<class T>
concept IsNode = std::is_base_of_v<NodeBase, T> && requires(int const v, NodePtrPair const c) {
    T(v, c);
};

template<class T>
concept IsParentListedNode = IsNode<T> && requires(T t) {
    { t.parentListIndex } -> std::same_as<size_t &>;
};

template<class T>
concept IsReferenceCountedNode = IsNode<T> && requires(T t) {
    { t.referenceCount } -> std::same_as<int &>;
};

template<class T>
concept IsDisconnectedCheckableNode = IsParentListedNode<T> || IsReferenceCountedNode<T>;

template<class T>
concept IsNodeState = std::regular<T> && requires(T const t) {
    { t.display() } -> std::same_as<std::string>;
    { t.hash() } -> std::same_as<size_t>;
};

template<class T>
concept IsNodeStateRef = IsNodeState<std::remove_reference_t<T>>;

template<class T>
concept IsStatefulNode = IsNode<T> && requires(T t, int n) {
    { t.state } -> IsNodeStateRef;
};

template<class T, class Node>
concept IsNodeStateManager = requires(
    T t,
    decltype(Node::state) const state,
    VertexList const E,
    T const constT, 
    int const v,
    bool const d,
    decltype(Node::state) const stateO,
    std::string const file,
    std::vector<std::pair<NodePtr, decltype(Node::state) *>> const states
) {
    { t.create0TerminalState() } -> std::same_as<decltype(Node::state)>;
    { t.create1TerminalState() } -> std::same_as<decltype(Node::state)>;
    { t.createRootState() } -> std::same_as<decltype(Node::state)>;
    { constT.isInfeasible(state) } -> std::same_as<bool>;
    { t.decide(state, v, d) } -> std::same_as<decltype(Node::state)>;
    { t.relaxes(state, stateO) } -> std::same_as<bool>;
    { t.writeState(file, states) } -> std::same_as<void>;
};

struct NodePtrPairHash;

/** pointer to a node in a BDD */
class NodePtr {

    template<IsNode Node>
    friend class BDDBase;
    template<IsNode Node>
    friend class UniqueTableBDD;
    template<IsStatefulNode Node, class NodeStateManager>
    requires IsNodeStateManager<NodeStateManager, Node>
    friend class StatefulBDD;
    template<IsNode Node>
    friend class ParentListIterator;
    friend class NodePtrHash;
    friend class ArcHash;
    friend class NodePtrPairHash;
    friend std::ostream & operator<<(std::ostream & out, NodePtr const & u) { out << u.index; return out; };
    friend NodePtrPair canonicallyOrderedPair(NodePtr const a, NodePtr const b) { return {std::min(a.index, b.index), std::max(a.index, b.index) }; }

private:

    size_t index;
    NodePtr(size_t const newIndex) : index{ newIndex } {}

public:

    NodePtr() : index{ NULL_INDEX } {}

    bool operator==(NodePtr const & o) const { return index == o.index; }

    bool operator<(NodePtr const & o) const { return index < o.index; }

};

/** hash function for NodePtr */
struct NodePtrHash {
    size_t operator()(NodePtr const & u) const { return u.index; }
};

/** an arc in a BDD given by its unique source and label */
struct Arc {
    NodePtr s;
    bool label;
    bool operator==(Arc const & o) const { return s == o.s && label == o.label; }
};

/** hash function for Arc */
struct ArcHash {
    size_t operator()(Arc const & a) const { return a.label ? -a.s.index : a.s.index; }
};

/** a path in a BDD given by it nodes and labels of connecting arcs */
struct Path {

    std::vector<NodePtr>        nodes;                          /** nodes[i - 1] is the i-th vertex on the path */
    std::vector<bool>           arcLabels;                      /** arcLabels[i] is the label of the arc connecting nodes[i] and nodes[i + 1] */

    /** the number of arcs on the path */
    size_t length() const { return arcLabels.size(); }

    /** returns the length (end - start) path formed by vertices nodes[start], nodes[start + 1], ..., nodes[end]
     *  and arcs arcLabels[start], arcLabels[start + 1], ..., arcLabels[end - 1]
     */
    Path subPath(size_t const start, size_t const end) const {
        return {
            {nodes.begin() + start, nodes.begin() + (end + 1)},
            {arcLabels.begin() + start, arcLabels.begin() + end}
        };
    }

};

/** represents a partial solution of a DD that conflictswith an edge */
struct ConflictSolution {
    std::vector<bool>           partialSolution;                /** partial binary vector solution indexed by layers that defines
                                                                 *  a conflicting path starting from the root
                                                                 *  it is guaranteed to end with a 1
                                                                 */
    std::vector<int>            edgeLayers;                     /** layers corresponding to the variables of the edge in ascending order */
};

/** represents a path in a DD that conflicts with an edge */
struct ConflictPath {
    Path                        path;                           /** path of which the arcs define a minimal path that contains the conflicting edge */
    std::vector<int>            edgeLayers;                     /** layers corresponding to the variables of the edge in ascending order */
};

/** hash function for NodePtrPair */
struct NodePtrPairHash {
    size_t operator()(NodePtrPair const & p) const {
        return 37 * p[0].index + p[1].index;
    }
};

/** base implementation for a BDD node */
class NodeBase {

    template<IsNode Node>
    friend class BDDBase;
    template<IsNode Node>
    friend class UniqueTableBDD;
    template<IsStatefulNode Node, class NodeStateManager>
    requires IsNodeStateManager<NodeStateManager, Node>
    friend class StatefulBDD;

private:

    int                         v;                              /** variable associated with this node
                                                                 *  this is n iff the node is a terminal,
                                                                 *  -1 iff the node is the super root
                                                                 *  and -2 iff the node is dead
                                                                 */
    NodePtrPair                 c;                              /** c[i] points to the i-child of this node
                                                                 *  these are NULL_INDEX if the node is a terminal
                                                                 */

    bool isDead() const { return v == -2; }

protected:

    void setDead() { v = -2; } // node states must be cleared when this is called

public:

    NodeBase(int const newV, NodePtrPair const newC) : v{ newV }, c{ newC } {}

};

/** internal linked list node used to store the in-neighbours of BDD nodes */
struct ParentListNode {
    NodePtr                     u;                              /** node represented by this linked list node */
    size_t                      next;                           /** location of next linked list node or NULL_INDEX if this is the last */
};

template<IsNode Node>
class BDDBase;

/** iterator used to iterate over the collection of parents of a DD node that is able to remove elements while iterating */
template<IsNode Node>
class ParentListIterator {
    
    friend class BDDBase<Node>;

private:

    BDDBase<Node> &             bdd;                            /** reference to the owning BDD */
    NodePtr                     u;                              /** the child node whose parents are being iterated over */
    size_t                      previousP;                      /** pointer to previous linked list node or NULL_INDEX if this is the first node */
    size_t                      p;                              /** pointer to the current linked list node */

    /** constructs an iterator over the parent list of a given node */
    ParentListIterator(BDDBase<Node> & newBdd, NodePtr const newU, size_t const newPreviousP, size_t const newP) :
        bdd{ newBdd }, u{ newU }, previousP{ newPreviousP }, p{ newP } {}

public:

    /** returns this iterator itself for compatibility with for loops */
    ParentListIterator<Node> begin() const { return { bdd, u, previousP, p }; }

    /** returns the iterator corresponding to the same node that points at the end of the linked list */
    ParentListIterator<Node> end() const { return { bdd, u, previousP, NULL_INDEX }; }

    /** move the iterator to the next element of the parent list */
    void operator++() {
        if (p == NULL_INDEX) {
            assert(p == previousP);
            p = bdd.nodes[u.index].parentListIndex;
        } else {
            previousP = p;
            p = bdd.parents[p].next;
        }
    }

    /** access the element stored at the current position of the list */
    NodePtr operator*() const { return bdd.parents[p].u; }

    /** compares this iterator against another iterator
     *  this should only be used to compare with the end of list iterator corresponding to the same DD node
     *  other comparisons are undefined
     */
    bool operator!=(ParentListIterator<Node> const & o) const { return p != o.p; }

    /** remove the current element that is being pointed towards
     *  after this, dereferencing this iterator is undefined
     *  instead one should increment this iterator
     */
    void remove() {
        ParentListNode &parent = bdd.parents[p];
        if (previousP == NULL_INDEX) bdd.nodes[u.index].parentListIndex = parent.next;
        else bdd.parents[previousP].next = parent.next;
        bdd.freeParents.push_back(p);
        p = previousP;
    }

};

/** base BDD implementation */
template<IsNode Node>
class BDDBase {

    friend class ParentListIterator<Node>;
    template<IsNode Node2>
    friend class UniqueTableBDD;
    template<IsStatefulNode Node2, class NodeStateManager>
    requires IsNodeStateManager<NodeStateManager, Node2>
    friend class StatefulBDD;

protected:

    int const                   n;                              /** number of variables of this BDD excluding the super root */
    NodePtr                     root{0};                        /** pointer to the current root of the diagram
                                                                 *  which is the unique node with v = -1 or the 0 terminal
                                                                 *  in which case the diagram consist of only the terminal nodes
                                                                 */
    std::vector<int>            pi;                             /** permutation on {0, 1, ..., n - 1} mapping layers to variables */
    std::vector<int>            piInv;                          /** permutation on {0, 1, ..., n - 1} mapping variables to layers */
    std::vector<Node>           nodes;                          /** container of all nodes of the diagram (including dead nodes)
                                                                 *  the 0 and 1 terminals are stored at indices 0 and 1 respectively
                                                                 */
    std::vector<NodePtr>        freeNodes;                      /** indices i such that nodes[i].isDead() */



    [[no_unique_address]] std::conditional_t<IsParentListedNode<Node>, std::vector<ParentListNode>, std::monostate>
                                parents;                        /** collection of linked list nodes forming all parent lists of all nodes */
    [[no_unique_address]] std::conditional_t<IsParentListedNode<Node>, std::vector<size_t>, std::monostate>
                                freeParents;                    /** indices pointing towards parent list nodes that are unused */



    /** construct a new BDD for the given number of variables
     *  together with the two terminal nodes
     *  the root is left undefined
     */
    BDDBase(int const newN, std::vector<int> const & newPi, std::vector<int> const & newPiInv) :
        n{ newN },
        pi{ newPi },
        piInv{ newPiInv }
    {
        nodes.push_back({n, {NULL_INDEX, NULL_INDEX}});
        nodes.push_back({n, {NULL_INDEX, NULL_INDEX}});
    }

    void populateTrivialDiagram() {
        // create root
        nodes.push_back({-1, {1, 0}});
        root = 2;
        addArc(getRoot(), getTerminal(0));
        addArc(getRoot(), getTerminal(1));
    }

public:

    typedef Node NodeType;



    /** construct a new BDD for the given number of variables
     *  together with the two terminal nodes and
     *  a super root at layer -1 with 0-child being the 1-terminal and 1-child being the 0-terminal
     */
    BDDBase(int const newN, std::vector<int> const & newPi) : BDDBase(newN, newPi, {})
    {
        // initialize identity variable order if unspecified
        if (pi.empty()) {
            pi.reserve(n);
            for (int v = 0; v < n; ++v)
                pi.push_back(v);
        }

        // derive inverse order
        piInv.resize(n);
        for (int l = 0; l < n; ++l)
            piInv[pi[l]] = l;

        populateTrivialDiagram();
    }

    

    /** returns the number of variables in the diagram, not counting the super root node */
    int getN() const { return n; }

    /** returns the root which is the unique vertex at layer -1 or the 0-terminal */
    NodePtr getRoot() const { return root; }

    /** returns the d-terminal */
    NodePtr getTerminal(const bool d) const { return d; }

    /** returns the number of alive nodes in this diagram */
    size_t nodeCount() const { return nodes.size() - freeNodes.size(); }

    /** get the variable corresponding to a layer */
    int getVariableByLayer(const int l) const { return (l != -1 && l != n) ? pi[l] : l; }

    /** get the layer corresponding to a variable */
    int getLayerByVariable(const int v) const { return (v != -1 && v != n) ? piInv[v] : v; }

    /** get the layer of a node */
    int getLayer(NodePtr const u) const { return getLayerByVariable(nodes[u.index].v); }

    /** get the variable of a node */
    int getVariable(NodePtr const u) const { return nodes[u.index].v; }

    /** get the children of a node */
    NodePtrPair getChildren(NodePtr const u) const { return nodes[u.index].c; }

    /** get the d-child of a node */
    NodePtr getChild(NodePtr const u, bool const d) const { return nodes[u.index].c[d]; }

    /** returns whether u is a terminal */
    bool isTerminal(NodePtr const u) const { return u.index == 0 || u.index == 1; }



    /** returns whether a node has zero indegree */
    bool hasZeroIndegree(NodePtr const t) const requires IsParentListedNode<Node> { return nodes[t.index].parentListIndex == NULL_INDEX; }

    /** returns whether a node has zero indegree */
    bool hasZeroIndegree(const NodePtr t) const requires IsReferenceCountedNode<Node> { return nodes[t.index].referenceCount == 0; }



    /** obtains an index in the nodes array and ensures that it becomes an alive node with the given variable and children
     *  this does not touch the unique table nor parent lists
     */
    NodePtr emplaceNode(int const v, NodePtrPair const c) {
        if (freeNodes.empty()) {
            // allocate new node
            NodePtr const u = nodes.size();
            nodes.emplace_back(v, c);
            assert( nodes[u.index].v != -1 || nodes[u.index].c != nodes[root.index].c );
            return u;
        } else {
            // reuse existing node
            NodePtr const u = freeNodes.back();
            freeNodes.pop_back();
            if constexpr (IsDisconnectedCheckableNode<Node>)
                assert( hasZeroIndegree(u) );
            assert( nodes[u.index].isDead() );
            nodes[u.index] = Node(v, c);
            assert( nodes[u.index].v != -1 || nodes[u.index].c != nodes[root.index].c );
            return u;
        }
    }

    /** marks the given node as dead and registers the position in the node array to be free */
    void removeNode(NodePtr const u) {
        if constexpr (IsDisconnectedCheckableNode<Node>)
            assert( getRoot() == getTerminal(0) || hasZeroIndegree(u) ); // do not check if root is being set to 0 due to setToZero clearing of diagram
        nodes[u.index].setDead();
        freeNodes.push_back(u);
    }

    /** in place modifies the d-child pointer of p to point to c
     *  this does not modify the unique table nor parent lists nor reference counts
     */
    void setChild(NodePtr const p, NodePtr const c, bool d) {
        assert( getLayer(p) < getLayer(c) );
        nodes[p.index].c[d] = c;
    }

    /** in place modifies the variable of u */
    void setVariable(NodePtr const u, int v) {
        nodes[u.index].v = v;
    }

    /** set the root of the diagram to a given node */
    void setRoot(NodePtr const u) {
        assert( getLayer(u) == -1 || u == getTerminal(0) );
        root = u;
    }



    void addArc(NodePtr const s, NodePtr const t) requires (!IsDisconnectedCheckableNode<Node>) {}

    void removeArc(NodePtr const s, NodePtr const t) requires (!IsDisconnectedCheckableNode<Node>) {}



    /** increase the reference count of u by 1 */
    void addArc(NodePtr const s, NodePtr const t) requires IsReferenceCountedNode<Node> { ++nodes[t.index].referenceCount; }

    /** decrease the reference count of u by 1 */
    void removeArc(NodePtr const s, NodePtr const t) requires IsReferenceCountedNode<Node> { --nodes[t.index].referenceCount; }



    /** obtain an iterator to iterate over the parents of t */
    ParentListIterator<Node> getParents(const NodePtr t) requires IsParentListedNode<Node> { return {*this, t, NULL_INDEX, nodes[t.index].parentListIndex}; }

    /** add s to the parent list of t */
    void addArc(NodePtr const s, NodePtr const t) requires IsParentListedNode<Node> {
        size_t newParent;
        if (freeParents.empty()) {
            // allocate new linked list node
            newParent = parents.size();
            parents.push_back({});
        } else {
            // reuse existing linked list node
            newParent = freeParents.back();
            freeParents.pop_back();
        }
        // insert into linked list
        ParentListNode & parent = parents[newParent];
        Node & node = nodes[t.index];
        parent.u = s;
        parent.next = node.parentListIndex;
        node.parentListIndex = newParent;
    }

    /** remove one occurence of s from the parent list of t */
    void removeArc(NodePtr const s, NodePtr const t) requires IsParentListedNode<Node> {
        for (auto it = getParents(t).begin(); it != getParents(t).end(); ++it) {
            if (*it == s) {
                it.remove();
                break;
            }
        }
    }

    /** remove all parents that are dead in the parent list of u
     *  returns whether the node has zero indegree after this operation
     *  no new nodes may be created between removing the intended parents
     *  and calling this function
     */
    void clearDeadInParentList(const NodePtr u) requires IsParentListedNode<Node> {
        for (auto it = getParents(u).begin(); it != getParents(u).end(); ++it)
            if (nodes[(*it).index].isDead())
                it.remove();
    }



    /** get a list of all nodes alive in the diagram */
    std::vector<NodePtr> getAliveNodes() const {
        std::vector<NodePtr> result;
        result.reserve(nodeCount());
        for (NodePtr u = 0; u.index < nodes.size(); u = u.index + 1)
            if (!nodes[u.index].isDead())
                result.push_back(u);
        return result;
    }

};

template<class T>
concept IsBDD = std::is_base_of_v<BDDBase<class T::NodeType>, T>;

/** base implementation of a BDD node with associated parent list */
struct SimpleParentListedNode : NodeBase {
    size_t parentListIndex = NULL_INDEX;
    SimpleParentListedNode(int const newV, NodePtrPair const newC) : NodeBase(newV, newC) {}
};

/** base implementation of a BDD reference counted node */
struct SimpleReferenceCountedNode : NodeBase {
    int referenceCount = 0;
    SimpleReferenceCountedNode(int const newV, NodePtrPair const newC) : NodeBase(newV, newC) {}
};