#pragma once

#include "bdd_base.hpp"

/** BDD that manages a unique table keeping track of all (v, c0, c1) combinations to ensure semantic uniqueness of its nodes */
template<IsNode Node>
class UniqueTableBDD : public BDDBase<Node> {

private:

    std::vector<std::unordered_map<NodePtrPair, NodePtr, NodePtrPairHash>>
                                uniqueTable;                    /** unique table of the diagram
                                                                 *  the v-th hash table contains (c0, c1) -> node mappings so that the nodes
                                                                 *  on variable v with children (c0, c1) are contained in that way in the table
                                                                 *  all alive nodes except the super root and terminals are represented in the table
                                                                 */

    /** construct a new BDD with a unique table for the given number of variables
     *  together with the two terminal nodes
     *  the root is left undefined
     */
    UniqueTableBDD(int const newN, std::vector<int> const & newPi, std::vector<int> const & newPiInv) :
        BDDBase<Node>(newN, newPi, newPiInv),
        uniqueTable(this->n, std::unordered_map<NodePtrPair, NodePtr, NodePtrPairHash>())
    {}

    /** removes all nonterminal nodes from the bdd and sets root to the 0-terminal */
    void setToZero() {
        this->setRoot(this->getTerminal(0));
        for (NodePtr u = 0; u.index < this->nodes.size(); u = u.index + 1) {
            if (!this->isTerminal(u) && !this->nodes[u.index].isDead()) {
                for (NodePtr const c : this->getChildren(u))
                    this->removeArc(u, c);
                removeFromUniqueTable(u);
                this->removeNode(u);
            }
        }
    }

    /** copies entire subtree of bdd rooted in o to this bdd and returns the corresponding root, reduces the diagram while doing so */
    template<IsBDD BDD>
    NodePtr copy(std::unordered_map<NodePtr, NodePtr, NodePtrHash> & copyMap, BDD const & bdd, NodePtr const o) {
        if (o == bdd.getTerminal(0))
            return this->getTerminal(0);
        if (o == bdd.getTerminal(1))
            return this->getTerminal(1);
        auto const result = copyMap.find(o);
        if (result == copyMap.end()) {
            NodePtr const c0 = copy(copyMap, bdd, bdd.getChild(o, 0));
            NodePtr const c1 = copy(copyMap, bdd, bdd.getChild(o, 1));
            NodePtr const u = uniqueE(bdd.getVariable(o), {c0, c1});
            copyMap.emplace(o, u);
            return u;
        } else {
            return (*result).second;
        }
    }

public:

    /** construct a new BDD with a unique table for the given number of variables
     *  together with the two terminal nodes and
     *  a super root at layer -1 with 0-child being the 1-terminal and 1-child being the 0-terminal
     */
    UniqueTableBDD(int const newN, std::vector<int> const & newPi) :
        BDDBase<Node>(newN, newPi),
        uniqueTable(this->n, std::unordered_map<NodePtrPair, NodePtr, NodePtrPairHash>())
    {}

    /** copy construct a BDD from another BDD and reduce it while doing so */
    template<IsBDD BDD>
    UniqueTableBDD(const BDD &bdd) : UniqueTableBDD(bdd.n, bdd.pi, bdd.piInv) { (*this) = bdd; }

    /** copy the given BDD and reduce it while doing so */
    template<IsBDD BDD>
    UniqueTableBDD<Node> & operator=(BDD const & bdd)
    {
        assert(this->n == bdd.n);
        setToZero();
        this->pi = bdd.pi;
        this->piInv = bdd.piInv;
        std::unordered_map<NodePtr, NodePtr, NodePtrHash> copyMap;
        this->setRoot(copy(copyMap, bdd, bdd.getRoot()));
        return *this;
    }

    /** insert a node into the unique table if it is not at layer -1 */
    void insertIntoUniqueTable(NodePtr const u) {
        int const v = this->getVariable(u);
        if (v == -1) return;
        assert(v < this->n);
        uniqueTable[v].emplace(this->getChildren(u), u);
    }

    /** remove an entry of the unique table corresponding to this node if it is not at layer -1 */
    void removeFromUniqueTable(NodePtr const u) {
        int const v = this->getVariable(u);
        if (v == -1) return;
        assert(v < this->n);
        uniqueTable[v].erase(this->getChildren(u));
    }

    /** returns the node stored at the given key in the unique table if it exists */
    std::optional<NodePtr> findInUniqueTable(int const v, NodePtrPair const c) {
        if (v == -1) return {};
        assert(v < this->n);
        auto const u = uniqueTable[v].find(c);
        if (u == uniqueTable[v].end()) return {};
        return (*u).second;
    }

    /** returns a unique node with a given semantic specification
     *  if both children are equal, it returns that child
     *  if a node with the given variable and children exists, that node is returned
     *  otherwise, this creates a new node with the given attributes and inserts the node into the unique table (if it is not at layer -1)
     *  and properly registers the new arcs
     */
    NodePtr uniqueE(int const v, NodePtrPair const c) {
        assert(this->getLayerByVariable(v) < this->getLayer(c[0]) && this->getLayerByVariable(v) < this->getLayer(c[1]));
        if (c[0] == c[1])
            return c[0];
        std::optional<NodePtr> const u = findInUniqueTable(v, c);
        if (u) {
            assert(this->getVariable(u.value()) == v && this->getChildren(u.value()) == c);
            return u.value();
        } else {
            NodePtr const result = this->emplaceNode(v, c);
            insertIntoUniqueTable(result);
            this->addArc(result, c[0]);
            this->addArc(result, c[1]);
            return result;
        }
    }

    /** locally swap the layers lLow and lLow + 1 and
     *  modifies the diagram to abide by the new variable order and maintain the same boolean function representation
     */
    void swapLayers(int const lLow) requires IsDisconnectedCheckableNode<Node> {
        assert( lLow >= 0 && lLow < this->n - 1 );

        // perform swap in pi and piInv
        int const vLow = this->getVariableByLayer(lLow);
        int const vHigh = this->getVariableByLayer(lLow + 1);
        std::swap(this->pi[lLow], this->pi[lLow + 1]);
        std::swap(this->piInv[vLow], this->piInv[vHigh]);

        std::vector<NodePtr> layer; // nodes associated with variable vLow
        layer.reserve(uniqueTable[vLow].size());
        for (auto const [_, u] : uniqueTable[vLow])
            layer.push_back(u);
        uniqueTable[vLow].clear();

        // readd nodes that do not depend on vHigh first (if they are not inserted first, equivalent nodes may have been added before them)
        for (NodePtr const u : layer)
            if (this->getVariable(this->getChild(u, 0)) != vHigh && this->getVariable(this->getChild(u, 1)) != vHigh)
                insertIntoUniqueTable(u);
        
        // move the remaining nodes on layer a layer downwards
        for (NodePtr const u : layer) {
            NodePtr const u0 = this->getChild(u, 0);
            NodePtr const u1 = this->getChild(u, 1);
            if (this->getVariable(u0) != vHigh && this->getVariable(u1) != vHigh)
                continue; // we already treated these nodes
            NodePtr const u00 = this->getVariable(u0) == vHigh ? this->getChild(u0, 0) : u0;
            NodePtr const u01 = this->getVariable(u0) == vHigh ? this->getChild(u0, 1) : u0;
            NodePtr const u10 = this->getVariable(u1) == vHigh ? this->getChild(u1, 0) : u1;
            NodePtr const u11 = this->getVariable(u1) == vHigh ? this->getChild(u1, 1) : u1;
            // garbage collect children with variable vHigh if possible
            for (NodePtr const c : {u0, u1}) {
                this->removeArc(u, c);
                if (this->hasZeroIndegree(c) && this->getVariable(c) == vHigh) { // if the child does not have variable vHigh, it will become a grandchild soon
                    removeFromUniqueTable(c);
                    this->removeNode(c);
                    // update grandchildren inneighbours
                    for (const NodePtr gc : this->getChildren(c)) {
                        this->removeArc(c, gc);
                        // it is guaranteed that grandchildren have references after finishing creating new children and relinking
                        // so we should not delete them here
                    }
                }
            }
            this->nodes[u.index].v = vHigh;
            NodePtr const newC0 = uniqueE(vLow, {u00, u10});
            NodePtr const newC1 = uniqueE(vLow, {u01, u11});
            this->setChild(u, newC0, 0);
            this->addArc(u, newC0);
            this->setChild(u, newC1, 1);
            this->addArc(u, newC1);
            insertIntoUniqueTable(u);
        }
    }

};

template<class T>
concept IsUniqueTableBDD = std::is_same_v<UniqueTableBDD<class T::NodeType>, T>;