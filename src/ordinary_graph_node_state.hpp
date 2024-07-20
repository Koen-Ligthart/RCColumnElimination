#pragma once

#include <sstream>
#include "bdd_stateful.hpp"
#include "ordinary_graph.h"

/** state related to a partial independent set in an ordinary graph
 *  states are represented by the set of vertices that can still be added to this partial independent set
 */
struct OrdinaryGraphNodeState {

    std::vector<bool>       eligibleVertices;   /** eligibleVertices[v] indicates whether v can still be added
                                                 *  to the partial independent set represented by this state
                                                 */

    /** create a trivial state and minimal memory usage */
    OrdinaryGraphNodeState() {};

    /** create a state given a set of eligible vertices */
    OrdinaryGraphNodeState(std::vector<bool> && newEligibleVertices) : eligibleVertices{ newEligibleVertices } {};

    /** creates a textual representation of this state */
    std::string display() const {
        if (eligibleVertices.empty())
            return "";
        std::ostringstream buf;
        buf << "{";
        bool first = true;
        for (size_t i = 0; i < eligibleVertices.size(); ++i) {
            if (eligibleVertices[i]) {
                if (first) {
                    first = false;
                } else {
                    buf << ",";
                }
                buf << i;
            }
        }
        buf << "}";
        return buf.str();
    }

    bool operator==(OrdinaryGraphNodeState const & state) const { return eligibleVertices == state.eligibleVertices; }

    size_t hash() const { return std::hash<std::vector<bool>>()(eligibleVertices); };

};

class OrdinaryGraphNodeStateManager {

private:

    OrdinaryGraph const &   G;                  /** associated graph */

    OrdinaryGraphNodeStateManager(OrdinaryGraphNodeStateManager const &) = delete;
    OrdinaryGraphNodeStateManager(OrdinaryGraphNodeStateManager &&) = delete;
    OrdinaryGraphNodeStateManager & operator=(OrdinaryGraphNodeStateManager const&) = delete;
    OrdinaryGraphNodeStateManager & operator=(OrdinaryGraphNodeStateManager const &&) = delete;

public:

    OrdinaryGraphNodeStateManager(OrdinaryGraph const & newG) : G{ newG } {};

    OrdinaryGraphNodeState create0TerminalState() { return {}; }
    OrdinaryGraphNodeState create1TerminalState() { return std::vector<bool>(G.n, false); }
    OrdinaryGraphNodeState createRootState() { return std::vector<bool>(G.n, true); }

    /** register the existence an edge */
    template<IsStatefulBDD BDD>
    requires std::is_same_v<class BDD::NodeStateManagerType, OrdinaryGraphNodeStateManager>
    void registerEdge(BDD & bdd, VertexList const & E) {
        // all edges are known up-front so this does not do anything
    }

    /** whether this state corresponds to a dependent set */
    bool isInfeasible(OrdinaryGraphNodeState const & state) const { return state.eligibleVertices.empty(); }

    /** computes the logically following state after deciding that v = d starting at the given state */
    OrdinaryGraphNodeState decide(OrdinaryGraphNodeState const & state, int const v, bool const d) {
        OrdinaryGraphNodeState result(state);
        if (v == -1)
            return result;
        if (d && !result.eligibleVertices[v])
            return {}; // results in infeasible state, i.e. cannot select this vertex
        result.eligibleVertices[v] = false;
        if (d) // if choosing to include v, remove adjacent vertices from eligible set
            for (int const neighbour : G.adj[v])
                result.eligibleVertices[neighbour] = false;
        return result;
    }

    /** returns whether relaxed relaxes strengthened */
    bool relaxes(OrdinaryGraphNodeState const & relaxed, OrdinaryGraphNodeState const & strengthened) {
        if (isInfeasible(strengthened))
            return true;
        if (isInfeasible(relaxed))
            return false;
        // checks whether relaxed is a super set of strengthened
        for (int i = 0; i < G.n; ++i)
            if (strengthened.eligibleVertices[i] && !relaxed.eligibleVertices[i])
                return false;
        return true;
    }

    /** to write a representation of the state of all nodes to a file */
    void writeState(std::string const & file, std::vector<std::pair<NodePtr, OrdinaryGraphNodeState *>> const & states) {
        // states are already visualized in display()
    }

};

/** BDD node that additionally has an associated state in terms of a list of eligible vertices */
struct OrdinaryGraphStatefulNode : public NodeBase {

    OrdinaryGraphNodeState  state;

    OrdinaryGraphStatefulNode(int const newV, NodePtrPair const newC) : NodeBase(newV, newC) {}

    void setDead() {
        NodeBase::setDead();
        state.eligibleVertices.resize(0);
    }

};