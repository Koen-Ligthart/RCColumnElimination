#pragma once

#include "primitives.h"

/** represents an explicit undirected hypergraph of which the edges are known up front */
struct Hypergraph {
    int                             n;              /** the number of vertices of the hypergraph */
    std::vector<VertexList>         E;              /** edges of the hypergraph */
};

/** very inefficient implementation of independent sets that compares against an explicitly enumerated set of edges */
class ExplicitIndependentSet : public IndependentSet {

private:

    std::vector<VertexList> const & E;              /** hypergraph */
    VertexList                      points;         /** current set of points */

    ExplicitIndependentSet(ExplicitIndependentSet const &) = delete;
    ExplicitIndependentSet(ExplicitIndependentSet &&) = delete;
    ExplicitIndependentSet & operator=(ExplicitIndependentSet const &) = delete;
    ExplicitIndependentSet & operator=(ExplicitIndependentSet &&) = delete;

public:

    /** construct a new independent set compared against the explicit list of edges newE
     *  all edges are assumed to be sorted on vertex index in increasing order
     *  the edges are assumed to form a clutter
     */
    ExplicitIndependentSet(std::vector<VertexList> const & newE);
    ~ExplicitIndependentSet();
    void add(int i);
    void removeLastAdded();
    bool isIndependent() const;
    VertexList findContainedCircuit() const;
    VertexList const & getElements() const;

};