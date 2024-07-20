#pragma once

#include "hypergraph.h"
#include <optional>

/** represents an ordinary graph of which the edges are known up front
 *  the graph is simple and undirected
 */
struct OrdinaryGraph : public Hypergraph {

    std::vector<VertexList>         adj;            /** adjacency lists */
    std::vector<std::vector<bool>>  m;              /** adjacency matrix */

    OrdinaryGraph(Hypergraph const && H);

};

/** reads an ordinary graph from a DIMACS formatted file
 *  every cardinality 2 edge is sorted on increasing vertex index
 */
OrdinaryGraph readDIMACSGraphFromFile(std::string const & file);

/** write an ordinary graph into DIMACS format */
void writeToDIMACSFile(OrdinaryGraph const & G, std::string const & file);

/** implementation of independent set that compares against an explicitly enumerated set of edges
 *  specialized for ordinary graphs for improved performance
 */
class OrdinaryGraphIndependentSet : public IndependentSet {

private:

    OrdinaryGraph const &           G;              /** graph */
    VertexList                      points;         /** current set of points */
    std::optional<VertexList>       conflict;       /** conflict if detected during an add invocation */

    OrdinaryGraphIndependentSet(OrdinaryGraphIndependentSet const &) = delete;
    OrdinaryGraphIndependentSet(OrdinaryGraphIndependentSet &&) = delete;
    OrdinaryGraphIndependentSet & operator=(OrdinaryGraphIndependentSet const &) = delete;
    OrdinaryGraphIndependentSet & operator=(OrdinaryGraphIndependentSet &&) = delete;

public:

    /** construct a new independent set compared against the the edges of newG */
    OrdinaryGraphIndependentSet(OrdinaryGraph const & newG);
    ~OrdinaryGraphIndependentSet();
    void add(int i);
    void removeLastAdded();
    bool isIndependent() const;
    VertexList findContainedCircuit() const;
    VertexList const & getElements() const;

};

/** find a maximal clique that contains the given vertex */
VertexList findMaximalCliqueContainingVertex(OrdinaryGraph const & G, int v);

/** find a maximal clique for every vertex in G and return a vector containing size[v]
 *  where size[v] is the size of a maximal clique in G that contains v
 */
std::vector<size_t> findMaximalCliqueSizes(OrdinaryGraph const & G);