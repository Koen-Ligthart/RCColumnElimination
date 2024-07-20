#include <cassert>
#include <fstream>
#include <ranges>
#include <set>
#include "ordinary_graph.h"
#include "profiler.h"

OrdinaryGraph::OrdinaryGraph(Hypergraph const && H) : Hypergraph(H), adj(n), m(n, std::vector<bool>(n, false)) {
    for (VertexList const & e : H.E) {
        assert( e.size() == 2 );
        int const a = e[0];
        int const b = e[1];
        // update adjacency matrix
        m[a][b] = true;
        m[b][a] = true;
        // update adjacency list
        adj[a].push_back(b);
        adj[b].push_back(a);
    }
}

OrdinaryGraph readDIMACSGraphFromFile(std::string const & file) {
    std::ifstream in(file);
    if (!in)
        throw std::invalid_argument("file could not be found");
    
    std::string token;
    int n;
    size_t m = std::numeric_limits<size_t>::max();
    size_t edgesRead = 0;
    std::vector<VertexList> E;
    std::set<std::pair<int, int>> edgeSet;

    while (edgesRead < m) {
        in >> token;
        if (token == "c") { // skip a comment
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else if (token == "p") { // graph size information
            in >> token;
            assert( token == "edge" );
            in >> n >> m;
        } else { // edge
            assert (token == "e");
            int v1;
            int v2;
            in >> v1 >> v2;
            int const l = std::min(v1, v2) - 1;
            int const u = std::max(v1, v2) - 1;
            if (!edgeSet.contains({l, u})) { // prevent duplicate edges
                edgeSet.insert({l, u});
                E.push_back({l, u});
            }
            ++edgesRead;
        }
    }
    return OrdinaryGraph({n, E});
}

void writeToDIMACSFile(OrdinaryGraph const & G, std::string const & file) {
    std::ofstream out(file);
    out << "p edge " << G.n << " " << G.E.size() << "\n";
    for (const VertexList& E : G.E)
        out << "e " << (std::min(E[0], E[1]) + 1) << " " << (std::max(E[0], E[1]) + 1) << "\n";
    out.flush();
}

OrdinaryGraphIndependentSet::OrdinaryGraphIndependentSet(OrdinaryGraph const & newG) : G{ newG } {}

OrdinaryGraphIndependentSet::~OrdinaryGraphIndependentSet() {}

void OrdinaryGraphIndependentSet::add(int const i) {
    PROFILE(PROFILE_IndependenceOracle);
    for (int const j : std::ranges::views::reverse(points)) {
        if (G.m[j][i]) {
            conflict = {j, i};
            break;
        }
    }
    points.push_back(i);
}

void OrdinaryGraphIndependentSet::removeLastAdded() {
    PROFILE(PROFILE_IndependenceOracle);
    points.pop_back();
    conflict = {};
}

bool OrdinaryGraphIndependentSet::isIndependent() const {
    PROFILE(PROFILE_IndependenceOracle);
    return !conflict;
}

VertexList OrdinaryGraphIndependentSet::findContainedCircuit() const {
    PROFILE(PROFILE_IndependenceOracle);
    return conflict.value();
}

const VertexList& OrdinaryGraphIndependentSet::getElements() const {
    return points;
}

VertexList findMaximalCliqueContainingVertex(OrdinaryGraph const & G, int const v) {
    VertexList clique{v};
    for (int w = 0; w < G.n; ++w) {
        // check if w can greedily be added to the clique, v is not added twice because G has no self-loops
        bool canBeAdded = true;
        for (int const z : clique) {
            if (!G.m[w][z]) { // w is not incident to z
                canBeAdded = false;
                break;
            }
        }
        if (canBeAdded)
            clique.push_back(w);
    }
    return clique;
}

std::vector<size_t> findMaximalCliqueSizes(OrdinaryGraph const & G) {
    std::vector<size_t> cliqueSizes(G.n, 1);
    for (int v = 0; v < G.n; ++v) {
        VertexList clique = findMaximalCliqueContainingVertex(G, v);
        // update sizes for all elements of the clique if this clique is larger than the known clique size
        for (int const w : clique)
            cliqueSizes[w] = std::max(cliqueSizes[w], clique.size());
    }
    return cliqueSizes;
}