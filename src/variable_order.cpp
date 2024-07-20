#include "variable_order.h"
#include <queue>

std::vector<int> computeMinWidthOrder(OrdinaryGraph const & G) {
    std::vector<int> variableOrder(G.n);
    std::vector<int> degree;
    degree.reserve(G.n);

    // maintain a priority queue of tuples (v, degree in G\set of vertices considered so far)
    // top node is the node with least degree in the graph of remaining vertices
    // ties are broken to prefer lower indexed vertices to match the implementation of A. Karahalios and W.-J. van Hoeve
    auto const minDegree = [](std::pair<int, size_t> const & a, std::pair<int, size_t> const & b){ return a.second == b.second ? a.first > b.first : a.second > b.second; };
    std::priority_queue<std::pair<int, size_t>, std::vector<std::pair<int, size_t>>, decltype(minDegree)> Q(minDegree);
    for (int v = 0; v < G.n; ++v) {
        degree.push_back(G.adj[v].size());
        Q.emplace(v, degree.back());
    }

    for (int i = G.n - 1; i >= 0; --i) {
        // find vertex with minimum degree in the remaining graph
        int v;
        while (degree[v = Q.top().first] < 0)
            Q.pop();
        Q.pop();

        variableOrder[i] = v;
        degree[v] = -1; // mark as already placed in order

        // update degrees of adjacent vertices
        for (int a : G.adj[v])
            if (degree[a] >= 0)
                Q.emplace(a, --degree[a]);
    }

    return variableOrder;
}

std::vector<int> computeMaxDegreeOrder(OrdinaryGraph const & G) {
    std::vector<int> variableOrder(G.n);
    std::vector<int> degree(G.n, 0);

    // maintain a priority queue of tuples (v, degree in G, degree w.r.t. vertices already added)
    // top node is node with largest degree w.r.t. nodes considered so far
    // ties are broken towards larger degree in G and smaller vertex index
    auto const maxDegree = [](std::tuple<int, size_t, size_t> const & a, std::tuple<int, size_t, size_t> const & b){
        if (std::get<2>(a) == std::get<2>(b)) {
            if (std::get<1>(a) == std::get<1>(b)) {
                return std::get<0>(a) > std::get<0>(b);
            } else {
                return std::get<1>(a) < std::get<1>(b);
            }
        } else {
            return std::get<2>(a) < std::get<2>(b);
        }
    };
    std::priority_queue<std::tuple<int, size_t, size_t>, std::vector<std::tuple<int, size_t, size_t>>, decltype(maxDegree)> Q(maxDegree);
    for (int v = 0; v < G.n; ++v)
        Q.emplace(v, G.adj[v].size(), 0);
    
    for (int i = 0; i < G.n; ++i) {
        // obtain next top node
        int v;
        while (degree[v = std::get<0>(Q.top())] < 0)
            Q.pop();
        Q.pop();

        variableOrder[i] = v;
        degree[v] = -1; // mark as already placed in order

        // update degrees of adjacent vertices
        for (int a : G.adj[v])
            if (degree[a] >= 0)
                Q.emplace(a, G.adj[a].size(), ++degree[a]);
    }

    return variableOrder;
}