#pragma once

#include "ordinary_graph.h"

/** compute a variable order for G based on the min width heuristic by A. Karahalios and W.-J. van Hoeve */
std::vector<int> computeMinWidthOrder(OrdinaryGraph const & G);

/** compute a variable order for G based on the max degree heuristic by W.-J. van Hoeve */
std::vector<int> computeMaxDegreeOrder(OrdinaryGraph const & G);