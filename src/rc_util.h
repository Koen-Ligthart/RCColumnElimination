#pragma once

#include "ordinary_graph.h"
#include "rc_instance.h"
#include "settings.h"

/** compute the hiding graph of an relaxation complexity instance */
OrdinaryGraph computeHidingGraph(RCInstance const & instance, AlgorithmSettings const & settings);

/** compute a new relaxation complexity instance of which all non-observer points of Y are removed */
RCInstance filterNonobservers(RCInstance const & instance, AlgorithmSettings const & settings);