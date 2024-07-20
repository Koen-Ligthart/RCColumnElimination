#pragma once

#include <functional>
#include <limits>
#include <optional>
#include <string>

enum ConflictSeparator {
    CONFLICT_SEPARATOR_Global,
    CONFLICT_SEPARATOR_Local,
    CONFLICT_SEPARATOR_LocalStateful
};

enum ArcRedirectionPolicy {
    ARC_REDIRECTION_POLICY_None,
    ARC_REDIRECTION_POLICY_DuringSplitDirectToExactState,
    ARC_REDIRECTION_POLICY_DuringSplitDirectToRelaxedState,
    ARC_REDIRECTION_POLICY_BottomUpDirectToExactState,
    ARC_REDIRECTION_POLICY_BottomUpDirectToRelaxedState
};

enum FlowChooser {
    FLOW_CHOOSER_PREFER_0ARC,               /** prefers 0-arc if it has flow */
    FLOW_CHOOSER_PREFER_1ARC,               /** prefers 1-arc if it has flow */
    FLOW_CHOOSER_PREFER_HIGHEST_THEN_0ARC,  /** prefers arc with highest flow and prefers 0-arc in case of a tie */
    FLOW_CHOOSER_PREFER_HIGHEST_THEN_1ARC   /** prefers arc with highest flow and prefers 1-arc in case of a tie */
};

enum VariableOrder {
    VARIABLE_ORDER_Identity,
    VARIABLE_ORDER_MinWidth,
    VARIABLE_ORDER_MaxDegree
};

/** container for settings for the algorithms */
class AlgorithmSettings {

private:

    /** series of tie breaking functions to use when decomposing a flow into paths */
    std::function<bool(double, double)> const
                                    FLOW_CHOOSER_prefer0Arc;                            /** prefers 0-arc if it has flow */
    std::function<bool(double, double)> const
                                    FLOW_CHOOSER_prefer1Arc;                            /** prefers 1-arc if it has flow */
    std::function<bool(double, double)> const
                                    FLOW_CHOOSER_preferHighestThen0Arc;                 /** prefers arc with highest flow and prefers 0-arc in case of a tie */
    std::function<bool(double, double)> const
                                    FLOW_CHOOSER_preferHighestThen1Arc;                 /** prefers arc with highest flow and prefers 1-arc in case of a tie */
    
public:

    ConflictSeparator               separator = CONFLICT_SEPARATOR_Global;              /** way of resolving conflicts in the BDD */
    std::optional<std::string>      fileOut;                                            /** target location where bounds and other result metadata are reported */
    bool                            generateDotFiles = false;                           /** whether to generate .dot files displaying the decision diagrams */
    bool                            generateLPFiles = false;                            /** whether to generate .lp files that describe the flow model */
    bool                            verboseSCIP = false;                                /** whether to output SCIP solving information and statistics */                                                                                     
    size_t                          maxConflictCount = std::numeric_limits<size_t>::max();
                                                                                        /** the maximum number of conflicts that may be resolved within a single iteration */
    bool                            heuristicClassGrowth = true;                        /** whether we should greedily add vertices to
                                                                                         *  the recently decomposed path color class in the primal heuristic 
                                                                                         */
    double                          eps;                                                /** numerical tolerance
                                                                                         *  this is default initialized by the default SCIP settings
                                                                                         */
    std::optional<std::string>      scipSet;                                            /** file that provides SCIP settings */
    bool                            startRelaxed = true;                                /** whether to first use the flow model LP relaxation to guide BDD construction
                                                                                         *  the algorithm will swap over to the IP once the optimal LP solution does not conflict */
    FlowChooser                     flowChooserConflictFinder = FLOW_CHOOSER_PREFER_1ARC;
                                                                                        /** heuristic decision maker for identifying paths in flow decomposition
                                                                                         *  during the procedure that identifies conflicts
                                                                                         */
    FlowChooser                     flowChooserPrimalHeuristic = FLOW_CHOOSER_PREFER_HIGHEST_THEN_1ARC;
                                                                                        /** heuristic decision maker for identifying paths in flow decomposition
                                                                                         *  during the procedure that finds a heuristic feasible solution
                                                                                         */
    bool                            sift = false;                                       /** whether to dynamically reorder variables of the diagram after
                                                                                         *  every separated conflict
                                                                                         *  this uses the sifting heuristic by R. Rudell
                                                                                         */
    int                             longestPathIterations = 0;                          /** the number of iterations in which the longest path
                                                                                         *  is used to detect conflicts for refining the diagram
                                                                                         *  this is done in a first phase before solving the flow model
                                                                                         */
    bool                            flowObjCut = false;                                 /** whether to use dual objective bounds in the flow model formulation as constraints */
    bool                            warmStart = false;                                  /** whether to supply SCIP with primal feasible flow model solutions obtained from the primal heuristic */
    bool                            cliqueCuts = false;                                 /** whether to use admissible cuts based on greedy maximal cliques to
                                                                                         *  restrict the flow model solution
                                                                                         */
    ArcRedirectionPolicy            arcRedirectionPolicy = ARC_REDIRECTION_POLICY_DuringSplitDirectToExactState;
    double                          minConflictFlow = 0.5;                              /** the minimum flow that needs to be on a path
                                                                                         *  for it to be considered for conflict detection
                                                                                         */
    VariableOrder                   variableOrder = VARIABLE_ORDER_MinWidth;            /** heuristic used for initial static variable ordering */
    bool                            removeNonObservers = false;                         /** whether to remove non-observers from Y */
    int                             maxConvexHullDimension = 8;                         /** maximum dimension of the instance for which the convex hull
                                                                                         *  is computed to filter observers and compute the hiding graph
                                                                                         */
    double                          heuristicTimeBudget = 0.2;                          /** maximum fraction of time that may be spent on the primal heuristic
                                                                                         *  when this fraction is exceeded, heuristic computation is skipped in an interation
                                                                                         */
    double                          arcRedirectionTimeBudget = 0.2;                     /** maximum fraction of time that may be spent on bottom-up arc redirection
                                                                                         *  when this fraction is exceeded, redirection is skipped in an interation
                                                                                         */

    AlgorithmSettings();

    /** |x| < eps */
    bool isZero(double x) const;

    /** x > eps */
    bool isPositive(double x) const;

    /** ceil(x - eps) */
    double ceilReal(double x) const;

    std::function<bool(double, double)> const & getFlowChooserByType(const FlowChooser flowChooser) const;
    
};