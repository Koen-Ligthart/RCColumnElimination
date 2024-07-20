#pragma once

#include <iomanip>
#include <filesystem>
#include "flow_model.hpp"
#include "separation.hpp"

/** iteratively refine the given BDD and attempt to find the chromatic number associated to the given independence system */
template<IsBDD BDD>
std::tuple<std::optional<std::vector<VertexList>>, size_t, int> iterativelyRefine(
    BDD &                           bdd,                    /** the initial BDD that is to be refined */
    IndependenceSystem const &      independenceSystem,     /** oracle providing access to whether subsets of vertices are independent */
    std::optional<int> const        previousIterationCount, /** number of iterations performed in a previous stage */
    size_t const                    newLb,                  /** known lower bound on the chromatic number from a previous stage */
    std::optional<std::vector<size_t>> const &
                                    cliqueSizes,            /** size of a maximal clique for each vertex */
    AlgorithmSettings const &       settings                /** settings to use in the algorithm */
) {

    assert( settings.separator != CONFLICT_SEPARATOR_LocalStateful || !settings.sift );
    assert( settings.minConflictFlow <= 1.0 );
    assert( cliqueSizes.has_value() == settings.cliqueCuts );

    static_assert(IsStatefulBDD<BDD> || (IsUniqueTableBDD<BDD> && IsDisconnectedCheckableNode<class BDD::NodeType>));
    if constexpr(IsStatefulBDD<BDD>)
        assert(settings.separator == CONFLICT_SEPARATOR_LocalStateful);
    if constexpr(IsReferenceCountedNode<class BDD::NodeType>)
        assert(settings.separator == CONFLICT_SEPARATOR_Global);
    if constexpr(IsParentListedNode<class BDD::NodeType>)
        assert(settings.separator == CONFLICT_SEPARATOR_Local || settings.separator == CONFLICT_SEPARATOR_LocalStateful);

    std::ostringstream namebuf{std::ios_base::ate};

    int iteration = 0;
    std::chrono::system_clock::time_point lastLog = std::chrono::system_clock::now(); // last time since logging the profiler information
    std::optional<std::vector<VertexList>> bestSol;
    size_t lb = newLb;
    size_t ub = std::numeric_limits<size_t>::max() - 1;

    bool relaxed = settings.startRelaxed;
    bool longestPath = settings.longestPathIterations > 0;

    while (true) {
        ++iteration;

        if (longestPath && iteration >= settings.longestPathIterations) {
            std::cout << "stop separation on longest paths and swap to flow model solving" << std::endl;
            longestPath = false;
        }

        // potentially log profiler information
        std::chrono::system_clock::time_point newTime = std::chrono::system_clock::now();
        if ((newTime - lastLog).count() > 10'000'000'000) {
            PROFILE(PROFILE_Output);
            lastLog = newTime;
            getProfiler().print();
        }

        // sift
        int const ddSizePreSift = bdd.nodeCount();
        if (settings.sift) {
            if constexpr (IsDisconnectedCheckableNode<class BDD::NodeType> && IsUniqueTableBDD<BDD>) {
                PROFILE(PROFILE_VariableReordering);
                sift(bdd);
            } else assert(false);
        }

        bool foundBetterHeuristicSol = false;
        std::vector<ConflictSolution> conflicts;

        double flowValue = 0;
        std::string activeNodeCount = "?";

        if (longestPath) {
            // attempt to find conflict on the longest path in the BDD
            
            Path path = findLongestPath(bdd);
            auto const conflict = findConflictOnPath(bdd, path, independenceSystem);
            if (conflict)
                conflicts.push_back(conflict.value());
            
            // output BDD and longest path to .dot file
            if (settings.generateDotFiles) {
                namebuf.str("out/decision-diagram-iteration-");
                namebuf << iteration << ".dot";
                std::unordered_map<Arc, double, ArcHash> flow;
                for (size_t i = 0; i < path.length(); ++i)
                    flow[{path.nodes[i], path.arcLabels[i]}] = 1;
                writeToDot(bdd, namebuf.str(), false, false, true, flow);
                if constexpr(IsStatefulBDD<BDD>) {
                    namebuf.str("out/decision-diagram-iteration-");
                    namebuf << iteration << "-state-manager.dot";
                    bdd.writeStateManagerState(namebuf.str());
                }
            }

        } else {
            // construct and solve the flow model to identify conflicts

            FlowModel<BDD> model(bdd, relaxed, independenceSystem, lb, cliqueSizes, settings, bestSol);
            activeNodeCount = std::to_string(model.activeNodeCount());

            std::optional<double> objVal = model.solve();
            assert(objVal); // cannot become infeasible without yielding a conflict of cardinality at most 1

            // output BDD and flow solution to .dot file
            if (settings.generateDotFiles) {
                PROFILE(PROFILE_Output);
                namebuf.str("out/decision-diagram-iteration-");
                namebuf << iteration << ".dot";
                writeToDot(bdd, namebuf.str(), false, false, true, model.getFlow());
                if constexpr(IsStatefulBDD<BDD>) {
                    namebuf.str("out/decision-diagram-iteration-");
                    namebuf << iteration << "-state-manager.dot";
                    bdd.writeStateManagerState(namebuf.str());
                }
            }

            // output flow model to .lp file
            if (settings.generateLPFiles) {
                PROFILE(PROFILE_Output);
                namebuf.str("out/flow-model-iteration-");
                namebuf << iteration << ".lp";
                model.writeToLP(namebuf.str());
            }

            flowValue = objVal.value();
            lb = std::max(lb, (size_t) settings.ceilReal(flowValue));
            FlowDecompositionResult result = model.findConflicts();
            conflicts = std::move(result.conflicts);

            if (conflicts.empty() && !relaxed) { // if an IP solution has no conflicts, we have found the chromatic number
                std::cout << "found optimal flow solution without conflicts" << std::endl;
                assert(lb == result.partition.value().size());
                // result.partition may be a covering where vertices are colored with multiple colors, therefore we transform it into a proper coloring partition
                std::vector<bool> colored(bdd.getN());
                bestSol = std::vector<VertexList>(lb);
                for (size_t c = 0; c < lb; ++c) {
                    for (int const v : result.partition.value()[c]) {
                        if (!colored[v]) {
                            colored[v] = true;
                            bestSol.value()[c].push_back(v);
                        }
                    }
                }
                ub = lb;
            } else {
                // attempt to find an improved heuristic solution
                std::optional<std::vector<VertexList>> heuristicSol = model.primalHeuristic(ub);
                if (heuristicSol) {
                    assert(ub > heuristicSol.value().size());
                    foundBetterHeuristicSol = true;
                    bestSol = heuristicSol;
                    ub = bestSol.value().size();
                }
            }

        }

        // check whether problem is infeasible by checking if conflicts of cardinality at most 1 exist
        for (ConflictSolution const & conflict : conflicts) {
            if (conflict.edgeLayers.size() <= 1) {
                lb = std::numeric_limits<size_t>::max(); // represent infeasibility by lb > ub
                assert(!bestSol);
                break;
            }
        }

        // output information
        {
            PROFILE(PROFILE_Output);
            std::cout << "iteration " << iteration << " dd size " << bdd.nodeCount() << " active " << activeNodeCount;
            if (settings.sift)
                std::cout << " (" << ddSizePreSift << " pre sift)";
            if (longestPath) {
                std::cout << " longest path";
            } else {
                std::cout << " " << (relaxed ? "LP" : "IP") << " flow " << flowValue << " lb " << lb << " ub " << ub << " heur " << (foundBetterHeuristicSol ? "new" : "old");
            }
            std::cout << " conflicts " << conflicts.size() << " count sat " << countSat(bdd) << std::endl;

            // output results to output file
            if (settings.fileOut) {
                // rename old file to backup file
                if (std::filesystem::is_regular_file(settings.fileOut.value()))
                    std::filesystem::rename(settings.fileOut.value(), settings.fileOut.value() + ".old");
                
                std::ofstream out(settings.fileOut.value());

                // write contents
                if (previousIterationCount)
                    out << std::setw(30) << std::left << "stage_one_iterations" << " = " << std::setw(30) << std::right << previousIterationCount.value() << "\n";
                out << std::setw(30) << std::left << "iterations" << " = " << std::setw(30) << std::right << iteration << "\n";
                out << std::setw(30) << std::left << "lb" << " = " << std::setw(30) << std::right << lb << "\n";
                out << std::setw(30) << std::left << "ub" << " = " << std::setw(30) << std::right << ub << "\n";
                out << std::setw(30) << std::left << "dd_size_total" << " = " << std::setw(30) << std::right << bdd.nodeCount() << "\n";
                if (activeNodeCount != "?")
                    out << std::setw(30) << std::left << "dd_size_active" << " = " << std::setw(30) << std::right << activeNodeCount << "\n";
                out << std::setw(30) << std::left << "relaxed" << " = " << std::setw(30) << std::right << (relaxed ? "true" : "false") << "\n";
                getProfiler().logTotal(out);
                if (bestSol) {
                    out << std::setw(30) << std::left << "best_sol" << " = [\n";
                    bool firstClass = true;
                    for (VertexList const & clazz : bestSol.value()) {
                        if (firstClass) {
                            firstClass = false;
                        } else {
                            out << ",\n";
                        }
                        out << "  [ ";
                        bool firstVertex = true;
                        for (int const v : clazz) {
                            if (firstVertex) {
                                firstVertex = false;
                            } else {
                                out << ", ";
                            }
                            out << v;
                        }
                        out << " ]";
                    }
                    out << "\n]\n";
                }
                out.flush();

                // remove backup file
                std::filesystem::remove(settings.fileOut.value() + ".old");
            }
        }

        // check for termination
        if (lb == ub) {
            std::cout << "found optimal solution" << std::endl;
            break;
        } else if (lb > ub) {
            std::cout << "encountered infeasible problem" << std::endl;
            break;
        }

        // check for absence of conflicts and swap between phases
        if (conflicts.empty()) {
            if (longestPath) {
                std::cout << "no more conflicts on longest path, swapping to flow model solving" << std::endl;
                longestPath = false;
            } else { // relaxed holds
                std::cout << "no more conflicts in relaxation solution, swapping to integral flow model optimization" << std::endl;
                relaxed = false;
            }
        }

        // register edges of all conflicts
        if (settings.separator == CONFLICT_SEPARATOR_LocalStateful) {
            PROFILE(PROFILE_EdgeRegistration)
            if constexpr(IsStatefulBDD<BDD>) {
                for (ConflictSolution const & conflict : conflicts) {
                    VertexList E; // sorted on layers
                    for (int const l : conflict.edgeLayers)
                        E.push_back(bdd.getVariableByLayer(l));
                    bdd.registerEdge(E);
                }
            } else assert(false);
        }

        // resolve conflicts
        {
            PROFILE(PROFILE_Separation);
            for (ConflictSolution const & conflict : conflicts) {
                std::optional<ConflictPath> const path = extractConflictPath(bdd, conflict);
                if (path) {
                    if (settings.separator == CONFLICT_SEPARATOR_Global) {
                        if constexpr(IsReferenceCountedNode<class BDD::NodeType>) {
                            separateGlobal(bdd, path.value());
                        } else assert(false);
                    } else if (settings.separator == CONFLICT_SEPARATOR_Local) {
                        if constexpr(IsParentListedNode<class BDD::NodeType> && !IsStatefulBDD<BDD>) {
                            separateLocal(bdd, path.value());
                        } else assert(false);
                    } else { // settings.separator == CONFLICT_SEPARATOR_LocalStateful
                        if constexpr(IsStatefulBDD<BDD>) {
                            separateLocalStateful(bdd, path.value(), settings);
                        } else assert(false);
                    }
                }
            }
        }

        // redirect arcs
        if constexpr(IsStatefulBDD<BDD>) {
            if (settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_BottomUpDirectToExactState || settings.arcRedirectionPolicy == ARC_REDIRECTION_POLICY_BottomUpDirectToRelaxedState) {
                PROFILE(PROFILE_RedirectArcs);
                redirectArcs(bdd, settings);
            }
        }
    }

    getProfiler().print();

    return {bestSol, lb, iteration};
}