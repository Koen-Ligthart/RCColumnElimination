#include "CLI11.hpp"
#include "algorithm.h"
#include "profiler.h"
#include "rc_util.h"
#include "scip_exception.hpp"
#include "variable_order.h"

int main(int const nArgs, char * * args) {
    try {
        // parse CLI arguments
        CLI::App app{ "Computing the relaxation complexity of two finite sets using decision diagrams"};
        AlgorithmSettings settings;
        args = app.ensure_utf8(args);

        app.set_config("--config")->check(CLI::ExistingFile);

        std::string fileX;
        app.add_option("fileX,-x,--fileX", fileX, "inner set of points X")->required()->check(CLI::ExistingFile);
        std::string fileY;
        app.add_option("fileY,-y,--fileY", fileY, "outer set of points Y")->required()->check(CLI::ExistingFile);
        std::string strategy;
        app.add_option("-s,--strategy", strategy, "solution strategy")->check(CLI::IsMember({
            "localSeparation", "globalSeparation", "statefulHypergraph", "statefulTwoStep"}
        ))->default_val("statefulHypergraph");

        app.add_option("fileOut,-o,--fileOut", settings.fileOut, "file containing result report");
        app.add_flag("-d,--debugDot", settings.generateDotFiles, "generate .dot files")->default_val(settings.generateDotFiles);
        app.add_flag("-l,--debugLp", settings.generateLPFiles, "generate .lp files for the diagram flow model")->default_val(settings.generateLPFiles);
        app.add_flag("-v,--verbose", settings.verboseSCIP, "whether SCIP should output solving information to stdout")->default_val(settings.verboseSCIP);
        app.add_option("-c,--conflictCount", settings.maxConflictCount, "maximum number of conflicts separated simultaneously")->default_val(settings.maxConflictCount)->check(CLI::PositiveNumber);
        app.add_flag("--heuristicClassGrowth", settings.heuristicClassGrowth, "whether to use heuristic color class growth in the primal heuristic")->default_val(settings.heuristicClassGrowth)->capture_default_str();
        app.add_option("--epsilon", settings.eps, "numerical tolerance")->default_val(settings.eps)->check(CLI::PositiveNumber);
        app.add_option("--scipSettings", settings.scipSet, "SCIP settings file")->check(CLI::ExistingFile);
        app.add_flag("--startRelaxed", settings.startRelaxed, "whether to start solving LP relaxations of the flow model rather than the IPs")->default_val(settings.startRelaxed);
        std::map<std::string, FlowChooser> const flowChooserMap{
            {"prefer0Arc",              FLOW_CHOOSER_PREFER_0ARC},
            {"prefer1Arc",              FLOW_CHOOSER_PREFER_1ARC},
            {"preferHighestThen0Arc",   FLOW_CHOOSER_PREFER_HIGHEST_THEN_0ARC},
            {"preferHighestThen1Arc",   FLOW_CHOOSER_PREFER_HIGHEST_THEN_1ARC}
        };
        app.add_option("--flowChooserConflictFinder", settings.flowChooserConflictFinder, "flow chooser to use for finding conflicts")->default_val(settings.flowChooserConflictFinder)->transform(CLI::CheckedTransformer(flowChooserMap));
        app.add_option("--flowChooserPrimalHeuristic", settings.flowChooserPrimalHeuristic, "flow chooser to use for constructing primal heuristic solutions")->default_val(settings.flowChooserPrimalHeuristic)->transform(CLI::CheckedTransformer(flowChooserMap));
        app.add_flag("--sift", settings.sift, "whether to reorder the decision diagram after every iteration")->default_val(settings.sift);
        app.add_option("--longestPathIterations", settings.longestPathIterations, "number of starting iterations in which conflicts are found on longest paths instead of flow solutions")->default_val(settings.longestPathIterations)->check(CLI::NonNegativeNumber);
        app.add_flag("--flowObjCut", settings.flowObjCut, "whether to use objective bounds in the flow model formulation as constraints")->default_val(settings.flowObjCut);
        app.add_flag("--warmStart", settings.warmStart, "whether to supply SCIP with primal feasible flow model solutions obtained from the primal heuristic")->default_val(settings.warmStart);
        app.add_flag("--cliqueCuts", settings.cliqueCuts, "whether to use greedy maximal clique cuts in the flow model")->default_val(settings.cliqueCuts);
        app.add_option("--arcRedirection", settings.arcRedirectionPolicy, "way in which arcs are heuristically redirected")->default_val(settings.arcRedirectionPolicy)->transform(CLI::CheckedTransformer(std::map<std::string, ArcRedirectionPolicy>{
            {"none",                    ARC_REDIRECTION_POLICY_None},
            {"splitExact",              ARC_REDIRECTION_POLICY_DuringSplitDirectToExactState},
            {"splitRelaxed",            ARC_REDIRECTION_POLICY_DuringSplitDirectToRelaxedState},
            {"alwaysExact",             ARC_REDIRECTION_POLICY_BottomUpDirectToExactState},
            {"alwaysRelaxed",           ARC_REDIRECTION_POLICY_BottomUpDirectToRelaxedState}
        }));
        app.add_option("--minConflictFlow", settings.minConflictFlow, "minimum amount of flow required to report a conflict")->default_val(settings.minConflictFlow)->check(CLI::Range(0.0, 1.0) & CLI::PositiveNumber);
        app.add_option("--variableOrder", settings.variableOrder, "static initial variable ordering")->default_val(settings.variableOrder)->transform(CLI::CheckedTransformer(std::map<std::string, VariableOrder>{
            {"identity",                VARIABLE_ORDER_Identity},
            {"minWidth",                VARIABLE_ORDER_MinWidth},
            {"maxDegree",               VARIABLE_ORDER_MaxDegree}
        }));
        app.add_flag("--removeNonObservers", settings.removeNonObservers, "whether to preprocess an instance and remove non observer points from Y")->default_val(settings.removeNonObservers);
        app.add_option("--maxConvexHullDimension", settings.maxConvexHullDimension, "maximum dimension of instance for which to use convex hull based approaches")->default_val(settings.maxConvexHullDimension);
        app.add_option("--heuristicTimeBudget", settings.heuristicTimeBudget, "fraction of total time at which heuristic computation is paused")->default_val(settings.heuristicTimeBudget)->check(CLI::Range(0.0, 1.0));
        app.add_option("--arcRedirectionTimeBudget", settings.arcRedirectionTimeBudget, "fraction of total time at which arc redirection is paused")->default_val(settings.arcRedirectionTimeBudget)->check(CLI::Range(0.0, 1.0));
        CLI11_PARSE(app, nArgs, args);

        if (settings.fileOut)
            std::filesystem::create_directories(std::filesystem::path(settings.fileOut.value()).parent_path());

        getProfiler(); // initialize timings

        RCInstance instance(fileX, fileY);
        instance.print(false);

        // remove non-observers
        if (settings.removeNonObservers) {
            std::cout << "removing non observers from Y" << std::endl;
            PROFILE(PROFILE_RemoveNonObservers);
            instance = filterNonobservers(instance, settings);
            instance.print(false);
        }
    
        // compute hiding graph
        std::optional<OrdinaryGraph> hidingGraph;
        if (strategy == "statefulHypergraph" || strategy == "statefulTwoStep" || settings.variableOrder != VARIABLE_ORDER_Identity || settings.cliqueCuts) {
            std::cout << "computing hiding graph" << std::endl;
            hidingGraph = computeHidingGraph(instance, settings);
        }

        // compute heuristic variable order
        std::vector<int> variableOrder;
        switch (settings.variableOrder) {
            case VARIABLE_ORDER_Identity:
                for (size_t i = 0; i < instance.Y.size(); ++i)
                    variableOrder.push_back(i);
                break;
            case VARIABLE_ORDER_MinWidth:
                std::cout << "computing min width variable order" << std::endl;
                variableOrder = computeMinWidthOrder(hidingGraph.value());
                break;
            case VARIABLE_ORDER_MaxDegree:
                std::cout << "computing max degree variable order" << std::endl;
                variableOrder = computeMaxDegreeOrder(hidingGraph.value());
                break;
        }

        // compute maximal cliques for each vertex
        std::optional<std::vector<std::size_t>> cliqueSizes;
        if (settings.cliqueCuts)
            cliqueSizes = findMaximalCliqueSizes(hidingGraph.value());

        if (strategy == "localSeparation" || strategy == "globalSeparation") {
            // state oblivious separation
            if (settings.arcRedirectionPolicy != ARC_REDIRECTION_POLICY_None) {
                std::cout << "arc redirection is incompatible with state oblivious diagrams" << std::endl;
                return 105;
            }
            settings.separator = strategy == "localSeparation" ? CONFLICT_SEPARATOR_Local : CONFLICT_SEPARATOR_Global;
            computeRC(instance, variableOrder, cliqueSizes, settings);
        } else {
            // stateful separation
            settings.separator = CONFLICT_SEPARATOR_LocalStateful;
            if (settings.sift && strategy == "statefulHypergraph") {
                std::cout << "sifting is incompatible with stateful diagrams" << std::endl;
                return 105;
            }
            if (strategy == "statefulHypergraph") {
                computeRCHypergraph(instance, variableOrder, cliqueSizes, hidingGraph.value(), settings);
            } else { // strategy == "statefulTwoStep"
                computeRCTwoStep(instance, variableOrder, cliqueSizes, hidingGraph.value(), settings);
            }
        }
    } catch (SCIPException & e) {
        std::cerr << "ERROR " << e.what() << std::endl;
        return -1;
    }
    return 0;
}