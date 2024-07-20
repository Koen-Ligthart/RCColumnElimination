#pragma once

#include <array>
#include <chrono>
#include <string>
#include <vector>

/** possible options for labelling timing information */
enum ProfileCategory {
    PROFILE_ComputeHidingGraph,
    PROFILE_IndependenceOracle,
    PROFILE_FlowModelConstruction,
    PROFILE_FlowModelDestruction,
    PROFILE_FlowModelSolving,
    PROFILE_LongestPath,
    PROFILE_FindConflicts,
    PROFILE_PrimalHeuristic,
    PROFILE_Separation,
    PROFILE_BDDGarbageCollection,
    PROFILE_BDDConjunction,
    PROFILE_VariableReordering,
    PROFILE_RedirectArcs,
    PROFILE_RemoveNonObservers,
    PROFILE_EdgeRegistration,
    PROFILE_Output
};

/** names of profile categories */
constexpr auto profileCategoryNames = std::to_array<std::pair<ProfileCategory, char const *>>({
    {PROFILE_ComputeHidingGraph,    "hidingGraph"},
    {PROFILE_IndependenceOracle,    "oracle"},
    {PROFILE_FlowModelConstruction, "modelConstr"},
    {PROFILE_FlowModelDestruction,  "modelDestr"},
    {PROFILE_FlowModelSolving,      "modelSolve"},
    {PROFILE_LongestPath,           "longestPath"},
    {PROFILE_FindConflicts,         "findConfl"},
    {PROFILE_PrimalHeuristic,       "heur"},
    {PROFILE_Separation,            "sep"},
    {PROFILE_BDDGarbageCollection,  "bddGc"},
    {PROFILE_BDDConjunction,        "bddConj"},
    {PROFILE_VariableReordering,    "varOrder"},
    {PROFILE_RedirectArcs,          "redirectArcs"},
    {PROFILE_RemoveNonObservers,    "removeNonObservers"},
    {PROFILE_EdgeRegistration,      "edgeRegistration"},
    {PROFILE_Output,                "output"}
});

/** node in the profiler trie */
struct ProfilerTrieNode {

    std::chrono::nanoseconds                totalTime{};                            /** total time elapsed for this node since profiler creation */
    std::chrono::nanoseconds                newTime{};                              /** time elapsed for this node since last profiler print */
    int                                     parent;                                 /** index of parent node in the trie */
    int                                     adj[std::size(profileCategoryNames)];   /** indices of child nodes indexed by child profile category */

    /** create a node with the given parent */
    ProfilerTrieNode(int newParent);

};

/** object containing profiling information for different (nested) profiling categories */
class Profiler {

private:

    std::vector<ProfilerTrieNode>           trie;                                   /** nodes in the profiling trie */
    int                                     current;                                /** index of the node in the trie to which
                                                                                     *  the current operations are being charged */
    std::chrono::system_clock::time_point   lastStart;                              /** time since last current node modification */

    /** change current node to a child of the current node indexed by the given category
     *  if that child does not exist yet, it will create the node and initialize it with 0 total time
     */
    void pushCategory(ProfileCategory category);

    /** change current node to the parent of the current node */
    void popCategory();

    /** using the given prefix, print to stdout and update timing information of the node and all its children */
    void print(int node, std::string const & prefix);

    /** using the given prefix, output the total timing information of the node and all its children */
    void logTotal(int node, std::string const & prefix, std::ostream & out) const;

    /** returns the amount of time spent in this node and all its children recursively */
    std::chrono::nanoseconds sumTime(int node);

    Profiler(Profiler const &) = delete;
    Profiler(Profiler &&) = delete;
    Profiler & operator=(Profiler const &) = delete;
    Profiler & operator=(Profiler &&) = delete;

public:

    friend class ProfilerMeasurement;

    /** create a profiler with empty profiling trie */
    Profiler();
    
    /** print all timing information gathered so far to stdout and update timing information */
    void print();

    /** output all total timing results */
    void logTotal(std::ostream & out) const;

    /** reports whether the total time spend in the current profiler node is at most budget times the total time spent in the application */
    bool currentNodeWithinTimeBudget(double budget);

};

/** utility struct to push and pop profiling measurements in a scope through the use of C++ destructors */
struct ProfilerMeasurement {

    /** push the category to the global profiler */
    ProfilerMeasurement(ProfileCategory category);

    /** pop the current node in the global profiler */
    ~ProfilerMeasurement();

};

/** allow for automatic variable naming to prevent name shadowing problems */
#define CONCATENATE(a, b) CONCATENATE_INNER(a, b)
#define CONCATENATE_INNER(a, b) a ## b

/** register elapsed time in this scope under the given category in the global profiler */
#define PROFILE(category) ProfilerMeasurement CONCATENATE(profilerMeasurementObject, __COUNTER__)(category);

/** obtain a reference to the global profiler */
Profiler & getProfiler();