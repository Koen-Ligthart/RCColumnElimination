#include <iomanip>
#include <iostream>
#include "profiler.h"

ProfilerTrieNode::ProfilerTrieNode(int const newParent) : parent{ newParent } {
    std::fill(adj, adj + std::size(profileCategoryNames), -1);
}

Profiler::Profiler() {
    // create root of trie
    trie.push_back({-1});

    // set current to root
    current = 0;

    lastStart = std::chrono::system_clock::now();
}

std::chrono::nanoseconds Profiler::sumTime(int const node) {
    std::chrono::nanoseconds sum = trie[node].totalTime + trie[node].newTime;
    for (auto & [category, name] : profileCategoryNames) {
        int const child = trie[node].adj[category];
        if (child != -1)
            sum += sumTime(child);
    }
    return sum;
}

bool Profiler::currentNodeWithinTimeBudget(double const budget) {
    return sumTime(current).count() <= budget * sumTime(0).count();
}

void Profiler::pushCategory(ProfileCategory const category) {
    // end current counting
    std::chrono::nanoseconds const duration = std::chrono::system_clock::now() - lastStart;
    trie[current].newTime += duration;
    // move deeper into trie
    int newCurrent = trie[current].adj[category];
    if (newCurrent == -1) {
        // create new trie node with parent being current
        newCurrent = trie[current].adj[category] = trie.size();
        trie.emplace_back(current);
    }
    current = newCurrent;
    // start new counting
    lastStart = std::chrono::system_clock::now();
}

void Profiler::popCategory() {
    // end current counting
    std::chrono::nanoseconds const duration = std::chrono::system_clock::now() - lastStart;
    trie[current].newTime += duration;
    // move back to parent
    current = trie[current].parent;
    // start new counting
    lastStart = std::chrono::system_clock::now();
}

void printNanoseconds(std::chrono::system_clock::rep duration) {
    // hh:mm:ss:ms0:us0:ns0
    int const ns = duration % 1000; duration /= 1000;
    int const us = duration % 1000; duration /= 1000;
    int const ms = duration % 1000; duration /= 1000;
    int const s  = duration % 60;   duration /= 60;
    int const m  = duration % 60;   duration /= 60;
    int const h  = duration;
    std::cout << std::setw(2) << std::right; if (h > 0)                                       std::cout << h  << ":" << std::setfill('0'); else std::cout << "" << " ";
    std::cout << std::setw(2) << std::right; if (h > 0 || m > 0)                              std::cout << m  << ":" << std::setfill('0'); else std::cout << "" << " ";
    std::cout << std::setw(2) << std::right; if (h > 0 || m > 0 || s > 0)                     std::cout << s  << ":" << std::setfill('0'); else std::cout << "" << " ";
    std::cout << std::setw(3) << std::right; if (h > 0 || m > 0 || s > 0 || ms > 0)           std::cout << ms << ":" << std::setfill('0'); else std::cout << "" << " ";
    std::cout << std::setw(3) << std::right; if (h > 0 || m > 0 || s > 0 || ms > 0 || us > 0) std::cout << us << ":" << std::setfill('0'); else std::cout << "" << " ";
    std::cout << std::setw(3) << std::right; std::cout << ns;
    std::cout << std::setfill(' ');
}

void Profiler::print(int const node, std::string const & prefix) {
    trie[node].totalTime += trie[node].newTime;
    std::cout << std::setw(50) << std::left << prefix;
    std::cout << " | ";
    printNanoseconds(trie[node].newTime.count());
    std::cout << " | ";
    printNanoseconds(trie[node].totalTime.count());
    std::cout << "\n";
    trie[node].newTime = std::chrono::nanoseconds{};
    for (auto & [category, name] : profileCategoryNames) {
        int const child = trie[node].adj[category];
        if (child != -1)
            print(child, node == 0 ? (std::string{"    "} + name) : (prefix + "." + name));
    }
}

void Profiler::print() {
    std::cout << "profiler:\n";
    std::cout << std::setw(50) << std::left << "    category" << " | " << std::setw(20) << std::right << "new" << " | " << std::setw(20) << std::right << "cumulative" << "\n";
    std::cout << std::setw(50) << std::left << "" << " | " << std::setw(20) << std::right << "hh:mm:ss:ms0:us0:ns0" << " | " << std::setw(20) << std::right << "hh:mm:ss:ms0:us0:ns0" << "\n";
    print(0, "    untracked");
    std::cout.flush();
}

void Profiler::logTotal(int const node, std::string const & prefix, std::ostream & out) const {
    out << std::setw(30) << std::left << prefix << " = " << std::setw(30) << std::right << (trie[node].totalTime + trie[node].newTime).count() << "\n";
    for (auto &[category, name] : profileCategoryNames) {
        const int child = trie[node].adj[category];
        if (child != -1)
            logTotal(child, node == 0 ? (std::string{"profiler_"} + name) : (prefix + "_" + name), out);
    }
}

void Profiler::logTotal(std::ostream & out) const {
    logTotal(0, "profiler_untracked", out);
}

static Profiler profiler;

ProfilerMeasurement::ProfilerMeasurement(ProfileCategory const category) {
    profiler.pushCategory(category);
}

ProfilerMeasurement::~ProfilerMeasurement() {
    profiler.popCategory();
}

Profiler & getProfiler() {
    return profiler;
}