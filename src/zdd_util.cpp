#include <cassert>
#include <cstring>
#include <ranges>
#include <sstream>
#include "cuddInt.h"
#include "zdd_util.h"

void zddVerifyNonNull(DdManager * const zdd, DdNode * const z) {
    if (z == NULL) {
        assert( Cudd_ReadErrorCode(zdd) == CUDD_MEMORY_OUT );
        throw std::runtime_error("CUDD out of memory");
    }
}

// by D. Knuth
DdNode * zddNonSup(DdManager * const zdd, DdNode * const w, DdNode * const z) {
    if (z == Cudd_ReadZero(zdd)) // if z = emptyset then return w
        return w;
    if (w == Cudd_ReadZero(zdd) || z == DD_ONE(zdd) || w == z)
        return Cudd_ReadZero(zdd);
    DdNode * result = cuddCacheLookup2Zdd(zdd, zddNonSup, w, z);
    if (result != NULL)
        return result;
    int const lW = Cudd_ReadPermZdd(zdd, Cudd_NodeReadIndex(w));
    int const lZ = Cudd_ReadPermZdd(zdd, Cudd_NodeReadIndex(z));
    if (lW > lZ) { // all sets of z that include the element of the variable associated to z are not contained in elements of w
        result = zddNonSup(zdd, w, Cudd_E(z));
        Cudd_Ref(result);
    } else {
        DdNode * rE;
        DdNode * rT;
        if (lW < lZ) {
            rE = zddNonSup(zdd, Cudd_E(w), z);
            Cudd_Ref(rE);
            rT = zddNonSup(zdd, Cudd_T(w), z);
            Cudd_Ref(rT);
        } else {
            rE = zddNonSup(zdd, Cudd_E(w), Cudd_E(z));
            Cudd_Ref(rE);
            DdNode* a = zddNonSup(zdd, Cudd_T(w), Cudd_T(z));
            Cudd_Ref(a);
            DdNode* b = zddNonSup(zdd, Cudd_T(w), Cudd_E(z));
            Cudd_Ref(b);
            rT = cuddZddIntersect(zdd, a, b);
            zddVerifyNonNull(zdd, rT);
            Cudd_Ref(rT);
            Cudd_RecursiveDerefZdd(zdd, a);
            Cudd_RecursiveDerefZdd(zdd, b);
        }
        result = cuddZddGetNode(zdd, Cudd_NodeReadIndex(w), rT, rE);
        zddVerifyNonNull(zdd, result);
        Cudd_Ref(result);
        Cudd_RecursiveDerefZdd(zdd, rE);
        Cudd_RecursiveDerefZdd(zdd, rT);
    }
    cuddCacheInsert2(zdd, zddNonSup, w, z, result);
    Cudd_Deref(result);
    return result;
}

DdNode * zddWeakInclNode(DdManager * const zdd, DdNode * const w, DdNode * const z) {
    if (z == Cudd_ReadZero(zdd))
        return (w == Cudd_ReadZero(zdd)) ? DD_ONE(zdd) : Cudd_ReadZero(zdd);
    if (w == Cudd_ReadZero(zdd) || z == DD_ONE(zdd) || w == z)
        return DD_ONE(zdd);
    DdNode * result = cuddCacheLookup2Zdd(zdd, zddWeakInclNode, w, z);
    if (result != NULL)
        return result;
    result = cuddCacheLookup2Zdd(zdd, zddNonSup, w, z);
    if (result != NULL) {
        Cudd_Ref(result);
        DdNode* actualResult = (result == Cudd_ReadZero(zdd)) ? DD_ONE(zdd) : Cudd_ReadZero(zdd);
        Cudd_RecursiveDerefZdd(zdd, result);
        return actualResult;
    }
    int const lW = Cudd_ReadPermZdd(zdd, Cudd_NodeReadIndex(w));
    int const lZ = Cudd_ReadPermZdd(zdd, Cudd_NodeReadIndex(z));
    if (lW > lZ) {
        result = zddWeakInclNode(zdd, w, Cudd_E(z));
    } else if (lW < lZ) {
        result = ((zddWeakInclNode(zdd, Cudd_E(w), z) == DD_ONE(zdd)) && (zddWeakInclNode(zdd, Cudd_T(w), z) == DD_ONE(zdd))) ? DD_ONE(zdd) : Cudd_ReadZero(zdd);
    } else if (zddWeakInclNode(zdd, Cudd_E(w), Cudd_E(z)) != DD_ONE(zdd)) {
        // check if rE != 0, in which case we can immediately conclude that weak inclusion does not hold
        result = Cudd_ReadZero(zdd);
    } else {
        // lazily evaluate as much as possible to come to a quicker conclusion
        // need to return Intersection(NonSup, NonSup) == 0
        // if any of the NonSups are 0, we can immediately return true
        DdNode * const a = zddNonSup(zdd, Cudd_T(w), Cudd_T(z));
        Cudd_Ref(a);
        if (a == Cudd_ReadZero(zdd)) {
            Cudd_RecursiveDerefZdd(zdd, a);
            result = DD_ONE(zdd);
        } else {
            DdNode * const b = zddNonSup(zdd, Cudd_T(w), Cudd_E(z));
            Cudd_Ref(b);
            if (b == Cudd_ReadZero(zdd)) {
                Cudd_RecursiveDerefZdd(zdd, b);
                result = DD_ONE(zdd);
            } else {
                DdNode * const intersection = cuddZddIntersect(zdd, a, b);
                zddVerifyNonNull(zdd, intersection);
                Cudd_Ref(intersection);
                result = (intersection == Cudd_ReadZero(zdd)) ? DD_ONE(zdd) : Cudd_ReadZero(zdd);
                Cudd_RecursiveDerefZdd(zdd, intersection);
            }
        }
    }
    cuddCacheInsert2(zdd, zddWeakInclNode, w, z, result);
    return result;
}

bool zddWeakIncl(DdManager * const zdd, DdNode * const w, DdNode * const z) {
    return zddWeakInclNode(zdd, w, z) == DD_ONE(zdd);
}

// by D. Knuth
DdNode * zddMinMal(DdManager * const zdd, DdNode * const w) {
    if (w == Cudd_ReadZero(zdd) || w == DD_ONE(zdd))
        return w;
    DdNode * result = cuddCacheLookup1Zdd(zdd, zddMinMal, w);
    if (result != NULL)
        return result;
    DdNode * const minMalE = zddMinMal(zdd, Cudd_E(w));
    Cudd_Ref(minMalE);
    DdNode * const minMalT = zddMinMal(zdd, Cudd_T(w));
    Cudd_Ref(minMalT);
    DdNode * const childT = zddNonSup(zdd, minMalT, minMalE);
    Cudd_Ref(childT);
    Cudd_RecursiveDerefZdd(zdd, minMalT);
    result = cuddZddGetNode(zdd, Cudd_NodeReadIndex(w), childT, minMalE);
    zddVerifyNonNull(zdd, result);
    Cudd_Ref(result);
    Cudd_RecursiveDerefZdd(zdd, minMalE);
    Cudd_RecursiveDerefZdd(zdd, childT);
    cuddCacheInsert1(zdd, zddMinMal, w, result);
    Cudd_Deref(result);
    return result;
}

DdNode * zddSingleton(DdManager * const zdd, VertexList const & element) {
    DdNode * current = DD_ONE(zdd);
    Cudd_Ref(current);
    for (int const v : std::ranges::views::reverse(element)) {
        DdNode * const parent = cuddUniqueInterZdd(zdd, v, current, Cudd_ReadZero(zdd));
        zddVerifyNonNull(zdd, parent);
        Cudd_Ref(parent);
        Cudd_RecursiveDerefZdd(zdd, current);
        current = parent;
    }
    Cudd_Deref(current);
    return current;
}

DdNode * zddSingletons(DdManager * const zdd, VertexList const & singletons) {
    DdNode * current = Cudd_ReadZero(zdd);
    Cudd_Ref(current);
    for (int const v : std::ranges::views::reverse(singletons)) {
        DdNode * const parent = cuddUniqueInterZdd(zdd, v, DD_ONE(zdd), current);
        zddVerifyNonNull(zdd, parent);
        Cudd_Ref(parent);
        Cudd_RecursiveDerefZdd(zdd, current);
        current = parent;
    }
    Cudd_Deref(current);
    return current;
}

void writeZDDToDot(int const n, DdManager * const dd, std::string const & file, std::vector<std::pair<std::string, DdNode *>> const & labelledNodes) {
    int const m = labelledNodes.size();

    char * * const varnames = new char *[n];

    std::ostringstream namebuf;
    for (int i = 0; i < n; ++i) {
        namebuf.str("");
        namebuf << "v" << i;
        varnames[i] = strdup(namebuf.str().c_str());
    }

    char const * * const nodenames = new char const *[labelledNodes.size()];
    DdNode * * const nodes = new DdNode *[m];
    for (int j = 0; j < m; ++j) {
        nodenames[j] = labelledNodes[j].first.c_str();
        nodes[j] = labelledNodes[j].second;
    }
    FILE * const cFile = fopen(file.c_str(), "w");
    if (cFile == NULL)
        throw std::invalid_argument("file could not be found");

    Cudd_zddDumpDot(dd, m, nodes, varnames, nodenames, cFile);

    fclose(cFile);
    delete[] nodenames;
    for (int i = 0; i < n; ++i)
        free(varnames[i]);
    delete[] varnames;
    delete[] nodes;
}