#pragma once

#include "cudd.h"
#include "primitives.h"

/** throws an exception when z is NULL */
void zddVerifyNonNull(DdManager * zdd, DdNode * z);

/** recursively computes the collection of all subsets of w that are not a superset of an element in z */
DdNode * zddNonSup(DdManager * zdd, DdNode * w, DdNode * z);

/** recursively computes whether for each element in w there is a subset in z
 *  and returns a terminal indicating the boolean value
 */
DdNode * zddWeakInclNode(DdManager * zdd, DdNode * w, DdNode * z);

/** recursively computes whether for each element in w there is a subset in z */
bool zddWeakIncl(DdManager * zdd, DdNode * w, DdNode * z);

/** recursively computes the collection of all minimal subsets of w */
DdNode * zddMinMal(DdManager * zdd, DdNode * w);

/** constructs a ZDD node that represents the collection containing only the specified subset element
 *  assumes that element is sorted on layer increasingly
 */
DdNode * zddSingleton(DdManager * zdd, VertexList const & element);

/** constructs a ZDD node that represents a collection of singletons for each of the elements in singletons
 *  assumes that singletons is sorted on layer increasingly
 */
DdNode * zddSingletons(DdManager * zdd, VertexList const & singletons);

/** writes a ZDD representation of a diagram to a .dot file */
void writeZDDToDot(
    int                                                     n,             /** the number of variables of the ZDD */
    DdManager *                                             dd,            /** the ZDD to represent */
    std::string const &                                     file,          /** name of the file */
    std::vector<std::pair<std::string, DdNode *>> const &   labelledNodes  /** the labels of the nodes of which the subdiagrams should be shown */
);