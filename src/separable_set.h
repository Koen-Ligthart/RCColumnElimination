#pragma once

#include "scip/scip.h"
#include "primitives.h"
#include "rc_instance.h"
#include "settings.h"

/** represents a subset of Y so that X and Y can be strictly linearly separated by a single hyperplane */
class SeparableSet : public IndependentSet {

private:

    RCInstance const &                  instance;       /** underlying rc problem instance */
    SCIP_LPI *                          scipLpi;        /** pointer to SCIP LPI model */
    VertexList                          points;         /** current set of points */
    AlgorithmSettings const &           settings;       /** settings */

    SeparableSet(SeparableSet const &) = delete;
    SeparableSet(SeparableSet &&) = delete;
    SeparableSet & operator=(SeparableSet const &) = delete;
    SeparableSet & operator=(SeparableSet &&) = delete;

public:

    /** constructs an empty separable subset of Y which initializes the SCIP LPI model with columns for X */
    SeparableSet(RCInstance const & newInstance, AlgorithmSettings const & newSettings);
    ~SeparableSet();
    void add(int i);
    void removeLastAdded();
    bool isIndependent() const;
    VertexList findContainedCircuit() const;
    VertexList const & getElements() const;

};