#pragma once

#include <functional>
#include <memory>

typedef std::vector<int> VertexList;

/** stateful oracle that a candidate independent set of an independence system
 *  it represents access to the underlying independence system
 *  the initial state of an independent set is empty, which is assumed to be independent
 *  no more elements are assumed to be added from the moment this set became dependent
 */
class IndependentSet {

public:

    virtual ~IndependentSet() = default;

    /** updates the independent set to include i */
    virtual void add(int i) = 0;

    /** restores the state of the set by removing the last added element
     *  this is only allowed to be invoked when the set was independent before addition of this last element
     *  hence, isIndependent holds after invoking this
     */
    virtual void removeLastAdded() = 0;

    /** checks whether the stored set is independent */
    virtual bool isIndependent() const = 0;

    /** finds a circuit contained in this set which is minimal subset of this set that is not dependent
     *  this may be an arbitrary circuit
     *  this function may only be invoked if the last added element to this set caused the set to become dependent
     *  vertices are returned in the order in which they were inserted in this set
     */
    virtual VertexList findContainedCircuit() const = 0;

    /** provides read only view of indices contained in this set in order of insertion */
    virtual VertexList const & getElements() const = 0;
    
};

/** function that returns a stateful oracle maintaining partial independent sets */
typedef std::function<std::unique_ptr<IndependentSet>()> IndependenceSystem;