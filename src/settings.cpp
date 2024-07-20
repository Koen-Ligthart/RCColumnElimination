#include <stdexcept>
#include "scip_exception.hpp"
#include "settings.h"

AlgorithmSettings::AlgorithmSettings() : 
    FLOW_CHOOSER_prefer0Arc             { [this](double const f0, double const f1){ return isPositive(f0     ) ? 0 : 1; } },
    FLOW_CHOOSER_prefer1Arc             { [this](double const f0, double const f1){ return isPositive(f1     ) ? 1 : 0; } },
    FLOW_CHOOSER_preferHighestThen0Arc  { [this](double const f0, double const f1){ return isPositive(f1 - f0) ? 1 : 0; } },
    FLOW_CHOOSER_preferHighestThen1Arc  { [this](double const f0, double const f1){ return isPositive(f0 - f1) ? 0 : 1; } }
{
    // get default value of SCIP eps
    SCIP * scip;
    SCIP_CALL_EXC( SCIPcreate(&scip) );
    SCIP_CALL_EXC( SCIPgetRealParam(scip, "numerics/epsilon", &eps) );
    SCIP_CALL_EXC( SCIPfree(&scip) );
}

bool AlgorithmSettings::isZero(double const x) const {
    return EPSZ(x, eps); // SCIPisZero copy
}

bool AlgorithmSettings::isPositive(double const x) const {
    return EPSP(x, eps); // SCIPisPositive copy
}

double AlgorithmSettings::ceilReal(double const x) const {
    return EPSCEIL(x, eps); // SCIPceil copy
}

std::function<bool(double, double)> const & AlgorithmSettings::getFlowChooserByType(FlowChooser const flowChooser) const {
    switch (flowChooser) {
        case FLOW_CHOOSER_PREFER_0ARC:              return FLOW_CHOOSER_prefer0Arc;
        case FLOW_CHOOSER_PREFER_1ARC:              return FLOW_CHOOSER_prefer1Arc;
        case FLOW_CHOOSER_PREFER_HIGHEST_THEN_0ARC: return FLOW_CHOOSER_preferHighestThen0Arc;
        case FLOW_CHOOSER_PREFER_HIGHEST_THEN_1ARC: return FLOW_CHOOSER_preferHighestThen1Arc;
    }
    throw std::runtime_error("");
}