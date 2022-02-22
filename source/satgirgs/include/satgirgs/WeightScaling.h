#pragma once

#include <vector>

#include <satgirgs/satgirgs_api.h>


namespace satgirgs {

SATGIRGS_API double estimateWeightScaling(const std::vector<double> &weights, double desiredAvgDegree, int dimension, double alpha);

SATGIRGS_API double estimateWeightScalingThreshold(const std::vector<double>& weights, double desiredAvgDegree, int dimension);

} // namespace girgs

