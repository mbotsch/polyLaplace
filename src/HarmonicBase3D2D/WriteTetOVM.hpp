#pragma once

#include <Eigen/Dense>
#include <vector>

void writeTetOVM(const std::vector<Eigen::Vector3d>& pts, const std::vector<std::vector<int>>& tets, const std::string& fname);
