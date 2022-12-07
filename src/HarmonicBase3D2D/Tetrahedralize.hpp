#pragma once

#include <Eigen/Dense>
#include <vector>

void
tetrahedralize(const std::vector<Eigen::Vector3d>& points,
               const std::vector<std::vector<int>>& faces,
               std::vector<Eigen::Vector3d>& outPts,
               std::vector<std::vector<int>>& outTets);

void
tetrahedralize(const std::vector<Eigen::Vector3d>& points,
               const std::vector<std::vector<int>>& faces,
               std::vector<Eigen::MatrixXd>& tets);
