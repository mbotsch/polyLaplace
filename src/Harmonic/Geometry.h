#pragma once

#include <Eigen/Dense>

bool triangleRayIntersection(const Eigen::Vector3d& v0,
                             const Eigen::Vector3d& v1,
                             const Eigen::Vector3d& v2,
                             const Eigen::Vector3d& orig,
                             const Eigen::Vector3d& dir, double& dist);

bool trianglePointOrientation(const Eigen::Vector3d& v0,
                              const Eigen::Vector3d& v1,
                              const Eigen::Vector3d& v2,
                              const Eigen::Vector3d& orig);
