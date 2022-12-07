#include "Geometry.hpp"

bool triangleRayIntersection(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& orig, const                            Eigen::Vector3d& dir, double& dist) {
    
    const double eps = 1e-4;
    
    Eigen::Vector3d v0v1 = v1 - v0;
    Eigen::Vector3d v0v2 = v2 - v0;
    Eigen::Vector3d pvec = dir.cross(v0v2);
    double det = v0v1.dot(pvec);
    
    if(abs(det) < eps)  return false;

    double invDet = 1. / det;
    Eigen::Vector3d tvec = orig - v0;
    double u = tvec.dot(pvec) * invDet;
    if (u < 0 || u > 1) return false;
    
    Eigen::Vector3d qvec = tvec.cross(v0v1);
    double v = dir.dot(qvec) * invDet;
    if (v < 0 || u + v > 1) return false;
    
    dist = v0v2.dot(qvec) * invDet;
    return true;
}

bool trianglePointOrientation(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& orig) {
    Eigen::Matrix3d A;
    A.row(0) = v0 - orig;
    A.row(1) = v1 - orig;
    A.row(2) = v2 - orig;
    
    return A.determinant() < 0;
}
