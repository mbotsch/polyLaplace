/*=====================================================================================*/
/*! \file		GaussQuadrature.h
	\author		peterkau
	\brief		Declaration of class GaussQuadrature
 */
/*=====================================================================================*/
#pragma once

#include <functional>
#include <Eigen/Dense>
#include <vector>

class GaussQuadrature
{
public:
	/*! Gets the \c i-th point of the \c totalEvalPoints - point gauss quadrature rule. \c i must be in
		0 .. \c totalEvalPoints - 1. \c position is in the range [-1, 1] and the weights sum up to 2. */
	static void GetEvalPoint(int totalEvalPoints, int i, double &weight, double &position);

	static void GetEvalPoint(int totalEvalPoints, int i, float &weight, float &position) {
		double dWeight;
		double dPosition;

		GetEvalPoint(totalEvalPoints, i, dWeight, dPosition);

		weight = (float)dWeight;
		position = (float)dPosition;
	}
};


double gridQuadrature(const Eigen::Vector3d& min, const Eigen::Vector3d& max,
                  //    std::function<bool(Eigen::Vector3d)> inside,
                      std::function<double(Eigen::Vector3d)> &fun, int pts = 10);


class PolyhedralQuadrature {
    
    std::vector<Eigen::Vector3d> p;
    std::vector<double> w;
    
    void add(const Eigen::Vector3d& pi, double wi);
    
public:
    
    double apply(const std::function<double(Eigen::Vector3d)>& fun);
    
    PolyhedralQuadrature(const Eigen::Vector3d& min, const Eigen::Vector3d& max,
                         std::function<bool(Eigen::Vector3d)> &inside,
                         int pts, int ptsg);
};



class TetrahedralQuadrature {
    
    static int N;
    static double quadPoints[8][5];
    
    std::vector<Eigen::MatrixXd> tets;
    std::vector<double> volumes;
public:
    
    double apply(const std::function<double(Eigen::Vector3d)>& fun);
    
    TetrahedralQuadrature();
    
    explicit TetrahedralQuadrature(const std::vector<Eigen::MatrixXd>& tets);
};


double gaussQuadratureTriangle(const Eigen::MatrixXd& tri, const std::function<double(double, double)>& g, int n);
