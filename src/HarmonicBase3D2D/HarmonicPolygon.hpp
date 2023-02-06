#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>



class HarmonicPolygon {

    int n = 0;
    int nnKernels;
    int nnProbes;

    Eigen::MatrixXd kernels;
    Eigen::MatrixXd probes;
//
//    const int nKernels =  20;// 5; // kernels per edge;
//    const int nProbes = 60;//15; // probes per edge;
//    const double eps = 1e-1;
//    const int nKernels = 2 * 10;// 5; // kernels per edge;
//    const int nProbes = 6 * 30;//15; // probes per edge;
    int nKernels = 20; // kernels per edge;
    int nProbes =80; // probes per edge;
//    int nKernels = 4 * 10;// 5; // kernels per edge;
//    int nProbes = 8 * 30;//15; // probes per edge;
    const double eps = 1e-3;


    Eigen::MatrixXd coefficients; // for each of the 'n' functions: nKernels + 3 coefficients. One for each kernel and three for the linear part.
    
    Eigen::SparseMatrix<double> Pr;
    
    double area = .0;
    double area2d = .0;
    double scale = .0;
    
    Eigen::Vector3d mean;
   
    Eigen::Vector3d a;
    double dPlane;
    
    Eigen::MatrixXd P;
    
    std::vector<Eigen::Vector2d> samples;
    
    std::vector<Eigen::MatrixXd> quadratureTriangles;
    
    void initQuadrature();
    
public:
    
    Eigen::Matrix<double, -1, 2> poly2d;

    Eigen::Vector3d getScaledNormal() const;
    
    double getD() const;
    
    const Eigen::Vector3d& getNormal() const;
    
    Eigen::Vector2d projectToPlane(const Eigen::Vector3d& p);
    
    Eigen::Vector3d unproject(const Eigen::Vector2d& p);
    
    std::vector<Eigen::Vector3d> unproject(const std::vector<Eigen::Vector2d>& p);
    
    std::vector<Eigen::Vector2d> generateRandomSamples(int nSamples) const;
  
    double evaluate(int id, const Eigen::Vector2d& p);
    
    Eigen::Vector2d evaluateGrad(int id, const Eigen::Vector2d& p);
      
    std::vector<double> evaluateAtPoints(int id, const std::vector<Eigen::Vector2d>& points);
    
    void dump2d(const std::string& fname);
    
    explicit HarmonicPolygon(const Eigen::MatrixXd& pts_);

    explicit HarmonicPolygon(const Eigen::MatrixXd& pts_,int nKernel, int nProbes);

    double quadrature(const std::function<double(double, double)>& g);
    
    void stiffnessMatrix(Eigen::MatrixXd& K);
    
    void massMatrix(Eigen::MatrixXd& M);

};
