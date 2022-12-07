#pragma once
#include <Eigen/Dense>
#include <iostream>
#include <vector>

class PolygonSampler {
    const int N = 1024;
    
    int n;
    
    bool pointInPolygon(const Eigen::Vector2d& p);
    
    std::vector<Eigen::Vector2d> edgeVectors;
    
    
public:
    static double samples[1024][2];
    
    std::vector<Eigen::Vector2d> polySamples;
    
    PolygonSampler(const Eigen::Matrix<double, -1, 2>& poly, const double area, const double avgNumSamples);
};

