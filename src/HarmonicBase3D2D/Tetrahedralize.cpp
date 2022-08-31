#include "Tetrahedralize.hpp"
#include "tetgen.h"



void
tetrahedralize(const std::vector<Eigen::Vector3d>& points,
               const std::vector<std::vector<int>>& faces,
               std::vector<Eigen::MatrixXd>& tets, const double area)
{
    
    std::vector<Eigen::Vector3d> outPts;
    std::vector<std::vector<int>> outTets;
    tetrahedralize(points, faces, outPts, outTets, area);

    for(auto& t : outTets) {
        Eigen::MatrixXd t2(4, 3);
        
        for(int i = 0; i < 4; ++i) {
            t2.row(i) = outPts[t[i]].transpose();
        }
        
        tets.push_back(t2);
    }
}

void
tetrahedralize(const std::vector<Eigen::Vector3d>& points,
               const std::vector<std::vector<int>>& faces,
               std::vector<Eigen::Vector3d>& outPts,
               std::vector<std::vector<int>>& outTets, const double area)
{
    using namespace std;

    tetgenio in, out;
    
    in.firstnumber = 0;
    in.numberofpoints = (int)points.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    

    for(int i = 0; i < in.numberofpoints; ++i)
        for(int j = 0; j < 3; ++j)
            in.pointlist[3 * i + j] = points[i](j);
    
    in.numberoffacets = (int)faces.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    fill_n(in.facetmarkerlist, in.numberoffacets, 1);
    
    
    for(int i = 0; i < in.numberoffacets; ++i)
    {
        auto f = &in.facetlist[i];
        f->numberofpolygons = 1;
        f->numberofholes = 0;
        f->polygonlist = new tetgenio::polygon;
        f->holelist = nullptr;
        auto p = f->polygonlist;
        p->numberofvertices = (int)faces[i].size();
        p->vertexlist = new int[p->numberofvertices];
        
        for(int j = 0; j < p->numberofvertices; ++j) {
            p->vertexlist[j] = faces[i][j];
        }
    }
    
    tetgenbehavior settings;
    string opts = string("pq1.2a") + to_string(area);
    
    settings.parse_commandline((char*)opts.c_str());
    settings.quiet = 1;
    
    tetrahedralize(&settings, &in, &out);
    
    for(int i = 0; i < out.numberoftetrahedra; ++i) {
        std::vector<int> tet(4);
        
        for(int j = 0; j < 4; ++j) {
            tet[j] = out.tetrahedronlist[4 * i + j];
        }
        
        outTets.push_back(tet);
    }
    
    for(int i = 0; i < out.numberofpoints; ++i) {
        Eigen::Vector3d p;
        
        for(int k = 0; k < 3; ++k) {
            p(k) = out.pointlist[3 * i + k];
        }
        
        outPts.push_back(p);
    }
}
