#include "Tetrahedralize.h"
#include "../Volume/diffgeo_3D.h"
#include "../Surface/diffgeo.h"

void tetrahedralize(const std::vector<Eigen::Vector3d>& points,
                    const std::vector<std::vector<int>>& faces,
                    std::vector<Eigen::MatrixXd>& tets)
{
    std::vector<Eigen::Vector3d> outPts;
    std::vector<std::vector<int>> outTets;
    tetrahedralize(points, faces, outPts, outTets);

    for (auto& t : outTets)
    {
        Eigen::MatrixXd t2(4, 3);
        for (int i = 0; i < 4; ++i)
        {
            t2.row(i) = outPts[t[i]].transpose();
        }
        tets.push_back(t2);
    }
}

double vol(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
           const Eigen::Vector3d& c, const Eigen::Vector3d& d)
{
    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;
    return A.determinant() / 6.;
}

Eigen::Vector3d cross(const Eigen::Vector3d& v0, Eigen::Vector3d& v1)
{
    Eigen::Vector3d cross = Eigen::Vector3d(v0[1] * v1[2] - v0[2] * v1[1],
                                            v0[2] * v1[0] - v0[0] * v1[2],
                                            v0[0] * v1[1] - v0[1] * v1[0]);
    return cross;
}

void tetrahedralize(const std::vector<Eigen::Vector3d>& points,
                    const std::vector<std::vector<int>>& faces,
                    std::vector<Eigen::Vector3d>& outPts,
                    std::vector<std::vector<int>>& outTets)
{
    using namespace std;
    outTets.clear();
    outPts.clear();
    int np = (int)points.size();
    outPts = points;
    std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>
        triangles;
    int np_face;
    for (const auto& face : faces)
    {
        //individual polygon faces of polyhedron
        np_face = (int)face.size();
        Eigen::MatrixXd poly;
        poly.resize(np_face, 3);
        for (int i = 0; i < np_face; ++i)
        {
            poly.row(i) = points[face[i]];
        }
        Eigen::VectorXd w;
        find_area_minimizer_weights(poly, w);
        Eigen::Vector3d areamin_point = poly.transpose() * w;
        outPts.emplace_back(areamin_point);
        Eigen::Matrix3d triangle_pos;
        for (int i = 0; i < np_face; i++)
        {
            Eigen::Vector3d p0, p1, p2;
            p0 = points[face[i]];
            p1 = points[face[(i + 1) % (np_face)]];
            p2 = areamin_point;
            triangle_pos.row(0) = p2;
            triangle_pos.row(1) = p0;
            triangle_pos.row(2) = p1;
            triangles.emplace_back(triangle_pos);
        }
    }
    Eigen::Vector3d vol_minimizer;
    find_volume_minimizing_point(triangles, vol_minimizer);
    outPts.emplace_back(vol_minimizer);
    int tet_np = (int)outPts.size();
    int i = 0;
    for (const auto& face : faces)
    {
        np_face = (int)face.size();
        for (int j = 0; j < (int)face.size(); j++)
        {
            std::vector<int> tet;
            if (vol(outPts[face[j]], outPts[face[(j + 1) % (np_face)]],
                    outPts[np + i], outPts[tet_np - 1]) > 0)
            {
                tet.emplace_back(face[j]);
                tet.emplace_back(face[(j + 1) % (np_face)]);
                tet.emplace_back(np + i);
                tet.emplace_back(tet_np - 1);
                outTets.emplace_back(tet);
            }
            else
            {
                tet.emplace_back(face[(j + 1) % (np_face)]);
                tet.emplace_back(face[j]);
                tet.emplace_back(np + i);
                tet.emplace_back(tet_np - 1);
                outTets.emplace_back(tet);
                if (vol(outPts[tet_np - 1], outPts[face[(j + 1) % (np_face)]],
                        outPts[np + i], outPts[face[j]]) < 0)
                {
                    std::cout << "negative volume!" << std::endl;
                }
            }
        }
        i++;
    }
}
