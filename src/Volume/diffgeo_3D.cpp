

#include "diffgeo_3D.h"
#include "../Surface/diffgeo.h"
#include <Eigen/StdVector>
#include <VolumeMesh.h>

//-----------------------------------------------------------------------------
enum VolumePoints
{
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

enum AreaPoints
{
    Quadratic_Areas_ = 0,
    Face_Centroid = 1
};

double volume_tetrahedron_signed(VolumeMesh::PointT i, VolumeMesh::PointT j,
                                 VolumeMesh::PointT k, VolumeMesh::PointT l)
{
    Eigen::Vector3d a, b, c, d;
    a << i[0], i[1], i[2];
    b << j[0], j[1], j[2];
    c << k[0], k[1], k[2];
    d << l[0], l[1], l[2];

    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;

    return A.determinant() / 6.;
}

//-----------------------------------------------------------------------------

void compute_3D_virtual_points_and_weights(VolumeMesh& mesh, int face_point,
                                           int cell_point)
{
    // save barycenter vertex for faces and cells
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");

    // we precompute them so they are ordered
    for (auto c_it = mesh.c_iter(); c_it.valid(); ++c_it)
    {
        //        std::vector<Eigen::Matrix3d> triangles; // Container for triangle vertex positions, needed for volume minimizer
        //        std::vector<Eigen::Vector3d> cell_v; // Container for cell vertices

        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d>>
            triangles;
        std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>>
            cell_v;

        auto c_f_it_pair = mesh.cell_faces(*c_it);
        for (auto c_f_it = c_f_it_pair.first; c_f_it != c_f_it_pair.second;
             ++c_f_it)
        {
            int val = (int)mesh.face(*c_f_it)
                          .halfedges()
                          .size();    //valence of the face
            Eigen::VectorXd w_f(val); // weight vector for virtual face vertex
            Eigen::MatrixXd poly(val, 3); // Polygon vertex positions
            auto f_v_it_pair = mesh.face_vertices(*c_f_it);
            int i = 0;
            for (auto f_v_it = f_v_it_pair.first; f_v_it != f_v_it_pair.second;
                 ++f_v_it)
            {
                VolumeMesh::PointT p = mesh.vertex(*f_v_it);
                for (int h = 0; h < 3; h++)
                {
                    poly.row(i)(h) = p[h];
                }
                i++;
            }
            Eigen::Vector3d min;
            if (face_point == Quadratic_Areas_)
            {
                find_area_minimizer_weights(poly, w_f);
                min = poly.transpose() *
                      w_f; //compute point obtained through weights
                VolumeMesh::PointT virtual_point(min(0), min(1), min(2));
                f_prop[*c_f_it] = virtual_point;
                f_w_prop[*c_f_it] = w_f;
            }
            else if (face_point == Face_Centroid)
            {
                VolumeMesh::PointT bary = mesh.barycenter(*c_f_it);
                f_prop[*c_f_it] = bary;
                min << bary[0], bary[1], bary[2];
                w_f = Eigen::MatrixXd::Ones(val, 1);
                w_f /= double(val);
                f_w_prop[*c_f_it] = w_f;
            }

            cell_v.emplace_back(
                min); // Virtual face vertex is part of the affine combination for the volume minimizer

            for (int k = 0; k < val; k++)
            {
                // Triangle always contains face minimizer and vertices sharing an edge of polygon
                Eigen::Matrix3d triangle;
                triangle.row(0) = min;
                triangle.row(1) = poly.row(k);
                if (k + 1 < val)
                {
                    triangle.row(2) = poly.row(k + 1);
                }
                else
                {
                    triangle.row(2) = poly.row(0);
                }
                triangles.emplace_back(triangle);
            }
        }
        auto c_v_it_pair =
            mesh.cell_vertices(*c_it); // iterate over cell vertices
        for (auto c_v_it = c_v_it_pair.first; c_v_it != c_v_it_pair.second;
             ++c_v_it)
        {
            VolumeMesh::PointT p = mesh.vertex(*c_v_it);
            cell_v.emplace_back(Eigen::Vector3d(p[0], p[1], p[2]));
        }

        Eigen::Vector3d vol_minimizer;
        Eigen::VectorXd cell_weights;
        if (cell_point == Quadratic_Volume_)
        {
            find_volume_minimizing_point(triangles, vol_minimizer);
            find_weights_for_point_3d(cell_v, vol_minimizer, cell_weights);
            VolumeMesh::PointT virtual_point(vol_minimizer(0), vol_minimizer(1),
                                             vol_minimizer(2));
            c_prop[*c_it] = virtual_point;
            c_w_prop[*c_it] = cell_weights;
        }
        else if (cell_point == Cell_Centroid_)
        {
            int cell_val = 0;
            VolumeMesh::PointT bary = VolumeMesh::PointT(0.0, 0.0, 0.0);
            for (auto c_v : cell_v)
            {
                bary += VolumeMesh::PointT(c_v[0], c_v[1], c_v[2]);
                cell_val++;
            }
            bary /= (double)cell_val;
            c_prop[*c_it] = bary;
            cell_weights = Eigen::MatrixXd::Ones(cell_val, 1);
            cell_weights /= double(cell_val);
            c_w_prop[*c_it] = cell_weights;
        }
    }
}

void find_volume_minimizing_point(
    const std::vector<Eigen::Matrix3d,
                      Eigen::aligned_allocator<Eigen::Matrix3d>>& tetrahedrons,
    Eigen::Vector3d& point)
{
    // 288*V² = cayley menger determinant  https://en.wikipedia.org/wiki/Cayley%E2%80%93Menger_determinant
    // find minimizer for sum of squared tetrahedron volumes per cell.
    //  std::vector<Eigen::MatrixXd> &tetrahedrons contains 3x3 matrices with the positions of the faces triangle vertices. (site of tetrahedron)

    Eigen::Matrix3d J = Eigen::MatrixXd::Zero(3, 3);
    Eigen::Vector3d b_ = Eigen::MatrixXd::Zero(3, 1);

    for (const auto& tetrahedron : tetrahedrons)
    {
        Eigen::Vector3d b, c, d;
        //Points of the triangle
        b = tetrahedron.row(0);
        c = tetrahedron.row(1);
        d = tetrahedron.row(2);

        // Squared Distances from b to c (23), b to d (24) and c to d (34)
        double d23, d24, d34;
        d23 = (b - c).squaredNorm();
        d24 = (b - d).squaredNorm();
        d34 = (c - d).squaredNorm();

        // Entries in Jacobian matrix
        for (int k = 0; k < 3; k++)
        {
            for (int l = 0; l < 3; l++)
            {
                //construct Linear system J*p + b = 0;
                if (k == l)
                {
                    J(k, k) += 2.0 * d23 * d23 -
                               d23 * (d24 * 4.0 + d34 * 4.0 +
                                      8.0 * (c(k) - d(k)) * (d(k) - b(k))) +
                               2.0 * d24 * d24 -
                               d24 * (d34 * 4.0 +
                                      8.0 * (d(k) - c(k)) * (c(k) - b(k))) +
                               2.0 * d34 * d34 +
                               d34 * 8.0 * (d(k) - b(k)) * (c(k) - b(k));
                }
                else
                {
                    J(k, l) += -d23 * (4.0 * (d(k) - b(k)) * (c(l) - d(l)) +
                                       4.0 * (c(k) - d(k)) * (d(l) - b(l))) -
                               d24 * (4.0 * (c(k) - b(k)) * (d(l) - c(l)) +
                                      4.0 * (d(k) - c(k)) * (c(l) - b(l))) +
                               d34 * (4.0 * (c(k) - b(k)) * (d(l) - b(l)) +
                                      4.0 * (d(k) - b(k)) * (c(l) - b(l)));
                }
            }
            b_(k) +=
                -2.0 * d23 * d23 * d(k) -
                d23 *
                    (d24 * (-2.0 * (c(k) + d(k))) + d34 * -2.0 * (b(k) + d(k)) +
                     2.0 * (d(k) - b(k)) * (d.squaredNorm() - c.squaredNorm()) +
                     2.0 * (c(k) - d(k)) *
                         (b.squaredNorm() - d.squaredNorm())) -
                2.0 * d24 * d24 * c(k) -
                d24 *
                    (-d34 * 2.0 * (b(k) + c(k)) +
                     2.0 * (c(k) - b(k)) * (c.squaredNorm() - d.squaredNorm()) +
                     2.0 * (d(k) - c(k)) *
                         (b.squaredNorm() - c.squaredNorm())) -
                2.0 * d34 * d34 * b(k) +
                d34 *
                    (2.0 * (c(k) - b(k)) * (b.squaredNorm() - d.squaredNorm()) +
                     2.0 * (d(k) - b(k)) * (b.squaredNorm() - c.squaredNorm()));
        }
    }
    //solve linear system (there are probably better solver for this)
    point = (J).completeOrthogonalDecomposition().solve(-b_);
}
//-----------------------------------------------------------------------------

//More or less the same code as before, but without Eigen Matrix for polygon.
//Should be cleaned up
void find_weights_for_point_3d(
    const std::vector<Eigen::Vector3d,
                      Eigen::aligned_allocator<Eigen::Vector3d>>& poly,
    const Eigen::Vector3d& point, Eigen::VectorXd& weights)
{
    const int n = (int)poly.size();
    Eigen::MatrixXd M(3, n - 1), M_(3, n);
    const Eigen::Vector3d& pn = poly[n - 1];
    M_.col(n - 1) = poly[n - 1];
    for (int i = 0; i < n - 1; i++)
    {
        M.col(i) = poly[i] - pn;
        M_.col(i) = poly[i];
    }
    Eigen::VectorXd w;

    weights.resize(n);
    Eigen::VectorXd b(3);
    b = point - pn;

    w = M.completeOrthogonalDecomposition().solve(b);
    weights.topRows(n - 1) = w;
    weights(n - 1) = 1.0 - w.sum();
}

//-----------------------------------------------------------------------------

VolumeMesh::PointT compute_triangle_normal(VolumeMesh::PointT a,
                                           VolumeMesh::PointT b,
                                           VolumeMesh::PointT c)
{
    VolumeMesh::PointT n = (b - a) % (c - a);
    return n.normalize();
}
