//=============================================================================

#include "AQAPoly_Laplacian.h"
#include <Misha/Exceptions.h>
#include <Misha/Geometry.h>
#include<Misha/Miscellany.h>
#include <Misha/SimplexRefinableMesh.h>
#include <Misha/Meshes.h>
#include <Misha/MGSolver.h>
#include <igl/slice.h>
#include <Eigen/CholmodSupport>
#include "basis.h"
#include "diffgeo.h"
#include "unsupported/Eigen/SparseExtra"
#include "Spectra/MatOp/SymShiftInvert.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "Spectra/SymGEigsShiftSolver.h"
#include "Spectra/SymGEigsSolver.h"
#include "Spectra/MatOp/SparseCholesky.h"
//=============================================================================

enum InsertedPoint
{
    Centroid = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2,
};


//=============================================================================


void setup_AQAPoly_matrices(pmp::SurfaceMesh &mesh,
                            Eigen::SparseMatrix<double> &M, MatrixType matrix,
                            CoarseDimension dof, int minpoint,
                            int degree)
{
    const int nv = mesh.n_vertices();
    const int ne = mesh.n_edges();
    std::vector<Eigen::Triplet<double>> trip;
    int coarseDim;
    if (dof == Vertices){
        coarseDim = 0;
    }else if(dof == Edges){
        coarseDim = 1;
    }
    if(degree == 1)
    {
        setup_local_AQAPoly_matrices<1>(mesh,trip, matrix, dof,minpoint);
    }else if(degree == 2)
    {
        setup_local_AQAPoly_matrices<2>(mesh,trip, matrix, dof,minpoint);
    }
    else if(degree == 3)
    {
        setup_local_AQAPoly_matrices<3>(mesh,trip, matrix, dof,minpoint);
    }
    if (degree == 1 || dof == Vertices)
    {
        M.resize(nv, nv);
        M.setFromTriplets(trip.begin(), trip.end());
    }
    else
    {
        M.resize(nv + ne, nv + ne);
        M.setFromTriplets(trip.begin(), trip.end());
    }
    if (matrix == Stiffness)
    {
        M *= -1.0;
    }
}

//----------------------------------------------------------------------------
//
template<unsigned int Degree>
void setup_local_AQAPoly_matrices(pmp::SurfaceMesh &mesh,std::vector<Eigen::Triplet<double>> &triplet, MatrixType matrix, CoarseDimension dof,int minpoint ){
    bool squaredAreaMin = true;
    if (minpoint != AreaMinimizer)
    {
        squaredAreaMin = false;
    }
    typedef typename SimplexRefinableElements<2,Degree>::NodeMultiIndex NodeMultiIndex;

    for (pmp::Face f : mesh.faces())
    {
        std::vector<::Point<double, 3>> vertices;
        std::vector<std::vector<unsigned int>> polygons;
        std::vector<unsigned int> face;
        int i = 0;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {

            pmp::Vertex v = mesh.from_vertex(he);
            vertices.emplace_back(::Point<double, 3>(mesh.position(v)[0],
                                                     mesh.position(v)[1],
                                                     mesh.position(v)[2]));
            face.emplace_back(i);
            i++;
        }
        polygons.emplace_back(face);
        Meshes::PolygonMesh<unsigned int>  polyMesh = Meshes::PolygonMesh<unsigned int>(polygons);

        SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
        eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

        std::function< Point< double ,3> ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
        {
            if( idx<vertices.size() ) return vertices[idx];
            ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
            return Point< double ,3 >();
        };
        std::function< Point< double , 3 > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction ,  squaredAreaMin  );

        int coarse_dim;
        switch (dof)
        {
            case Vertices:
                coarse_dim = 0;
                break;
            case Edges:
                coarse_dim = 1;
                break;
            case Refined_mesh:
                coarse_dim = 2;
                break;
            default:
                std::cout << "should not happen";
        }
        // solve system directly
        SimplexRefinableCellMesh< 2 , Degree > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim , false );
        Eigen::MatrixXd S_;
        if (matrix == Stiffness)
        {
            S_ = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
        }
        else
        {
            S_ =simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
        }
        if( dof == Vertices){
            local_matrix_to_global(mesh, f, S_, triplet, Vertices);
        }  else if(dof == Edges){
            local_matrix_to_global(mesh, f, S_, triplet, Edges);
        }
    }
}

//=============================================================================

void local_matrix_to_global(pmp::SurfaceMesh &mesh, pmp::Face &f,
                            Eigen::MatrixXd &M,
                            std::vector<Eigen::Triplet<double>> &triplet,
                            int degree)
{
    int i, j;
    if (degree == 1)
    {
        i = 0;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            j = 0;
            for (pmp::Halfedge hhe : mesh.halfedges(f))
            {
                pmp::Vertex vv = mesh.from_vertex(hhe);
                triplet.emplace_back(vv.idx(), v.idx(), M(i, j));
                j++;
            }
            i++;
        }
    }
    else
    {
        int nv = mesh.n_vertices();

        int v_idx, vv_idx, e_idx, ee_idx;

        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            pmp::Edge e = mesh.edge(he);

            v_idx = vertex_node_to_idx(mesh, f, he);
            e_idx = edge_node_to_idx(mesh, f, he);
            for (pmp::Halfedge hhe : mesh.halfedges(f))
            {
                pmp::Vertex vv = mesh.from_vertex(hhe);
                pmp::Edge ee = mesh.edge(hhe);

                vv_idx = vertex_node_to_idx(mesh, f, hhe);
                ee_idx = edge_node_to_idx(mesh, f, hhe);

                // The enumeration is circular (v_i = 2*i, e_i = 2*i+1, v_i+1 = 2*(i+1) ... etc.) until the last two edge elements.
                // Reason: the first vertex node starts with (0,0), the last vertex is (val-1,val-1).
                // The edge between those vertices is (val-1,0), the edge node leading to the last vertex is (val-1,val-2). Mishas ordering goes from
                // smallest to tallest, therefore the edge indices  in this special case are swapped and not ciruclar as before. the last vertex gets the highest idx.

                triplet.emplace_back(v.idx(), vv.idx(), M(v_idx, vv_idx));
                triplet.emplace_back(nv + e.idx(), vv.idx(), M(e_idx, vv_idx));

                triplet.emplace_back(v.idx(), nv + ee.idx(), M(v_idx, ee_idx));
                triplet.emplace_back(nv + e.idx(), nv + ee.idx(),
                                     M(e_idx, ee_idx));
            }
        }
    }
}
void local_matrix_to_global(pmp::SurfaceMesh &mesh, pmp::Face &f,
                            Eigen::MatrixXd &M,
                            std::vector<Eigen::Triplet<double>> &triplet,
                            CoarseDimension dof)
{
    int i, j;
    if (dof == Vertices)
    {
        i = 0;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            j = 0;
            for (pmp::Halfedge hhe : mesh.halfedges(f))
            {
                pmp::Vertex vv = mesh.from_vertex(hhe);
                triplet.emplace_back(vv.idx(), v.idx(), M(i, j));
                j++;
            }
            i++;
        }
    }
    else
    {
        int nv = mesh.n_vertices();

        int v_idx, vv_idx, e_idx, ee_idx;

        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            pmp::Edge e = mesh.edge(he);

            v_idx = vertex_node_to_idx(mesh, f, he);
            e_idx = edge_node_to_idx(mesh, f, he);
            for (pmp::Halfedge hhe : mesh.halfedges(f))
            {
                pmp::Vertex vv = mesh.from_vertex(hhe);
                pmp::Edge ee = mesh.edge(hhe);

                vv_idx = vertex_node_to_idx(mesh, f, hhe);
                ee_idx = edge_node_to_idx(mesh, f, hhe);

                // The enumeration is circular (v_i = 2*i, e_i = 2*i+1, v_i+1 = 2*(i+1) ... etc.) until the last two edge elements.
                // Reason: the first vertex node starts with (0,0), the last vertex is (val-1,val-1).
                // The edge between those vertices is (val-1,0), the edge node leading to the last vertex is (val-1,val-2). Mishas ordering goes from
                // smallest to tallest, therefore the edge indices  in this special case are swapped and not ciruclar as before. the last vertex gets the highest idx.

                triplet.emplace_back(v.idx(), vv.idx(), M(v_idx, vv_idx));
                triplet.emplace_back(nv + e.idx(), vv.idx(), M(e_idx, vv_idx));

                triplet.emplace_back(v.idx(), nv + ee.idx(), M(v_idx, ee_idx));
                triplet.emplace_back(nv + e.idx(), nv + ee.idx(),
                                     M(e_idx, ee_idx));
            }
        }
    }
}

int edge_node_to_idx(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Halfedge &he)
{
    int i = 0;
    pmp::Vertex v = mesh.from_vertex(he);
    for (auto hhe : mesh.halfedges(f))
    {
        pmp::Vertex vv = mesh.from_vertex(hhe);
        if ((i == 2 * mesh.valence(f) - 4) && (v.idx() == vv.idx()))
        {
            return 2 * mesh.valence(f) - 2;
        }
        else if (i == 2 * mesh.valence(f) - 2)
        {
            return 2 * mesh.valence(f) - 3;
        }
        else if (v.idx() == vv.idx())
        {
            return i + 1;
        }
        i += 2;
    }
}

int vertex_node_to_idx(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Halfedge &he)
{
    int i = 0;
    pmp::Vertex v = mesh.from_vertex(he);
    for (auto hhe : mesh.halfedges(f))
    {
        pmp::Vertex vv = mesh.from_vertex(hhe);
        if ((i == 2 * mesh.valence(f) - 2) && (v.idx() == vv.idx()))
        {
            return 2 * mesh.valence(f) - 1;
        }
        else if (v.idx() == vv.idx())
        {
            return i;
        }
        i += 2;
    }
}

//================ Explicit quadratic triangle stiffness/mass matrices=====================================

void setup_triangle_FEM_stiffness_matrix(pmp::SurfaceMesh &mesh,
                                         Eigen::SparseMatrix<double> &S,
                                         bool linear)
{
    const int nv = mesh.n_vertices();
    const int ne = mesh.n_edges();

    std::vector<Eigen::Triplet<double>> trip;

    for (pmp::Face f : mesh.faces())
    {
        local_triangle_fem_matrix(mesh, f, trip, linear);
    }

    if (linear)
        S.resize(nv, nv);
    else
    {
        S.resize(nv + ne, nv + ne);
    }
    S.setFromTriplets(trip.begin(), trip.end());
    S *= -1.0;
}

void local_triangle_fem_matrix(pmp::SurfaceMesh &mesh, pmp::Face f,
                               std::vector<Eigen::Triplet<double>> &trip,
                               bool linear)
{
    int nb = 6;
    pmp::vec2 z0, z1, z2;
    Eigen::MatrixXd Jacobian;
    Jacobian.resize(2, 2);
    std::vector<pmp::Point> points;
    std::vector<int> global_idx(6);
    Eigen::MatrixXd basis; // rows = nr basis, cols = dim
    if (linear)
    {
        nb = 3;
        global_idx.resize(3);
        linear_basis_gradient(basis);
    }
    int vidx = 0;
    for (auto he : mesh.halfedges(f))
    {
        pmp::Vertex v = mesh.from_vertex(he);
        points.push_back(mesh.position(v));
        global_idx[vidx] = v.idx();
        if (!linear)
            global_idx[3 + vidx] = mesh.n_vertices() + mesh.edge(he).idx();
        vidx++;
    }
    project_triangle(points[0], points[1], points[2], z0, z1, z2);

    Jacobian(0, 0) = z2[1] - z0[1];
    Jacobian(0, 1) = z0[0] - z2[0];
    Jacobian(1, 0) = z0[1] - z1[1];
    Jacobian(1, 1) = z1[0] - z0[0];

    double det =
        (z2[1] - z0[1]) * (z1[0] - z0[0]) - (z0[0] - z2[0]) * (z0[1] - z1[1]);
    Jacobian /= det;

    for (int j = 0; j < nb; j++)
    {
        for (int k = 0; k < nb; k++)
        {
            double val;
            if (!linear)
                val = gauss_quadrature_gradient(Jacobian, j, k);
            else
            {
                val = 0.5 * (basis.row(j)) * Jacobian * Jacobian.transpose() *
                      (basis.row(k).transpose());
            }
            int row, col;
            row = global_idx[j];
            col = global_idx[k];

            trip.emplace_back(row, col, val * det);
        }
    }
}

double TriangleJacobian(Eigen::Matrix3d Triangle, Eigen::Matrix2d &Jacobian){
    Eigen::Matrix2d Triangle2d;
    project_triangle(Triangle, Triangle2d);
    Eigen::Vector2d z0 = Triangle2d.row(0);
    Eigen::Vector2d z1 = Triangle2d.row(1);
    Eigen::Vector2d z2 = Triangle2d.row(2);

    Jacobian(0, 0) = z2(1) - z0(1);
    Jacobian(0, 1) = z0(0) - z2(0);
    Jacobian(1, 0) = z0(1)  - z1(1) ;
    Jacobian(1, 1) = z1(0) - z0(0);

    double det =
        (z2(1) - z0(1)) * (z1(0) - z0(0)) - (z0(0) - z2(0)) * (z0(1) - z1(1));
    Jacobian /= det;
    return det;
}
void setup_triangle_FEM_mass_matrix(pmp::SurfaceMesh &mesh,
                                    Eigen::SparseMatrix<double> &M)
{
    const int nv = mesh.n_vertices();
    const int ne = mesh.n_edges();

    std::vector<Eigen::Triplet<double>> trip;

    for (pmp::Face f : mesh.faces())
    {
        local_triangle_mass_fem_matrix(mesh, f, trip);
    }

    M.resize(nv + ne, nv + ne);

    M.setFromTriplets(trip.begin(), trip.end());
}

void local_triangle_mass_fem_matrix(pmp::SurfaceMesh &mesh, pmp::Face f,
                                    std::vector<Eigen::Triplet<double>> &trip)
{
    int nb = 6;
    pmp::vec2 z0, z1, z2;
    Eigen::MatrixXd Jacobian;
    Jacobian.resize(2, 2);
    std::vector<pmp::Point> points;
    std::vector<int> global_idx(6);
    Eigen::MatrixXd basis; // rows = nr basis, cols = dim

    int vidx = 0;
    for (auto he : mesh.halfedges(f))
    {

        pmp::Vertex v = mesh.from_vertex(he);
        points.push_back(mesh.position(v));
        global_idx[vidx] = v.idx();
        global_idx[3 + vidx] = mesh.n_vertices() + mesh.edge(he).idx();
        vidx++;
    }

    project_triangle(points[0], points[1], points[2], z0, z1, z2);

    // For the mass only the determinant of the Jacobian is needed
    double det =
        (z2[1] - z0[1]) * (z1[0] - z0[0]) - (z0[0] - z2[0]) * (z0[1] - z1[1]);

    for (int j = 0; j < nb; j++)
    {
        for (int k = 0; k < nb; k++)
        {
            double val = gauss_quadrature_mass(j, k);
            int row, col;
            row = global_idx[j];
            col = global_idx[k];
            trip.emplace_back(row, col, val * det);
        }
    }
}

/**
     * \brief Computes the coordinates of the vertices of a triangle
     * in a local 2D orthonormal basis of the triangle's plane.
     * \param[in] p0 , p1 , p2 the 3D coordinates of the vertices of
     *   the triangle
     * \param[out] z0 , z1 , z2 the 2D coordinates of the vertices of
     *   the triangle
     */
void project_triangle(const pmp::Point &p0, const pmp::Point &p1,
                      const pmp::Point &p2, pmp::vec2 &z0, pmp::vec2 &z1,
                      pmp::vec2 &z2)
{
    pmp::Point X = p1 - p0;
    X.normalize(); // normalized by dividing x,y,z with length of the vector
    pmp::Point Z = cross(X, (p2 - p0));
    Z.normalize();
    pmp::Point Y = cross(Z, X); //cross product
    const pmp::Point &O = p0;

    double x0 = 0;
    double y0 = 0;
    double x1 = norm(p1 - O);
    double y1 = 0;
    double x2 = dot((p2 - O), X);
    double y2 = dot((p2 - O), Y);

    z0 = pmp::vec2(x0, y0);
    z1 = pmp::vec2(x1, y1);
    z2 = pmp::vec2(x2, y2);
}

void project_triangle(Eigen::Matrix3d &Triangle3, Eigen::Matrix2d &Triangle2)
{
    Eigen::Vector3d p0 = Triangle3.row(0);
    Eigen::Vector3d p1 = Triangle3.row(1);
    Eigen::Vector3d p2 = Triangle3.row(2);

    Eigen::Vector3d X = p1 - p0;
    X.normalize(); // normalized by dividing x,y,z with length of the vector
    Eigen::Vector3d Z = (X.cross(p2 - p0));
    Z.normalize();
    Eigen::Vector3d Y = Z.cross(X); //cross product
    Eigen::Vector3d &O = p0;

    double x0 = 0;
    double y0 = 0;
    double x1 = (p1 - O).norm();
    double y1 = 0;
    double x2 = (p2 - O).dot(X);
    double y2 = (p2 - O).dot(Y);

    Triangle2.row(0) = Eigen::Vector2d(x0, y0);
    Triangle2.row(1) = Eigen::Vector2d(x1, y1);
    Triangle2.row(2) = Eigen::Vector2d(x2, y2);
}

void problem_faces(pmp::SurfaceMesh &mesh, int minpoint)
{

}
//====================================new multigrid code====================================================
double Franke( Point< double , 2 > p )
{
    double x = p[0] , y = p[1];
    double cx2 = (9. * x - 2.) * (9. * x - 2.);
    double cy2 = (9. * y - 2.) * (9. * y - 2.);

    double cx1 = (9. * x + 1.) * (9. * x + 1.);
    double cx7 = (9. * x - 7.) * (9. * x - 7.);

    double cy3 = (9. * y - 3.) * (9. * y - 3.);
    double cx4 = (9. * x - 4.) * (9. * x - 4.);

    double cy7 = (9. * y - 7.) * (9. * y - 7.);

    return (3. / 4.) * exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2) +
           (3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10.) +
           (1. / 2.) * exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3) -
           (1. / 5.) * exp(-cx4 - cy7);
}

double FrankeLaplacian( Point< double , 2 > p )
{
    double x = p[0] , y = p[1];
    double mathematica =
        64.8 * exp(-pow(-4. + 9. * x, 2.0) - pow(-7. + 9. * y, 2.0)) -
        40.5 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) -
        60.75 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) -
        1.8720918367346937 * exp(-0.02040816326530612 * pow(1. + 9. * x, 2) -
                                 0.1 * (1. + 9. * y)) +
        10.125 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) *
            pow(-7. + 9. * x, 2) -
        64.8 * exp(-pow(-4. + 9. * x, 2) - pow(-7. + 9. * y, 2)) *
            pow(-4. + 9. * x, 2) +
        15.1875 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) *
            pow(-2. + 9. * x, 2) +
        0.1012078300708038 *
            exp(-0.02040816326530612 * pow(1. + 9. * x, 2) -
                0.1 * (1. + 9. * y)) *
            pow(1. + 9. * x, 2) -
        64.8 * exp(-pow(-4. + 9. * x, 2) - pow(-7. + 9. * y, 2)) *
            pow(-7. + 9. * y, 2) +
        10.125 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) *
            pow(-3. + 9. * y, 2) +
        15.1875 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) *
            pow(-2. + 9. * y, 2);

    return mathematica;
}
template< unsigned int Degree>
double Execute_Direct
    (
        std::function< Point< double , 2 > ( unsigned int ) > &fullVertexPositionFunction,
        const SimplexRefinableCellMesh< 2 , Degree > &simplexRefinableCellMesh ,
        const std::function< Point< double ,2 > ( typename SimplexMesh< 2 , Degree >::NodeMultiIndex ) > &NodePosition ,
        const std::function< bool ( typename SimplexMesh< 2 , Degree >::NodeMultiIndex ) > &IsBoundaryNode
    )
{
    typedef typename SimplexMesh<2 , Degree >::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S , P;
    Timer timer;

    // Get the system matrices
    P = simplexRefinableCellMesh.P();
    Eigen::saveMarket(P,"P.mtx");

    timer.reset();
    M = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
    S = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
    std::cout << "Got system matrices: " << timer.elapsed() << std::endl;
    Eigen::saveMarket(S,"S.mtx");
    Eigen::saveMarket(simplexRefinableCellMesh.simplexMesh().stiffness(),"Sf.mtx");

    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

    // Identify the nodes that are locked
    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

    Eigen::VectorXd b( simplexRefinableCellMesh.nodes() );

    // Initiailize the constraints to the evaluation of the Laplacian of the Franke function
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) b[i] = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

    // Compute the dual representation of the constraints
    b = M * b;

    // Set the constraints at the locked vertices to the evluation of the Franke function
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) if( lockedNodes[i] ) b[i] = Franke( NodePosition( nodeMultiIndices[i] ) );

    // Adjust the right-hand-side to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b[ iter.row() ] -= b[ iter.col() ] * iter.value();

    // Adjust the system matrix to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
            else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

//    Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > solver;
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver;

        // Compute the Cholesky factorization
    timer.reset();
    solver.compute( S );
    switch( solver.info() )
    {
        case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
        case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
        case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
    }

    std::cout << "Performed Cholesky (LDLt) factorization: " << timer.elapsed() << std::endl;
    std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;

    timer.reset();
    Eigen::VectorXd x = solver.solve( b );

    std::cout << "Solved direct system: " << timer.elapsed() << std::endl;
    double eIn = 0 , eOut = 0;
    Eigen::VectorXd r = b - S* x;
    for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
    std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;

    double rms = 0;
    //=============Gauss Quadrature ==============
    Eigen::VectorXd PSol = P*x;
    Eigen::MatrixXd QuadWeigths {{0.33333333333333, 0.33333333333333, -0.56250000000000},
                                      {0.2, 0.2, 0.52083333333333},
                                      {0.2, 0.6, 0.52083333333333},
                                      {0.6, 0.2, 0.52083333333333}};
    Eigen::MatrixXd Tri;

    Tri.resize(2,3);
    double l2_error = 0;
    double area_sum = 0;
    for(unsigned  int i=0 ; i< simplexRefinableCellMesh.simplexMesh().simplices();i++){
       SimplexIndex< 2 , unsigned int > s =  simplexRefinableCellMesh.simplexMesh().simplex(i);
       for( unsigned int n=0 ; n<3 ; n++)
       {
           Point<double, 2> p = fullVertexPositionFunction(s.idx[n]);
           Tri.col(n) = Eigen::Vector2d(p[0], p[1]);
       }
       Eigen::Vector2d a,b,c;
       a = Tri.col(0);
       b = Tri.col(2);
       c = Tri.col(1);
       double det = 0.5 * (-a(1) * b(0) + a(0) * b(1) + a(1) * c(0) - b(1) * c(0) - a(0) * c(1) + b(0) * c(1));
       area_sum+=det;
       //Iterate over Quadrature points
       for(unsigned int q = 0; q<QuadWeigths.rows(); q++){
           // initialize Sample
           Eigen::Vector3d bcCoordinates, globalPos;
           Point< double , 2 > p(QuadWeigths.row(q)(0),QuadWeigths.row(q)(1)) ;
           typename SimplexMesh< 2 >::Sample sample;

           sample.sIdx = i;
           sample.bcCoordinates[0] = 1. - p[0] - p[1];
           sample.bcCoordinates[1] = p[0];
           sample.bcCoordinates[2] = p[1];
           bcCoordinates <<1.-p[0]-p[1],p[0],p[1];
           globalPos = Tri*bcCoordinates;
           Point< double , 2 > globalP(globalPos[0],globalPos[1]);
           double v_exact =  Franke(globalP);
           double v_approx = simplexRefinableCellMesh.simplexMesh().evaluate(PSol,sample);
//           std::cout << "v_exact: " << v_exact << " v_approx: " << v_approx << std::endl;
           l2_error+=(v_exact-v_approx)*(v_exact-v_approx)*det*QuadWeigths.row(q)(2);
       }
    }
    std::cout << "area sum : "<<area_sum << std::endl;
    l2_error = sqrt(fabs(l2_error));

     //===========================
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
    std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes() ) << std::endl;
    std::cout << "L2 Error: " << l2_error << std::endl;
//    return sqrt( rms / simplexRefinableCellMesh.nodes() );
    return l2_error;
}

template< unsigned int Degree , typename RelaxerType >
double Execute_MG
    (
        const HierarchicalSimplexRefinableCellMesh< 2, Degree > &simplexRefinableCellMesh ,
        int CoarseNodeDimension,
        const std::function< Point< double , 2 > ( typename SimplexMesh< 2 , Degree >::NodeMultiIndex ) > &NodePosition ,
        const std::function< bool ( typename SimplexMesh< 2, Degree >::NodeMultiIndex ) > &IsBoundaryNode ,
        unsigned int vCycles , unsigned int gsIters
    )
{
    typedef typename SimplexMesh< 2 , Degree >::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S;
    std::vector< Eigen::SparseMatrix< double > > P( CoarseNodeDimension);
    Timer timer;

    // Get the system matrices
    timer.reset();
    for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension ; d++ ) P[d] = simplexRefinableCellMesh.P( d+1 , d );
    {
        Eigen::SparseMatrix< double > _P , _Pt;
        _P = simplexRefinableCellMesh.P( simplexRefinableCellMesh.maxLevel() , CoarseNodeDimension );
        _Pt = _P.transpose();
        M = _Pt * simplexRefinableCellMesh.simplexMesh().mass() * _P;
        S = _Pt * simplexRefinableCellMesh.simplexMesh().stiffness() * _P;
    }
     std::cout << "Got system matrices: " << timer.elapsed() << std::endl;


    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes( CoarseNodeDimension ) );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap( CoarseNodeDimension ) ) nodeMultiIndices[idx] = nmi;

    // Identify the nodes that are locked
    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes( CoarseNodeDimension ) );
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

    Eigen::VectorXd b( simplexRefinableCellMesh.nodes( CoarseNodeDimension) );

    // Initiailize the constraints to the evaluation of the Laplacian of the Franke function
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension ) ; i++ ) b[i] = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

    // Compute the dual representation of the constraints
    b = M * b;

    // Set the constraints at the locked vertices to the evluation of the Franke function
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) if( lockedNodes[i] ) b[i] = Franke( NodePosition( nodeMultiIndices[i] ) );

    // Adjust the right-hand-side to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b[ iter.row() ] -= b[ iter.col() ] * iter.value();

    // Adjust the system matrix to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
            else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

    // Compute the hierarchy of systems
    timer.reset();
    MGSolver::Solver< RelaxerType > mgSolver( S , P , true);
    std::cout << "Constructed multigrid system: " << timer.elapsed() << std::endl;

    std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;

        timer.reset();
        Eigen::VectorXd x = mgSolver.solve( b , vCycles , gsIters , true);

        std::cout << "Solved multigrid system: " << timer.elapsed() << std::endl;
        double eIn = 0 , eOut = 0;
        Eigen::VectorXd r = b - S* x;
        for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
        std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;


        double rms = 0;
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension ) ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
        std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension ) ) << std::endl;
    return  sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension ) ) ;
}

template<unsigned int Degree>
double solve_AQAPoly_Poisson_mg(pmp::SurfaceMesh &mesh, CoarseDimension dof, bool direct, MG_Solver solver, int vcycles, int iterations)
{
    std::vector<::Point<double, 2>> vertices;
    std::vector<std::vector<unsigned int>> polygons;
    for (pmp::Vertex v : mesh.vertices()){
        vertices.emplace_back(::Point<double, 2>(
            mesh.position(v)[0], mesh.position(v)[1]));
    }
    for (pmp::Face f : mesh.faces())
    {
        std::vector<unsigned int> face;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            face.emplace_back(v.idx());
        }
        polygons.emplace_back(face);
    }
    Meshes::PolygonMesh<unsigned int>  polyMesh = Meshes::PolygonMesh<unsigned int>(polygons);
    typedef typename SimplexMesh< 2 , Degree >::NodeMultiIndex NodeMultiIndex;

    std::set< unsigned int > boundaryVertices;
    {
        typedef MultiIndex< 2 , unsigned int , false > EdgeIndex;
        std::map< EdgeIndex , unsigned int > edgeCount;

        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ )
        {
            Meshes::Polygon< unsigned int > poly = polyMesh.polygon(i);
            for( unsigned int j=0 ; j<poly.size() ; j++ )
            {
                unsigned int idx[] = { poly[j] , poly[(j+1)%poly.size()] };
                edgeCount[ EdgeIndex(idx) ]++;
            }
        }
        for( auto const & [ edge , count ] : edgeCount ) if( count==1 ) boundaryVertices.insert( edge[0] ) , boundaryVertices.insert( edge[1] );
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , 2> ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double ,2 >();
    };
    std::function< Point< double , 2 > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , 2 > p;
        for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
        return p/(double)Degree;
    };
    int coarse_dim;
    switch (dof)
    {
        case Vertices:
            coarse_dim = 0;
            break;
        case Edges:
            coarse_dim = 1;
            break;
        case Refined_mesh:
            coarse_dim = 2;
            break;
        default:
            std::cout << "should not happen";
    }
    double error;
    // solve system directly
    if(direct)
    {
        SimplexRefinableCellMesh< 2 , Degree > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim , true );
        error = Execute_Direct< Degree >(fullVertexPositionFunction, simplexRefinableCellMesh , NodePosition , IsBoundaryNode );
    }else // solve via multigrid
    {
        HierarchicalSimplexRefinableCellMesh<2, Degree > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim);
        switch( solver )
        {
            case MG_Solver::RELAXER_JACOBI:
                error = Execute_MG< Degree , MGSolver::JacobiRelaxer>( simplexRefinableCellMesh ,coarse_dim, NodePosition , IsBoundaryNode , vcycles , iterations );
                break;
            case MG_Solver::RELAXER_GAUSS_SEIDEL:
                error = Execute_MG< Degree , MGSolver::GaussSeidelRelaxer >( simplexRefinableCellMesh,coarse_dim, NodePosition , IsBoundaryNode , vcycles , iterations );
                break;
            case MG_Solver::RELAXER_PARALLEL_GAUSS_SEIDEL:
                error = Execute_MG< Degree , MGSolver::ParallelGaussSeidelRelaxer< 20 > >( simplexRefinableCellMesh ,coarse_dim, NodePosition , IsBoundaryNode ,vcycles , iterations  );
                break;
            default: ERROR_OUT( "Unrecognized relaxer type: " , solver );
        }
    }
    return error;
}

double solve_AQAPoly_Poisson_mg(pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree, bool direct, MG_Solver solver , int vcycles , int iterations){
    if (degree == 1)
    {
        return solve_AQAPoly_Poisson_mg<1>(mesh,dof, direct, solver,vcycles, iterations);
    }
    else if (degree == 2)
    {
        return solve_AQAPoly_Poisson_mg<2>(mesh,dof, direct, solver,vcycles, iterations);
    }
    else if (degree == 3)
    {
        return solve_AQAPoly_Poisson_mg<3>(mesh,dof, direct, solver,vcycles, iterations);
    }
}

template<unsigned int Degree>
double solve_2D_AQAPoly_Poisson(pmp::SurfaceMesh &mesh, CoarseDimension dof){
    std::vector<::Point<double, 2>> vertices;
    std::vector<std::vector<unsigned int>> polygons;
    for (pmp::Vertex v : mesh.vertices()){
        vertices.emplace_back(::Point<double, 2>(
            mesh.position(v)[0], mesh.position(v)[1]));
    }
    for (pmp::Face f : mesh.faces())
    {
        std::vector<unsigned int> face;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            face.emplace_back(v.idx());
        }
        polygons.emplace_back(face);
    }
    Meshes::PolygonMesh<unsigned int>  polyMesh = Meshes::PolygonMesh<unsigned int>(polygons);
    typedef typename SimplexMesh< 2 , Degree >::NodeMultiIndex NodeMultiIndex;

    std::set< unsigned int > boundaryVertices;
    {
        typedef MultiIndex< 2 , unsigned int , false > EdgeIndex;
        std::map< EdgeIndex , unsigned int > edgeCount;

        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ )
        {
            Meshes::Polygon< unsigned int > poly = polyMesh.polygon(i);
            for( unsigned int j=0 ; j<poly.size() ; j++ )
            {
                unsigned int idx[] = { poly[j] , poly[(j+1)%poly.size()] };
                edgeCount[ EdgeIndex(idx) ]++;
            }
        }
        for( auto const & [ edge , count ] : edgeCount ) if( count==1 ) boundaryVertices.insert( edge[0] ) , boundaryVertices.insert( edge[1] );
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , 2> ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double ,2 >();
    };
    std::function< Point< double , 2 > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , 2 > p;
        for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
        return p/(double)Degree;
    };
    int coarse_dim;
    switch (dof)
    {
        case Vertices:
            coarse_dim = 0;
            break;
        case Edges:
            coarse_dim = 1;
            break;
        case Refined_mesh:
            coarse_dim = 2;
            break;
        default:
            std::cout << "should not happen";
    }

    SimplexRefinableCellMesh< 2 , Degree > simplexRefinableCellMesh;
    simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim , true );
    typedef typename SimplexMesh<2 , Degree >::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S , P;
    Timer timer;

    // Get the system matrices
    P = simplexRefinableCellMesh.P();
    Eigen::saveMarket(P,"P.mtx");
    timer.reset();
    M = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
    S = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
    std::cout << "Got system matrices: " << timer.elapsed() << std::endl;

    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

    // Identify the nodes that are locked
    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

    Eigen::VectorXd b( simplexRefinableCellMesh.nodes() );

    // Initiailize the constraints to the evaluation of the Laplacian of the Franke function
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) b[i] = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

    // Compute the dual representation of the constraints
    b = M * b;

    // Set the constraints at the locked vertices to the evluation of the Franke function
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) if( lockedNodes[i] ) b[i] = Franke( NodePosition( nodeMultiIndices[i] ) );

    // Adjust the right-hand-side to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b[ iter.row() ] -= b[ iter.col() ] * iter.value();

    // Adjust the system matrix to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
            else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

    //    Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > solver;
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver;

    // Compute the Cholesky factorization
    timer.reset();
    solver.compute( S );
    switch( solver.info() )
    {
        case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
        case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
        case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
    }

    std::cout << "Performed Cholesky (LDLt) factorization: " << timer.elapsed() << std::endl;
    std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;
    timer.reset();
    Eigen::VectorXd x = solver.solve( b );

    std::cout << "Solved direct system: " << timer.elapsed() << std::endl;
    double eIn = 0 , eOut = 0;
    Eigen::VectorXd r = b - S* x;
    for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
    std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;

    double rms = 0;
    //=============Gauss Quadrature ==============
//    Eigen::VectorXd PSol = P*x;
//    Eigen::MatrixXd QuadWeigths {{0.33333333333333, 0.33333333333333, -0.56250000000000},
//                                {0.2, 0.2, 0.52083333333333},
//                                {0.2, 0.6, 0.52083333333333},
//                                {0.6, 0.2, 0.52083333333333}};
//    Eigen::MatrixXd Tri;
//
//    Tri.resize(2,3);
//    double l2_error = 0;
//    double area_sum = 0;
//    for(unsigned  int i=0 ; i< simplexRefinableCellMesh.simplexMesh().simplices();i++){
//        SimplexIndex< 2 , unsigned int > s =  simplexRefinableCellMesh.simplexMesh().simplex(i);
//        for( unsigned int n=0 ; n<3 ; n++)
//        {
//            Point<double, 2> p = fullVertexPositionFunction(s.idx[n]);
//            Tri.col(n) = Eigen::Vector2d(p[0], p[1]);
//        }
//        Eigen::Vector2d a,b,c;
//        a = Tri.col(0);
//        b = Tri.col(2);
//        c = Tri.col(1);
//        double det = 0.5 * (-a(1) * b(0) + a(0) * b(1) + a(1) * c(0) - b(1) * c(0) - a(0) * c(1) + b(0) * c(1));
//        area_sum+=det;
//        //Iterate over Quadrature points
//        for(unsigned int q = 0; q<QuadWeigths.rows(); q++){
//            // initialize Sample
//            Eigen::Vector3d bcCoordinates, globalPos;
//            Point< double , 2 > p(QuadWeigths.row(q)(0),QuadWeigths.row(q)(1)) ;
//            typename SimplexMesh< 2 >::Sample sample;
//
//            sample.sIdx = i;
//            sample.bcCoordinates[0] = 1. - p[0] - p[1];
//            sample.bcCoordinates[1] = p[0];
//            sample.bcCoordinates[2] = p[1];
//            bcCoordinates <<1.-p[0]-p[1],p[0],p[1];
//            globalPos = Tri*bcCoordinates;
//            Point< double , 2 > globalP(globalPos[0],globalPos[1]);
//            double v_exact =  Franke(globalP);
//            double v_approx = simplexRefinableCellMesh.simplexMesh().evaluate(PSol,sample);
//            //           std::cout << "v_exact: " << v_exact << " v_approx: " << v_approx << std::endl;
//            l2_error+=(v_exact-v_approx)*(v_exact-v_approx)*det*QuadWeigths.row(q)(2);
//        }
//    }
//    std::cout << "area sum : "<<area_sum << std::endl;
//    l2_error = sqrt(fabs(l2_error));
//
//    //===========================
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
    std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes() ) << std::endl;
//    std::cout << "L2 Error: " << l2_error << std::endl;
    //    return sqrt( rms / simplexRefinableCellMesh.nodes() );
    std::cout << "i was here" << std::endl;
    return sqrt( rms / simplexRefinableCellMesh.nodes());
}

double solve_2D_AQAPoly_Poisson( pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree){
        if (degree == 1)
        {
            return solve_2D_AQAPoly_Poisson< 1 >(mesh, dof);
        }
        else if (degree == 2)
        {
            return solve_2D_AQAPoly_Poisson<2>(mesh, dof);
        }
        else if (degree == 3)
        {
            return solve_2D_AQAPoly_Poisson<3>(mesh, dof);
        }
}

template<unsigned int Degree>
double AQAPoly_condition_nr(pmp::SurfaceMesh &mesh, CoarseDimension dof, std::vector<double> &condnr){
    std::vector<::Point<double, 2>> vertices;
    std::vector<std::vector<unsigned int>> polygons;
    for (pmp::Vertex v : mesh.vertices()){
        vertices.emplace_back(::Point<double, 2>(
            mesh.position(v)[0], mesh.position(v)[1]));
    }
    for (pmp::Face f : mesh.faces())
    {
        std::vector<unsigned int> face;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            face.emplace_back(v.idx());
        }
        polygons.emplace_back(face);
    }
    Meshes::PolygonMesh<unsigned int>  polyMesh = Meshes::PolygonMesh<unsigned int>(polygons);
    typedef typename SimplexMesh< 2 , Degree >::NodeMultiIndex NodeMultiIndex;

    std::set< unsigned int > boundaryVertices;
    {
        typedef MultiIndex< 2 , unsigned int , false > EdgeIndex;
        std::map< EdgeIndex , unsigned int > edgeCount;

        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ )
        {
            Meshes::Polygon< unsigned int > poly = polyMesh.polygon(i);
            for( unsigned int j=0 ; j<poly.size() ; j++ )
            {
                unsigned int idx[] = { poly[j] , poly[(j+1)%poly.size()] };
                edgeCount[ EdgeIndex(idx) ]++;
            }
        }
        for( auto const & [ edge , count ] : edgeCount ) if( count==1 ) boundaryVertices.insert( edge[0] ) , boundaryVertices.insert( edge[1] );
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , 2> ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double ,2 >();
    };
    std::function< Point< double , 2 > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , 2 > p;
        for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
        return p/(double)Degree;
    };
    int coarse_dim;
    switch (dof)
    {
        case Vertices:
            coarse_dim = 0;
            break;
        case Edges:
            coarse_dim = 1;
            break;
        case Refined_mesh:
            coarse_dim = 2;
            break;
        default:
            std::cout << "should not happen";
    }
    SimplexRefinableCellMesh< 2 , Degree > simplexRefinableCellMesh;
    simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim , true );

    typedef typename SimplexMesh<2 , Degree >::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S;
    Timer timer;

    // Get the system matrices
    timer.reset();
    M = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
    S = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
    std::cout << "Got system matrices: " << timer.elapsed() << std::endl;

    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

////    // Identify the nodes that are locked
//    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
//    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );
//
//    // Adjust the system matrix to account for the locked nodes
//    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
//            if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
//            else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;



    using OpType =  Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    // Offset the stiffness matrix so that it becomes positive definite
//    double offset = 1e-04;
//    Eigen::SparseMatrix<double> S_ = S + M * offset;
//
    OpType op(S, M);
    BOpType Bop(M);
//
//    // Construct generalized eigen solver object, seeking three generalized
//    // eigenvalues that are closest to zero. This is equivalent to specifying
//    // a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule
    int num_eval = 2;
    int converge_speed = 15 * num_eval;
    Spectra::SymGEigsShiftSolver<OpType, BOpType, Spectra::GEigsMode::ShiftInvert> geigs(op, Bop, num_eval, converge_speed, -0.1);
    geigs.init();
    geigs.compute(Spectra::SortRule::LargestMagn);

 //---------------------------------------------
    // Construct matrix operation objects using the wrapper classes
    Spectra::SparseSymMatProd<double> op2(S);
    Spectra::SparseCholesky<double>  Bop2(M);

    // Construct generalized eigen solver object, requesting the largest three generalized eigenvalues
    Spectra::SymGEigsSolver< Spectra::SparseSymMatProd<double> , Spectra::SparseCholesky<double>, Spectra::GEigsMode::Cholesky> geigs2(op2, Bop2, 1, converge_speed);

    // Initialize and compute
    geigs2.init();
    int nconv = geigs2.compute(Spectra::SortRule::LargestAlge);
    // Retrieve results
    Eigen::VectorXd evalues;
    if (geigs2.info() == Spectra::CompInfo::Successful)
    {
        evalues = geigs2.eigenvalues();
    } else if (geigs2.info() == Spectra::CompInfo::NotComputed)
        fprintf(stderr, "[ERROR] Not computed\n"), exit(0);
    else if (geigs2.info() == Spectra::CompInfo::NotConverging)
        fprintf(stderr, "[ERROR] Not converging\n"), exit(0);
    else if (geigs2.info() == Spectra::CompInfo::NumericalIssue)
        fprintf(stderr, "[ERROR] Numerical issue\n"), exit(0);
    else
        fprintf(stderr, "[ERROR] Failed\n"), exit(0);
 //---------------------------------------------

    std::cout << geigs.eigenvalues()<< std::endl;
    std::cout << geigs2.eigenvalues()<< std::endl;
    std::cout << "DoF: " << simplexRefinableCellMesh.nodes()  << std::endl;
    double l_max =  geigs2.eigenvalues()(0);
    std::cout << "Max evalue: " << l_max << std::endl;
//    geigs.compute(Spectra::SortRule::SmallestMagn);
    double l_min = geigs.eigenvalues()(0);
    std::cout << "Min evalue: "  << l_min << std::endl;
    double cond = l_max/l_min;
    std::cout << "condition nr: " << cond << std::endl;
    condnr.resize(2);
    condnr[0] = cond;
    condnr[1] = cond/simplexRefinableCellMesh.nodes() ;

    return cond;
}


double AQAPoly_condition_nr(pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree, std::vector<double> &condnr)
{
    if (degree == 1)
    {
        return AQAPoly_condition_nr<1>(mesh,dof,condnr);
    }
    else if (degree == 2)
    {
        return AQAPoly_condition_nr<2>(mesh,dof, condnr);
    }
    else if (degree == 3)
    {
        return AQAPoly_condition_nr<3>(mesh,dof,condnr);
    }
    else if (degree == 4)
    {
        return AQAPoly_condition_nr<4>(mesh,dof,condnr);
    }
}


//=================================================================
double naive_Franke(pmp::SurfaceMesh &mesh, std::string &filename){
    std::vector<::Point<double, 2>> vertices;
    std::vector<std::vector<unsigned int>> polygons;

     for (pmp::Vertex v : mesh.vertices()){
        vertices.emplace_back(::Point<double, 2>(
            mesh.position(v)[0], mesh.position(v)[1]));
    }
    for (pmp::Face f : mesh.faces())
    {
        std::vector<unsigned int> face;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            face.emplace_back(v.idx());
        }
        polygons.emplace_back(face);
    }
    Meshes::PolygonMesh<unsigned int>  polyMesh = Meshes::PolygonMesh<unsigned int>(polygons);
    typedef typename SimplexMesh< 2 , 2 >::NodeMultiIndex NodeMultiIndex;

    std::set< unsigned int > boundaryVertices;
    {
        typedef MultiIndex< 2 , unsigned int , false > EdgeIndex;
        std::map< EdgeIndex , unsigned int > edgeCount;

        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ )
        {
            Meshes::Polygon< unsigned int > poly = polyMesh.polygon(i);
            for( unsigned int j=0 ; j<poly.size() ; j++ )
            {
                unsigned int idx[] = { poly[j] , poly[(j+1)%poly.size()] };
                edgeCount[ EdgeIndex(idx) ]++;
            }
        }
        for( auto const & [ edge , count ] : edgeCount ) if( count==1 ) boundaryVertices.insert( edge[0] ) , boundaryVertices.insert( edge[1] );
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<2 ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , 2> ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double ,2 >();
    };
    std::function< Point< double , 2 > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , 2 > p;
        for( unsigned int d=0 ; d<2 ; d++ ) p += fullVertexPositionFunction( nmi[d] );
        return p/(double)2;
    };
    int coarse_dim = 1;

    double error;
    // solve system directly

        SimplexRefinableCellMesh< 2 , 2 > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< 2 >( fullVertexPositionFunction , eWeights , coarse_dim , true );

    //-----------Execute direct----------------------------------------
        typedef typename SimplexMesh<2 , 2>::NodeMultiIndex NodeMultiIndex;
        Eigen::SparseMatrix< double > M , S ,M2,S2, Sf,P, Pt;
        Timer timer;
        Eigen::loadMarket(P,filename);
//       std::ifstream infile(filename);
//       std::vector<Eigen::Triplet<double>> tripP;
//       int ro, co, row,col;
//       double  v;
//       int c = 0;
//       while (infile >>ro >> co >> v )
//       {
//
//           if(c == 0){
//               std::cout << ro <<", " << co << ", "<<v <<std::endl;
//               row = ro;
//               col = co;
//           }else{
//               tripP.emplace_back(ro,co,v);
//           }
//           c++;
//       }
//       P.resize(row,col);
//       P.setFromTriplets(tripP.begin(),tripP.end());
       Pt = P.transpose();
        // Get the system matrices
//        P = simplexRefinableCellMesh.P();
//        timer.reset();
        std::cout << filename << std::endl;
        std::cout << "P :" << P.rows() << "," << P.cols() << std::endl;
        Sf = simplexRefinableCellMesh.simplexMesh().stiffness();
        std::cout << "Sf :" << Sf.rows() << "," << Sf.cols() << std::endl;

        M =  Pt  * simplexRefinableCellMesh.simplexMesh().mass() * P ;
        S =  Pt  * simplexRefinableCellMesh.simplexMesh().stiffness() * P;
//
        M2 = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
        S2 = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();

        std::cout << "Got system matrices: " << timer.elapsed() << std::endl;

        // Get the list of node multi-indices
        std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
        for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

        // Identify the nodes that are locked
        std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
        for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

        Eigen::VectorXd b( simplexRefinableCellMesh.nodes() );

        // Initiailize the constraints to the evaluation of the Laplacian of the Franke function
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) b[i] = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

        // Compute the dual representation of the constraints
        b = M * b;

        // Set the constraints at the locked vertices to the evluation of the Franke function
        for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) if( lockedNodes[i] ) b[i] = Franke( NodePosition( nodeMultiIndices[i] ) );

        // Adjust the right-hand-side to account for the locked nodes
        for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
                if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b[ iter.row() ] -= b[ iter.col() ] * iter.value();

        // Adjust the system matrix to account for the locked nodes
        for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
                if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
                else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

//
        // Adjust the system matrix to account for the locked nodes
        for( unsigned int i=0 ; i<S2.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S2,i) ; iter ; ++iter )
                if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
                else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

//            Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver;
        Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver;

        Eigen::saveMarket(S,"System_naive.mtx");
        Eigen::saveMarket(S2,"System_misha.mtx");
        Eigen::loadMarket(S,"./System_naive.mtx");
        // Compute the Cholesky factorization
        timer.reset();
        solver.compute( S );
        switch( solver.info() )
        {
            case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
            case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
            case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
        }

        std::cout << "Performed Cholesky (LDLt) factorization: " << timer.elapsed() << std::endl;
        std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;

        timer.reset();
        Eigen::VectorXd x = solver.solve( b );

        std::cout << "Solved direct system: " << timer.elapsed() << std::endl;
        double eIn = 0 , eOut = 0;
        Eigen::VectorXd r = b - S* x;
        for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
        std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;

        double rms = 0;
        //=============Gauss Quadrature ==============
        Eigen::VectorXd PSol = P*x;
        Eigen::MatrixXd QuadWeigths {{0.33333333333333, 0.33333333333333, -0.56250000000000},
                                    {0.2, 0.2, 0.52083333333333},
                                    {0.2, 0.6, 0.52083333333333},
                                    {0.6, 0.2, 0.52083333333333}};
        Eigen::MatrixXd Tri;

        Tri.resize(2,3);
        double l2_error = 0;
        double area_sum = 0;
        for(unsigned  int i=0 ; i< simplexRefinableCellMesh.simplexMesh().simplices();i++){
            SimplexIndex< 2 , unsigned int > s =  simplexRefinableCellMesh.simplexMesh().simplex(i);
            for( unsigned int n=0 ; n<3 ; n++)
            {
                Point<double, 2> p = fullVertexPositionFunction(s.idx[n]);
                Tri.col(n) = Eigen::Vector2d(p[0], p[1]);
            }
            Eigen::Vector2d a,b,c;
            a = Tri.col(0);
            b = Tri.col(2);
            c = Tri.col(1);
            double det = 0.5 * (-a(1) * b(0) + a(0) * b(1) + a(1) * c(0) - b(1) * c(0) - a(0) * c(1) + b(0) * c(1));
            area_sum+=det;
            //Iterate over Quadrature points
            for(unsigned int q = 0; q<QuadWeigths.rows(); q++){
                // initialize Sample
                Eigen::Vector3d bcCoordinates, globalPos;
                Point< double , 2 > p(QuadWeigths.row(q)(0),QuadWeigths.row(q)(1)) ;
                typename SimplexMesh< 2 >::Sample sample;

                sample.sIdx = i;
                sample.bcCoordinates[0] = 1. - p[0] - p[1];
                sample.bcCoordinates[1] = p[0];
                sample.bcCoordinates[2] = p[1];
                bcCoordinates <<1.-p[0]-p[1],p[0],p[1];
                globalPos = Tri*bcCoordinates;
                Point< double , 2 > globalP(globalPos[0],globalPos[1]);
                double v_exact =  Franke(globalP);
                double v_approx = simplexRefinableCellMesh.simplexMesh().evaluate(PSol,sample);
//                std::cout << "v_exact: " << v_exact << " v_approx: " << v_approx << std::endl;
                l2_error+=(v_exact-v_approx)*(v_exact-v_approx)*det*QuadWeigths.row(q)(2);
            }
        }
        std::cout << "area sum : "<<area_sum << std::endl;
        l2_error = sqrt(fabs(l2_error));

        //===========================
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
        error = sqrt( rms / simplexRefinableCellMesh.nodes() ) ;
        std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes() ) << std::endl;
//        std::cout << "L2 Error: " << l2_error << std::endl;
        //    return sqrt( rms / simplexRefinableCellMesh.nodes() );

    //-----------------------------------------------------------------

    return error;
}


//=================================================================
// Quadratic reproduction

double q_function( Point< double , 2 > p )
{
    //Assume planar face
    double x = p[0] , y = p[1];
//    return 2*x*x+2*y*y+x*y+3*x+3*y+1;

    return -2*x*x+2*y*y+x*y+3*x-3*y+1;

//    return 3*x-4*y+1;

}

double test_quadratic_reproduction(pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree)
{
    if (degree == 1)
    {
        return test_quadratic_reproduction<1>(mesh,dof);
    }
    else if (degree == 2)
    {
        return test_quadratic_reproduction<2>(mesh,dof);
    }
    else if (degree == 3)
    {
        return test_quadratic_reproduction<3>(mesh,dof);
    }
    else if (degree == 4)
    {
        return test_quadratic_reproduction<4>(mesh,dof);
    }
}

template<unsigned int Degree>
double test_quadratic_reproduction(pmp::SurfaceMesh &mesh, CoarseDimension dof){
    std::vector<::Point<double, 2>> vertices;
    std::vector<std::vector<unsigned int>> polygons;
    for (pmp::Vertex v : mesh.vertices()){
        vertices.emplace_back(::Point<double, 2>(
            mesh.position(v)[0], mesh.position(v)[1]));
    }
    for (pmp::Face f : mesh.faces())
    {
        std::vector<unsigned int> face;
        for (pmp::Halfedge he : mesh.halfedges(f))
        {
            pmp::Vertex v = mesh.from_vertex(he);
            face.emplace_back(v.idx());
        }
        polygons.emplace_back(face);
    }
    Meshes::PolygonMesh<unsigned int>  polyMesh = Meshes::PolygonMesh<unsigned int>(polygons);
    typedef typename SimplexMesh< 2 , Degree >::NodeMultiIndex NodeMultiIndex;
    std::set< unsigned int > boundaryVertices;
    {
        typedef MultiIndex< 2 , unsigned int , false > EdgeIndex;
        std::map< EdgeIndex , unsigned int > edgeCount;

        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ )
        {
            Meshes::Polygon< unsigned int > poly = polyMesh.polygon(i);
            for( unsigned int j=0 ; j<poly.size() ; j++ )
            {
                unsigned int idx[] = { poly[j] , poly[(j+1)%poly.size()] };
                edgeCount[ EdgeIndex(idx) ]++;
            }
        }
        for( auto const & [ edge , count ] : edgeCount ) if( count==1 ) boundaryVertices.insert( edge[0] ) , boundaryVertices.insert( edge[1] );
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , 2> ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double ,2 >();
    };
    std::function< Point< double , 2 > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , 2 > p;
        for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
        return p/(double)Degree;
    };
    int coarse_dim;
    switch (dof)
    {
        case Vertices:
            coarse_dim = 0;
            break;
        case Edges:
            coarse_dim = 1;
            break;
        case Refined_mesh:
            coarse_dim = 2;
            break;
        default:
            std::cout << "should not happen";
    }

    SimplexRefinableCellMesh< 2 ,  Degree> simplexRefinableCellMesh;
    simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh<Degree  >( fullVertexPositionFunction , eWeights , coarse_dim , true );
    typedef typename SimplexMesh<2 , Degree>::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S , P;
    Timer timer;

    // Get the system matrices
    P = simplexRefinableCellMesh.P();
    timer.reset();

    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

    // Identify the nodes that are locked
    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

    Eigen::VectorXd x( simplexRefinableCellMesh.nodes() );

    // Initiailize the quadratic Function at the nodes
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) x[i] = q_function( NodePosition( nodeMultiIndices[i] ) );

//    // Compute the dual representation of the constraints
//    b = M * b;

    //=============Gauss Quadrature ==============
        Eigen::VectorXd PSol = P*x;
        Eigen::MatrixXd QuadWeigths {{0.33333333333333, 0.33333333333333, -0.56250000000000},
                                    {0.2, 0.2, 0.52083333333333},
                                    {0.2, 0.6, 0.52083333333333},
                                    {0.6, 0.2, 0.52083333333333}};
        Eigen::MatrixXd Tri;

        Tri.resize(2,3);
        double l2_error = 0;
        double area_sum = 0;
        for(unsigned  int i=0 ; i< simplexRefinableCellMesh.simplexMesh().simplices();i++){
            SimplexIndex< 2 , unsigned int > s =  simplexRefinableCellMesh.simplexMesh().simplex(i);
            for( unsigned int n=0 ; n<3 ; n++)
            {
                Point<double, 2> p = fullVertexPositionFunction(s.idx[n]);
                Tri.col(n) = Eigen::Vector2d(p[0], p[1]);
            }
            Eigen::Vector2d a,b,c;
            a = Tri.col(0);
            b = Tri.col(2);
            c = Tri.col(1);
            double det = 0.5 * (-a(1) * b(0) + a(0) * b(1) + a(1) * c(0) - b(1) * c(0) - a(0) * c(1) + b(0) * c(1));
            area_sum+=det;
            //Iterate over Quadrature points
            for(unsigned int q = 0; q<QuadWeigths.rows(); q++){
                // initialize Sample
                Eigen::Vector3d bcCoordinates, globalPos;
                Point< double , 2 > p(QuadWeigths.row(q)(0),QuadWeigths.row(q)(1)) ;
                typename SimplexMesh< 2 >::Sample sample;

                sample.sIdx = i;
                sample.bcCoordinates[0] = 1. - p[0] - p[1];
                sample.bcCoordinates[1] = p[0];
                sample.bcCoordinates[2] = p[1];
                bcCoordinates <<1.-p[0]-p[1],p[0],p[1];
                globalPos = Tri*bcCoordinates;
                Point< double , 2 > globalP(globalPos[0],globalPos[1]);
                double v_exact =  q_function(globalP);
                double v_approx = simplexRefinableCellMesh.simplexMesh().evaluate(PSol,sample);
                           std::cout << "v_exact: " << v_exact << " v_approx: " << v_approx << std::endl;
                l2_error+=(v_exact-v_approx)*(v_exact-v_approx)*det*QuadWeigths.row(q)(2);
            }
        }
        std::cout << "area sum : "<<area_sum << std::endl;
        l2_error = sqrt(fabs(l2_error));

    //===========================
    std::cout << "L2 Error: " << l2_error << std::endl;
    return l2_error;
}
