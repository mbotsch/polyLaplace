//=============================================================================

#include "[dGBD20]Laplace.h"
#include "[AW11]Laplace.h"
#include "diffgeo.h"

//=============================================================================

using SurfaceMesh = pmp::SurfaceMesh;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

// recommended in paper
float disney_laplace_lambda_ = 1.0;

//=============================================================================

void setup_disney_gradient_operator(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G) {

// ------------------------------------ Gradient in Paper -----------------------------------------
    G.resize(3 * mesh.n_faces(), mesh.n_vertices());

    std::vector<Eigen::Triplet<double>> triplets;
    int cnt = 0;
//
    for (auto f: mesh.faces()) {
        const unsigned int n = mesh.valence(f);
        Eigen::VectorXi F(n);

        int i = 0;
        for (auto v : mesh.vertices(f)) {
            F(i) = v.idx();
            ++i;
        }


        Eigen::MatrixXd E, B, Avg, EtAvg;
        Avg = Eigen::MatrixXd::Zero(n, n);
        for (int i = 0; i < n; i++) {
            if (i == n - 1) { Avg(i, 0) = 0.5; }
            else { Avg(i, i + 1) = 0.5; }
            Avg(i, i) = 0.5;
        }
        setup_E_and_B_perFace(mesh, f, E, B);
        EtAvg.resize(3, n);

        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
        Eigen::Vector3d nf = af;
        EtAvg = E.transpose() * Avg;
        nf.normalize();


        Eigen::MatrixXd d(n, n);
        Eigen::MatrixXd U_f, G_f;
        setup_disney_sharp_operator_per_face(mesh, f, U_f);
        d.setZero();
        G_f.resize(3,n);
        for (unsigned int i = 0; i < n; ++i) {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }

        double area = -1.0 / af.norm();
        int j = 0;
        for (pmp::Vertex v : mesh.vertices(f)) {
            Eigen::Vector3d x = EtAvg.col(j);
            Eigen::Vector3d gradient = nf.cross(x);
            G_f.col(j) = area*gradient;
            triplets.emplace_back(3 * f.idx(), v.idx(), area * gradient(0));
            triplets.emplace_back(3 * f.idx() + 1, v.idx(), area * gradient(1));
            triplets.emplace_back(3 * f.idx() + 2, v.idx(), area * gradient(2));
            j++;
        }
//
//        std::cout << "Gradient G_f: \n" << G_f << std::endl;
//        std::cout << "Gradient U_f*D_f: \n" << U_f*d << std::endl;
//        std::cout << "--------------------------------" << std::endl;

//        std::cout << "Norm Gradients: \n" << (G_f-U_f*d).norm() << std::endl;
//        std::cout << "Gradient U_f*D_f: \n" << U_f*d << std::endl;
//        std::cout << "--------------------------------" << std::endl;
    }

    G.setFromTriplets(triplets.begin(), triplets.end());


// ------------------------------------ Gradient as Difference Matrix -----------------------------------------

//
//    std::vector<Eigen::Triplet<double>> triplets;
//    int cnt = 0;
//
//    for (auto f : mesh.faces()) {
//        const unsigned int n = mesh.valence(f);
//        Eigen::VectorXi F(n);
//
//        int i = 0;
//        for (auto v : mesh.vertices(f)) {
//            F(i) = v.idx();
//            ++i;
//        }
//
//        for (unsigned int i = 0; i < n; ++i) {
//            triplets.emplace_back(cnt, F(i), -1);
//            triplets.emplace_back(cnt, F((i + 1) % n), 1);
//            ++cnt;
//        }
//    }
//
//    G.resize(cnt, mesh.n_vertices());
//    G.setFromTriplets(triplets.begin(), triplets.end());
//

// ------------------------------------ Gradient from Vertex gradients -----------------------------------------


//    std::vector<Eigen::Triplet<double>> triplets;
//
//    int cnt = 0;
//    for (auto f : mesh.faces()) {
//        const unsigned int n = mesh.valence(f);
//        Eigen::VectorXi F(n);
//
//        int i = 0;
//        for (Vertex v : mesh.vertices(f)) {
//            F(i) = v.idx();
//            ++i;
//        }
//        Eigen::MatrixXd E, B, Gf;
//
//        setup_E_and_B_perFace(mesh, f, E, B);
//
//
//        // compute vector area
//        const Eigen::Matrix3d Af = E.transpose() * B;
//        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
//        Eigen::Vector3d nf = af;
//        nf.normalize();
////        area_triplets.emplace_back(3*f.idx(), 3*f.idx(), af.norm());
////        area_triplets.emplace_back(3*f.idx()+1, 3*f.idx()+1, af.norm());
////        area_triplets.emplace_back(3*f.idx()+2, 3*f.idx()+2, af.norm());
//
//
//        // compute d
//        Eigen::MatrixXd d(n, n);
//        Eigen::MatrixXd U_f, G_f;
//        setup_disney_sharp_operator_per_face(mesh, f, U_f);
//        d.setZero();
//        G_f.resize(3,n);
//        for (unsigned int i = 0; i < n; ++i) {
//            d(i, i) = -1;
//            d(i, (i + 1) % n) = 1;
//        }
//
//        double area = 1.0 / (2.0*af.norm());
//
//        for (unsigned int i = 0; i < n; ++i) {
//            Point x_i = mesh.position(Vertex(F(i)));
//            Point x_ii = mesh.position(Vertex(F((i + 1) % n)));
//            x_i -= x_ii;
//            Eigen::Vector3d diff;
//            diff << x_i[0],x_i[1],x_i[2];
//            Eigen::Vector3d g_vf = area*nf.cross(diff);
//            G_f.col(F(i)) = g_vf;
//            triplets.emplace_back(3*f.idx(), F(i), g_vf(0));
//            triplets.emplace_back(3*f.idx()+1, F(i), g_vf(1));
//            triplets.emplace_back(3*f.idx()+2, F(i), g_vf(2));
//
//        }
//
//        std::cout << "Gradient G_f: \n" << G_f << std::endl;
//        std::cout << "Gradient U_f*D_f: \n" << U_f*d << std::endl;
//
//    }
//
//    G.resize(3*mesh.n_faces(), mesh.n_vertices());
//    G.setFromTriplets(triplets.begin(), triplets.end());
}

void setup_disney_divergence_operator(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &D) {

//---------------------------- Divergence as Gradient Transpose * area ------------------------------

    Eigen::SparseMatrix<double> G, A;
    setup_disney_gradient_operator(mesh,G);

    setup_disney_face_area_matrix(mesh,A);
    D.resize(mesh.n_vertices(),3*mesh.n_faces());
    D= -G.transpose()*A;
//    std::cout << "Gradient: " << G << std::endl;
//    std::cout << " Area: " << A << std::endl;
//    std::cout << " Div:" << -G.transpose()*A << std::endl;
//    std::cout << "Laplace: " << D*G << std::endl;


//---------------------------- Divergence in Paper------------------------------


//    std::vector<Eigen::Triplet<double>> triplets;
//    int colCnt = 0;
//
//    for (auto f : mesh.faces()) {
//        const unsigned int n = mesh.valence(f);
//        Eigen::VectorXi F(n);
//        Eigen::MatrixXd E, B;
//        setup_E_and_B_perFace(mesh, f, E, B);
//        int i = 0;
//        for (auto v : mesh.vertices(f)) {
//            F(i) = v.idx();
//            ++i;
//        }
//        // compute vector area
//        const Eigen::Matrix3d Af = E.transpose() * B;
//        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
//
//        const double area = af.norm();
//        af /= area;
//
//        // compute d
//        Eigen::MatrixXd d(n, n);
//        Eigen::MatrixXd M, Lf;
//        setup_disney_inner_product_per_face(mesh,f,M);
//        d.setZero();
//
//        for (unsigned int i = 0; i < n; ++i) {
//            d(i, i) = -1;
//            d(i, (i + 1) % n) = 1;
//        }
//
//        Lf = d.transpose() * M;
//
//        // add local laplacian to global matrix entries
//        for (unsigned int k = 0; k < n; ++k) {
//            for (unsigned int l = 0; l < n; ++l) {
//                triplets.push_back(
//                        Eigen::Triplet<double>(F(k), colCnt + l, - Lf(k, l)));
//            }
//        }
//
//        colCnt += n;
//    }
//
//    D.resize(mesh.n_vertices(), colCnt);
//    D.setFromTriplets(triplets.begin(), triplets.end());
}

void setup_disney_sharp_operator_per_face(pmp::SurfaceMesh &mesh, pmp::Face f, Eigen::MatrixXd &U) {

    const unsigned int n = mesh.valence(f);
    U.resize(3, n);
    Eigen::MatrixXd E, B, Bt;
    setup_E_and_B_perFace(mesh, f, E, B);
    Bt.resize(3, n);
    Bt = B.transpose();
    pmp::Point c = centroid(mesh, f);
    Eigen::Vector3d cf;
    cf << c[0], c[1], c[2];
    // compute vector area
    const Eigen::Matrix3d Af = E.transpose() * B;
    Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    Eigen::Vector3d nf = af;
    nf.normalize();

    double area = 1.0 / af.norm();
    for (int i = 0; i < n; i++) {
        Eigen::Vector3d x = Bt.col(i) - cf;
        Eigen::Vector3d cross_prod = nf.cross(x);
        U.col(i) = area * cross_prod;
    }
}

void setup_disney_discrete_flat_operator_per_face(pmp::SurfaceMesh &mesh, pmp::Face f, Eigen::MatrixXd &V) {

    const unsigned int n = mesh.valence(f);
    V.resize(n, 3);
    Eigen::MatrixXd E, B, NNT;
    setup_E_and_B_perFace(mesh, f, E, B);

    // compute vector area
    const Eigen::Matrix3d Af = E.transpose() * B;
    Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    Eigen::Vector3d nf = af;
    nf.normalize();
    NNT = Eigen::MatrixXd::Identity(3, 3);
    NNT -= nf * nf.transpose();
    V = E * NNT;
}

void setup_disney_projection_operator_per_face(pmp::SurfaceMesh &mesh, pmp::Face f, Eigen::MatrixXd &P){

    const unsigned int n = mesh.valence(f);
    P.resize(n,n);
    P.setIdentity();
    Eigen::MatrixXd U, V;
    setup_disney_discrete_flat_operator_per_face(mesh,f, V);
    setup_disney_sharp_operator_per_face(mesh,f, U);

    P -= V * U;
}

void setup_disney_inner_product_per_face(pmp::SurfaceMesh &mesh,pmp::Face f, Eigen::MatrixXd &M){
    const unsigned int n = mesh.valence(f);

    double lambda = disney_laplace_lambda_;
    M.resize(n, n);
    Eigen::MatrixXd E, B, U,P;
    setup_E_and_B_perFace(mesh, f, E, B);
    setup_disney_sharp_operator_per_face(mesh,f,U);
    setup_disney_projection_operator_per_face(mesh,f,P);
    const Eigen::Matrix3d Af = E.transpose() * B;
    Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    double area = af.norm();

    M = area*U.transpose()*U+lambda*P.transpose()*P;

}

void setup_disney_laplace_operator(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &L){

    L.resize(mesh.n_vertices(), mesh.n_vertices());

    std::vector<Eigen::Triplet<double>> triplets;
    int cnt = 0;

    for (auto f: mesh.faces()) {

        const unsigned int n = mesh.valence(f);

        Eigen::VectorXi F(n);
        int i = 0;
        for (auto v : mesh.vertices(f)) {
            F(i) = v.idx();
            ++i;
        }
        Eigen::MatrixXd M, Lf;
        setup_disney_inner_product_per_face(mesh,f,M);
        // compute d
        Eigen::MatrixXd d(n, n);
        d.setZero();
        for (unsigned int i = 0; i < n; ++i) {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }
        Lf = d.transpose()*M*d;

        // add local laplacian to global matrix entries
        for (unsigned int k = 0; k < n; ++k) {
            for (unsigned int l = 0; l < n; ++l) {
                triplets.push_back(
                        Eigen::Triplet<double>(F(k), F(l), - Lf(k, l)));
            }
        }
    }
    L.setFromTriplets(triplets.begin(),triplets.end());
    std::cout << "Disney Hyperparamter: " << disney_laplace_lambda_ << std::endl;

//    Eigen::SparseMatrix<double> Grad, Div;
//    setup_disney_gradient_operator(mesh,Grad);
//    setup_disney_divergence_operator(mesh,Div);
//    std::cout << "Laplace inner Prod: " << L << std::endl;
//    std::cout << "Norm: " << (L-(Div*Grad)).norm() << std::endl;

//    std::cout << "Div*Grad: \n" << Div*Grad << std::endl;



}
void setup_disney_mass_matrix(pmp::SurfaceMesh &mesh,
                            Eigen::SparseMatrix<double> &M) {
    M.resize(mesh.n_vertices(), mesh.n_vertices());

    std::vector<Triplet> tripletsM;
    double sum = 0.0;
    Eigen::MatrixXd E, B;
    for (auto v : mesh.vertices()) {
        for (auto f : mesh.faces(v)) {
            const unsigned int n = mesh.valence(f);
//            Eigen::MatrixXd X(n, 3);
//            int i = 0;
//            for (auto v : mesh.vertices(f)) {
//                const auto p = mesh.position(v);
//                X(i, 0) = p[0];
//                X(i, 1) = p[1];
//                X(i, 2) = p[2];
//                ++i;
//            }
//            Eigen::Vector3d af = Eigen::MatrixXd::Zero(3,1);
//            for (unsigned int i = 0; i < n; ++i) {
//                Eigen::Vector3d xi, xii;
//                xi = X.row(i);
//                xii = X.row((i + 1) % n);
//                af += xii.cross(xi);
//            }
//            af *= 0.5;

            setup_E_and_B_perFace(mesh, f, E, B);
//             compute vector area

            const Eigen::Matrix3d Af = E.transpose() * B;
            Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
            const double area = af.norm();
            sum += area / (double)(n);
        }
        tripletsM.emplace_back(v.idx(), v.idx(), sum);
        sum = 0.0;
    }
    // build sparse matrix from triplets
    M.setFromTriplets(tripletsM.begin(), tripletsM.end());
}

void setup_disney_face_area_matrix(pmp::SurfaceMesh &mesh,
                                   Eigen::SparseMatrix<double> &A){

    std::vector<Eigen::Triplet<double>> area_triplets;
    A.resize(3*mesh.n_faces(),3*mesh.n_faces());
    for (auto f : mesh.faces()) {

        Eigen::MatrixXd E, B, Gf;
        setup_E_and_B_perFace(mesh, f, E, B);


        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
        area_triplets.emplace_back(3*f.idx(), 3*f.idx(), af.norm());
        area_triplets.emplace_back(3*f.idx()+1, 3*f.idx()+1, af.norm());
        area_triplets.emplace_back(3*f.idx()+2, 3*f.idx()+2, af.norm());

    }
    A.setFromTriplets(area_triplets.begin(), area_triplets.end());

}
void normalize_disney_gradients(pmp::SurfaceMesh &mesh, Eigen::VectorXd &g,
                              const Eigen::VectorXd &h) {

    double lambda = disney_laplace_lambda_;

    assert(h.rows() == mesh.n_vertices());
    //  assert(mesh.n_halfedges() == g.rows());

    int hedgeCnt = 0;

    for (auto f : mesh.faces()) {
        const unsigned int n = mesh.valence(f);
        Eigen::VectorXi F(n);
        int i = 0;
        for (auto v : mesh.vertices(f)) {
            F(i) = v.idx();
            ++i;
        }

        Eigen::MatrixXd E, B;
        setup_E_and_B_perFace(mesh, f, E, B);

        //compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));

        const double area = af.norm();
        af /= area;

        Eigen::MatrixXd Mf, Lf, Vf;
        setup_disney_inner_product_per_face(mesh,f,Mf);
        setup_disney_discrete_flat_operator_per_face(mesh, f, Vf);
        // compute d
        Eigen::MatrixXd d(n, n);
        d.setZero();
        for (unsigned int i = 0; i < n; ++i) {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }
        Lf =  d.transpose()*Mf*d;


        // get heat values for face f

        Eigen::VectorXd hi(n);
        for (unsigned int i = 0; i < n; ++i)
            hi(i) = h(F(i));
        const double factor = sqrt(hi.transpose().dot(Lf * hi) / area);
//        for (unsigned int i = 0; i < n; ++i)
        for (unsigned int i = 0; i < n; ++i)
            g(hedgeCnt++) /= factor;
    }

    assert(g.rows() == hedgeCnt);
}