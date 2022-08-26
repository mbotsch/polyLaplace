//=============================================================================

#include <Misha/Geometry.h>
#include <Misha/SimplexRefinableBasis.h>
#include <Misha/Miscellany.h>
#include <Misha/Spectrum.h>
#include <Misha/Meshes.h>
#include <fstream>
#include <Eigen/src/SparseCore/SparseMatrix.h>
#include <igl/slice.h>
#include <Misha/MGSolver.h>
#include <Eigen/CholmodSupport>

#include "AQAPoly_Laplacian_3D.h"
#include "Eigenmodes.h"
#include "unsupported/Eigen/SparseExtra"

//=============================================================================

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

//=============================================================================

void ReadPolyhedra( const std::string &fileName , std::vector< std::vector< std::pair< unsigned int , bool > > > &polyhedra , std::vector< Meshes::Polygon< unsigned int > > &polygons , std::vector< Point< double , 3 > > &vertices )
{
    typedef std::pair< unsigned int , unsigned int > Edge;
    std::vector< Edge > edges;
    std::ifstream file( fileName );
    if( !file.is_open() ) ERROR_OUT( "Could not open file for reading: " , fileName );

    auto CheckString = [&]( const std::string& s )
    {
        std::string buf;
        file >> buf;
        if( s.compare(buf) )
        {
            std::cout << "reading ovm: string " << s << " not found. Returning." << std::endl;
            return false;
        }
        return true;
    };

    if( !CheckString( "OVM" ) ) ERROR_OUT( "Failed to parse 'OVM'" );
    if( !CheckString( "ASCII" ) ) ERROR_OUT( "Failed to parse 'ASCII'" );

    {
        if( !CheckString( "Vertices" ) ) ERROR_OUT( "Failed to parse 'Vertices'" );

        int nv;
        file >> nv;
        vertices.reserve( nv );

        double x, y, z;
        for( int i=0 ; i<nv ; i++ )
        {
            file >> x >> y >> z;
            vertices.push_back( Point< double , 3 >( x , y , z ) );
        }
    }

    {

        if( !CheckString( "Edges" ) ) ERROR_OUT( "Could not find 'Edges'" );
        int ne;
        file >> ne;
        edges.reserve( ne );

        Edge e;
        for( int i=0 ; i<ne ; i++ )
        {
            file >> e.first >> e.second;
            edges.push_back(e);
        }
    }

    {
        if( !CheckString( "Faces" ) ) ERROR_OUT( "Could not find 'Faces'" );
        int nf;
        file >> nf;
        polygons.resize( nf );

        for( int i=0 ; i<nf ; i++ )
        {
            int nei;
            file >> nei;
            polygons[i].reserve(nei);
            unsigned int id;

            for( int j=0 ; j<nei ; j++ )
            {
                file >> id;
                if( id&1 ) polygons[i].push_back( edges[ id>>1 ].second );
                else       polygons[i].push_back( edges[ id>>1 ].first );
            }
        }
    }

    {
        if( !CheckString( "Polyhedra" ) ) ERROR_OUT( "Could not find 'Polyhedra'" );
        int nc;
        file >> nc;
        polyhedra.resize( nc );

        for( int i=0 ; i<nc ; i++ )
        {
            int nfi;
            file >> nfi;
            polyhedra[i].reserve( nfi );

            unsigned int id;
            for( int j=0 ; j<nfi ; j++ )
            {
                file >> id;
                polyhedra[i].push_back( std::pair< unsigned int , bool >( id>>1 , !(id&1) ) );
            }
        }
    }

    file.close();
}
//---------------------------------------------------------------------------------------------------------
double Franke( Point< double , 3 > p )
{
    double x=p[0] , y=p[1] , z=p[2];

    double cx2 = (9. * x - 2.) * (9. * x - 2.);
    double cy2 = (9. * y - 2.) * (9. * y - 2.);
    double cz2 = (9. * z - 2.) * (9. * z - 2.);

    double cx1 = (9. * x + 1.) * (9. * x + 1.);
    double cx7 = (9. * x - 7.) * (9. * x - 7.);

    double cy3 = (9. * y - 3.) * (9. * y - 3.);
    double cx4 = (9. * x - 4.) * (9. * x - 4.);

    double cy7 = (9. * y - 7.) * (9. * y - 7.);
    double cz5 = (9. * z - 5.) * (9. * z - 5.);

    return (3. / 4.) *
               exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2 - (1. / 4.) * cz2) +
           (3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10. -
                           (9. / 10.) * z - 1. / 10.) +
           (1. / 2.) *
               exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3 - (1. / 4.) * cz5) -
           (1. / 5.) * exp(-cx4 - cy7 - cz5);
}

//---------------------------------------------------------------------------------------------------------

double FrankeLaplacian( Point< double , 3 > p )
{
    double x=p[0] , y=p[1] , z=p[2];

    return (243. * (-2299. + 1800. * x * (2. + 9. * x))) /
               (480200. *
                exp((9. * (12. + 10. * x * (2. + 9. * x) + 49. * y + 49. * z)) /
                    490.)) -
           (486. *
            exp(-pow(4. - 9. * x, 2) - pow(7. - 9. * y, 2) -
                pow(5. - 9. * z, 2)) *
            (59. + 6. * x * (-8. + 9. * x) + 6. * y * (-14. + 9. * y) +
             6. * z * (-10 + 9 * z))) /
               5. +
           (81. *
            exp((-pow(7. - 9. * x, 2) - 9 * pow(1. - 3. * y, 2) -
                 pow(5. - 9. * z, 2)) /
                4.) *
            (77. + 9. * x * (-14. + 9. * x) + 27. * y * (-2. + 3. * y) +
             9 * z * (-10. + 9. * z))) /
               8. +
           (729. * (2. + 3. * x * (-4. + 9. * x) + 3. * y * (-4. + 9. * y) +
                    3. * z * (-4. + 9. * z))) /
               (16. *
                exp((3. * (4. + 3. * x * (-4. + 9. * x) +
                           3. * y * (-4. + 9. * y) + 3. * z * (-4. + 9. * z))) /
                    4.0));
}

template<unsigned int Degree>
double solve_3D_AQAPoly_Poisson(std::string &meshname, CoarseDimension dof,int function){

    const unsigned int Dim = 3;
//    const unsigned int Degree = 2;
    std::vector<std::vector<std::pair<unsigned int, bool>>> polyhedra;
    std::vector<Meshes::Polygon<unsigned int>> polygons;
    std::vector<Point<double, 3>> vertices;
    std::vector<Eigen::Triplet<double>> trip;
    ReadPolyhedra(meshname,polyhedra,polygons,vertices);
    Meshes::PolyhedronMesh<unsigned int> polyMesh = Meshes::PolyhedronMesh<unsigned int>(polyhedra, polygons);

    typedef typename SimplexMesh<Dim, Degree>::NodeMultiIndex NodeMultiIndex;

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-08;
    SimplexRefinableElements< Dim >::EnergyType eType( eWeights );

    eType.isIntegrationFaceFunctor = [&]( SimplexRefinableElements<Dim>::FaceMultiIndex fmi )
    {
        // Integrate over any face connected to the cell's virtual center
        return fmi[0]>=( polyMesh.vertices() + polyMesh.polygons() );
    };

    std::function< Point< double , Dim > ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double , Dim >();
    };

   std::function< Point< double , Dim > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , false );

    //std::function< Point< double , Dim > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;

    switch(dof)
    {
        case Vertices:
            simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,0 );
            break;
        case Edges:
            simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,1 );
            break;
        case Virtual_Edges:
            simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,2 );
            break;
        case Refined_mesh:
            simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,3 );
            break;
        default:
            std::cout << "Coarsening Type not defined " << std::endl;
    }
    Eigen::SparseMatrix<double > M,S, S2,P;
    P = simplexRefinableCellMesh.P();
    S = simplexRefinableCellMesh.Pt()*simplexRefinableCellMesh.simplexMesh().stiffness()* simplexRefinableCellMesh.P();
    M = simplexRefinableCellMesh.Pt()*simplexRefinableCellMesh.simplexMesh().mass()*simplexRefinableCellMesh.P();

//    Eigen::saveMarket(P,"P.mtx");
//    std::cout << "-----------------------" << std::endl;
//    std::cout << "AQAPoly stiffness: \n" << S << std::endl;
//    std::cout << "-----------------------" << std::endl;
//    lump_matrix(M);
//    std::cout << "AQAPoly mass: \n" << M << std::endl;
    std::set< unsigned int > boundaryVertices;
    {
        std::vector< unsigned int > faceCount( polyMesh.polygons() , 0 );
        for( unsigned int i=0 ; i<polyMesh.polyhedra() ; i++ ) for( unsigned int j=0 ; j<polyMesh.polyhedronFaces(i) ; j++ )
                faceCount[ polyMesh.polyhedronFace(i,j).first ]++;
        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ ) if( faceCount[i]==1 )
            {
                Meshes::Polygon< unsigned int > polygon = polyMesh.polygon( i );
                for( unsigned int j=0 ; j<polygon.size() ; j++ ) boundaryVertices.insert( polygon[j] );
                boundaryVertices.insert( (unsigned int)vertices.size() + i );
            }
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , Dim > p;
        for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
        return p/(double)Degree;
    };


    //-------------------------------solve poisson---------------------------------------------

    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

    // Identify the nodes that are locked
    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

    Eigen::MatrixXd b( simplexRefinableCellMesh.nodes(),1);

    if(function == 0){
        std::cout << "setup franke right-hand side" <<std::endl;
        // Initiailize the constraints to the evaluation of the Laplacian of the Franke function
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) b(i,0) = - FrankeLaplacian( NodePosition( nodeMultiIndices[i] ) );

    }else{
        std::cout << "setup  right-hand side for linear precision" <<std::endl;
        b.resize( simplexRefinableCellMesh.nodes(),3 );
        b.setZero();
//        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ){
//            Point<double,3> pos =  NodePosition( nodeMultiIndices[i]);
//            for(int j = 0; j<3 ; j++){
//                b(i,j) =  pos[j];
//            }
//        }
    }

    // Compute the dual representation of the constraints
//    std::cout << M << std::endl;
    double count = 0.0;
    for (int k = 0; k < M.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it) {
            if (it.row() == it.col()) {
                double val = it.value();
            }
            count += it.value();
        }
    }
    std::cout << "Volume mass matrix: " << count << std::endl;
    if(function == 0){
        b = M * b;
    }

    // Set the constraints at the locked vertices to the evluation of the Franke function
    if(function == 0){
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) if( lockedNodes[i] ) b(i,0) = Franke( NodePosition( nodeMultiIndices[i] ) );

    // Adjust the right-hand-side to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) b(iter.row(),0) -= b( iter.col(),0) * iter.value();
    }else{
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ){
            if( lockedNodes[i] ){
                for(int j = 0; j<3;j++){
                    b(i,j) = NodePosition( nodeMultiIndices[i])[j] ;
                }
            }
        }

        // Adjust the right-hand-side to account for the locked nodes
        for( unsigned int i=0 ; i<S.outerSize() ; i++ ){
            for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            {
                if( !lockedNodes[ iter.row() ] && lockedNodes[ iter.col() ] ) {
                    for(int j = 0; j<3; j++){
                        b(iter.row(),j) -= b(iter.col(),j) * iter.value();
                    }
                }
            }
        }
    }
    // Adjust the system matrix to account for the locked nodes
    for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            if( lockedNodes[ iter.row() ] ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
            else if( lockedNodes[ iter.col() ] ) iter.valueRef() = 0;

    Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver;
     solver.compute( S );
     switch( solver.info() )
    {
        case Eigen::NumericalIssue: ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- numerical issue" );
        case Eigen::NoConvergence:  ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- no convergence" );
        case Eigen::InvalidInput:   ERROR_OUT( "Eigen::SimplicialLDLT failed to factorize matrix -- invalid input" );
    }
    double rms = 0;
    if(function == 0){
        Eigen::VectorXd x = solver.solve( b );
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ )
        {
            rms += pow(x(i, 0) - Franke(NodePosition(nodeMultiIndices[i])), 2.);
        }
    }else{
        Eigen::MatrixXd x = solver.solve( b );
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ){
            Eigen::Vector3d pos(NodePosition( nodeMultiIndices[i] )[0],NodePosition( nodeMultiIndices[i] )[1],NodePosition( nodeMultiIndices[i] )[2]);
            Eigen::Vector3d X(x(i,0),x(i,1),x(i,2));
                if((X-pos).norm()>0.0001){
                std::cout << "x: " << x.row(i) << " Orig. Pos : " <<  NodePosition( nodeMultiIndices[i] )<< " Difference: " << X-pos <<std::endl;
            }
           for(int j = 0; j < 3 ; j++){
               rms +=pow( x(i,j) - NodePosition( nodeMultiIndices[i] )[j] , 2. );
           }
        }
    }
//    std::cout << "Size x :" <<  x.size() <<std::endl;
//    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ){
////        std::cout << "x: " << x[i] << " Franke : " <<  Franke( NodePosition( nodeMultiIndices[i] ) )<<std::endl;
//        if(function == 0){
//            rms += pow( x(i,0) - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
//        }else{
//             std::cout << "x: " << x.row(i) << " Orig. Pos : " <<  NodePosition( nodeMultiIndices[i] )<<std::endl;
//            for(int j = 0; j < 3 ; j++){
//                rms +=pow( x(i,j) - NodePosition( nodeMultiIndices[i] )[j] , 2. );
//            }
//        }
//    }
        std::cout << "DoFs: " << simplexRefinableCellMesh.nodes() << std::endl;
        std::cout << "Matrix entries: " << S.nonZeros() << std::endl;

        if(function == 0){
            std::cout << "RMS: " << sqrt( rms / (simplexRefinableCellMesh.nodes())) << std::endl;
            return sqrt( rms / simplexRefinableCellMesh.nodes() );
        }else{
            std::cout << "RMS: " << sqrt( rms / (simplexRefinableCellMesh.nodes())*3.0) << std::endl;
            return sqrt( rms / simplexRefinableCellMesh.nodes()*3.0 );
        }

}
//---------------------------------------------------------------------------------------------------------

double solve_3D_AQAPoly_Poisson(std::string &meshname, CoarseDimension dof, int degree, int function){
    if (degree == 1)
    {
        return solve_3D_AQAPoly_Poisson< 1 >(meshname, dof, function);
    }
    else if (degree == 2)
    {
        return solve_3D_AQAPoly_Poisson<2>(meshname, dof,function);
    }
    else if (degree == 3)
    {
        return solve_3D_AQAPoly_Poisson<3>(meshname, dof,function);
    }
}

//---------------------------------------------------------------------------------------------------------

template<unsigned int Degree>
double solve_3D_AQAPoly_EigenModes(std::string &meshname, CoarseDimension dof,std::string &tesselation){

        const unsigned int Dim = 3;
        //    const unsigned int Degree = 2;
        std::vector<std::vector<std::pair<unsigned int, bool>>> polyhedra;
        std::vector<Meshes::Polygon<unsigned int>> polygons;
        std::vector<Point<double, 3>> vertices;
        std::vector<Eigen::Triplet<double>> trip;
        ReadPolyhedra(meshname,polyhedra,polygons,vertices);
        Meshes::PolyhedronMesh<unsigned int> polyMesh = Meshes::PolyhedronMesh<unsigned int>(polyhedra, polygons);
        typedef typename SimplexMesh<Dim, Degree>::NodeMultiIndex NodeMultiIndex;

        SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
        eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;
        SimplexRefinableElements< Dim >::EnergyType eType( eWeights );

        eType.isIntegrationFaceFunctor = [&]( SimplexRefinableElements<Dim>::FaceMultiIndex fmi )
        {
            // Integrate over any face connected to the cell's virtual center
            return fmi[0]>=( polyMesh.vertices() + polyMesh.polygons() );
        };

        std::function< Point< double , Dim > ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
        {
            if( idx<vertices.size() ) return vertices[idx];
            ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
            return Point< double , Dim >();
        };

        std::function< Point< double , Dim > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

        SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;

        switch(dof)
        {
            case Vertices:
                simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,0 );
                break;
            case Edges:
                simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,1 );
                break;
            case Virtual_Edges:
                simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights ,2 );
                break;
            default:
                std::cout << "Coarsening Type not defined " << std::endl;
        }
        Eigen::SparseMatrix<double > M,S;

//        S = -simplexRefinableCellMesh.Pt()*simplexRefinableCellMesh.simplexMesh().stiffness()* simplexRefinableCellMesh.P();
        M = simplexRefinableCellMesh.Pt()*simplexRefinableCellMesh.simplexMesh().mass()*simplexRefinableCellMesh.P();

        Eigen::loadMarket(S,"../../S.mtx");
//        Eigen::loadMarket(M,"../../M.mtx");
        S*=-1.0;

        std::set< unsigned int > boundaryVertices;
        {
            std::vector< unsigned int > faceCount( polyMesh.polygons() , 0 );
            for( unsigned int i=0 ; i<polyMesh.polyhedra() ; i++ ) for( unsigned int j=0 ; j<polyMesh.polyhedronFaces(i) ; j++ )
                    faceCount[ polyMesh.polyhedronFace(i,j).first ]++;
            for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ ) if( faceCount[i]==1 )
                {
                    Meshes::Polygon< unsigned int > polygon = polyMesh.polygon( i );
                    for( unsigned int j=0 ; j<polygon.size() ; j++ ) boundaryVertices.insert( polygon[j] );
                    boundaryVertices.insert( (unsigned int)vertices.size() + i );
                }
        }

        auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
        {
            for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
            return true;
        };
        auto NodePosition = [&]( NodeMultiIndex nmi )
        {
            Point< double , Dim > p;
            for( unsigned int d=0 ; d<Degree ; d++ ) p += fullVertexPositionFunction( nmi[d] );
            return p/(double)Degree;
        };

        // Get the list of node multi-indices
        std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes() );
        for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ) nodeMultiIndices[idx] = nmi;

        // Identify the nodes that are locked
        std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes() );
        for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i]);

        //-------------------------------------------------------------------
        std::string filename;
        switch (dof)
        {
            case Vertices:
                filename = "eigenmodes_[BHBK20]";
                break;
            case Edges:
                filename = "eigenmodes_AQAPoly";
                break;
            case Virtual_Edges:
                filename = "eigenmodes_AQAPoly-l2";
                break;
            case Refined_mesh:
                filename = "eigenmodes_AQAPoly-refined";
                break;
            default:
                std::cout << "should not happen";
        }

        filename +="_"+tesselation+".csv";
        std::cout << filename << std::endl;
        std::ofstream ev_file(filename);
        ev_file << "computed,analytic,offset" << std::endl;

        std::vector<int> indices;

        for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap() ){
            if(!IsBoundaryNode( nmi)){
                indices.push_back(idx);
            }
        }
        std::cout << "inner indices: " << indices.size() << std::endl;

        Eigen::VectorXi in(indices.size());
        //Rewrite indices to Eigen::Vector
        for (int i = 0; i < indices.size(); ++i)
        {
            in(i) = indices[i];
        }

        Eigen::SparseMatrix<double> S_in_in, M_in_in;

        //slice matrices so that only rows and cols for inner vertices remain

        igl::slice(S, in, in, S_in_in);
        igl::slice(M, in, in, M_in_in);

        int num_eval = 34;
        int converge_speed = 5 * num_eval;
        // S and M are sparse
        using OpType =  Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
        using BOpType = Spectra::SparseSymMatProd<double>; OpType op(S_in_in, M_in_in);
        BOpType Bop(M_in_in);

        // Construct generalized eigen solver object, seeking three generalized
        // eigenvalues that are closest to zero. This is equivalent to specifying
        // a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule
        Spectra::SymGEigsShiftSolver<OpType, BOpType,
                                     Spectra::GEigsMode::ShiftInvert>
            geigs(op, Bop, num_eval, converge_speed, 1e-06);
        geigs.init();
        geigs.compute(Spectra::SortRule::LargestMagn);

        Eigen::VectorXd evectors, analytic,evalues;
        // Retrieve results
        if (geigs.info() == Spectra::CompInfo::Successful)
        {
            evalues = geigs.eigenvalues();
        }
        analytic_eigenvalues_unitBall(analytic, num_eval);
        double error = 0.0;
        for (int i = 0; i < evalues.size(); i++)
        {
            ev_file << -evalues(i) << "," << analytic(i) << ","
                    << -evalues(i) - analytic(i) << std::endl;
            std::cout << "Computed evalue: " << -evalues(i)
                      << " analytical Bessel: " << analytic(i) << std::endl;

            error += pow(-evalues(i) - analytic(i), 2);
        }
        ev_file.close();
        error = sqrt(error / (double)evalues.size());
        std::cout << "Root mean squared error: " << error << std::endl;
        return error;
}

//---------------------------------------------------------------------------------------------------------

double solve_3D_AQAPoly_EigenModes( std::string &meshname, CoarseDimension dof, int degree,std::string &tesselation){
    std::cout << degree << " " << dof << std::endl;
    if (degree == 1)
    {
        return solve_3D_AQAPoly_EigenModes< 1 >(meshname, dof,tesselation);
    }
    else if (degree == 2)
    {
        return solve_3D_AQAPoly_EigenModes<2>(meshname, dof,tesselation);
    }
    else if (degree == 3)
    {
        return solve_3D_AQAPoly_EigenModes<3>(meshname, dof,tesselation);
    }
}





//=========================================================================================================================

template< unsigned int Degree , typename RelaxerType>
double Execute_MG
    (
        const HierarchicalSimplexRefinableCellMesh< 3 , Degree > &simplexRefinableCellMesh ,
        int CoarseNodeDimension,
        const std::function< Point< double , 3 > ( typename SimplexMesh< 3 , Degree >::NodeMultiIndex ) > &NodePosition ,
        const std::function< bool ( typename SimplexMesh< 3 , Degree >::NodeMultiIndex ) > &IsBoundaryNode ,
        unsigned int vCycles , unsigned int gsIters,
        std::ofstream &timings_file,
        Eigen::VectorXd &solution,
        bool write_mg = true
    )
{
    typedef typename SimplexMesh< 3 , Degree >::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S;
    std::vector< Eigen::SparseMatrix< double > > P( CoarseNodeDimension );
    Timer timer;
    double time;
//     Get the system matrices
    timer.reset();
    // Todo: Dimension P(d+1,d)
    for( unsigned int d=0 ; d<(unsigned int)CoarseNodeDimension; d++ ){P[d] = simplexRefinableCellMesh.P( d+1 , d );}
    {
        Eigen::SparseMatrix< double > _P , _Pt;
        _P = simplexRefinableCellMesh.P( simplexRefinableCellMesh.maxLevel() , CoarseNodeDimension );
        _Pt = _P.transpose();
        M = _Pt * simplexRefinableCellMesh.simplexMesh().mass() * _P;
        S = _Pt * simplexRefinableCellMesh.simplexMesh().stiffness() * _P;
    }
    time = timer.elapsed();
//    std::cout << "Got system matrices: " << timer.elapsed() << std::endl;
    if(write_mg){
        timings_file << "Vcycles " << vCycles << " it " << gsIters<< "," <<time << ",,," ;
    }

    // Get the list of node multi-indices
    std::vector< NodeMultiIndex > nodeMultiIndices( simplexRefinableCellMesh.nodes( CoarseNodeDimension ) );
    for( auto const &[ nmi , idx ] : simplexRefinableCellMesh.nodeMap( CoarseNodeDimension) ) nodeMultiIndices[idx] = nmi;

    // Identify the nodes that are locked
    std::vector< bool > lockedNodes( simplexRefinableCellMesh.nodes( CoarseNodeDimension ) );
    for( unsigned int i=0 ; i<lockedNodes.size() ; i++ ) lockedNodes[i] = IsBoundaryNode( nodeMultiIndices[i] );

    Eigen::VectorXd b( simplexRefinableCellMesh.nodes( CoarseNodeDimension ) );

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

    //with verbose
//    std::cout << "================Multigrid====================" << std::endl;
//    MGSolver::Solver< RelaxerType > mgSolver( S , P, true );

//    without verbose
    MGSolver::Solver< RelaxerType > mgSolver( S , P,false );
    time = timer.elapsed();
    std::cout << "Constructed multigrid system: " << time << std::endl;
    if(write_mg)
    {
        timings_file << time << ",";
    }
   std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;


   timer.reset();
   Eigen::VectorXd x = mgSolver.solve(b,vCycles , gsIters,false );
   time = timer.elapsed();

   std::cout << "Solved multigrid system: " << time  << std::endl;
   if(write_mg)
   {
       timings_file << time  << ",";
   }
   double eIn = 0 , eOut = 0;
   Eigen::VectorXd r = b - S* x;
   for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
   std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;
   double rms = 0;
   solution = x;
   for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes( CoarseNodeDimension) ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
   std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension) ) << std::endl;
   if(write_mg)
   {
       timings_file << sqrt(rms /simplexRefinableCellMesh.nodes(CoarseNodeDimension))<< ",";
   }
   return sqrt( rms / simplexRefinableCellMesh.nodes( CoarseNodeDimension) );
}
//-------------------------------------------------------------------------------------------------
template< unsigned int Degree>
double Execute_Direct
    (
        std::function< Point< double , 3 > ( unsigned int ) > &fullVertexPositionFunction,
        const SimplexRefinableCellMesh< 3, Degree > &simplexRefinableCellMesh ,
        const std::function< Point< double ,3> ( typename SimplexMesh< 3 , Degree >::NodeMultiIndex ) > &NodePosition ,
        const std::function< bool ( typename SimplexMesh< 3 , Degree >::NodeMultiIndex ) > &IsBoundaryNode,
        std::ofstream &timings_file,
        Eigen::VectorXd &solution,
        bool write = true
    )
{
    typedef typename SimplexMesh< 3 , Degree >::NodeMultiIndex NodeMultiIndex;
    Eigen::SparseMatrix< double > M , S, P;
    Timer timer;

    // Get the system matrices
    P = simplexRefinableCellMesh.P();
//    std::cout << P << std::endl;
//    std::cout << "==================================" << std::endl;
//    timer.reset();
    M = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().mass() * simplexRefinableCellMesh.P();
    timer.reset();
    S = simplexRefinableCellMesh.Pt() * simplexRefinableCellMesh.simplexMesh().stiffness() * simplexRefinableCellMesh.P();
//    Eigen::saveMarket(S, "S.mtx");
//    std::cout << S << std::endl;
//    std::cout << "------------------" << std::endl;
    std::cout << "Got stiffness matrices: " << timer.elapsed() << std::endl;
    if(write){
        timings_file << "AQAPoly direct," <<timer.elapsed() << "," ;
    }


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

    //====================================

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
        if(write)
        {
            timings_file << timer.elapsed() << ",";
        }
        std::cout << "DoFs / Non-zero matrix entries / Entries per row: " << S.rows() << " / " << S.nonZeros() << " / " << S.nonZeros()/S.rows() << std::endl;


    timer.reset();
    Eigen::VectorXd x = solver.solve( b );
    solution = x;

    std::cout << "Solved direct system: " << timer.elapsed() << std::endl;
    if(write)
    {
        timings_file << timer.elapsed() << ",,,";
    }
    //=============Gauss Quadrature ==============

    Eigen::VectorXd PSol = P*x;
    Eigen::MatrixXd QuadWeigths {{0.310885919263300609797345733763457832992617413439748379102751750,0.310885919263300609797345733763457832992617413439748379102751750,0.310885919263300609797345733763457832992617413439748379102751750,0.112687925718015850799185652333286333810467794083662916750111314},
                                {0.454496295874350350508119473720660560934367578380943975698540157,0.454496295874350350508119473720660560934367578380943975698540157,0.045503704125649649491880526279339439065632421619056024301459843,0.042546020777081466438069428120257441777627278117123448683280996},
                                {0.045503704125649649491880526279339439065632421619056024301459843,0.454496295874350350508119473720660560934367578380943975698540157,0.454496295874350350508119473720660560934367578380943975698540157,0.042546020777081466438069428120257441777627278117123448683280996},
                                {0.454496295874350350508119473720660560934367578380943975698540157,0.045503704125649649491880526279339439065632421619056024301459843,0.045503704125649649491880526279339439065632421619056024301459843,0.042546020777081466438069428120257441777627278117123448683280996},
                                {0.0927352503108912264023239137370306052452033537605457630284133352,0.0927352503108912264023239137370306052452033537605457630284133352,0.721794249067326320793028258788908184264389938718362710914759994,0.0734930431163619495437102054863275035230912887406519102249671927},
                                {0.310885919263300609797345733763457832992617413439748379102751750,0.310885919263300609797345733763457832992617413439748379102751750,0.0673422422100981706079627987096265010221477596807548626917447502,0.112687925718015850799185652333286333810467794083662916750111314},
                                {0.0927352503108912264023239137370306052452033537605457630284133352,0.0927352503108912264023239137370306052452033537605457630284133352,0.0927352503108912264023239137370306052452033537605457630284133352,0.0734930431163619495437102054863275035230912887406519102249671927},
                                {0.0927352503108912264023239137370306052452033537605457630284133352,0.721794249067326320793028258788908184264389938718362710914759994,0.0927352503108912264023239137370306052452033537605457630284133352,0.0734930431163619495437102054863275035230912887406519102249671927},
                                {0.310885919263300609797345733763457832992617413439748379102751750,0.0673422422100981706079627987096265010221477596807548626917447502,0.310885919263300609797345733763457832992617413439748379102751750,0.112687925718015850799185652333286333810467794083662916750111314},
                                {0.721794249067326320793028258788908184264389938718362710914759994,0.0927352503108912264023239137370306052452033537605457630284133352,0.0927352503108912264023239137370306052452033537605457630284133352,0.0734930431163619495437102054863275035230912887406519102249671927},
                                {0.0673422422100981706079627987096265010221477596807548626917447502,0.310885919263300609797345733763457832992617413439748379102751750,0.310885919263300609797345733763457832992617413439748379102751750,0.112687925718015850799185652333286333810467794083662916750111314},
                                {0.045503704125649649491880526279339439065632421619056024301459843,0.045503704125649649491880526279339439065632421619056024301459843,0.454496295874350350508119473720660560934367578380943975698540157,0.042546020777081466438069428120257441777627278117123448683280996},
                                {0.045503704125649649491880526279339439065632421619056024301459843,0.454496295874350350508119473720660560934367578380943975698540157,0.045503704125649649491880526279339439065632421619056024301459843,0.042546020777081466438069428120257441777627278117123448683280996},
                                {0.454496295874350350508119473720660560934367578380943975698540157,0.045503704125649649491880526279339439065632421619056024301459843,0.454496295874350350508119473720660560934367578380943975698540157,0.042546020777081466438069428120257441777627278117123448683280996}};
    Eigen::MatrixXd Tet;

    Tet.resize(3,4);
    double l2_error = 0;
    for(unsigned  int i=0 ; i< simplexRefinableCellMesh.simplexMesh().simplices();i++){
        SimplexIndex< 3 , unsigned int > s =  simplexRefinableCellMesh.simplexMesh().simplex(i);
        for( unsigned int n=0 ; n<4 ; n++)
        {
            Point<double, 3> p = fullVertexPositionFunction(s.idx[n]);
            Tet.col(n) = Eigen::Vector3d(p[0], p[1],p[2]);
        }
        Eigen::Matrix3d A;

        A.col(0) = Tet.col(1) - Tet.col(0);
        A.col(1) = Tet.col(2) - Tet.col(0);
        A.col(2) = Tet.col(3) - Tet.col(0);

        double det = A.determinant() / 6.;
        //Iterate over Quadrature points
        for(unsigned int q = 0; q<QuadWeigths.rows(); q++){
            // initialize Sample
            Eigen::Vector4d bcCoordinates;
            Eigen::Vector3d globalPos;
            Point< double , 3 > p(QuadWeigths.row(q)(0),QuadWeigths.row(q)(1),QuadWeigths.row(q)(2)) ;
            typename SimplexMesh< 3 >::Sample sample;

            sample.sIdx = i;
            sample.bcCoordinates[0] = 1. - p[0] - p[1] - p[2];
            sample.bcCoordinates[1] = p[0];
            sample.bcCoordinates[2] = p[1];
            sample.bcCoordinates[3] = p[2];
            bcCoordinates <<1.-p[0]-p[1]-p[2],p[0],p[1],p[2];
            globalPos = Tet*bcCoordinates;
            Point< double , 3 > globalP(globalPos[0],globalPos[1],globalPos[2]);
            double v_exact =  Franke(globalP);
            double v_approx = simplexRefinableCellMesh.simplexMesh().evaluate(PSol,sample);
            l2_error+=(v_exact-v_approx)*(v_exact-v_approx)*det*QuadWeigths.row(q)(3);
        }
    }
    l2_error = sqrt(fabs(l2_error));
    std::cout << "L2 error: " << l2_error << std::endl;
    //===========================

    double eIn = 0 , eOut = 0;
    Eigen::VectorXd r = b - S* x;
    for( unsigned int i=0 ; i<S.rows() ; i++ ) eIn += b[i] * b[i] , eOut += r[i] * r[i];
    std::cout << "Error: " << sqrt(eIn) << " -> " << sqrt(eOut) << std::endl;

    double rms = 0;
    for( unsigned int i=0 ; i<simplexRefinableCellMesh.nodes() ; i++ ) rms += pow( x[i] - Franke( NodePosition( nodeMultiIndices[i] ) ) , 2. );
    std::cout << "RMS: " << sqrt( rms / simplexRefinableCellMesh.nodes() ) << std::endl;
    if(write)
    {
        timings_file << sqrt(rms / simplexRefinableCellMesh.nodes()) << ",,"
                     << simplexRefinableCellMesh.nodes() << ",\n";
    }
    return sqrt( rms / simplexRefinableCellMesh.nodes() );
}
//=================================================================================
template<unsigned int Degree>
double solve_3D_AQAPoly_Poisson_mg(std::string &meshname,std::ofstream &timings_file, CoarseDimension dof, bool direct , MG_Solver solver , int vcycles , int iterations){
    const unsigned int Dim = 3;
    std::vector<std::vector<std::pair<unsigned int, bool>>> polyhedra;
    std::vector<Meshes::Polygon<unsigned int>> polygons;
    std::vector<Point<double, 3>> vertices;
    std::vector<Eigen::Triplet<double>> trip;
    ReadPolyhedra(meshname,polyhedra,polygons,vertices);

    Meshes::PolyhedronMesh<unsigned int> polyMesh = Meshes::PolyhedronMesh<unsigned int>(polyhedra, polygons);
    typedef typename SimplexMesh<Dim, Degree>::NodeMultiIndex NodeMultiIndex;

    typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;

    std::set< unsigned int > boundaryVertices;
    {
        std::vector< unsigned int > faceCount( polyMesh.polygons() , 0 );
        for( unsigned int i=0 ; i<polyMesh.polyhedra() ; i++ ) for( unsigned int j=0 ; j<polyMesh.polyhedronFaces(i) ; j++ )
                faceCount[ polyMesh.polyhedronFace(i,j).first ]++;
        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ ) if( faceCount[i]==1 )
            {
                Meshes::Polygon< unsigned int > polygon = polyMesh.polygon( i );
                for( unsigned int j=0 ; j<polygon.size() ; j++ ) boundaryVertices.insert( polygon[j] );
                boundaryVertices.insert( (unsigned int)vertices.size() + i );
            }
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , Dim > ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double , Dim >();
    };
    std::function< Point< double , Dim > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , Dim > p;
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
        case Virtual_Edges:
            coarse_dim = 2;
            break;
        case Refined_mesh:
            coarse_dim = 3;
            break;
        default:
            std::cout << "should not happen";
    }
//    timings_file << "Coarse Dim " << coarse_dim <<" ";
    double error;
    Eigen::VectorXd solution;
    // solve system directly
    if(direct)
    {
        SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim);
        error =  Execute_Direct< Degree >(fullVertexPositionFunction,simplexRefinableCellMesh , NodePosition , IsBoundaryNode,timings_file,solution );
    }else // solve via multigrid
    {
        HierarchicalSimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim);
        switch( solver )
        {
            case MG_Solver::RELAXER_JACOBI:
                error = Execute_MG< Degree , MGSolver::JacobiRelaxer>( simplexRefinableCellMesh ,coarse_dim, NodePosition , IsBoundaryNode , vcycles , iterations,timings_file,solution );
                break;
            case MG_Solver::RELAXER_GAUSS_SEIDEL:
                error = Execute_MG< Degree , MGSolver::GaussSeidelRelaxer >( simplexRefinableCellMesh,coarse_dim, NodePosition , IsBoundaryNode , vcycles , iterations,timings_file,solution  );
                break;
            case MG_Solver::RELAXER_PARALLEL_GAUSS_SEIDEL:
                error = Execute_MG< Degree , MGSolver::ParallelGaussSeidelRelaxer< 20 > >( simplexRefinableCellMesh ,coarse_dim, NodePosition , IsBoundaryNode ,vcycles , iterations ,timings_file,solution  );
                break;
            default: ERROR_OUT( "Unrecognized relaxer type: " , solver );
        }
        timings_file << "\n " ;
    }

    std::cout << "This error came out : " << error << std::endl;
    return error;
}

double solve_3D_AQAPoly_Poisson_mg(std::string &meshname, std::ofstream &timings_file,CoarseDimension dof, int degree, bool direct, MG_Solver solver, int vcycles, int iterations){
    if (degree == 1)
    {
        return solve_3D_AQAPoly_Poisson_mg<1>(meshname,timings_file,dof, direct, solver,vcycles, iterations);
    }
    else if (degree == 2)
    {
        return solve_3D_AQAPoly_Poisson_mg<2>(meshname,timings_file,dof, direct, solver,vcycles, iterations);
    }
    else if (degree == 3)
    {
        return solve_3D_AQAPoly_Poisson_mg<3>(meshname,timings_file,dof, direct, solver,vcycles, iterations);
    }
}

//================================================================================================================================

template<unsigned int Degree>
double compare_3D_AQAPoly_Poisson_mg(std::string &meshname,std::ofstream &timings_file, CoarseDimension dof, bool direct, MG_Solver solver, int vcycles, int iterations, bool write_direct){
    const unsigned int Dim = 3;
    std::vector<std::vector<std::pair<unsigned int, bool>>> polyhedra;
    std::vector<Meshes::Polygon<unsigned int>> polygons;
    std::vector<Point<double, 3>> vertices;
    std::vector<Eigen::Triplet<double>> trip;
    ReadPolyhedra(meshname,polyhedra,polygons,vertices);

    Meshes::PolyhedronMesh<unsigned int> polyMesh = Meshes::PolyhedronMesh<unsigned int>(polyhedra, polygons);
    typedef typename SimplexMesh<Dim, Degree>::NodeMultiIndex NodeMultiIndex;

    typedef typename SimplexMesh< Dim , Degree >::NodeMultiIndex NodeMultiIndex;

    std::set< unsigned int > boundaryVertices;
    {
        std::vector< unsigned int > faceCount( polyMesh.polygons() , 0 );
        for( unsigned int i=0 ; i<polyMesh.polyhedra() ; i++ ) for( unsigned int j=0 ; j<polyMesh.polyhedronFaces(i) ; j++ )
                faceCount[ polyMesh.polyhedronFace(i,j).first ]++;
        for( unsigned int i=0 ; i<polyMesh.polygons() ; i++ ) if( faceCount[i]==1 )
            {
                Meshes::Polygon< unsigned int > polygon = polyMesh.polygon( i );
                for( unsigned int j=0 ; j<polygon.size() ; j++ ) boundaryVertices.insert( polygon[j] );
                boundaryVertices.insert( (unsigned int)vertices.size() + i );
            }
    }

    auto IsBoundaryNode = [&]( NodeMultiIndex nmi )
    {
        for( unsigned int d=0 ; d<Degree ; d++ ) if( boundaryVertices.find( nmi[d] )==boundaryVertices.end() ) return false;
        return true;
    };

    SimplexRefinableElements<>::EnergyWeights eWeights( SimplexRefinableElements<>::EnergyWeights::CROSS_FACE_GRADIENT_DIFFERENCE );
    eWeights[ SimplexRefinableElements<>::EnergyWeights::GRADIENT_SQUARE_NORM ] = 1e-8;

    std::function< Point< double , Dim > ( unsigned int ) > vertexPositionFunction = [&]( unsigned int idx )
    {
        if( idx<vertices.size() ) return vertices[idx];
        ERROR_OUT( "Bad vertex index: " , idx , " / " ,  vertices.size() );
        return Point< double , Dim >();
    };
    std::function< Point< double , Dim > ( unsigned int ) > fullVertexPositionFunction = polyMesh.fullVertexPositionFunction( vertexPositionFunction , true );

    auto NodePosition = [&]( NodeMultiIndex nmi )
    {
        Point< double , Dim > p;
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
        case Virtual_Edges:
            coarse_dim = 2;
            break;
        case Refined_mesh:
            coarse_dim = 3;
            break;
        default:
            std::cout << "should not happen";
    }
//    timings_file << "Coarse Dim " << coarse_dim <<" ";
    double error_direct, error_mg;
    Eigen::VectorXd solution_direct, solution_mg;

    // solve system directly

        SimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMesh;
        simplexRefinableCellMesh = polyMesh.template simplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim);
        error_direct =  Execute_Direct< Degree >( fullVertexPositionFunction,simplexRefinableCellMesh , NodePosition , IsBoundaryNode,timings_file,solution_direct, write_direct);
    // solve via multigrid

        HierarchicalSimplexRefinableCellMesh< Dim , Degree > simplexRefinableCellMeshMG;
        simplexRefinableCellMeshMG = polyMesh.template hierarchicalSimplexRefinableCellMesh< Degree >( fullVertexPositionFunction , eWeights , coarse_dim);
        switch( solver )
        {
            case MG_Solver::RELAXER_JACOBI:
                error_mg = Execute_MG< Degree , MGSolver::JacobiRelaxer>( simplexRefinableCellMeshMG ,coarse_dim, NodePosition , IsBoundaryNode , vcycles , iterations,timings_file,solution_mg, !direct);
                break;
            case MG_Solver::RELAXER_GAUSS_SEIDEL:
                error_mg = Execute_MG< Degree , MGSolver::GaussSeidelRelaxer >( simplexRefinableCellMeshMG,coarse_dim, NodePosition , IsBoundaryNode , vcycles , iterations,timings_file,solution_mg,!direct);
                break;
            case MG_Solver::RELAXER_PARALLEL_GAUSS_SEIDEL:
                error_mg = Execute_MG< Degree , MGSolver::ParallelGaussSeidelRelaxer< 20 > >( simplexRefinableCellMeshMG ,coarse_dim, NodePosition , IsBoundaryNode ,vcycles , iterations ,timings_file,solution_mg,!direct);
                break;
            default: ERROR_OUT( "Unrecognized relaxer type: " , solver );
        }
    std::cout << "This error came out : " << error_mg << std::endl;
    double solution_rms=0.0;
    for( unsigned int i=0 ; i<simplexRefinableCellMeshMG.nodes(coarse_dim) ; i++ ) solution_rms += pow( solution_mg(i)- solution_direct(i), 2. );
    if(!direct){
        timings_file<< sqrt( solution_rms / simplexRefinableCellMeshMG.nodes( coarse_dim) ) << "," <<  simplexRefinableCellMeshMG.nodes( coarse_dim) << std::endl;
    }
    std::cout << "RMS: " << sqrt( solution_rms / simplexRefinableCellMeshMG.nodes(coarse_dim)) << std::endl;
    if(direct){
        return  error_direct;
    }
    else{
        return error_mg;
    }
}


double solve_3D_AQAPoly_Poisson_mg_extended(std::string &meshname, std::ofstream &timings_file,CoarseDimension dof, int degree, bool direct, MG_Solver solver, int vcycles, int iterations, bool write_direct){
    if (degree == 1)
    {
        return compare_3D_AQAPoly_Poisson_mg<1>(meshname,timings_file,dof, direct, solver,vcycles, iterations, write_direct);
    }
    else if (degree == 2)
    {
        return compare_3D_AQAPoly_Poisson_mg<2>(meshname,timings_file,dof, direct, solver,vcycles, iterations, write_direct);
    }
    else if (degree == 3)
    {
        return compare_3D_AQAPoly_Poisson_mg<3>(meshname,timings_file,dof, direct, solver,vcycles, iterations, write_direct);
    }
}