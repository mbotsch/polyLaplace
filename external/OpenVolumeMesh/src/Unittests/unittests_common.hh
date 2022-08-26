#ifndef INCLUDE_UNITTESTS_COMMON_HH
#define INCLUDE_UNITTESTS_COMMON_HH


#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Core/PropertyDefines.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

#ifdef __clang__
#  pragma GCC diagnostic ignored "-Weverything"
#  pragma GCC diagnostic ignored "-Wundef"
#  pragma GCC diagnostic ignored "-Wglobal-constructors"
#  pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#  pragma GCC diagnostic ignored "-Wmissing-noreturn"
#endif

#include <gtest/gtest.h>

#define EXPECT_HANDLE_EQ(a, b)  EXPECT_EQ((a).idx(), (b).idx())
#define EXPECT_HANDLE_NE(a, b)  EXPECT_NE((a).idx(), (b).idx())


/*
 * Simple test setting for polyhedral meshes
 */

typedef OpenVolumeMesh::GeometricPolyhedralMeshV3d PolyhedralMesh;

class PolyhedralMeshBase: public testing::Test {

protected:

  typedef OpenVolumeMesh::VertexHandle    VertexHandle;
  typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
  typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
  typedef OpenVolumeMesh::FaceHandle      FaceHandle;
  typedef OpenVolumeMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp() {

    // Do some initial stuff with the member data here...
    mesh_.enable_deferred_deletion(false);
    mesh_.enable_fast_deletion(false);
  }

  // This function is called after all tests are through
  virtual void TearDown() {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generatePolyhedralMesh(PolyhedralMesh& _mesh);

  // This member will be accessible in all tests
  PolyhedralMesh mesh_;
};

/*
 * Simple test setting for hexahedral meshes
 */

typedef OpenVolumeMesh::GeometricHexahedralMeshV3d HexahedralMesh;

class HexahedralMeshBase: public testing::Test {

protected:

  typedef OpenVolumeMesh::VertexHandle    VertexHandle;
  typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
  typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
  typedef OpenVolumeMesh::FaceHandle      FaceHandle;
  typedef OpenVolumeMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp() {

    // Do some initial stuff with the member data here...
    mesh_.enable_deferred_deletion(false);
    mesh_.enable_fast_deletion(false);
  }

  // This function is called after all tests are through
  virtual void TearDown() {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generateHexahedralMesh(HexahedralMesh& _mesh);

  // This member will be accessible in all tests
  HexahedralMesh mesh_;
};


/*
 * Simple test setting for tetrahedral meshes
 */

typedef OpenVolumeMesh::GeometricTetrahedralMeshV3d TetrahedralMesh;

class TetrahedralMeshBase: public testing::Test {

protected:

  typedef OpenVolumeMesh::VertexHandle    VertexHandle;
  typedef OpenVolumeMesh::HalfEdgeHandle  HalfEdgeHandle;
  typedef OpenVolumeMesh::EdgeHandle      EdgeHandle;
  typedef OpenVolumeMesh::HalfFaceHandle  HalfFaceHandle;
  typedef OpenVolumeMesh::FaceHandle      FaceHandle;
  typedef OpenVolumeMesh::CellHandle      CellHandle;

  // This function is called before each test is run
  virtual void SetUp() {

    // Do some initial stuff with the member data here...
    mesh_.enable_deferred_deletion(false);
    mesh_.enable_fast_deletion(false);
  }

  // This function is called after all tests are through
  virtual void TearDown() {

    // Do some final stuff with the member data here...
  }

  // Generate a basic hexahedral mesh
  void generateTetrahedralMesh(TetrahedralMesh& _mesh);

  // This member will be accessible in all tests
  TetrahedralMesh mesh_;
};


// Printer class (for STL compliance test)
class Print {
public:
  explicit Print(bool _mute = false) : mute_(_mute) {}
  void mute(bool _mute) { mute_ = _mute; }
  bool mute() const { return mute_; }
  void operator()(const OpenVolumeMesh::OpenVolumeMeshHandle& _h) const {
    if(!mute_) std::cerr << "Handle: " << _h.idx() << std::endl;
  }
private:
  bool mute_;
};

#endif // INCLUDE GUARD
