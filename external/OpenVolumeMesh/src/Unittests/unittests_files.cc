#include "unittests_common.hh"

#include <OpenVolumeMesh/FileManager/FileManager.hh>

using namespace OpenVolumeMesh;

using namespace Geometry;

TEST_F(PolyhedralMeshBase, LoadFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, LoadNonManifoldMesh) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("NonManifold.ovm", mesh_));

  EXPECT_EQ(12u, mesh_.n_vertices());
  EXPECT_EQ(20u, mesh_.n_edges());
  EXPECT_EQ(11u, mesh_.n_faces());
  EXPECT_EQ(2u, mesh_.n_cells());
}

TEST_F(HexahedralMeshBase, LoadFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.isHexahedralMesh("Cylinder.ovm"));
  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, SaveFile) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  // Write file
  ASSERT_TRUE(fileManager.writeFile("Cylinder.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cylinder.copy.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());
}

TEST_F(PolyhedralMeshBase, SaveFileWithProps) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  // Attach non-persistent properties
  HalfFacePropertyT<float> hfprop = mesh_.request_halfface_property<float>("MyHalfFaceProp");
  VertexPropertyT<unsigned int> vprop = mesh_.request_vertex_property<unsigned int>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      hfprop[HalfFaceHandle(i)] = (float)i/2.0f;
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      vprop[VertexHandle(i)] = i;
  }

  mesh_.set_persistent(hfprop);
  mesh_.set_persistent(vprop);

  // Write file
  ASSERT_TRUE(fileManager.writeFile("Cylinder.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cylinder.copy.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(1u, mesh_.n_vertex_props());

  HalfFacePropertyT<float> hfprop2 = mesh_.request_halfface_property<float>("MyHalfFaceProp");
  VertexPropertyT<unsigned int> vprop2 = mesh_.request_vertex_property<unsigned int>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      EXPECT_FLOAT_EQ((float)i/2.0f, hfprop2[HalfFaceHandle(i)]);
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      EXPECT_EQ(i, vprop2[VertexHandle(i)]);
  }
}

TEST_F(PolyhedralMeshBase, SaveFileWithVectorProps) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  // Attach non-persistent properties
  HalfFacePropertyT<Vec3d> hfprop = mesh_.request_halfface_property<Vec3d>("MyHalfFaceProp");
  VertexPropertyT<Vec2i> vprop = mesh_.request_vertex_property<Vec2i>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      hfprop[HalfFaceHandle(i)] = Vec3d((double)i/2.0, (double)i/2.0, (double)i/2.0);
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      vprop[VertexHandle(i)] = Vec2i(i, i);
  }

  mesh_.set_persistent(hfprop);
  mesh_.set_persistent(vprop);

  // Write file
  ASSERT_TRUE(fileManager.writeFile("Cylinder.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cylinder.copy.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(1u, mesh_.n_vertex_props());

  HalfFacePropertyT<Vec3d> hfprop2 = mesh_.request_halfface_property<Vec3d>("MyHalfFaceProp");
  VertexPropertyT<Vec2i> vprop2 = mesh_.request_vertex_property<Vec2i>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      HalfFaceHandle hfh(i);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop2[hfh][0]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop2[hfh][1]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop2[hfh][2]);
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      VertexHandle vh(i);
      EXPECT_EQ((int)i, vprop2[vh][0]);
      EXPECT_EQ((int)i, vprop2[vh][1]);
  }
}

TEST_F(PolyhedralMeshBase, SerializeVectorValuedProperties) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cylinder.ovm", mesh_));

  EXPECT_EQ(399u, mesh_.n_vertices());
  EXPECT_EQ(1070u, mesh_.n_edges());
  EXPECT_EQ(960u, mesh_.n_faces());
  EXPECT_EQ(288u, mesh_.n_cells());

  // Attach persistent properties
  HalfFacePropertyT<Vec3d> hfprop = mesh_.request_halfface_property<Vec3d>("MyHalfFaceProp");
  VertexPropertyT<Vec2i> vprop = mesh_.request_vertex_property<Vec2i>("MyVertexProp");

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      HalfFaceHandle hfh(i);
      hfprop[hfh] = Vec3d((double)i/2.0, (double)i/2.0, (double)i/2.0);
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      VertexHandle vh(i);
      vprop[vh] = Vec2i(i, i);
  }

  mesh_.set_persistent(hfprop);
  mesh_.set_persistent(vprop);

  std::ofstream ofs1("hfVecPropTest");
  std::ofstream ofs2("vVecPropTest");

  hfprop.serialize(ofs1);
  vprop.serialize(ofs2);

  ofs1.close();
  ofs2.close();

  /*
   * Change property values
   */
  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
	  hfprop[HalfFaceHandle(i)] = Vec3d((double)i/3.0, (double)i/3.0, (double)i/3.0);
  }

  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
	  vprop[VertexHandle(i)] = Vec2i(2*i, 2*i);
  }

  std::ifstream ifs1("hfVecPropTest");
  std::ifstream ifs2("vVecPropTest");

  hfprop.deserialize(ifs1);
  vprop.deserialize(ifs2);

  ifs1.close();
  ifs2.close();

  for(unsigned int i = 0; i < mesh_.n_halffaces(); ++i) {
      HalfFaceHandle hfh(i);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop[hfh][0]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop[hfh][1]);
      EXPECT_DOUBLE_EQ((double)i/2.0, hfprop[hfh][2]);
  }
  for(unsigned int i = 0; i < mesh_.n_vertices(); ++i) {
      VertexHandle vh(i);
      EXPECT_EQ((int)i, vprop[vh][0]);
      EXPECT_EQ((int)i, vprop[vh][1]);
  }
}

TEST_F(PolyhedralMeshBase, LoadFileWithProps) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_vertex_props());
  EXPECT_EQ(1u, mesh_.n_edge_props());
  EXPECT_EQ(0u, mesh_.n_halfedge_props());
  EXPECT_EQ(1u, mesh_.n_face_props());
  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(0u, mesh_.n_cell_props());
}

TEST_F(PolyhedralMeshBase, SaveFileWithProps2) {

  OpenVolumeMesh::IO::FileManager fileManager;

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_vertex_props());
  EXPECT_EQ(1u, mesh_.n_edge_props());
  EXPECT_EQ(0u, mesh_.n_halfedge_props());
  EXPECT_EQ(1u, mesh_.n_face_props());
  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(0u, mesh_.n_cell_props());

  ASSERT_TRUE(fileManager.writeFile("Cube_with_props.copy.ovm", mesh_));

  mesh_.clear();

  ASSERT_TRUE(fileManager.readFile("Cube_with_props.copy.ovm", mesh_));

  EXPECT_EQ(8u, mesh_.n_vertices());
  EXPECT_EQ(12u, mesh_.n_edges());
  EXPECT_EQ(6u, mesh_.n_faces());
  EXPECT_EQ(1u, mesh_.n_cells());

  EXPECT_EQ(1u, mesh_.n_vertex_props());
  EXPECT_EQ(1u, mesh_.n_edge_props());
  EXPECT_EQ(0u, mesh_.n_halfedge_props());
  EXPECT_EQ(1u, mesh_.n_face_props());
  EXPECT_EQ(1u, mesh_.n_halfface_props());
  EXPECT_EQ(0u, mesh_.n_cell_props());
}

