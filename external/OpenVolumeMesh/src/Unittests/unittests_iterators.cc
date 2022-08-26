#include "unittests_common.hh"


using namespace OpenVolumeMesh;

TEST_F(HexahedralMeshBase, HexVertexIterTest) {

    generateHexahedralMesh(mesh_);

    HexVertexIter hv_it = mesh_.hv_iter(CellHandle(0));

    EXPECT_TRUE(hv_it.valid());

    EXPECT_HANDLE_EQ(VertexHandle(0), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(1), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(2), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(3), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(4), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(7), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(6), *hv_it); ++hv_it;
    EXPECT_HANDLE_EQ(VertexHandle(5), *hv_it);
}

TEST_F(TetrahedralMeshBase, VertexVertexIteratorTest) {

    generateTetrahedralMesh(mesh_);

    {

        VertexVertexIter vv_it = mesh_.vv_iter(VertexHandle(0));

        EXPECT_TRUE(vv_it.valid());

        std::set<VertexHandle> onering;
        int valence = 0;

        while (vv_it.valid())
        {
          ++valence;
          onering.insert(*vv_it);
          ++vv_it;
        }

        // check that there have been three adjacent vertices
        EXPECT_EQ(3, valence);

        // check that no vertex was visited twice
        EXPECT_EQ(3u, onering.size());

        // check that no invalid vertex was adjacent
        EXPECT_EQ(onering.end(), std::find(onering.begin(), onering.end(), VertexHandle(-1)));

    }

#if __cplusplus >= 201103L || _MSC_VER >= 1800 // C++11
    {

      std::set<VertexHandle> onering;
      int valence = 0;

      for (auto vh : mesh_.vertex_vertices(VertexHandle(0)))
      {
        ++valence;
        onering.insert(vh);
      }

      // check that there have been three adjacent vertices
      EXPECT_EQ(3, valence);

      // check that no vertex was visited twice
      EXPECT_EQ(3u, onering.size());

      // check that no invalid vertex was adjacent
      EXPECT_EQ(onering.end(), std::find(onering.begin(), onering.end(), VertexHandle(-1)));

    }
#endif

}

TEST_F(TetrahedralMeshBase, VertexFaceIteratorTest) {

    generateTetrahedralMesh(mesh_);

    {

        VertexFaceIter vf_it = mesh_.vf_iter(VertexHandle(0));

        EXPECT_TRUE(vf_it.valid());

        std::set<FaceHandle> incident_faces;
        int valence = 0;

        while (vf_it.valid())
        {
          ++valence;
          incident_faces.insert(*vf_it);
          ++vf_it;
        }

        // check that there have been three adjacent vertices
        EXPECT_EQ(3, valence);

        // check that no vertex was visited twice
        EXPECT_EQ(3u, incident_faces.size());

        // check that no invalid vertex was adjacent
        EXPECT_EQ(incident_faces.end(), std::find(incident_faces.begin(), incident_faces.end(), FaceHandle(-1)));

    }

#if __cplusplus >= 201103L || _MSC_VER >= 1800 // C++11
    {

      std::set<VertexHandle> onering;
      int valence = 0;

      for (auto vh : mesh_.vertex_vertices(VertexHandle(0)))
      {
        ++valence;
        onering.insert(vh);
      }

      // check that there have been three adjacent vertices
      EXPECT_EQ(3, valence);

      // check that no vertex was visited twice
      EXPECT_EQ(3u, onering.size());

      // check that no invalid vertex was adjacent
      EXPECT_EQ(onering.end(), std::find(onering.begin(), onering.end(), VertexHandle(-1)));

    }
#endif

}

#if __cplusplus >= 201103L || _MSC_VER >= 1800 // C++11
TEST_F(HexahedralMeshBase, RangeForTest) {
    // no EXPECTs here, if it compiles, it'll work.
    generateHexahedralMesh(mesh_);
    VertexHandle _dummy; // use vh to avoid compiler warnings
    for (const auto& vh: mesh_.vertices()) { _dummy = vh;}
    const auto& constref = mesh_;
    for (const auto& vh: constref.vertices()) { _dummy = vh;}
}
#endif
