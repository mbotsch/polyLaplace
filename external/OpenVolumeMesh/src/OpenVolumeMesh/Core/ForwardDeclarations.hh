#pragma once
#include "Entities.hh"

namespace OpenVolumeMesh {

class BaseProperty;

template <class T>
class OpenVolumeMeshPropertyT;

template <class PropT, typename Entity>
class PropertyPtr;

template<typename T, typename Entity>
class PropertyTT;

template<typename T> using VertexPropertyT   = PropertyTT<T, Entity::Vertex>;
template<typename T> using EdgePropertyT     = PropertyTT<T, Entity::Edge>;
template<typename T> using HalfEdgePropertyT = PropertyTT<T, Entity::HalfEdge>;
template<typename T> using FacePropertyT     = PropertyTT<T, Entity::Face>;
template<typename T> using HalfFacePropertyT = PropertyTT<T, Entity::HalfFace>;
template<typename T> using CellPropertyT     = PropertyTT<T, Entity::Cell>;
template<typename T> using MeshPropertyT     = PropertyTT<T, Entity::Mesh>;

} // namespace OVM
