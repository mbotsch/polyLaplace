#pragma once

#include "OpenVolumeMesh/Config/Export.hh"
#include <type_traits>

namespace OpenVolumeMesh {

namespace Entity {
    struct OVM_EXPORT Vertex   { Vertex() = delete;};
    struct OVM_EXPORT Edge     {};
    struct OVM_EXPORT HalfEdge {};
    struct OVM_EXPORT Face     {};
    struct OVM_EXPORT HalfFace {};
    struct OVM_EXPORT Cell     {};
    struct OVM_EXPORT Mesh     {};
}

template<typename T>
struct is_entity : std::false_type {};

template<> struct OVM_EXPORT is_entity<Entity::Vertex>   : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Edge>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::HalfEdge> : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Face>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::HalfFace> : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Cell>     : std::true_type {};
template<> struct OVM_EXPORT is_entity<Entity::Mesh>     : std::true_type {};

} // namespace OpenVolumeMesh
