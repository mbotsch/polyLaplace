/*
 * Vector.cc
 *
 *  Created on: Aug 8, 2013
 *      Author: kremer
 */

#include <string>

#include "VectorT.hh"

namespace OpenVolumeMesh {

using namespace Geometry;

template <> const std::string typeName<Vec2f>() { return "vec2f"; }
template <> const std::string typeName<Vec2d>() { return "vec2d"; }
template <> const std::string typeName<Vec2i>() { return "vec2i"; }
template <> const std::string typeName<Vec2ui>() { return "vec2ui"; }

template <> const std::string typeName<Vec3f>() { return "vec3f"; }
template <> const std::string typeName<Vec3d>() { return "vec3d"; }
template <> const std::string typeName<Vec3i>() { return "vec3i"; }
template <> const std::string typeName<Vec3ui>() { return "vec3ui"; }

template <> const std::string typeName<Vec4f>() { return "vec4f"; }
template <> const std::string typeName<Vec4d>() { return "vec4d"; }
template <> const std::string typeName<Vec4i>() { return "vec4i"; }
template <> const std::string typeName<Vec4ui>() { return "vec4ui"; }

}
