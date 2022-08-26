/*===========================================================================*\
 *                                                                           *
 *                            OpenVolumeMesh                                 *
 *        Copyright (C) 2011 by Computer Graphics Group, RWTH Aachen         *
 *                        www.openvolumemesh.org                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenVolumeMesh.                                     *
 *                                                                           *
 *  OpenVolumeMesh is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenVolumeMesh is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenVolumeMesh.  If not,                              *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision$                                                         *
 *   $Date$                    *
 *   $LastChangedBy$                                                *
 *                                                                           *
\*===========================================================================*/

#ifndef TETRAHEDRALMESH_HH_
#define TETRAHEDRALMESH_HH_

#include "TetrahedralMeshTopologyKernel.hh"
#include "../Core/GeometryKernel.hh"

namespace OpenVolumeMesh {

/*
 * Predefines for most common mesh types
 */
typedef GeometryKernel<Geometry::Vec2i, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV2i;
typedef GeometryKernel<Geometry::Vec2ui, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV2ui;
typedef GeometryKernel<Geometry::Vec2f, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV2f;
typedef GeometryKernel<Geometry::Vec2d, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV2d;
typedef GeometryKernel<Geometry::Vec2c, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV2c;
typedef GeometryKernel<Geometry::Vec2uc, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV2uc;
typedef GeometryKernel<Geometry::Vec3i, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV3i;
typedef GeometryKernel<Geometry::Vec3ui, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV3ui;
typedef GeometryKernel<Geometry::Vec3f, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV3f;
typedef GeometryKernel<Geometry::Vec3d, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV3d;
typedef GeometryKernel<Geometry::Vec3c, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV3c;
typedef GeometryKernel<Geometry::Vec3uc, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV3uc;
typedef GeometryKernel<Geometry::Vec4i, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV4i;
typedef GeometryKernel<Geometry::Vec4ui, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV4ui;
typedef GeometryKernel<Geometry::Vec4f, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV4f;
typedef GeometryKernel<Geometry::Vec4d, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV4d;
typedef GeometryKernel<Geometry::Vec4c, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV4c;
typedef GeometryKernel<Geometry::Vec4uc, TetrahedralMeshTopologyKernel> GeometricTetrahedralMeshV4uc;

typedef TetrahedralMeshTopologyKernel TopologicTetrahedralMesh;

} // Namespace OpenVolumeMesh

#endif /* TETRAHEDRALMESH_HH_ */
