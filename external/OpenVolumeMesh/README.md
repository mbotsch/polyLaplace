-----
                            OpenVolumeMesh                                 
        Copyright (C) 2011-2018 by Computer Graphics Group, RWTH Aachen         
                         www.openvolumemesh.org                             
-----
                                                                           
  OpenVolumeMesh is free software: you can redistribute it and/or modify   
  it under the terms of the GNU Lesser General Public License as           
  published by the Free Software Foundation, either version 3 of           
  the License, or (at your option) any later version with the              
  following exceptions:                                                    
                                                                           
  OpenVolumeMesh is distributed in the hope that it will be useful,        
  but WITHOUT ANY WARRANTY; without even the implied warranty of           
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            
  GNU Lesser General Public License for more details.                      
                                                                           
  You should have received a copy of the GNU LesserGeneral Public          
  License along with OpenVolumeMesh. If not,                              
  see <http://www.gnu.org/licenses/>.                                      
                                                                           


# Overview:

0. Introduction
1. System Requirements
2. Building OpenVolumeMesh
3. License Information


# 0. Introduction

Thank you for downloading and using the OpenVolumeMesh library. OpenVolumeMesh
is a generic data structure for the comfortable handling of arbitrary
polyhedral meshes. Its concepts are closely related to OpenMesh
<http://www.openmesh.org>. In particular, OpenVolumeMesh carries the general
idea of storing edges as so-called (directed) half-edges over to the face
definitions. So, faces are split up into so-called half-faces having opposing
orientations. Furthermore, in OpenVolumeMesh the data is arranged in a top-down
hierarchy, meaning that each entity of dimension n is defined through a
(ordered) tuple of entities of dimension (n-1). These are the intrinsic
adjacency relations of the volumentric meshes. One can additionally compute
bottom-up adjacencies which means that for each entity of dimension n, we also
store its local adjacencies to entities of dimension (n+1). These adjacency
relations have to be computed explicitly which can be performed in linear time
complexity. Both adjacency relations, the top-down and the bottom-up
adjacencies, are used to provide a set of iterators and circulators that are
comfortable in use. As in OpenMesh, OpenVolumeMesh provides an entirely generic
underlying property system that allows attaching properties of any kind to the
entities. In order to learn more about the implementational details of
OpenVolumeMesh, please refer to the only documentation which is available
at <http://openvolumemesh.org/Documentation/OpenVolumeMesh-Doc-Latest>.

OpenVolumeMesh is entirely written in C++ making heavy use of the
standard template library as well as template programming paradigms.
Although OpenVolumeMesh has been developed to the best of my knowledge,
it does not claim to be free from defects nor does it raises the claim to
have inveterate underlying implemented concepts. So, any ambitious developer
is invited to participate in the development process to make OpenVolumeMesh
a well-working, reliable, and useful library. Please feel free to commit
suggestions, bug reports, or patches to the OpenVolumeMesh project management
system (redmine), which you find at

<https://www.graphics.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh>,

or you can send them directly to the following e-mail address which is

<moebius@cs.rwth-aachen.de>.


# 1. System Requirements

OpenVolumeMesh is shipped as source project and can be built on all common
architectures and operating systems including:

- Windows 7/8/10 with Visual Studio 2013 or higher

- MacOSX with latest XCode 

- Linux with GCC 6 and higher

Note that OpenVolumeMesh uses CMake as build system, so make sure that
the latest version of CMake is installed on your computer. Download
CMake under <http://www.cmake.org>.

Note also that, in order to build the documentation, you will need to
have Doxygen installed on your computer. Download Doxygen under
<http://www.doxygen.org>. The use of Doxygen is not mandatory though.


# 2. Building OpenVolumeMesh

OpenVolumeMesh is equipped with a CMake build system <http://www.cmake.org>.
Make sure that at least version 3.7 of CMake is installed on your computer.
Once CMake is installed, perform the following steps (on Linux, Mac OSX):

- Create a build directory somewhere outside the source directory
  of OpenVolumeMesh. Name it e.g. "OVM-build-release", or whatever
  name might be suitable.

- Change into the recently created directory and type
  "cmake /path/to/OpenVolumeMesh/sources".

- If you want to change the build configuration, say from debug to release,
  type "ccmake .". Note that unless not explicitly specified otherwise, CMake
  sets the build configuration such that it will be built with debug flags on.

- Once everything is configured to your satisfaction, type "make" followed
  by "sudo make install" in order to build and install the library.

On Windows, start the CMake gui tool and set the path to OpenVolumeMesh's
sources. Then select the target project type (Visual Studio 20xx)
and click on "Configure". Once everything is configured to your satisfaction,
click on "Generate". You will now find a Visual Studio project file
in the specified build folder (which is "Build" per default). Open this
file in Visual Studio and select "Build all".


# 3. License Information

OpenVolumeMesh is free software licensed under the terms of the
GNU Lesser General Public License Version 3 as published by the Free Software
Foundation. You can redistribute and/or modify it as stated in the
above mentioned license terms. A copy of the license can be found
in the license sub-folder of this source-tree or under
<http://www.gnu.org/licenses>. By downloading and using the OpenVolumeMesh
library you automatically agree to these terms.
