#!/bin/bash

COMPILER=$1
LANGUAGE=$2

# Exit script on any error
set -e 

OPTIONS=""
BUILDPATH=""

if [ "$COMPILER" == "gcc" ]; then
  echo "Building with GCC";
  BUILDPATH="gcc"
elif [ "$COMPILER" == "clang" ]; then

  OPTIONS="$OPTIONS -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_C_COMPILER=clang"
  echo "Building with CLANG";
  BUILDPATH="clang"  
fi  

if [ "$LANGUAGE" == "C++98" ]; then
  echo "Building with C++98";
  BUILDPATH="$BUILDPATH-cpp98"
elif [ "$LANGUAGE" == "C++11" ]; then
  echo "Building with C++11";
  OPTIONS="$OPTIONS -DCMAKE_CXX_FLAGS='-std=c++11' "
  BUILDPATH="$BUILDPATH-cpp11"  
fi  

#########################################

# Make release build folder
if [ ! -d build-release-$BUILDPATH ]; then
  mkdir build-release-$BUILDPATH
fi

cd build-release-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Release -DSTL_DEBUG=ON $OPTIONS ../

#build it
make

#build the unit tests
make unittests


#########################################
# Run Release Unittests
#########################################
cd Unittests

#execute tests
./unittests --gtest_color=yes --gtest_output=xml

cd ..
cd ..


#########################################
# Build Debug version and Unittests
#########################################

if [ ! -d build-debug-$BUILDPATH ]; then
  mkdir build-debug-$BUILDPATH
fi

cd build-debug-$BUILDPATH

cmake -DCMAKE_BUILD_TYPE=Debug -DSTL_DEBUG=ON $OPTIONS ../

#build the unit tests
make unittests


#########################################
# Run Debug Unittests
#########################################

cd Unittests

# Run the unittests
./unittests --gtest_color=yes --gtest_output=xml


