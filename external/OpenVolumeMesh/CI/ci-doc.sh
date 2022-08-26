#!/bin/bash

#########################################

# Make release build folder
if [ ! -d build-release ]; then
  mkdir build-release
fi

cd build-release

cmake -DCMAKE_BUILD_TYPE=Release ../

#build it
make doc
