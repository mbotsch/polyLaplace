cmake_minimum_required(VERSION 3.1)
project(Eigen)

add_library(eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS eigen)

target_include_directories(eigen INTERFACE include/)
