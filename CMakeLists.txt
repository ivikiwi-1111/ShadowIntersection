cmake_minimum_required(VERSION 3.23)
project(shadow_intersection.0)

set(CMAKE_CXX_STANDARD 14)

add_executable(shadow_intersection_modified main.cpp vec3.h quartic.h orbit_params.h calc_shadow_intersections.h)
