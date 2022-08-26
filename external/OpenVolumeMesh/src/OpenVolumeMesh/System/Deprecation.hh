#pragma once

#include "OpenVolumeMesh/Config/DeprecationConfig.hh"

#if defined(__cplusplus) && (__cplusplus >= 201402L)
#  define OVM_DEPRECATED(msg) [[deprecated(msg)]]
#elif defined(__GNUC__) || defined(__clang__)
#  define OVM_DEPRECATED(msg) __attribute__((deprecated))
#elif defined(_MSC_VER)
#  define OVM_DEPRECATED(msg) __declspec(deprecated)
#else
#  pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#  define OVM_DEPRECATED(msg)
#endif


