#include "OpenVolumeMesh/Config/Version.hh"

// C++ version when compiling (e.g. client code),
// be careful not to change ABI depending on these defines.

#if __cplusplus >= 201703L
    #define OVM_CXX_17 1
#else
    #define OVM_CXX_17 0
#endif
