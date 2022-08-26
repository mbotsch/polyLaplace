#pragma once

#include <string>
#include <typeinfo>

/// Get an internal name for a type. Important: this differs between
/// compilers and versions, do NOT use in file formats!
/// We need this in order to provide property type safety when
/// only limited RTTI support is available.
template<typename T>
inline std::string get_type_name()
{
#ifdef _MSC_VER
    // MSVC's type_name only returns a friendly name with .name(),
    // get the more unique mangled name using .raw_name():
    return typeid(T).raw_name();
#else
    // GCC and clang currently return the mangled name as .name(),
    // there is no .raw_name()
    return typeid(T).name();
#endif
}

