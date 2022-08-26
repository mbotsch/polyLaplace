#pragma once

#include <pmp/visualization/GL.h>
#include <pmp/MatVec.h>
#include <vector>

namespace pmp {

//=============================================================================

//! \addtogroup visualization visualization
//! @{

//=============================================================================

//! shader class for easy handling of the shader
class MyShader
{
public:
    //! default constructor
    MyShader();

    //! default destructor
    ~MyShader();

    //! is shader valid (ID != 0)
    bool is_valid() const { return pid_ != 0; }

    //! load (from file), compile, and link vertex and fragment shader,
    //! and optional geometry and tessellation shaders.
    //! unused shaders should be NULL.
    //! \param vshader string with the adress to the vertex shader
    //! \param fshader string with the adress to the fragment shader
    //! \param gshader string with the adress to the geometry shader
    //! \param tcshader string with the adress to the tessellation control shader
    //! \param teshader string with the adress to the tessellation evaluation shader
    bool source(const char* vshader,
              const char* fshader,
              const char* gshader=nullptr,
              const char* tcshader=nullptr,
              const char* teshader=nullptr);

    //! enable/bind this shader program
    void use();

    //! disable/unbind this shader program
    void disable();

    //! bind attribute to location
    void bind_attribute(const char* name, GLuint index);

    //! upload float uniform
    //! \param name string of the uniform name
    //! \param value the value for the uniform
    void set_uniform(const char* name, float value);

    //! upload int uniform
    //! \param name string of the uniform name
    //! \param value the value for the uniform
    void set_uniform(const char* name, int value);

    //! upload vec3 uniform
    //! \param name string of the uniform name
    //! \param vec the value for the uniform
    void set_uniform(const char* name, const vec3& vec);

    //! upload vec4 uniform
    //! \param name string of the uniform name
    //! \param vec the value for the uniform
    void set_uniform(const char* name, const vec4& vec);

    //! upload mat3 uniform
    //! \param name string of the uniform name
    //! \param mat the value for the uniform
    void set_uniform(const char* name, const mat3& mat);

    //! upload mat4 uniform
    //! \param name string of the uniform name
    //! \param mat the value for the uniform
    void set_uniform(const char* name, const mat4& mat);

private:
    //! deletes all shader and frees GPU shader capacities
    void cleanup();

    //! load shader from file, return as string
    bool load(const char* filename, std::string& source);

    //! compile a vertex/fragmend shader
    //! \param shader source
    //! \param type the type of the shader (vertex, fragment)
    GLint compile(const char* source, GLenum type);

    //! loads a vertex/fragmend shader from a file and compiles it
    //! \param filename the location and name of the shader
    //! \param type the type of the shader (vertex, geometry, fragment)
    GLint load_and_compile(const char* filename, GLenum type);

    //! relink: use this after setting/changing attrib location
    bool link();

private:
    //! id of the linked shader program
    GLint pid_;

    //! id of the vertex shader
    std::vector<GLint> shaders_;
};

//=============================================================================
//! @}
//=============================================================================
} // namespace pmp
//=============================================================================
