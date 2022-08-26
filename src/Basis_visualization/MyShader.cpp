#include "MyShader.h"
#include <iostream>
#include <fstream>

//=============================================================================

namespace pmp {

//=============================================================================

    MyShader::MyShader() : pid_(0)
    {
    }

//-----------------------------------------------------------------------------

    MyShader::~MyShader()
    {
        cleanup();
    }

//-----------------------------------------------------------------------------

    void MyShader::cleanup()
    {
        if (pid_)
        {
            glDeleteProgram(pid_);
            pid_ = 0;
        }

        for (GLint id : shaders_)
        {
            glDeleteShader(id);
        }
        shaders_.clear();
    }

//-----------------------------------------------------------------------------

    bool MyShader::source(const char* vshader, const char* fshader,
                      const char* gshader, const char* tcshader, const char* teshader)
    {
        GLint id;

        // cleanup existing shaders first
        cleanup();

        // create program
        pid_ = glCreateProgram();

        // vertex shader
        id = compile(vshader, GL_VERTEX_SHADER);
        if (!id)
        {
            std::cerr << "Cannot compile vertex shader!\n";
            return false;
        }
        glAttachShader(pid_, id);
        shaders_.push_back(id);

        // fragment shader
        id = compile(fshader, GL_FRAGMENT_SHADER);
        if (!id)
        {
            std::cerr << "Cannot compile fragment shader!\n";
            return false;
        }
        glAttachShader(pid_, id);
        shaders_.push_back(id);

        // tessellation control shader
        if (tcshader)
        {
            id = compile(tcshader, GL_TESS_CONTROL_SHADER);
            if (!id)
            {
                std::cerr << "Cannot compile tessellation control shader!\n";
                return false;
            }
            glAttachShader(pid_, id);
            shaders_.push_back(id);
        }

        // tessellation evaluation shader
        if (teshader)
        {
            id = compile(teshader, GL_TESS_EVALUATION_SHADER);
            if (!id)
            {
                std::cerr << "Cannot compile tessellation evaluation shader!\n";
                return false;
            }
            glAttachShader(pid_, id);
            shaders_.push_back(id);
        }

        // geometry shader
        if (gshader)
        {
            id = compile(gshader, GL_GEOMETRY_SHADER);
            if (!id)
            {
                std::cerr << "Cannot compile geometry shader!\n";
                return false;
            }
            glAttachShader(pid_, id);
            shaders_.push_back(id);
        }

        // link program
        if (!link())
        {
            std::cerr << "Cannot link program!\n";
            return false;
        }

        return true;
    }

//-----------------------------------------------------------------------------

    bool MyShader::link()
    {
        glLinkProgram(pid_);
        GLint status;
        glGetProgramiv(pid_, GL_LINK_STATUS, &status);
        if (status == GL_FALSE)
        {
            GLint length;
            glGetProgramiv(pid_, GL_INFO_LOG_LENGTH, &length);

            auto* info = new GLchar[length + 1];
            glGetProgramInfoLog(pid_, length, nullptr, info);
            std::cerr << "Shader: Cannot link program:\n" << info << std::endl;
            delete[] info;

            cleanup();

            return false;
        }

        return true;
    }

//-----------------------------------------------------------------------------

    bool MyShader::load(const char* filename, std::string& source)
    {
        std::ifstream ifs(filename);
        if (!ifs)
        {
            std::cerr << "Shader: Cannot open file \"" << filename << "\"\n";
            return false;
        }

        std::stringstream ss;
        ss << ifs.rdbuf();
        source = ss.str();

        ifs.close();
        return true;
    }

//-----------------------------------------------------------------------------

    GLint MyShader::compile(const char* source, GLenum type)
    {
        // create shader
        GLint id = glCreateShader(type);
        if (!id)
        {
            std::cerr << "Shader: Cannot create shader object\n";
            return 0;
        }

        // compile vertex shader
        glShaderSource(id, 1, &source, nullptr);
        glCompileShader(id);

        // check compile status
        GLint status;
        glGetShaderiv(id, GL_COMPILE_STATUS, &status);
        if (status == GL_FALSE)
        {
            GLint length;
            glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);

            auto* info = new GLchar[length + 1];
            glGetShaderInfoLog(id, length, nullptr, info);
            std::cerr << "Shader: Cannot compile shader\n" << info << std::endl;
            delete[] info;

            glDeleteShader(id);

            return 0;
        }

        return id;
    }

//-----------------------------------------------------------------------------

    GLint MyShader::load_and_compile(const char* filename, GLenum type)
    {
        std::string source;
        if (!load(filename, source))
        {
            std::cerr << "Shader: Cannot open file \"" << filename << "\"\n";
            return 0;
        }

        return compile(source.c_str(), type);
    }

//-----------------------------------------------------------------------------

    void MyShader::use()
    {
        if (pid_)
            glUseProgram(pid_);
    }

//-----------------------------------------------------------------------------

    void MyShader::disable()
    {
        glUseProgram(0);
    }

//-----------------------------------------------------------------------------

    void MyShader::bind_attribute(const char* name, GLuint index)
    {
        if (!pid_)
            return;
        glBindAttribLocation(pid_, index, name);
        link(); // have to re-link now!
    }

//-----------------------------------------------------------------------------

    void MyShader::set_uniform(const char* name, float value)
    {
        if (!pid_)
            return;
        int location = glGetUniformLocation(pid_, name);
        if (location == -1)
        {
            std::cerr << "Invalid uniform location for: " << name << std::endl;
            return;
        }
        glUniform1f(location, value);
    }

//-----------------------------------------------------------------------------

    void MyShader::set_uniform(const char* name, int value)
    {
        if (!pid_)
            return;
        int location = glGetUniformLocation(pid_, name);
        if (location == -1)
        {
            std::cerr << "Invalid uniform location for: " << name << std::endl;
            return;
        }
        glUniform1i(location, value);
    }

//-----------------------------------------------------------------------------

    void MyShader::set_uniform(const char* name, const vec3& vec)
    {
        if (!pid_)
            return;
        int location = glGetUniformLocation(pid_, name);
        if (location == -1)
        {
            std::cerr << "Invalid uniform location for: " << name << std::endl;
            return;
        };
        glUniform3f(location, vec[0], vec[1], vec[2]);
    }

//-----------------------------------------------------------------------------

    void MyShader::set_uniform(const char* name, const vec4& vec)
    {
        if (!pid_)
            return;
        int location = glGetUniformLocation(pid_, name);
        if (location == -1)
        {
            std::cerr << "Invalid uniform location for: " << name << std::endl;
            return;
        }
        glUniform4f(location, vec[0], vec[1], vec[2], vec[3]);
    }

//-----------------------------------------------------------------------------

    void MyShader::set_uniform(const char* name, const mat3& mat)
    {
        if (!pid_)
            return;
        int location = glGetUniformLocation(pid_, name);
        if (location == -1)
        {
            std::cerr << "Invalid uniform location for: " << name << std::endl;
            return;
        }
        glUniformMatrix3fv(location, 1, false, mat.data());
    }

//-----------------------------------------------------------------------------

    void MyShader::set_uniform(const char* name, const mat4& mat)
    {
        if (!pid_)
            return;
        int location = glGetUniformLocation(pid_, name);
        if (location == -1)
        {
            std::cerr << "Invalid uniform location for: " << name << std::endl;
            return;
        }
        glUniformMatrix4fv(location, 1, false, mat.data());
    }

//=============================================================================
} // namespace pmp
//=============================================================================
