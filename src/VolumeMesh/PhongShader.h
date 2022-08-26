//=============================================================================
// Copyright (C) 2011-2020 The pmp-library developers
//
// This file is part of the Polygon Mesh Processing Library.
// Distributed under a MIT-style license, see LICENSE.txt for details.
//
// SPDX-License-Identifier: MIT-with-employer-disclaimer
//=============================================================================

// clang-format off

static const char* phong_vshader = 
#ifndef __EMSCRIPTEN__
    "#version 330"
#else
    "#version 300 es"
#endif
R"glsl(
layout (location=0) in vec4 v_position;
layout (location=1) in vec3 v_f_normal;
layout (location=2) in vec3 v_embedding;
layout (location=3) in vec2 v_texture;

out vec3 v2f_normal;
out vec2 v2f_tex;
out vec3 v2f_view;

uniform mat4   modelview_projection_matrix;
uniform mat4   modelview_matrix;
uniform mat3   normal_matrix;
uniform float  point_size;
uniform bool   embedded;
uniform float  embedding_size;

void main()
{
    v2f_normal   = normal_matrix * v_f_normal;
    v2f_tex      = v_texture;
    vec4 pos     = v_position;
    if (embedded) {
        pos = pos + embedding_size * vec4(v_embedding, 0);
    }
    v2f_view     = -(modelview_matrix * pos).xyz;
    gl_PointSize = point_size;
    gl_Position  = modelview_projection_matrix * pos;
}
)glsl";


static const char* phong_fshader = 
#ifndef __EMSCRIPTEN__
    "#version 330"
#else
    "#version 300 es"
#endif
R"glsl(
precision mediump float;

in vec3  v2f_normal;
in vec2  v2f_tex;
in vec3  v2f_view;

uniform bool   use_lighting;
uniform bool   use_texture;
uniform vec3   front_color;
uniform vec3   back_color;
uniform float  ambient;
uniform float  diffuse;
uniform float  specular;
uniform float  shininess;
uniform float  alpha;
uniform vec3   light1;
uniform vec3   light2;

uniform sampler2D mytexture;

out vec4 f_color;

void main()
{
    vec3 color = front_color;
    vec3 rgb;

    if (use_lighting)
    {
        vec3 L1 = normalize(light1);
        vec3 L2 = normalize(light2);
        vec3 V  = normalize(v2f_view);
        vec3 N  = normalize(v2f_normal);
        vec3 R;
        float NL, RV;

        // front-facing or back-facing?
        // (gl_FrontFacing does not work with Apple's shitty OpenGL drivers)
        if (!gl_FrontFacing)
        {
            N = -N;
            color = back_color;
        }

        rgb = ambient * 0.1 * color;

        NL = dot(N, L1);
        if (NL > 0.0)
        {
            rgb += diffuse * NL * color;
            R  = normalize(-reflect(L1, N));
            RV = dot(R, V);
            if (RV > 0.0) 
            {
                rgb += vec3( specular * pow(RV, shininess) );
            }
        }

        NL = dot(N, L2);
        if (NL > 0.0)
        {
            rgb += diffuse * NL * color;
            R  = normalize(-reflect(L2, N));
            RV = dot(R, V);
            if (RV > 0.0) 
            {
                rgb += vec3( specular * pow(RV, shininess) );
            }
        }
    }

    // do not use lighting
    else
    {
        rgb = color;
    }

    if (use_texture) rgb *= texture(mytexture, v2f_tex).xyz;

    f_color = vec4(rgb, alpha);
}
)glsl";


//=============================================================================
// clang-format on
//=============================================================================
