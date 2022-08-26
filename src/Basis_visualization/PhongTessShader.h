// clang-format off

static const char* phong_tess_vshader =
R"glsl(
#version 410 core
layout (location=0) in vec4 v_position;
layout (location=1) in vec3 v_normal;
layout (location=2) in vec2 v_tex;
layout (location=3) in float v_coeffs;

out vec3 v2tc_normal;
out vec2 v2tc_tex;
out vec3 v2tc_view;
out float v2tc_coeffs;

uniform mat4 modelview_matrix;
uniform mat3 normal_matrix;
uniform float point_size;
uniform bool show_texture_layout;

void main()
{
  v2tc_normal  = normal_matrix * v_normal;
  v2tc_tex     = v_tex;
  vec4 pos     = show_texture_layout ? vec4(v_tex, 0.0, 1.0) : v_position;
  v2tc_view    = -(modelview_matrix * pos).xyz;
  v2tc_coeffs  = v_coeffs;
  gl_PointSize = point_size;
  gl_Position  = modelview_matrix * pos;
}
)glsl";

static const char* phong_tess_tcshader =
    R"glsl(
#version 410 core

layout(vertices = 6) out;

in vec3 v2tc_normal[];
in vec2 v2tc_tex[];
in vec3 v2tc_view[];
in float v2tc_coeffs[];

out vec3 tc2te_normal[];
out vec2 tc2te_tex[];
out vec3 tc2te_view[];
out float tc2te_coeffs[];

uniform float tess_level;

void main()
{
    gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
    gl_out[gl_InvocationID].gl_PointSize = gl_in[gl_InvocationID].gl_PointSize;

    tc2te_normal[gl_InvocationID] = v2tc_normal[gl_InvocationID];
    tc2te_tex[gl_InvocationID] = v2tc_tex[gl_InvocationID];
    tc2te_view[gl_InvocationID] = v2tc_view[gl_InvocationID];
    tc2te_coeffs[gl_InvocationID] = v2tc_coeffs[gl_InvocationID];

    if(gl_InvocationID == 0) {
        gl_TessLevelInner[0] = tess_level;
        gl_TessLevelOuter[0] = tess_level;
        gl_TessLevelOuter[1] = tess_level;
        gl_TessLevelOuter[2] = tess_level;
    }
}
)glsl";

static const char* phong_tess_teshader =
    R"glsl(
#version 410 core

layout(triangles) in;

in vec3 tc2te_normal[];
in vec2 tc2te_tex[];
in vec3 tc2te_view[];
in float tc2te_coeffs[];

out vec3 te2g_normal;
out vec2 te2g_tex;
out vec3 te2g_view;

uniform mat4 projection_matrix;
uniform bool use_basis_function_texture;
uniform bool use_basis_function_elevation;
uniform float basis_function_elevation;

//#define QUADRATIC_INTERPOLATION(type) \
//    type quadratic_interpolation(float u, float v, type a, type b, type c, type d, type e, type f) { \
//        return tc2te_coeffs[0] * (1-u-v)*(1-2*u-2*v) * a \
//            +  tc2te_coeffs[1] * u*(2*u-1)           * b \
//            +  tc2te_coeffs[2] * v*(2*v-1)           * c \
//            +  tc2te_coeffs[3] * 4*u*(1-u-v)         * d \
//            +  tc2te_coeffs[4] * 4*u*v               * e \
//            +  tc2te_coeffs[5] * 4*v*(1-u-v)         * f;\
//    }

#define QUADRATIC_INTERPOLATION(type) \
    type quadratic_interpolation(float u, float v, type a, type b, type c, type d, type e, type f) { \
        return tc2te_coeffs[0] * (u+v-1)*(2*u+2*v-1) * a \
            +  tc2te_coeffs[1] * u*(2*u-1)           * b \
            +  tc2te_coeffs[2] * v*(2*v-1)           * c \
            +  tc2te_coeffs[3] * -4*u*(u+v-1)        * d \
            +  tc2te_coeffs[4] * 4*u*v               * e \
            +  tc2te_coeffs[5] * -4*v*(u+v-1)        * f;\
    }

QUADRATIC_INTERPOLATION(vec2)
QUADRATIC_INTERPOLATION(vec4)

void main()
{
    float u = gl_TessCoord.x;
    float v = gl_TessCoord.y;
    float w = gl_TessCoord.z;

    //vec3 n = quadratic_interpolation(v, w, tc2te_normal[0], tc2te_normal[1], tc2te_normal[2], tc2te_normal[3], tc2te_normal[4], tc2te_normal[5]);
    vec3 n = u * tc2te_normal[0] + v * tc2te_normal[1] + w * tc2te_normal[2];
    te2g_normal = n;

    vec2 tex;
    if (!use_basis_function_texture) {
        //tex = quadratic_interpolation(v, w, tc2te_tex[0], tc2te_tex[1], tc2te_tex[2], tc2te_tex[3], tc2te_tex[4], tc2te_tex[5]);
        tex = u * tc2te_tex[0] + v * tc2te_tex[1] + w * tc2te_tex[2];
    } else {
        tex = (tc2te_coeffs[0] * (1-v-w)*(1-2*v-2*w)
                +  tc2te_coeffs[1] * v*(2*v-1)
                +  tc2te_coeffs[2] * w*(2*w-1)
                +  tc2te_coeffs[3] * 4*v*(1-v-w)
                +  tc2te_coeffs[4] * 4*v*w
                +  tc2te_coeffs[5] * 4*w*(1-v-w)) * vec2(1.0, 0.0);

        tex[0] += 1;

        tex[0] /= 2.0;
    }
    te2g_tex = tex;

    //vec3 view = quadratic_interpolation(v, w, tc2te_view[0], tc2te_view[1], tc2te_view[2], tc2te_view[3], tc2te_view[4], tc2te_view[5]);
    vec3 view = u * tc2te_view[0] + v * tc2te_view[1] + w * tc2te_view[2];
    te2g_view = view;

    //float pointSize = quadratic_interpolation(v, w, gl_in[0].gl_PointSize, gl_in[1].gl_PointSize, gl_in[2].gl_PointSize, gl_in[3].gl_PointSize, gl_in[4].gl_PointSize, gl_in[5].gl_PointSize);
    float pointSize = u * gl_in[0].gl_PointSize + v * gl_in[1].gl_PointSize + w * gl_in[2].gl_PointSize;
    gl_PointSize = pointSize;

    //vec4 p = quadratic_interpolation(v, w, gl_in[0].gl_Position, gl_in[1].gl_Position, gl_in[2].gl_Position, gl_in[3].gl_Position, gl_in[4].gl_Position, gl_in[5].gl_Position);
    vec4 p = u * gl_in[0].gl_Position + v * gl_in[1].gl_Position + w * gl_in[2].gl_Position;
    if(use_basis_function_elevation) {
        p = p + (tc2te_coeffs[0] * (1-v-w)*(1-2*v-2*w)
                +  tc2te_coeffs[1] * v*(2*v-1)
                +  tc2te_coeffs[2] * w*(2*w-1)
                +  tc2te_coeffs[3] * 4*v*(1-v-w)
                +  tc2te_coeffs[4] * 4*v*w
                +  tc2te_coeffs[5] * 4*w*(1-v-w)) * basis_function_elevation * vec4(n, 1.0);
    }
    gl_Position = projection_matrix * p;
}
)glsl";

static const char* phong_tess_gshader =
    R"glsl(
#version 410 core

layout(triangles) in;
layout(triangle_strip, max_vertices = 3) out;

in vec3 te2g_normal[];
in vec2 te2g_tex[];
in vec3 te2g_view[];

out vec3 g2f_normal;
out vec2 g2f_tex;
out vec3 g2f_view;
out vec3 g2f_bary;

uniform bool draw_edges = false;

void main() {

    g2f_bary = vec3(1.0,1.0,1.0);

    if (draw_edges)
        g2f_bary    = vec3(1.0, 0.0, 0.0);
    g2f_normal  = te2g_normal[0];
    g2f_tex     = te2g_tex[0];
    g2f_view    = te2g_view[0];
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    if (draw_edges)
        g2f_bary    = vec3(0.0, 1.0, 0.0);
    g2f_normal  = te2g_normal[1];
    g2f_tex     = te2g_tex[1];
    g2f_view    = te2g_view[1];
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    if (draw_edges)
        g2f_bary    = vec3(0.0, 0.0, 1.0);
    g2f_normal  = te2g_normal[2];
    g2f_tex     = te2g_tex[2];
    g2f_view    = te2g_view[2];
    gl_Position = gl_in[2].gl_Position;
    EmitVertex();
}
)glsl";

static const char* phong_tess_fshader =
    R"glsl(
#version 410 core

precision mediump float;

in vec3 g2f_normal;
in vec2 g2f_tex;
in vec3 g2f_view;
in vec3 g2f_bary;

uniform bool   use_lighting;
uniform bool   use_texture;
uniform bool   use_srgb;
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

float edgeFactor()
{
    vec3 dx = dFdx(g2f_bary);
    vec3 dy = dFdy(g2f_bary);
    vec3 d  = sqrt(dx*dx + dy*dy);
    vec3 a3 = smoothstep(vec3(0.0), d*1.5, g2f_bary);
    return min(min(a3.x, a3.y), a3.z);
}

void main()
{
    vec3 color = gl_FrontFacing ? front_color : back_color;

    vec3 rgb;

    if (use_lighting)
    {
        vec3 L1 = normalize(light1);
        vec3 L2 = normalize(light2);
        vec3 V  = normalize(g2f_view);
        vec3 N  = gl_FrontFacing ? normalize(g2f_normal) : -normalize(g2f_normal);
        vec3 R;
        float NL, RV;

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

    if (use_texture) rgb *= texture(mytexture, g2f_tex).xyz;
    if (use_srgb)    rgb  = pow(clamp(rgb, 0.0, 1.0), vec3(0.45));

    rgb *= edgeFactor();

    f_color = vec4(rgb, alpha);
}
)glsl";