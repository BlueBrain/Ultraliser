#version 450 core

layout(triangles) in;
layout(line_strip, max_vertices=3) out;

uniform mat4 modelViewProj;

#define EPSILON_TEST 0.0001

void main()
{
  const vec4 a = gl_in[0].gl_Position;
  const vec4 b = gl_in[1].gl_Position;
  const vec4 c = gl_in[2].gl_Position;

  const float ablen = length(a.xyz - b.xyz);
  const float aclen = length(a.xyz - c.xyz);

  vec4 pivot, first, last;

  if(abs(ablen - aclen) > EPSILON_TEST)
  {
    pivot = a;
    first = ablen < aclen? b : c;
    last = ablen < aclen? c : b;
  }
  else
  {
    pivot = ablen < aclen? b : c;
    first = a;
    last = c;
  }

  gl_Position = modelViewProj * pivot;
  EmitVertex();
  gl_Position = modelViewProj * first;
  EmitVertex();
  gl_Position = modelViewProj * last;
  EmitVertex();
  
  EndPrimitive();
  
}