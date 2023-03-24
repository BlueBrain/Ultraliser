#version 450 core

layout (location=0) in vec4 inPosition;

layout (location=0) out vec3 outLocalPos;
layout (location=1) out vec3 outCamSpacePos;

uniform mat4 modelView;
uniform mat4 modelViewProj;

void main( void )
{
  outLocalPos = inPosition.xyz;
  outCamSpacePos = (modelView * inPosition).xyz;
  gl_Position = modelViewProj * vec4(inPosition.xyz, 1.0);
}