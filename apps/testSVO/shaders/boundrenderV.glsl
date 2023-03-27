#version 450 core

layout (location=0) in vec4 inPosition;

layout (std140, binding=0) buffer NodeData
{ 
  vec4 nodes[];
};

void main( void )
{
  vec3 scaled = inPosition.xyz * nodes[gl_InstanceID].w;
  vec3 translated = scaled + nodes[gl_InstanceID].xyz;
  gl_Position = vec4(translated, 1);
}