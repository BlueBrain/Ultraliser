#version 450 core

layout (location=0) in vec4 inPosition;

void main( void )
{
  gl_Position = inPosition;
}