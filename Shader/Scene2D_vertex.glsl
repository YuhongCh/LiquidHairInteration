
#version 450 core

layout(location = 0) in vec2 aPos;
uniform vec2 minCoord;
uniform vec2 maxCoord;
uniform float pointSize;

void main()
{
    vec2 point = (aPos - minCoord) / (maxCoord - minCoord) * 2.0 - 1.0;
    gl_Position = vec4(point, 0.0, 1.0);
    gl_PointSize = pointSize;
}
