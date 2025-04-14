
#version 450 core

layout(location = 0) in vec3 aPos;
uniform vec3 minCoord;
uniform vec3 maxCoord;
uniform float pointSize;

void main() {
    vec3 point = (aPos - minCoord) / (maxCoord - minCoord) * 2.0 - 1.0;
    gl_Position = vec4(point, 1.0);
    gl_PointSize = pointSize;
}
