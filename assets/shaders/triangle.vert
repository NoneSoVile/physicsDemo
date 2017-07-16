//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
#version 100
attribute highp vec2 aPosition;

uniform mediump mat4 uMVP;

varying lowp vec4 vColor;

void main(void)
{
    gl_Position = uMVP * vec4(aPosition.x, aPosition.y, 0.0, 1.0);
    vColor = vec4(aPosition.x * 0.5 + 0.5, aPosition.y * 0.5 + 0.5, aPosition.x * 0.5 + 0.5, 1.0);
}
