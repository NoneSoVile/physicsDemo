//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
#version 100
attribute highp vec4 aPosition;
attribute highp vec3 aColor;

uniform mediump mat4 uMVP;

varying lowp vec4 vColor;

void main(void)
{
    gl_Position = uMVP * vec4(aPosition.x, aPosition.y, aPosition.z, 1.0);
    vColor = vec4(pow(aColor.x, 1.0), pow(aColor.y, 1.0), aColor.z,1.0);
}
