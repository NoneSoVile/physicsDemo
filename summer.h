//----------------------------------------------------------------------------------

//
//----------------------------------------------------------------------------------
#include "NvAppBase/gl/NvSampleAppGL.h"

#include "KHR/khrplatform.h"
#include "NvGamepad/NvGamepad.h"
#include "NV/NvMath.h"


class NvStopWatch;
class NvGLSLProgram;
class NvFramerateCounter;

class LineSegment;
class Triangle;
class CollisionDetection;
#define TRIANGLE_NUMBER 30
class Summer : public NvSampleAppGL
{
public:
    Summer();
    ~Summer();
    
    void initRendering(void);
	void initScene();
    void shutdownRendering(void);
    void initUI(void);
    void draw(void);
    void reshape(int32_t width, int32_t height);

    void configurationCallback(NvGLConfiguration& config);

private:
    NvGLSLProgram* mProgram;
	NvGLSLProgram* mTriangleProgram;
	NvGLSLProgram* mSegmentProgram;
	LineSegment* mTopWall, *mBottomWall, *mLeftWall, *mRightWall;
	Triangle* mTriangles[TRIANGLE_NUMBER];
	CollisionDetection* mCollisionDetection;
};
