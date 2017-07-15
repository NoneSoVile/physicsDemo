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

class Summer : public NvSampleAppGL
{
public:
    Summer();
    ~Summer();
    
    void initRendering(void);
    void shutdownRendering(void);
    void initUI(void);
    void draw(void);
    void reshape(int32_t width, int32_t height);

    void configurationCallback(NvGLConfiguration& config);

private:
    NvGLSLProgram* mProgram;
};
