//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
#include "Summer.h"
#include "NvAppBase/NvInputTransformer.h"
#include "NvAssetLoader/NvAssetLoader.h"
#include "NvGLUtils/NvGLSLProgram.h"
#include "NvGLUtils/NvShapesGL.h"
#include "NvUI/NvTweakBar.h"
#include "NV/NvLogs.h"
#include "NV/NvVector.h"
#include "NV/NvMatrix.h"
using nv::vec2f;
using nv::vec3f;
using nv::vec4f;
using nv::matrix4f;
using nv::quaternionf;
using nv::normalize;
using nv::cross;
using nv::dot;
using nv::length;
#define PI NV_PI


static const float QUAD_COORDS[] = {
	-1.0f, -1.0f, -1.0f, 1.0f,
	1.0f, -1.0f, -1.0f, 1.0f,
	-1.0f,  1.0f, -1.0f, 1.0f,
	1.0f,  1.0f, -1.0f, 1.0f };

static const float TRIANGLE_COORDS[] = {
	-1.0f, -1.0f, -1.0f, 1.0f,        1.0f, 0, 0,
	1.0f, -1.0f, -1.0f, 1.0f,         1.0f, 0, 0,
	-1.0f,  1.0f, -1.0f, 1.0f,        1.0f, 0, 0, };

static const float gravity = 9.8f;
#define epsilon 0.01
#define epsilon_vec4f (vec4f(0.0001))
#define epsilon_vec3f (vec3f(0.0001))

class BaseShader : public NvGLSLProgram
{
public:
	BaseShader(const char *vertexProgramPath,
		const char *fragmentProgramPath) :
		positionAHandle(-1)
	{
		setSourceFromFiles(vertexProgramPath, fragmentProgramPath);
		positionAHandle = getAttribLocation("a_vPosition");
		colorAHandle = getAttribLocation("a_vColor");

		CHECK_GL_ERROR();
	}

	GLint positionAHandle;
	GLint colorAHandle;
};

class Triangle
{
public:
	vec3f vVelocity;
	vec3f vAngularVelocity;
	vec3f vAcceleration;
	vec3f vAngularAcceleration;

	vec3f vForce;
	float collideForce;
	vec4f local_vertices[3];
	vec4f vertices[3];   //CC Clock
	vec4f vertices_colors[3];   //CC Clock
	vec4f edges[3];
	vec4f edges_N[3];
	vec4f vCenterOfMass;  //cener of mass,  object origin
	quaternionf orientation;
	matrix4f mTranslateMatrix;
	matrix4f mRotationMatrix;
	matrix4f mScaleMatrix;
	matrix4f mModelMatrix;
	vec3f mTraslateVec;
	vec4f mRotateVec;
	vec4f mScaleVec;
	float thetaX, thetaY, thetaZ;
	float fArea;
	float fDensity;
	float fMass;
	float fInertia;
	float fInertiaInverse;
	float fSpeed;  //vector length of velocity;

public:
	int getNormalCount() { return ARRAYSIZE(edges_N); }
	int getVerticesCount() { return ARRAYSIZE(vertices); }
	int getEdgeCount() { return ARRAYSIZE(edges); }

	Triangle(vec4f translatePos, float angle, float scale) {
		init(translatePos, angle, scale);
	}
	void init(vec4f translatePos, float angle, float scale) {
		int stride = 7;
		for (int i = 0; i < ARRAYSIZE(TRIANGLE_COORDS); i = i + stride) {
			matrix4f rotation;
			nv::rotationZ(rotation, angle);
			vertices[i / stride].set_value(TRIANGLE_COORDS + i);
			local_vertices[i / stride].set_value(TRIANGLE_COORDS + i);
			vertices[i / stride] *= scale;
			vertices[i / stride].z = -1;
			vertices[i / stride] = rotation * vertices[i / stride];
			vertices[i / stride] += translatePos;
			vertices_colors[i / stride] = vec4f(.0f);
		}
		

		vCenterOfMass = (vertices[0] + vertices[1] + vertices[2]) / 3;
		edges[0] = (vertices[1] - vertices[0]);
		edges[1] = (vertices[2] - vertices[1]);
		edges[2] = (vertices[0] - vertices[2]);
		vec3f area_vec = cross((vec3f)edges[0], (vec3f)edges[1]);
		fArea = 0.5*nv::length(area_vec);
		fDensity = 1.0f;
		fMass = fArea * fDensity;
		matrix4f r_half_pi;
		nv::rotationZ(r_half_pi, -0.5f*PI);
		for (int i = 0; i < 3; i++)
			edges_N[i] = normalize(edges[i])*r_half_pi;

		vVelocity = 0.0;
		vVelocity.x = 0.0f;
		vVelocity.y = -0.0f;
		vAngularVelocity = 0.0;
		vAcceleration = 0.0;
		vAngularAcceleration = 0.0;
		fSpeed = length(vVelocity);
		thetaX = 0; thetaY = 0; thetaZ = 0;
		vForce = vec3f(0, 0, 1.0f);
		collideForce = 0.0f;
	}

	void translateAndRotate(vec3f translatePos) {
		mTraslateVec += translatePos;
		nv::translation(mTranslateMatrix, mTraslateVec.x, mTraslateVec.y, mTraslateVec.z);
		for (int i = 0; i < ARRAYSIZE(vertices); i++) {
			vertices[i] += vec4f(translatePos, 0);
		}

		vCenterOfMass = (vertices[0] + vertices[1] + vertices[2]) / 3;
		//rotation
		float theta = vAngularVelocity.z;
				
		for (int i = 0; i < ARRAYSIZE(vertices); i++) {
			vertices[i] -= vCenterOfMass;
			float x = vertices[i].x, y = vertices[i].y;
			vertices[i].x = cosf(theta)*x - sinf(theta)*y;
			vertices[i].y = sinf(theta)*x + cosf(theta)*y;
			vertices[i] += vCenterOfMass;
		}
		edges[0] = (vertices[1] - vertices[0]);
		edges[1] = (vertices[2] - vertices[1]);
		edges[2] = (vertices[1] - vertices[2]);
		vec3f area_vec = cross((vec3f)edges[0], (vec3f)edges[1]);
		fArea = 0.5*length(area_vec);
		fDensity = 1.0f;
		fMass = fArea * fDensity;
		matrix4f r_half_pi;
		nv::rotationZ(r_half_pi, -0.5f*PI);
		for (int i = 0; i < 3; i++)
			edges_N[i] = normalize(edges[i])*r_half_pi;
	}

	void rAngle(float theta) {

	}

	matrix4f rotateAtCenter(const float radians, vec4f center) const
	{
		nv::matrix4f tmp;
		nv::matrix4f transBackM = nv::translation(tmp,
			center.x,
			center.y,
			center.z);
		nv::matrix4f transOrigM = nv::translation(tmp,
			-center.x,
			-center.y,
			-center.z);

		nv::matrix4f rotation = nv::rotationZ(tmp, radians);

		return transBackM * rotation * transOrigM;
	}

	void draw(const unsigned int attrHandle) {
			glVertexAttribPointer(attrHandle, 4, GL_FLOAT, GL_FALSE, 4* sizeof(float),
				(float*)vertices);
			glEnableVertexAttribArray(attrHandle);
			glDrawArrays(GL_TRIANGLE_STRIP, 0, 3);
			glDisableVertexAttribArray(attrHandle);
	}

	void draw(const unsigned int posHandle, const unsigned int colorHandle) {
		glVertexAttribPointer(posHandle, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
			(float*)vertices);
		glVertexAttribPointer(colorHandle, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
			(float*)vertices_colors);
		glEnableVertexAttribArray(posHandle);
		glEnableVertexAttribArray(colorHandle);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 3);
		glDisableVertexAttribArray(colorHandle);
		glDisableVertexAttribArray(posHandle);
	}

	void processCollision(int vertexIndex, vec3f collisionNormal) {
		if (vertexIndex >= 0) {
			vertices_colors[vertexIndex].x = 1.0f;
			float projLen = dot(vVelocity, collisionNormal);
			vec3f half = vVelocity - collisionNormal * projLen;
			vVelocity = 2.0*half - vVelocity;
			vAngularVelocity = -vAngularVelocity;
			collideForce = 3.0f;
		}
	}

	void processCollision(int edgeIndex, int edgeT) {
		vertices_colors[edgeIndex].x = 1.0f;
		vertices_colors[(edgeIndex + 1) % 3].y = 1.0f;
		float projLen = dot(vVelocity, (vec3f)edges_N[edgeIndex]);
		vec3f half = vVelocity - (vec3f)edges_N[edgeIndex] * projLen;
		vVelocity = 2.0*half - vVelocity;
	}

	void processNoCollision() {
		for(int i = 0; i < ARRAYSIZE(vertices_colors);i++)
			vertices_colors[i] = vec4f(0.0f);
	}

	float distanceFromEdge(vec3f point, int segmentIndex, float& tCloset) {
		if (segmentIndex >= 3) {
			return -1;
		}
		vec3f startPoint(vertices[segmentIndex]);
		vec3f diff = point - startPoint;
		vec3f edge = vec3f(edges[segmentIndex]);
		vec3f direction = normalize(edge);
		tCloset = dot(diff, edge) / length(edge);
		if (tCloset > length(edge) || tCloset < 0) {
			return -1;
		}
		diff -= tCloset*direction;
		float distance = length(diff);
		return distance;
	}

	bool isPointInside(vec4f point) {
		for (int i = 0; i < 3; i++) {
			vec3f v1 = point - vertices[i];
			float determin = dot((vec3f)edges_N[i], v1);
			if (determin < 0.0f)
				return false;

		}
		return true;
		
	}

	void updateDynamics(float deltaTime) {
		vec3f travel = vVelocity * deltaTime;
		
			
		translateAndRotate((travel));
	}
};

class rectangle
{};


class circle
{};

float LINE_SEGMENT_COORDS[] = {
	-1.0f, 0.0f, -1.0f, 1.0f,
	1.0f, 0.0f, -1.0f, 1.0f,
};
class LineSegment {
public:
	vec4f vertices[2];
	vec4f local_vertices[2];
	vec4f vertices_colors[2];
	vec4f centerOfMass;
	vec4f segment;
	vec4f direction;
	vec4f normal;
	vec4f mTranslateVec;
	vec4f mScale;
	float angleZ;

	LineSegment(vec4f translatePos, float angle, float scale) {
		init(translatePos, angle, scale);
	}

	void init(vec4f translatePos, float angle, float scale) {
		mTranslateVec = translatePos;
		for (int i = 0; i < ARRAYSIZE(LINE_SEGMENT_COORDS); i = i + 4) {
			matrix4f rotation;
			vertices[i / 4].set_value(LINE_SEGMENT_COORDS + i);
			local_vertices[i / 4].set_value(LINE_SEGMENT_COORDS + i);
			nv::rotationZ(rotation, angle);
			vertices[i / 4] *= scale;
			vertices[i / 4].z = -1;
			vertices[i / 4] = rotation * vertices[i / 4];
			vertices[i / 4] += translatePos;
			vertices_colors[i / 4] = vec4f(.0f, 1.0f, 0.0f,1.0f);
		}
		centerOfMass = (vertices[0] + vertices[1]) / 2;
		segment = direction = vertices[1] - vertices[0];
		direction = normalize(direction);
		matrix4f r_half_pi;
		nv::rotationZ(r_half_pi, -0.5f*PI);
		normal = r_half_pi *  direction;
	}

	float distanceFrom(vec3f point, float& tCloset) {
		vec3f diff = point - (vec3f)vertices[0];
		float segmentLength = length((vec3f)segment);
		tCloset = dot(diff, (vec3f)segment) / segmentLength;
		
		if (tCloset > segmentLength || tCloset < 0) {
			return -1;
		}
		diff -= tCloset*(vec3f)direction;
		float distance = length(diff);
		LOGI("segment distance = %f", distance);
		return distance;
	}

	void processCollision() {
		vertices_colors[0] = vec4f(1.0f, 0, 0, 0);
		vertices_colors[1] = vec4f(1.0f, 0, 0, 0);
	}

	void processNoCollision() {
		vertices_colors[0] = vec4f(0.0f, 1, 0, 1);
		vertices_colors[1] = vec4f(0.0f, 1, 0, 1);
	}

	void draw(const unsigned int posHandle, const unsigned int colorHandle) {
		glVertexAttribPointer(posHandle, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
			(float*)vertices);
		glVertexAttribPointer(colorHandle, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float),
			(float*)vertices_colors);
		glEnableVertexAttribArray(posHandle);
		glEnableVertexAttribArray(colorHandle);
		glDrawArrays(GL_LINES, 0, 2);
		glDisableVertexAttribArray(colorHandle);
		glDisableVertexAttribArray(posHandle);
	}
};

class CollisionDetection {
public:
	bool testIntersection(Triangle* T0, Triangle* T1) {
		Triangle* triangles[2];
		triangles[0] = T0;
		triangles[1] = T1;
		for (int t = 0; t < 2; t++) {
			Triangle* curTriangle = triangles[t];
			for (int i = 0; i < curTriangle->getNormalCount(); i++) {
				float lambda0_min = 99999.0f, lambda0_max = -lambda0_min;
				float lambda1_min = 99999.0f, lambda1_max = -lambda1_min;

				vec4f dir = curTriangle->edges_N[i];
				for (int j = 0; j < T0->getVerticesCount(); j++) {
					float lambda = dot(dir, T0->vertices[j]);
					lambda0_min = lambda0_min > lambda ? lambda : lambda0_min;
					lambda0_max = lambda0_max < lambda ? lambda : lambda0_max;


				}
				for (int k = 0; k < T1->getVerticesCount(); k++) {
					float lambda = dot(dir, T1->vertices[k]);
					lambda1_min = lambda1_min > lambda ? lambda : lambda1_min;
					lambda1_max = lambda1_max < lambda ? lambda : lambda1_max;
				}
				if (lambda0_min > lambda1_max) {
					return false;
				}
				if (lambda1_min > lambda0_max) {
					return false;
				}
			}
		}

		return true;
	}

	bool findIntersection(Triangle* T0, Triangle* T1) {
		//edge-vertex
		bool result = false;
		for (int i = 0; i < T0->getVerticesCount(); i++) {
			vec3f curVertex = T0->vertices[i];
			for (int j = 0; j < T1->getEdgeCount(); j++) {
				float closestT;
				float distance = T1->distanceFromEdge(curVertex, j, closestT);
				
				if ((distance >= 0 && distance < epsilon) ||
					(distance >= -epsilon && distance <= 0) 
					) {
					T0->processCollision(i, (vec3f)T1->edges_N[j]);
					T1->processCollision(j, closestT);
					result = true;
					return result;
				}

			}

			bool isPtInside = T1->isPointInside(vec4f(curVertex, 1));
			if (isPtInside) {
				T0->vVelocity = vec3f(-T0->vVelocity.x, -T0->vVelocity.y, 0);
				T0->vAngularVelocity = -T0->vAngularVelocity;
				result = true;
				return result;
			}
		}

		for (int i = 0; i < T1->getVerticesCount(); i++) {
			vec3f curVertex = T1->vertices[i];
			for (int j = 0; j < T0->getEdgeCount(); j++) {
				float closestT;
				float distance = T0->distanceFromEdge(curVertex, j, closestT);
				
				if ((distance >= 0 && distance < epsilon) ||
					(distance >= -epsilon && distance <= 0)
					) {
					T1->processCollision(i, (vec3f)T0->edges_N[j]);
					T0->processCollision(j, closestT);
					result = true;
					return result;
				}

			}

			bool isPtInside = T0->isPointInside(vec4f(curVertex, 1));
			if (isPtInside) {
				T1->vVelocity = vec3f(-T1->vVelocity.x, -T1->vVelocity.y, 0);
				T1->vAngularVelocity = -T1->vAngularVelocity;
				//T0->vVelocity = -2.0f*-T0->vVelocity;
				result = true;
				return result;
			}
		}
		//vertex-vertex
		//edge-edge

		//process no collison
		//T0->processNoCollision();
		//T1->processNoCollision();
		return result;
	}

	bool findIntersection(LineSegment* lineSeg, Triangle* tri) {
		for (int i = 0; i < tri->getVerticesCount(); i++) {
			vec3f curVertex = tri->vertices[i];
			float closestT;
			float distance = lineSeg->distanceFrom(curVertex, closestT);
			if (distance >= 0 && distance < epsilon) {
				tri->processCollision(i, (vec3f)lineSeg->normal);
				lineSeg->processCollision();
				return true;
			}
		}

		//process no collison
		//tri->processNoCollision();
		//lineSeg->processNoCollision();
		return false;
	}

	bool testIntersection(LineSegment* lineSeg, Triangle* tri) {
		return true;
	}
};


Summer::Summer()
{
    m_transformer->setTranslationVec(nv::vec3f(0.0f, 0.0f, -2.2f));
    m_transformer->setRotationVec(nv::vec3f(0.0f, 0.0f, 0.0f));

    forceLinkHack();
}

Summer::~Summer()
{
    LOGI("Basic: destroyed\n");
}

void Summer::configurationCallback(NvGLConfiguration& config)
{ 
    config.depthBits = 24; 
    config.stencilBits = 0; 
    config.apiVer = NvGLAPIVersionGL4();
}

void Summer::initRendering(void) {
	NV_APP_BASE_SHARED_INIT();
	NvAssetLoaderAddSearchPath("summer/");

	mProgram = new BaseShader("shaders/plain.vert", "shaders/plain.frag");//NvGLSLProgram::createFromFiles("shaders/plain.vert", "shaders/plain.frag");
	mTriangleProgram = new BaseShader("shaders/triangle.vert", "shaders/triangle.frag");
	mSegmentProgram = new BaseShader("shaders/segment.vert", "shaders/segment.frag");
	initScene();
}

#define SCALE_WALL_V 1.3f
#define SCALE_WALL_H 2.2f
#define TRANSLATE_WALL_H 2.2f
#define TRANSLATE_WALL_V 1.3f
void Summer::initScene() {
	mTopWall = new LineSegment(vec4f(0, TRANSLATE_WALL_V, 0.0f, 0), 0, SCALE_WALL_H);
	mBottomWall = new LineSegment(vec4f(0, -TRANSLATE_WALL_V, 0.0f, 0), 0, SCALE_WALL_H);
	mLeftWall = new LineSegment(vec4f(-TRANSLATE_WALL_H, 0.0f, 0.0f, 0), 0.5f*PI, SCALE_WALL_V);
	mRightWall = new LineSegment(vec4f(TRANSLATE_WALL_H, .0f, 0.0f, 0), 0.5f*PI, SCALE_WALL_V);
	for (int i = 0; i < ARRAYSIZE(mTriangles); i++) {
	
		mTriangles[i] = new Triangle(vec4f(0.01*(i + 1)*cosf(i*1.0f), 0.01*(i + 1)*sinf(i*1.0f), 0.0f, 0), 0.15*PI*i, i>=1 ? 0.10f : 0.2f - 0.03f*i);
		mTriangles[i]->vVelocity = 2.0*vec3f(0.001*(i + 1)*cosf(i*.1f), 0.001*(i + 1)*sinf(i*.1f), 0);
		mTriangles[i]->vAngularVelocity =  vec3f(0, 0, i >= 5 ? 0.003f*i : 0.01*(i + 1));
	}
	

	mCollisionDetection = new CollisionDetection();
}

void Summer::shutdownRendering(void) {

    // destroy other resources here
}

void Summer::initUI(void) {
    if (mTweakBar) {
//        NvTweakVarBase *var;

        mTweakBar->syncValues();
    }
}


void Summer::reshape(int32_t width, int32_t height)
{
    glViewport( 0, 0, (GLint) width, (GLint) height );
}


void Summer::draw(void)
{
	CHECK_GL_ERROR();
	glClear(GL_COLOR_BUFFER_BIT);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	CHECK_GL_ERROR();


	int triangleCount = ARRAYSIZE(mTriangles);
	const int wallCount = 4;
	LineSegment* walls[wallCount] = { mTopWall,
								mBottomWall,
								mLeftWall,
								mRightWall
	};

	
		
	for (int j = 0; j < wallCount; j++) {
		walls[j]->processNoCollision();
	}

	for (int i = 0; i < triangleCount; i++) {
		mTriangles[i]->processNoCollision();
		//mTriangles[i]->updateDynamics(1.0f);
	}
	

	for (int i = 0; i < triangleCount; i++)
		for (int j = 0; j < wallCount; j++) {
			if (mCollisionDetection->testIntersection(walls[j], mTriangles[i])) {
				if (mCollisionDetection->findIntersection(walls[j], mTriangles[i])) {

				}
			}
		}


	for (int i = 0; i < triangleCount; i++) {
		
		for (int j = i + 1;j < triangleCount; j++) {
			if (i !=j && mCollisionDetection->testIntersection(mTriangles[j], mTriangles[i])) {
				mCollisionDetection->findIntersection(mTriangles[j], mTriangles[i]);
			}
		}
	}

	for (int i = 0; i < triangleCount; i++) {
		//mTriangles[i]->processNoCollision();
		mTriangles[i]->updateDynamics(1.0f);
	}
	
	
    nv::matrix4f projection_matrix;
    nv::perspective(projection_matrix, 3.14f * 0.25f, m_width/(float)m_height, 0.1f, 30.0f);
    nv::matrix4f camera_matrix = projection_matrix * m_transformer->getModelViewMat();
    CHECK_GL_ERROR();

	mTriangleProgram->enable();
	mTriangleProgram->setUniformMatrix4fv("uMVP", camera_matrix._array, 1, GL_FALSE);
	for(int i = 0; i < ARRAYSIZE(mTriangles);i++)
		mTriangles[i]->draw(mTriangleProgram->getAttribLocation("aPosition"),
			mTriangleProgram->getAttribLocation("aColor"));
	CHECK_GL_ERROR();
	mTopWall->draw(mTriangleProgram->getAttribLocation("aPosition"),
		mTriangleProgram->getAttribLocation("aColor"));
	mBottomWall->draw(mTriangleProgram->getAttribLocation("aPosition"),
		mTriangleProgram->getAttribLocation("aColor"));
	mLeftWall->draw(mTriangleProgram->getAttribLocation("aPosition"),
		mTriangleProgram->getAttribLocation("aColor"));
	mRightWall->draw(mTriangleProgram->getAttribLocation("aPosition"),
		mTriangleProgram->getAttribLocation("aColor"));
	CHECK_GL_ERROR();
	mTriangleProgram->disable();
}


NvAppBase* NvAppFactory() {
    return new Summer();
}


