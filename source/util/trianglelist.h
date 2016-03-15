#ifndef TRIANGLELIST_H
#define TRIANGLELIST_H

#include <cl/cl.h>

extern "C" {
#include "opencl/raytrace.h"
}

cl_float GetPointToTriangleSqLen(cl_float3 a, cl_float3 b, cl_float3 c, cl_float3 point);

// These classes are used to optimize the scene and not have to go through every triangle to check for collisions

// CameraTriangleList is the simple class: it lists all triangles that could be in the rays cast inside each camera pixel.
//  width x height unsorted lists of triangles are kept in (triangleListStart and triangleList)
//  the reason the lists are unsorted is because triangle proximity is not a partial order, even if we don't accept crossing triangles
class CameraTriangleList {
private:
	cl_uint* triangleListStart;
	cl_uint* triangleListEnd;
	cl_uint* triangleList;
	ptrdiff_t triangleListSize;

public:
	CameraTriangleList(cl_uint pixelCount, ptrdiff_t pixelTriangleCount);
	~CameraTriangleList();

	static CameraTriangleList* New(cl_uint2 cameraImageDimension, cl_float3 cameraEye, cl_float3 cameraEyeToTopLeftVector, cl_float3 cameraLeftToRightPixelSizeVector, cl_float3 cameraTopToBottomPixelSizeVector, cl_float cameraPixelSizeInv, cl_uint vertexCount, cl_uint triangleCount, cl_float3* vertex, cl_int3* triangleVertexIndex);

	cl_uint* GetTriangleListStart();
	cl_uint* GetTriangleListEnd();
	cl_uint* GetTriangleList();
	ptrdiff_t GetTriangleListSize();
};

// SceneTriangleList is a more complex class: is splits up the scene into AXES_DIVISION * AXES_DIVISION * AXES_DIVISION boxes (AXES_DIVISION per dimension)
//  the boxes are sized to spread as much as possible the triangles over each axis
// To check for collision with a triangle, one must rasterize a vector segment in a 3D grid where the boxes are of different sizes and then check every box
//  independently (see RayIntersectsTriangles() in raytrace_opencl.h and FillCube() in trianglelist.cpp)
//  AXES_DIVISION x AXES_DIVISION x AXES_DIVISION unsorted lists of triangles are kept in (triangleListStart and triangleList)
class SceneTriangleList {
public:
	enum {
		AXES_DIVISION = 256,
	};
private:
	cl_float3 boxMin[AXES_DIVISION + 1];
	cl_uint* triangleListStart;
	cl_uint* triangleList;

public:
	SceneTriangleList(cl_uint pixelCountInc, ptrdiff_t pixelTriangleCount);
	~SceneTriangleList();

	static SceneTriangleList* New(cl_uint vertexCount, cl_uint triangleCount, cl_float3* vertex, cl_int3* triangleVertexIndex);

	cl_float3* GetBoxMin();
	cl_uint* GetTriangleListStart();
	cl_uint* GetTriangleList();
};

#endif // TRIANGLELIST_H
