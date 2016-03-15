#ifndef RAYTRACE_H
#define RAYTRACE_H

#include <cl/cl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// OPENCL-TO-C EXPRESSIONS
#define __global
#define __private
#define __kernel

// OPENCL-TO-C TYPES
#define bool   cl_bool
#define double cl_double
#define float  cl_float
#define float2 cl_float2
#define float3 cl_float3
#define int    cl_int
#define int2   cl_int2
#define int3   cl_int3
#define uchar3 cl_uchar3
#define uint   cl_uint
#define uint2  cl_uint2
#define ulong  cl_ulong
#define ushort cl_ushort

#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))

#define M_PI 3.14159265f

#include "raytrace_opencl.h"

cl_float dot(cl_float3 a, cl_float3 b);
cl_float3 cross(cl_float3 a, cl_float3 b);
cl_float3 normalize(cl_float3 v);
cl_float3 vector(cl_float3 a, cl_float3 b);
cl_float bindf(cl_float value, cl_float a, cl_float b);
cl_float GetPointToLineSqLen(cl_float3 origin, cl_float3 destination, cl_float3 point);
cl_bool RayIntersectsTriangle(cl_float3 origin, cl_float3 ray, cl_float minDistance, cl_float maxDistance, cl_float3 a, cl_float3 b, cl_float3 c, cl_float* outRayMult, cl_float* outABL, cl_float* outACL);
cl_int3 GetBoxAddress(cl_int axesDivCount, cl_float3* boxMin, cl_float3 position);

void InitOpenCL();
void ResetComputationType();
cl_bool GetIsComputationTypeUpdated();
size_t GetComputationTypeCount();
cl_bool GetComputationTypeName(size_t id, size_t strLen, cl_char* str);

cl_float GetProgress();
void SetProgress(cl_float p);
clock_t GetStartTime();
clock_t GetEndTime();
void ResetTime();

cl_bool
RaytraceAll(cl_uint computationType,

			cl_uint2 cameraImageDimension,
			cl_float3 cameraEye,
			cl_float3 cameraEyeToTopLeftVector,
			cl_float3 cameraLeftToRightPixelSizeVector,
			cl_float3 cameraTopToBottomPixelSizeVector,
			cl_float cameraPixelSizeInv,

			cl_uint* cameraPixelTriangleListStart,
			cl_uint* cameraPixelTriangleListEnd,
			cl_uint* cameraPixelTriangleList,
			ptrdiff_t cameraPixelTriangleListSize,

			cl_uint sampleCount,

			cl_uint vertexCount,
			cl_float3* vertex,

			cl_uint triangleCount,
			cl_int3* triangleVertexIndex,
			cl_int* triangleMaterialId,
			cl_float2* triangleUv, // 3x triangleCount
			cl_float3* triangleNormal, // 3x triangleCount

			cl_int axesDivCount,
			cl_float3* sceneBoxMin,
			cl_uint* scenePixelTriangleListStart,
			cl_uint* scenePixelTriangleList,

			cl_uint materialCount,
			cl_uint2* materialImageSize, // materialCount * _MATERIAL_CHANNEL_COUNT
			cl_int* materialImageStart, // materialCount * _MATERIAL_CHANNEL_COUNT

			cl_uint texturesSize,
			cl_uchar3* textures,

			cl_uint lightCount,
			cl_int* lightType,
			cl_float3* lightPosition,
			cl_float3* lightDirection,
			cl_float3* lightColour,
			cl_float* lightRadius,
			cl_float* lightHalfAttenuationDistance,

			cl_ushort* outputRed,
			cl_ushort* outputGreen,
			cl_ushort* outputBlue);

#endif // RAYTRACE_H
