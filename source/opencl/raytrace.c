#pragma warning(disable : 4756)

#ifdef _WIN32
	#include <windows.h>
#endif // WIN32

#include "opencl/raytrace_opencl.bin.h"
#include "opencl/raytrace_opencl.bin.c"

#include "raytrace.h"
#include <stdio.h>
#include <string.h>

// OPENCL-TO-C FUNCTIONS - functions that are provided by OpenCL standard but not C that we used
int get_global_id(uint dim) {
	return -1;
}
float dot(float3 a, float3 b) {
	return a.s[0] * b.s[0] + a.s[1] * b.s[1] + a.s[2] * b.s[2];
}
float3 cross(float3 a, float3 b) {
	float3 c;
	c.s[0] = a.s[1] * b.s[2] - a.s[2] * b.s[1];
	c.s[1] = a.s[2] * b.s[0] - a.s[0] * b.s[2];
	c.s[2] = a.s[0] * b.s[1] - a.s[1] * b.s[0];
	return c;
}
float3 normalize(float3 v) {
	float len = (float)sqrt(dot(v, v));
	float3 r;
	r.s[0] = v.s[0] / len;
	r.s[1] = v.s[1] / len;
	r.s[2] = v.s[2] / len;
	return r;
}
float3 vector(float3 a, float3 b) {
	float3 ab;
	ab.s[0] = b.s[0] - a.s[0];
	ab.s[1] = b.s[1] - a.s[1];
	ab.s[2] = b.s[2] - a.s[2];
	return ab;
}
float bindf(float value, float a, float b) {
	return min(max(value, a), b);
}

// OPENCL-TO-C HACKS FOR COMPATIBILITY
//
// Note: This will not make every OpenCL code compatible with C, only the current one.
//       Some syntax is really hard to convert to C and we avoided it for easier debug.
//
// TODO test OpenCL optimizations to make sure we are not shooting ourselves in the foot
//      by not using them.
//
// Examples of unsupported syntax with our current scheme:
// int3 a, b, c;
// a = b + c;
// a.xy = b.xy + c.xy;
// a.xy = b.yz + c.zx;
//
// Example of supported syntax:
// a.x = b.y + c.z;
#define false CL_FALSE
#define true CL_TRUE
#define x s[0]
#define y s[1]
#define z s[2]

// We include the OpenCL code with our functions and defines to try to get it to compile in regular C
#include "raytrace_opencl.c"


// Initialize OpenCL platforms and devices and get their names for selection
cl_bool isPlatformDeviceUpdated = CL_FALSE;
cl_int platformDeviceCount = 0;
cl_int2 platformDevice[256] = { 0 };
cl_char platformDeviceName[256 * 256] = { 0 };
void InitOpenCL() {
	cl_int err = 0;
	size_t _platformDeviceCount = 0;
	cl_int2 _platformDevice[256];
	cl_char _platformDeviceName[256 * 256];
	cl_platform_id platforms[64] = { 0 };
	cl_uint i, j, platformsCount = 0;
	err = clGetPlatformIDs(sizeof(platforms) / sizeof(cl_platform_id), platforms, &platformsCount);
	if (CL_SUCCESS == err) {
		cl_context_properties context_properties[] = { CL_CONTEXT_PLATFORM, 0, 0 };
		for (i = 0; i < platformsCount; ++i) {
			size_t len = 0;
			cl_context deviceContext;
			clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 256, &_platformDeviceName[256 * _platformDeviceCount], NULL);
			strcat(&_platformDeviceName[256 * _platformDeviceCount], " ");
			len = strlen(&_platformDeviceName[256 * _platformDeviceCount]);
			context_properties[1] = (cl_context_properties)platforms[i];
			deviceContext = clCreateContextFromType(context_properties, CL_DEVICE_TYPE_ALL, NULL, NULL, &err);
			if (CL_SUCCESS == err && deviceContext) {
				cl_device_id devices[64] = { 0 };
				size_t devicesSize;
				err = clGetContextInfo(deviceContext, CL_CONTEXT_DEVICES, sizeof(devices), devices, &devicesSize);
				if (CL_SUCCESS == err) {
					size_t deviceCount = devicesSize / sizeof(cl_device_id);
					for (j = 0; j < deviceCount; ++j) {
						_platformDevice[_platformDeviceCount].s[0] = i;
						_platformDevice[_platformDeviceCount].s[1] = j;
						if (0 < j) memcpy(&_platformDeviceName[256 * _platformDeviceCount], &_platformDeviceName[256 * (_platformDeviceCount - 1)], len);
						clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 256 - len, &_platformDeviceName[256 * _platformDeviceCount++ + len], NULL);
						clReleaseDevice(devices[j]); devices[j] = 0;
						if (256 <= _platformDeviceCount) break;
					}
					if (256 <= _platformDeviceCount) break;
				}
				clReleaseContext(deviceContext);
				deviceContext = NULL;
			}
		}
	}
	memcpy(platformDeviceName, _platformDeviceName, 256 * _platformDeviceCount);
	memcpy(platformDevice, _platformDevice, sizeof(cl_int2) * _platformDeviceCount);
	platformDeviceCount = (cl_int)_platformDeviceCount;
	isPlatformDeviceUpdated = CL_TRUE;
}

// Functions to control the initialization of OpenCL
void ResetComputationType() {
	if (isPlatformDeviceUpdated) {
		isPlatformDeviceUpdated = CL_FALSE;
		platformDeviceCount = 0;
	}
}
cl_bool GetIsComputationTypeUpdated() {
	return isPlatformDeviceUpdated;
}
size_t GetComputationTypeCount() {
	return 1 + platformDeviceCount;
}
cl_bool GetComputationTypeName(size_t id, size_t strLen, cl_char* str) {
	if (0 == id--) {
		cl_char t[] = "Local CPU single thread";
		if (strlen(t) <= strLen) {
			strcpy(str, t);
			return CL_TRUE;
		}
	}
	else {
		if (id < platformDeviceCount) {
			if (strlen(&platformDeviceName[256 * id]) <= strLen) {
				strcpy(str, &platformDeviceName[256 * id]);
				return CL_TRUE;
			}
		}
	}
	return CL_FALSE;
}

// Functions to time the OpenCL kernel processing
cl_float progress;
clock_t startTime, endTime;
cl_float GetProgress() {
	return progress;
}
void SetProgress(cl_float p) {
	progress = p;
}
clock_t GetStartTime() {
	return startTime;
}
clock_t GetEndTime() {
	return endTime;
}
void ResetTime() {
	startTime = 0;
	endTime = 0;
}

// This either runs the OpenCL code on this thread (computationType = 0)
// or it compiles the OpenCL kernel and splits up the data to be run in parallel
// It is optimized for usage locally and on rentmypc.net network
// Visit the rentmypc.net forums for information on how to optimize network OpenCL code

// input:
//   computationType: the machine code will run on. 0 means CPU single thread (this one), x = 1 to n means the OpenCL platform/device # (x - 1)
//   cameraImageDimension: size of the image to be filled
//   cameraEye: location of the camera in the scene
//   cameraEyeToTopLeftVector: vector pointing from the camera location to the top left screen location in the scene
//   cameraLeftToRightPixelSizeVector: vector orthogonal to the camera (look at middle of screen) of the size of one pixel at the distance of cameraEyeToTopLeftVector
//   cameraTopToBottomPixelSizeVector: vector orthogonal to the camera (look at middle of screen) and of cameraLeftToRightPixelSizeVector of the size of one pixel at the distance of cameraEyeToTopLeftVector
//   cameraPixelSizeInv: inverse of the size of the pixel at unit vector level

//   cameraPixelTriangleListStart: this is an array of size cameraImageDimension+1, each pixel will point to the start of a list of triangles (cameraPixelTriangleList) that is possibly in this pixel
//   cameraPixelTriangleList: array of size cameraPixelTriangleListStart[cameraImageDimension] of triangle pointers that is pointed to by cameraPixelTriangleListStart
//       for example, a 2x2 image of a 3 triangle scene could have something like this:
//         listStart = {0, 2, 2, 5, 6}
//         list = {0, 1, 0, 1, 2, 2}
//       this means pixel 0,0 has triangle 0 and 1
//                  pixel 0,1 has no triangle
//                  pixel 1,0 has triangle 0,1,2
//                  pixel 1,1 has triangle 2

//   sampleCount: number of rays sampling the scene per pixel
//   vertexCount: number of vertices in the vertex array
//   vertex: array of vertexCount vertices
//   triangleCount: number of triangles
//   triangleVertexIndex: array of 3 vertex pointers of size triangleCount
//   triangleMaterialId: array of material pointers of size triangleCount
//   triangleUv: array of uv of size 3x triangleCount
//   triangleNormal: array of normals of size 3x triangleCount

//       the scene is divided into arbitrarily sized boxes and the triangles are placed in the boxes if they intersect
//   axesDivCount: number of divisions in all 3 axes. There are axesDivCount*axesDivCount*axesDivCount boxes
//   sceneBoxMin: array of size axesDivCount+1. Limits for boxes in x, y and z dimension
//   scenePixelTriangleListStart: array of size axesDivCount*axesDivCount*axesDivCount+1. Pointer to the start of the triangles that intersect the box
//   scenePixelTriangleList: array of size scenePixelTriangleListStart[axesDivCount*axesDivCount*axesDivCount] of triangle pointers pointed to by scenePixelTriangleListStart. Similar to the way it works with camera except in 3d instead of 2d

//   materialCount: number of materials in the scene
//   materialImageSize: array size _MATERIAL_CHANNEL_COUNT * materialCount. Texture sizes for the material. There are _MATERIAL_CHANNEL_COUNT textures for every material.
//   materialImageStart: array size _MATERIAL_CHANNEL_COUNT * materialCount. Pointer to the start of the texture in textures
//   texturesSize: number of pixels in the array textures
//   textures: array of size texturesSize pointed to by materialImageStart

//   lightCount: number of lights in the scene
//   lightType: array of size lightCount, type of the light e.g. _LIGHT_TYPE_OMNI
//   lightPosition: array of size lightCount, position of the light
//   lightDirection: array of size lightCount, direction of the light
//   lightColour: array of size lightCount, colour of the light
//   lightRadius: array of size lightCount, radius of the light, 0.f means a pin light
//   lightHalfAttenuationDistance: array of size lightCount, distance from the light where it looses half it's intensity

// output:
//   outputImage: a preallocated image of size cameraImageDimension
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
			cl_uint2* materialImageSize, // _MATERIAL_CHANNEL_COUNT x materialCount
			cl_int* materialImageStart, // _MATERIAL_CHANNEL_COUNT x materialCount

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
			cl_ushort* outputBlue) {
	cl_uint nextPixelId = 0;
	cl_uint i, j;
	cl_int err = 0;
	cl_uint sampleId = 0;
	if (0 < computationType--) {
		#define ARGUMENT_COUNT 35
		char* source[2] = { NULL };
		size_t sourceLen[2] = { 0 };
		cl_int2 platDev = platformDevice[computationType];
		cl_context context = NULL;
		cl_device_id device = NULL;
		cl_command_queue commandQueue = NULL;
		cl_program program = NULL;
		cl_kernel kernelRaytraceAll = NULL;
		cl_mem memobjs[ARGUMENT_COUNT] = { NULL };
		// Getting the specified platform and device, creating context
		{
			cl_platform_id platforms[64] = { 0 };
			cl_uint platformsCount = 0;
			err |= clGetPlatformIDs(sizeof(platforms) / sizeof(cl_platform_id), platforms, &platformsCount);
			if (CL_SUCCESS == err) {
				cl_context_properties context_properties[] = { CL_CONTEXT_PLATFORM, 0, 0 };
				context_properties[1] = (cl_context_properties)platforms[platDev.s[0]];
				cl_context deviceContext = clCreateContextFromType(context_properties, CL_DEVICE_TYPE_ALL, NULL, NULL, &err);
				if (CL_SUCCESS == err && deviceContext) {
					cl_device_id devices[64] = { 0 };
					size_t devicesSize;
					err |= clGetContextInfo(deviceContext, CL_CONTEXT_DEVICES, sizeof(devices), devices, &devicesSize);
					if (CL_SUCCESS == err) {
						size_t deviceCount = devicesSize / sizeof(cl_device_id);
						if (platDev.s[1] < deviceCount) {
							context = deviceContext;
							device = devices[platDev.s[1]];
						}
						for (j = 0; j < deviceCount; ++j) { if (devices[j] != device) clReleaseDevice(devices[j]); devices[j] = 0; }
					}
					if (context != deviceContext) clReleaseContext(deviceContext); deviceContext = NULL;
				}
			}
		}
		if (!context || !device) goto cleanup;

		// fetching the sources before compiling
		source[0] = (char*)source_opencl_raytrace_opencl_h;
		sourceLen[0] = source_opencl_raytrace_opencl_h_len;
		source[1] = (char*)source_opencl_raytrace_opencl_c;
		sourceLen[1] = source_opencl_raytrace_opencl_c_len;

		// Command queue creation and program compiling
		commandQueue = clCreateCommandQueue(context, device, CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE /*| CL_QUEUE_PROFILING_ENABLE*/, &err);
		if (err || NULL == commandQueue) goto cleanup;
		program = clCreateProgramWithSource(context, 2, source, sourceLen, &err);
		if (NULL == program) goto cleanup;
		if (CL_SUCCESS != (err = clBuildProgram(program, 1, &device, "", NULL, NULL))) {
			char* outputLog = (char*)malloc(128 * 1024 * sizeof(char));
			if (outputLog) {
				err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(outputLog), outputLog, NULL);
				printf(outputLog);
				free(outputLog);
			}
			goto cleanup;
		}

		// Creating paramters and setting up all parameters for kernel Raytrace
		// FILE* f;
		kernelRaytraceAll = clCreateKernel(program, "Raytrace", &err);
		memobjs[0] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(nextPixelId), &nextPixelId, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 0, sizeof(cl_mem), (void*)&memobjs[0]);
		// f = fopen("memobj00", "w+b"); fwrite(&nextPixelId, sizeof(nextPixelId), 1, f); fclose(f);
		memobjs[1] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cameraImageDimension), &cameraImageDimension, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 1, sizeof(cl_mem), (void*)&memobjs[1]);
		// f = fopen("memobj01", "w+b"); fwrite(&cameraImageDimension, sizeof(cameraImageDimension), 1, f); fclose(f);
		memobjs[2] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cameraEye), &cameraEye, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 2, sizeof(cl_mem), (void*)&memobjs[2]);
		// f = fopen("memobj02", "w+b"); fwrite(&cameraEye, sizeof(cameraEye), 1, f); fclose(f);
		memobjs[3] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cameraEyeToTopLeftVector), &cameraEyeToTopLeftVector, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 3, sizeof(cl_mem), (void*)&memobjs[3]);
		// f = fopen("memobj03", "w+b"); fwrite(&cameraEyeToTopLeftVector, sizeof(cameraEyeToTopLeftVector), 1, f); fclose(f);
		memobjs[4] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cameraLeftToRightPixelSizeVector), &cameraLeftToRightPixelSizeVector, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 4, sizeof(cl_mem), (void*)&memobjs[4]);
		// f = fopen("memobj04", "w+b"); fwrite(&cameraLeftToRightPixelSizeVector, sizeof(cameraLeftToRightPixelSizeVector), 1, f); fclose(f);
		memobjs[5] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cameraTopToBottomPixelSizeVector), &cameraTopToBottomPixelSizeVector, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 5, sizeof(cl_mem), (void*)&memobjs[5]);
		// f = fopen("memobj05", "w+b"); fwrite(&cameraTopToBottomPixelSizeVector, sizeof(cameraTopToBottomPixelSizeVector), 1, f); fclose(f);
		memobjs[6] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(cameraPixelSizeInv), &cameraPixelSizeInv, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 6, sizeof(cl_mem), (void*)&memobjs[6]);
		// f = fopen("memobj06", "w+b"); fwrite(&cameraPixelSizeInv, sizeof(cameraPixelSizeInv), 1, f); fclose(f);
		memobjs[7] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, (cameraImageDimension.s[0] * cameraImageDimension.s[1]) * sizeof(cl_uint), cameraPixelTriangleListStart, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 7, sizeof(cl_mem), (void*)&memobjs[7]);
		// f = fopen("memobj07", "w+b"); fwrite(cameraPixelTriangleListStart, (cameraImageDimension.s[0] * cameraImageDimension.s[1]) * sizeof(cl_uint), 1, f); fclose(f);
		memobjs[8] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, (cameraImageDimension.s[0] * cameraImageDimension.s[1]) * sizeof(cl_uint), cameraPixelTriangleListEnd, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 8, sizeof(cl_mem), (void*)&memobjs[8]);
		// f = fopen("memobj08", "w+b"); fwrite(cameraPixelTriangleListStart, (cameraImageDimension.s[0] * cameraImageDimension.s[1]) * sizeof(cl_uint), 1, f); fclose(f);
		memobjs[9] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, (size_t)cameraPixelTriangleListSize * sizeof(cl_uint), cameraPixelTriangleList, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 9, sizeof(cl_mem), (void*)&memobjs[9]);
		// f = fopen("memobj09", "w+b"); fwrite(cameraPixelTriangleList, cameraPixelTriangleListStart[cameraImageDimension.s[0] * cameraImageDimension.s[1]] * sizeof(cl_uint), 1, f); fclose(f);
		memobjs[10] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, sizeof(sampleId), &sampleId, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 10, sizeof(cl_mem), (void*)&memobjs[10]);
		// f = fopen("memobj10", "w+b"); fwrite(&sampleId, sizeof(sampleId), 1, f); fclose(f);
		memobjs[11] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(sampleCount), &sampleCount, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 11, sizeof(cl_mem), (void*)&memobjs[11]);
		// f = fopen("memobj11", "w+b"); fwrite(&sampleCount, sizeof(sampleCount), 1, f); fclose(f);
		memobjs[12] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, vertexCount * sizeof(cl_float3), vertex, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 12, sizeof(cl_mem), (void*)&memobjs[12]);
		// f = fopen("memobj12", "w+b"); fwrite(vertex, vertexCount * sizeof(cl_float3), 1, f); fclose(f);
		memobjs[13] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(triangleCount), &triangleCount, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 13, sizeof(cl_mem), (void*)&memobjs[13]);
		// f = fopen("memobj13", "w+b"); fwrite(&triangleCount, sizeof(triangleCount), 1, f); fclose(f);
		memobjs[14] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, triangleCount * sizeof(cl_int3), triangleVertexIndex, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 14, sizeof(cl_mem), (void*)&memobjs[14]);
		// f = fopen("memobj14", "w+b"); fwrite(triangleVertexIndex, triangleCount * sizeof(cl_int3), 1, f); fclose(f);
		memobjs[15] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, triangleCount * sizeof(cl_int), triangleMaterialId, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 15, sizeof(cl_mem), (void*)&memobjs[15]);
		// f = fopen("memobj15", "w+b"); fwrite(triangleMaterialId, triangleCount * sizeof(cl_int), 1, f); fclose(f);
		memobjs[16] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, 3 * triangleCount * sizeof(cl_float2), triangleUv, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 16, sizeof(cl_mem), (void*)&memobjs[16]);
		// f = fopen("memobj16", "w+b"); fwrite(triangleUv, 3 * triangleCount * sizeof(cl_float2), 1, f); fclose(f);
		memobjs[17] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, 3 * triangleCount * sizeof(cl_float3), triangleNormal, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 17, sizeof(cl_mem), (void*)&memobjs[17]);
		// f = fopen("memobj17", "w+b"); fwrite(triangleNormal, 3 * triangleCount * sizeof(cl_float3), 1, f); fclose(f);
		memobjs[18] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(axesDivCount), &axesDivCount, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 18, sizeof(cl_mem), (void*)&memobjs[18]);
		// f = fopen("memobj18", "w+b"); fwrite(&axesDivCount, sizeof(axesDivCount), 1, f); fclose(f);
		memobjs[19] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, (axesDivCount + 1) * sizeof(cl_float3), sceneBoxMin, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 19, sizeof(cl_mem), (void*)&memobjs[19]);
		// f = fopen("memobj19", "w+b"); fwrite(sceneBoxMin, (axesDivCount + 1) * sizeof(cl_float3), 1, f); fclose(f);
		memobjs[20] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, (axesDivCount * axesDivCount * axesDivCount + 1) * sizeof(cl_uint), scenePixelTriangleListStart, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 20, sizeof(cl_mem), (void*)&memobjs[20]);
		// f = fopen("memobj20", "w+b"); fwrite(scenePixelTriangleListStart, (axesDivCount * axesDivCount * axesDivCount + 1) * sizeof(cl_uint), 1, f); fclose(f);
		memobjs[21] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, scenePixelTriangleListStart[axesDivCount * axesDivCount * axesDivCount] * sizeof(cl_uint), scenePixelTriangleList, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 21, sizeof(cl_mem), (void*)&memobjs[21]);
		// f = fopen("memobj21", "w+b"); fwrite(scenePixelTriangleList, scenePixelTriangleListStart[axesDivCount * axesDivCount * axesDivCount] * sizeof(cl_uint), 1, f); fclose(f);
		memobjs[22] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, materialCount * _MATERIAL_CHANNEL_COUNT * sizeof(cl_uint2), materialImageSize, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 22, sizeof(cl_mem), (void*)&memobjs[22]);
		// f = fopen("memobj22", "w+b"); fwrite(materialImageSize, materialCount * _MATERIAL_CHANNEL_COUNT * sizeof(cl_uint2), 1, f); fclose(f);
		memobjs[23] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, materialCount * _MATERIAL_CHANNEL_COUNT * sizeof(cl_int), materialImageStart, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 23, sizeof(cl_mem), (void*)&memobjs[23]);
		// f = fopen("memobj23", "w+b"); fwrite(materialImageStart, materialCount * _MATERIAL_CHANNEL_COUNT * sizeof(cl_int), 1, f); fclose(f);
		memobjs[24] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, materialImageStart[materialCount * _MATERIAL_CHANNEL_COUNT] * sizeof(cl_uchar3), textures, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 24, sizeof(cl_mem), (void*)&memobjs[24]);
		// f = fopen("memobj24", "w+b"); fwrite(textures, materialImageStart[materialCount * _MATERIAL_CHANNEL_COUNT] * sizeof(cl_uchar3), 1, f); fclose(f);
		memobjs[25] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(lightCount), &lightCount, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 25, sizeof(cl_mem), (void*)&memobjs[25]);
		// f = fopen("memobj25", "w+b"); fwrite(&lightCount, sizeof(lightCount), 1, f); fclose(f);
		memobjs[26] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, lightCount * sizeof(cl_int), lightType, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 26, sizeof(cl_mem), (void*)&memobjs[26]);
		// f = fopen("memobj26", "w+b"); fwrite(lightType, lightCount * sizeof(cl_int), 1, f); fclose(f);
		memobjs[27] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, lightCount * sizeof(cl_float3), lightPosition, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 27, sizeof(cl_mem), (void*)&memobjs[27]);
		// f = fopen("memobj27", "w+b"); fwrite(lightPosition, lightCount * sizeof(cl_float3), 1, f); fclose(f);
		memobjs[28] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, lightCount * sizeof(cl_float3), lightDirection, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 28, sizeof(cl_mem), (void*)&memobjs[28]);
		// f = fopen("memobj28", "w+b"); fwrite(lightDirection, lightCount * sizeof(cl_float3), 1, f); fclose(f);
		memobjs[29] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, lightCount * sizeof(cl_float3), lightColour, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 29, sizeof(cl_mem), (void*)&memobjs[29]);
		// f = fopen("memobj29", "w+b"); fwrite(lightColour, lightCount * sizeof(cl_float3), 1, f); fclose(f);
		memobjs[30] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, lightCount * sizeof(cl_float), lightRadius, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 30, sizeof(cl_mem), (void*)&memobjs[30]);
		// f = fopen("memobj30", "w+b"); fwrite(lightRadius, lightCount * sizeof(cl_float), 1, f); fclose(f);
		memobjs[31] = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, lightCount * sizeof(cl_float), lightHalfAttenuationDistance, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 31, sizeof(cl_mem), (void*)&memobjs[31]);
		// f = fopen("memobj31", "w+b"); fwrite(lightHalfAttenuationDistance, lightCount * sizeof(cl_float), 1, f); fclose(f);
		// I will force the write memory to 0 here for two reasons:
		// 1) Every buffer is sent to the client because we cannot detect if a kernel has written a bit if the bit is already at the boolean value. We need to send the buffer to be in sync with client. 0 is easier to compress than noise.
		// 2) Our sampling progressively increases the value from 0. We could clear the buffer on the client, but this would slightly increase the number of kernel calls
		memset(outputRed, 0, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort));
		memobjs[32] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort), outputRed, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 32, sizeof(cl_mem), (void*)&memobjs[32]);
		// f = fopen("memobj32", "w+b"); fwrite(outputImage, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort), 1, f); fclose(f);
		memset(outputGreen, 0, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort));
		memobjs[33] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort), outputGreen, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 33, sizeof(cl_mem), (void*)&memobjs[33]);
		// f = fopen("memobj33", "w+b"); fwrite(outputImage, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort), 1, f); fclose(f);
		memset(outputBlue, 0, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort));
		memobjs[34] = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort), outputBlue, NULL);
		err |= clSetKernelArg(kernelRaytraceAll, 34, sizeof(cl_mem), (void*)&memobjs[34]);
		// f = fopen("memobj34", "w+b"); fwrite(outputImage, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_ushort), 1, f); fclose(f);

		{
			// We will split up the work into chunks that are worth it to run on a network
			// We have to upload the scene and download the image, plus the network agent
			// has to compile the code. We need a good uptime: we are paying for the total
			// time. On the other hand, we don't want chunks that run for hours: it might
			// crash and we don't want to wait hours for the result. We will try to get a
			// 90% uptime. Assuming that the scene is rather large, upload and compile
			// could take 5-10 minutes. We will try to aim for a 45 minutes on average
			// computation on both CPU and GPU.
			#define NUMBER_OF_RAYS_PER_FORTYFIVE_MIN_CPU (60 * 1024 * 1024)
			#define NUMBER_OF_RAYS_PER_FORTYFIVE_MIN_GPU (128 * 1024 * 1024)

			// MAX_WORK_GROUP_SIZE_SQRT x MAX_WORK_GROUP_SIZE_SQRT must be big enough to
			// fill all parallel processors with some to spare. With modern GPU having
			// 4096 and more processors, we want this to be relatively large.
			#define MAX_WORK_GROUP_SIZE_SQRT 128

			size_t global_work_size[2] = { MAX_WORK_GROUP_SIZE_SQRT, MAX_WORK_GROUP_SIZE_SQRT };
			size_t global_work_offset[2] = { 0, 0 };
			size_t rayCount = 0;
			size_t numberOfRaysPerTwentyFiveMinGpu = NUMBER_OF_RAYS_PER_FORTYFIVE_MIN_CPU;
			size_t workId, workDone = 0, workCount = 0, maxWorkCount;
			cl_device_type deviceType;
			cl_event* partEndEvent;

			err |= clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(deviceType), &deviceType, NULL);
			if (err) goto cleanup;
			if (CL_DEVICE_TYPE_GPU & deviceType) {
				numberOfRaysPerTwentyFiveMinGpu = NUMBER_OF_RAYS_PER_FORTYFIVE_MIN_GPU;
			}

			maxWorkCount = (((cameraImageDimension.s[0] + global_work_size[0] - 1) / global_work_size[0]) *
							((cameraImageDimension.s[1] + global_work_size[1] - 1) / global_work_size[1]));
			partEndEvent = (cl_event*)malloc(maxWorkCount * sizeof(cl_event));
			memset(partEndEvent, 0, maxWorkCount * sizeof(cl_event));
			// We will enqueue kernels to be sent to the compute device. As long as we don't clFlush, clFinish
			// or enqueue an IO operation (clEnqueueReadBuffer, clEnqueueWriteBuffer, ...), the kernel calls are
			// accumulated to be sent as a batch to the device. Since we are thinking about network computation,
			// we need to make sure we send a batch that is reasonable in size.
			// TODO we might want to test the computation time before sending it out: create a small batch of
			// ray traces and (CL_QUEUE_PROFILING_ENABLE) profile them on the local CPU to have an idea of how
			// big the scene is to make sure we don't send batches that take less than 5 min to complete
			// (expensive downtime) or more than an hour (lack of parallelism, crash danger, waiting for results)
			while (global_work_offset[1] < cameraImageDimension.s[1]) {
				global_work_offset[0] = 0;
				while (global_work_offset[0] < cameraImageDimension.s[0]) {
					cl_event e0 = NULL, e1 = NULL;
					for (cl_uint sId = 0; sId < sampleCount; ++sId) {
						err |= clEnqueueNDRangeKernel(commandQueue, kernelRaytraceAll, 2, global_work_offset, global_work_size, NULL, e0 ? 1 : 0, e0 ? &e0 : NULL, &e1);
						if (e0) clReleaseEvent(e0); e0 = e1;
						e1 = NULL;
					}
					global_work_offset[0] += global_work_size[0];
					rayCount += (size_t)sampleCount * global_work_size[0] * global_work_size[1];
					if (partEndEvent[workCount]) clReleaseEvent(partEndEvent[workCount]);
					partEndEvent[workCount] = e0;
					if (numberOfRaysPerTwentyFiveMinGpu <= rayCount) {
						// The accumulated batch is large enough, flush it to a computer
						++workCount;
						clFlush(commandQueue);
						rayCount = 0;
					}
				}
				global_work_offset[1] += global_work_size[1];
			}
			if (0 < rayCount) {
				// If the batch is non null, flush it to a computer
				++workCount;
				clFlush(commandQueue);
				rayCount = 0;
			}
			progress = 0.f;
			startTime = clock();
			endTime = startTime;
			while (CL_TRUE) {
				for (workId = 0; workId < workCount; ++workId) {
					// Check for done batches to give the completion percentage
					if (partEndEvent[workId]) {
						cl_int eventStatus;
						err |= clGetEventInfo(partEndEvent[workId], CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(eventStatus), &eventStatus, NULL);
						if (err || (CL_QUEUED != eventStatus && CL_SUBMITTED != eventStatus && CL_RUNNING != eventStatus)) {
							clReleaseEvent(partEndEvent[workId]);
							partEndEvent[workId] = NULL;
							++workDone;
						}
					}
				}
				// Progress will be stuck to 99.9% until the bitmap is copied to local memory
				progress = 0.999f * (cl_float)workDone / (cl_float)workCount;
				if (workCount <= workDone) break;
				#ifdef _WIN32
					Sleep(1000);
				#else
					usleep(uSEC_PER_MSEC * 1000);
				#endif
			}
			//err |= clEnqueueReadBuffer(commandQueue, memobjs[32], CL_TRUE, 0, cameraImageDimension.s[0] * cameraImageDimension.s[1] * sizeof(cl_float3), outputImage, 0, NULL, NULL);
			err |= clFinish(commandQueue);
			if (err) goto cleanup;
			endTime = clock();
		}

	cleanup:
		if (context) clReleaseContext(context); context = NULL;
		if (device) clReleaseDevice(device); device = NULL;
		if (commandQueue) clReleaseCommandQueue(commandQueue); commandQueue = NULL;
		if (program) clReleaseProgram(program); program = NULL;
		if (kernelRaytraceAll) clReleaseKernel(kernelRaytraceAll); kernelRaytraceAll = NULL;
		for (i = 0; i < ARGUMENT_COUNT; ++i) {
			if (memobjs[i]) clReleaseMemObject(memobjs[i]); memobjs[i] = 0;
		}
	}
	else {
		// Local single thread simulation of the OpenCL code. Mainly used for debugging purposes.
		cl_uint imageSize = cameraImageDimension.s[0] * cameraImageDimension.s[1];
		cl_float incProgress = 0.999f / ((cl_float)imageSize * (cl_float)sampleCount);
		nextPixelId = 0;
		progress = 0.f;
		startTime = clock();
		endTime = startTime;
		for (; nextPixelId < imageSize; ++nextPixelId) {
			cl_uint sampleId = 0;
			progress = 0.999f * (cl_float)nextPixelId / (cl_float)imageSize;
			while (sampleId < sampleCount) {
				Raytrace(	&nextPixelId,
							&cameraImageDimension,
							&cameraEye,
							&cameraEyeToTopLeftVector,
							&cameraLeftToRightPixelSizeVector,
							&cameraTopToBottomPixelSizeVector,
							&cameraPixelSizeInv,
							cameraPixelTriangleListStart,
							cameraPixelTriangleListEnd,
							cameraPixelTriangleList,
							&sampleId,
							&sampleCount,
							vertex,
							&triangleCount,
							triangleVertexIndex,
							triangleMaterialId,
							triangleUv,
							triangleNormal,
							&axesDivCount,
							sceneBoxMin,
							scenePixelTriangleListStart,
							scenePixelTriangleList,
							materialImageSize,
							materialImageStart,
							textures,
							&lightCount,
							lightType,
							lightPosition,
							lightDirection,
							lightColour,
							lightRadius,
							lightHalfAttenuationDistance,
							outputRed,
							outputGreen,
							outputBlue);
				progress += incProgress;
			}
		}
		endTime = clock();
	}
	return 0 == err;
}
