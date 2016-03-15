ulong rotl64(ulong v, int n) {
	return (v << n) | (v >> (64 - n));
}

ulong xorshift64star(ulong v) {
	v ^= v >> 12; // a
	v ^= v << 25; // b
	v ^= v >> 27; // c
	return v * 2685821657736338717ul;
}

float randF(__private ulong* rndData, float min, float max) {
	*rndData ^= xorshift64star((rotl64(*rndData, 55) ^ rotl64(*rndData, 3)) * 0xc23f3c0ad9da6357ul);
	*rndData ^= xorshift64star((rotl64(*rndData, 35) ^ rotl64(*rndData, 3)) ^ 0xce84d6af03c16b89ul);
	*rndData ^= xorshift64star((rotl64(*rndData, 63) ^ rotl64(*rndData, 35)) * 0xf097ef8bbe03ddccul);
	*rndData ^= xorshift64star((rotl64(*rndData, 41) ^ rotl64(*rndData, 12)) ^ 0x48302294fbfe30bful);
	*rndData ^= xorshift64star((rotl64(*rndData, 1) ^ rotl64(*rndData, 62)) * 0x79e7425e3f4f147dul);
	*rndData ^= xorshift64star((rotl64(*rndData, 42) ^ rotl64(*rndData, 29)) ^ 0x14d1d30856e5be9aul);
	*rndData ^= xorshift64star((rotl64(*rndData, 47) ^ rotl64(*rndData, 45)) * 0x24289d47a66617c3ul);
	*rndData ^= xorshift64star((rotl64(*rndData, 39) ^ rotl64(*rndData, 6)) ^ 0x5576fb2f80a05d14ul);

	return min + (max - min) * (float)((double)*rndData / (double)0xfffffffffffffffful);
}

float positive_modf(float i) {
	double intpart;
	return (float)modf(modf((double)i, &intpart) + 1., &intpart);
}

float3 GetSpherePoint(__private ulong* rndData, float radius) {
	float3 spherePt;
	float len;
	float tmp;
	do {
		spherePt.x = randF(rndData, -1, 1);
		spherePt.y = randF(rndData, -1, 1);
		spherePt.z = randF(rndData, -1, 1);
		len = (float)sqrt(dot(spherePt, spherePt));
	} while (len <= 0.f);
	tmp = (float)sqrt(randF(rndData, 0, 1)) * radius / len;
	spherePt.x = tmp * spherePt.x;
	spherePt.y = tmp * spherePt.y;
	spherePt.z = tmp * spherePt.z;
	return spherePt;
}

/* Not used anymore, but keeping it just in case
cl_float GetPointToSegmentSqLen(cl_float3 point, cl_float3 segmentOrigin, cl_float3 segment, cl_float segmentMult) {
	cl_float distanceSegmentToLightCenterSq;
	cl_float3 eyeToLight;
	eyeToLight.x = point.x - segmentOrigin.x;
	eyeToLight.y = point.y - segmentOrigin.y;
	eyeToLight.z = point.z - segmentOrigin.z;
	cl_float proj1 = dot(eyeToLight, segment);
	if (proj1 <= 0) {
		cl_float3 rayToLight;
		rayToLight.x = point.x - segmentOrigin.x;
		rayToLight.y = point.y - segmentOrigin.y;
		rayToLight.z = point.z - segmentOrigin.z;
		distanceSegmentToLightCenterSq = dot(rayToLight, rayToLight);
	}
	else {
		cl_float proj2 = dot(segment, segment);
		if (segmentMult < INFINITY && proj2 <= proj1) {
			cl_float3 lightToRayEnd;
			lightToRayEnd.x = point.x - (segmentOrigin.x + segmentMult * segment.x);
			lightToRayEnd.y = point.y - (segmentOrigin.y + segmentMult * segment.y);
			lightToRayEnd.z = point.z - (segmentOrigin.z + segmentMult * segment.z);
			distanceSegmentToLightCenterSq = dot(lightToRayEnd, lightToRayEnd);
		}
		else {
			cl_float closestPointPct = proj1 / proj2;
			cl_float3 lightToRay;
			lightToRay.x = segmentOrigin.x + closestPointPct * segment.x;
			lightToRay.y = segmentOrigin.y + closestPointPct * segment.y;
			lightToRay.z = segmentOrigin.z + closestPointPct * segment.z;
			distanceSegmentToLightCenterSq = dot(lightToRay, lightToRay);
		}
	}
	return distanceSegmentToLightCenterSq;
}*/

float GetPointToLineSqLen(float3 origin, float3 destination, float3 point) {
	float3 od, op, pointProj, distVect;
	float odSqLen, opProjLen;
	od.x = destination.x - origin.x;
	od.y = destination.y - origin.y;
	od.z = destination.z - origin.z;
	odSqLen = dot(od, od);
	op.x = point.x - origin.x;
	op.y = point.y - origin.y;
	op.z = point.z - origin.z;
	opProjLen = dot(op, od) / odSqLen;
	pointProj.x = origin.x + opProjLen * od.x;
	pointProj.y = origin.y + opProjLen * od.y;
	pointProj.z = origin.z + opProjLen * od.z;
	distVect.x = pointProj.x - point.x;
	distVect.y = pointProj.y - point.y;
	distVect.z = pointProj.z - point.z;
	return dot(distVect, distVect);
}

float3 Get2dTableValue3(__global uchar3* table, uint2 tableSize, float2 aUV, float2 bUV, float2 cUV, float abL, float acL) {
	float2 texPct;
	float2 texLoc;
	int2 texLocFloor;
	int texLocMap;
	uchar3 tableValue;
	float3 returnValue;
	texPct.x = positive_modf(aUV.x + (bUV.x - aUV.x) * abL + (cUV.x - aUV.x) * acL);
	texPct.y = positive_modf(aUV.y + (bUV.y - aUV.y) * abL + (cUV.y - aUV.y) * acL);
	texLoc.x = texPct.x * (float)(tableSize.x - 1);
	texLoc.y = texPct.y * (float)(tableSize.y - 1);
	texLocFloor.x = (int)floor(texLoc.x);
	texLocFloor.y = (int)floor(texLoc.y);
	texLocMap = texLocFloor.x + texLocFloor.y * tableSize.x;
	tableValue = table[texLocMap];
	returnValue.x = (tableValue.x) / 255.f;
	returnValue.y = (tableValue.y) / 255.f;
	returnValue.z = (tableValue.z) / 255.f;
	return returnValue;
}

bool RayIntersectsTriangle(float3 origin, float3 ray, float minDistance, float maxDistance, float3 a, float3 b, float3 c, __private float* outRayMult, __private float* outABL, __private float* outACL) {
	// http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()
	bool intersects = false;
	float3 ab;
	float3 ac;
	float3 ao;
	float3 normal;
	ab.x = b.x - a.x;
	ab.y = b.y - a.y;
	ab.z = b.z - a.z;
	ac.x = c.x - a.x;
	ac.y = c.y - a.y;
	ac.z = c.z - a.z;
	ao.x = origin.x - a.x;
	ao.y = origin.y - a.y;
	ao.z = origin.z - a.z;
	normal = cross(ac, ab);
	*outRayMult = -dot(normal, ao) / dot(normal, ray);

	if (minDistance < *outRayMult && *outRayMult < maxDistance) {
		float3 proj;

		float abab = dot(ab, ab);
		float abac = dot(ab, ac);
		float acac = dot(ac, ac);
		float D = 1.f / (abac * abac - abab * acac);

		float3 aproj;
		float aprojab;
		float aprojac;

		proj.x = origin.x + *outRayMult * ray.x;
		proj.y = origin.y + *outRayMult * ray.y;
		proj.z = origin.z + *outRayMult * ray.z;

		aproj.x = proj.x - a.x;
		aproj.y = proj.y - a.y;
		aproj.z = proj.z - a.z;

		aprojab = dot(aproj, ab);
		aprojac = dot(aproj, ac);

		*outABL = (abac * aprojac - acac * aprojab) * D;
		*outACL = (abac * aprojab - abab * aprojac) * D;

		intersects = (0 <= *outABL && 0 <= *outACL && *outABL + *outACL <= 1.f);
	}
	return intersects;
}

int3 GetBoxAddress(int axesDivCount, __global float3* boxMin, float3 position) {
	int3 axisMin = { 0, 0, 0 };
	while (1 < axesDivCount) {
		int midway;
		axesDivCount /= 2;
		midway = axisMin.x + axesDivCount;
		if (boxMin[midway].x < position.x) {
			axisMin.x = midway;
		}
		midway = axisMin.y + axesDivCount;
		if (boxMin[midway].y < position.y) {
			axisMin.y = midway;
		}
		midway = axisMin.z + axesDivCount;
		if (boxMin[midway].z < position.z) {
			axisMin.z = midway;
		}
	}
	return axisMin;
}

float3 GetTriangleNormal(	float3 location,
							float3 rayOrigin,
							float3 rayVector,
							uint triangleId,
							float triangleABL,
							float triangleACL,
							__global int3* triangleVertexIndex,
							__global float3* triangleNormal, // triangleCount * 3
							__global int* triangleMaterialId,
							__global float2* triangleUv, // triangleCount * 3
							__global float3* vertex,
							__global uint2* materialImageSize, // materialCount * _MATERIAL_CHANNEL_COUNT
							__global int* materialImageStart, // materialCount * _MATERIAL_CHANNEL_COUNT
							__global uchar3* textures,
							float3 cameraTopToBottomPixelSizeVector,
							float3 cameraLeftToRightPixelSizeVector,
							float cameraPixelSizeInv) {
	int3 vi = triangleVertexIndex[triangleId];
	float3 aNorm = triangleNormal[3 * triangleId + 0];
	float3 bNorm = triangleNormal[3 * triangleId + 1];
	float3 cNorm = triangleNormal[3 * triangleId + 2];
	int materialId = triangleMaterialId[triangleId];
	// Along the edge AB, C must not affect the normal, otherwise the 2 triangle will not be smooth on edge AB
	float3 a = vertex[vi.x];
	float3 b = vertex[vi.y];
	float3 c = vertex[vi.z];
	float dab = (float)sqrt(GetPointToLineSqLen(a, b, location));
	float dbc = (float)sqrt(GetPointToLineSqLen(b, c, location));
	float dca = (float)sqrt(GetPointToLineSqLen(c, a, location));
	float dabcInv = 1.f / (dab + dbc + dca);
	float3 localNormal;
	uint2 bumpSize = materialImageSize[_MATERIAL_CHANNEL_COUNT * materialId + _MATERIAL_CHANNEL_BUMP];
	localNormal.x = (dab * cNorm.x + dbc * aNorm.x + dca * bNorm.x) * dabcInv;
	localNormal.y = (dab * cNorm.y + dbc * aNorm.y + dca * bNorm.y) * dabcInv;
	localNormal.z = (dab * cNorm.z + dbc * aNorm.z + dca * bNorm.z) * dabcInv;

	if (0 <= materialId && 0 < bumpSize.x) {
		float2 uv0 = triangleUv[3 * triangleId + 0];
		float2 uv1 = triangleUv[3 * triangleId + 1];
		float2 uv2 = triangleUv[3 * triangleId + 2];
		float mult, abL, acL;
		float3 rayS, rayE;
		float3 height, heightS, heightE;
		float xPart, yPart, normalPart, lenInv;
		__global uchar3* bump = &textures[materialImageStart[_MATERIAL_CHANNEL_COUNT * materialId + _MATERIAL_CHANNEL_BUMP]];
		height = Get2dTableValue3(bump, bumpSize, uv0, uv1, uv2, triangleABL, triangleACL);
		rayS.x = rayVector.x + cameraTopToBottomPixelSizeVector.x;
		rayS.y = rayVector.y + cameraTopToBottomPixelSizeVector.y;
		rayS.z = rayVector.z + cameraTopToBottomPixelSizeVector.z;
		RayIntersectsTriangle(rayOrigin, rayS, 0.f, INFINITY, vertex[vi.x], vertex[vi.y], vertex[vi.z], &mult, &abL, &acL);
		heightS = Get2dTableValue3(bump, bumpSize, uv0, uv1, uv2, abL, acL);
		rayE.x = rayVector.x + cameraLeftToRightPixelSizeVector.x;
		rayE.y = rayVector.y + cameraLeftToRightPixelSizeVector.y;
		rayE.z = rayVector.z + cameraLeftToRightPixelSizeVector.z;
		RayIntersectsTriangle(rayOrigin, rayE, 0.f, INFINITY, vertex[vi.x], vertex[vi.y], vertex[vi.z], &mult, &abL, &acL);
		heightE = Get2dTableValue3(bump, bumpSize, uv0, uv1, uv2, abL, acL);
		xPart = (float)sin((heightE.x - height.x) * M_PI / 2.f);
		yPart = (float)sin((heightS.x - height.x) * M_PI / 2.f);
		normalPart = (float)cos((heightE.x - height.x) * M_PI / 2.f) * (float)cos((heightS.x - height.x) * M_PI / 2.f);
		localNormal.x = normalPart * localNormal.x / cameraPixelSizeInv + xPart * cameraLeftToRightPixelSizeVector.x + yPart * cameraTopToBottomPixelSizeVector.x;
		localNormal.y = normalPart * localNormal.y / cameraPixelSizeInv + xPart * cameraLeftToRightPixelSizeVector.y + yPart * cameraTopToBottomPixelSizeVector.y;
		localNormal.z = normalPart * localNormal.z / cameraPixelSizeInv + xPart * cameraLeftToRightPixelSizeVector.z + yPart * cameraTopToBottomPixelSizeVector.z;
		lenInv = 1.f / (float)sqrt(dot(localNormal, localNormal));
		localNormal.x *= lenInv;
		localNormal.y *= lenInv;
		localNormal.z *= lenInv;
	}
	return localNormal;
}

bool BindInCube(float3* origin, float3 ray, float3 cubeMin, float3 cubeMax) {
	float t;
	if (origin[0].x < cubeMin.x) {
		if (ray.x <= 0) {
			return false;
		}
		t = (cubeMin.x - origin[0].x) / ray.x;
		origin[0].x += t * ray.x;
		origin[0].y += t * ray.y;
		origin[0].z += t * ray.z;
	}
	if (cubeMax.x < origin[0].x) {
		if (0 <= ray.x) {
			return false;
		}
		t = (cubeMax.x - origin[0].x) / ray.x;
		origin[0].x += t * ray.x;
		origin[0].y += t * ray.y;
		origin[0].z += t * ray.z;
	}
	if (origin[0].y < cubeMin.y) {
		if (ray.y <= 0) {
			return false;
		}
		t = (cubeMin.y - origin[0].y) / ray.y;
		origin[0].x += t * ray.x;
		origin[0].y += t * ray.y;
		origin[0].z += t * ray.z;
	}
	if (cubeMax.y < origin[0].y) {
		if (0 <= ray.y) {
			return false;
		}
		t = (cubeMax.y - origin[0].y) / ray.y;
		origin[0].x += t * ray.x;
		origin[0].y += t * ray.y;
		origin[0].z += t * ray.z;
	}
	if (origin[0].z < cubeMin.z) {
		if (ray.z <= 0) {
			return false;
		}
		t = (cubeMin.z - origin[0].z) / ray.z;
		origin[0].x += t * ray.x;
		origin[0].y += t * ray.y;
		origin[0].z += t * ray.z;
	}
	if (cubeMax.z < origin[0].z) {
		if (0 <= ray.z) {
			return false;
		}
		t = (cubeMax.z - origin[0].z) / ray.z;
		origin[0].x += t * ray.x;
		origin[0].y += t * ray.y;
		origin[0].z += t * ray.z;
	}
	return true;
}

uint RayIntersectsTriangles(float3 origin,
							float3 ray,
							float minDistance,
							float maxDistance,
							uint excludedTriangleIndex,
							__global int3* triangleVertexIndex,
							__global float3* triangleNormal, // triangleCount * 3
							__global int* triangleMaterialId,
							__global float2* triangleUv, // triangleCount * 3
							__global float3* vertex,
							__global uint2* materialImageSize, // materialCount * _MATERIAL_CHANNEL_COUNT
							__global int* materialImageStart, // materialCount * _MATERIAL_CHANNEL_COUNT
							__global uchar3* textures,
							float3 cameraTopToBottomPixelSizeVector,
							float3 cameraLeftToRightPixelSizeVector,
							float cameraPixelSizeInv,
							int axesDivCount,
							__global float3* sceneBoxMin,
							__global uint* scenePixelTriangleListStart,
							__global uint* scenePixelTriangleList,
							__private float* outRayMult,
							__private float* outABL,
							__private float* outACL) {
	uint closestTriangleIndex = (uint)-1;
	int3 cubeId, endCubeId = { -1, -1, -1 };
	float3 distance;
	float3 startLocation, endLocation;
	startLocation.x = origin.x + minDistance * ray.x;
	startLocation.y = origin.y + minDistance * ray.y;
	startLocation.z = origin.z + minDistance * ray.z;
	BindInCube(&startLocation, ray, sceneBoxMin[0], sceneBoxMin[axesDivCount]);
	cubeId = GetBoxAddress(axesDivCount, sceneBoxMin, startLocation);
	if (maxDistance < INFINITY) {
		endLocation.x = origin.x + maxDistance * ray.x;
		endLocation.y = origin.y + maxDistance * ray.y;
		endLocation.z = origin.z + maxDistance * ray.z;
		BindInCube(&endLocation, ray, sceneBoxMin[0], sceneBoxMin[axesDivCount]);
		endCubeId = GetBoxAddress(axesDivCount, sceneBoxMin, endLocation);
	}
	// Rasterising on cubes of different dimensions
	while (true) {
		uint i, longId = cubeId.x + axesDivCount * cubeId.y + axesDivCount * axesDivCount * cubeId.z;
		*outRayMult = maxDistance;
		for (i = scenePixelTriangleListStart[longId]; i < scenePixelTriangleListStart[longId + 1]; ++i) {
			uint triangleIndex = scenePixelTriangleList[i];
			if (excludedTriangleIndex != triangleIndex) {
				int3 t = triangleVertexIndex[triangleIndex];
				float mult, abL, acL;
				if (RayIntersectsTriangle(origin, ray, minDistance, *outRayMult, vertex[t.x], vertex[t.y], vertex[t.z], &mult, &abL, &acL)) {
					closestTriangleIndex = scenePixelTriangleList[i];
					*outRayMult = mult;
					*outABL = abL;
					*outACL = acL;
				}
			}
		}
		if ((closestTriangleIndex != (uint)-1) ||
			(cubeId.x == endCubeId.x && cubeId.y == endCubeId.y && cubeId.z == endCubeId.z)) break;

		distance.x = (sceneBoxMin[cubeId.x + (0 <= ray.x)].x - origin.x) / ray.x;
		distance.y = (sceneBoxMin[cubeId.y + (0 <= ray.y)].y - origin.y) / ray.y;
		distance.z = (sceneBoxMin[cubeId.z + (0 <= ray.z)].z - origin.z) / ray.z;

		if ((distance.x < distance.y) & (distance.x < distance.z)) {
			cubeId.x += (0 <= ray.x) ? 1 : -1;
			if (cubeId.x < 0 || axesDivCount <= cubeId.x) break;
		}
		else if (distance.y < distance.z) {
			cubeId.y += (0 <= ray.y) ? 1 : -1;
			if (cubeId.y < 0 || axesDivCount <= cubeId.y) break;
		}
		else {
			cubeId.z += (0 <= ray.z) ? 1 : -1;
			if (cubeId.z < 0 || axesDivCount <= cubeId.z) break;
		}
	};
	return closestTriangleIndex;
}


#define MAX_RAY_BUFFER_SIZE 12

__kernel void
Raytrace(	__global uint* nextPixelId, 

			__global uint2* cameraImageDimensionP,
			__global float3* cameraEyeP,
			__global float3* cameraEyeToTopLeftVectorP,
			__global float3* cameraLeftToRightPixelSizeVectorP,
			__global float3* cameraTopToBottomPixelSizeVectorP,
			__global float* cameraPixelSizeInvP,

			__global uint* cameraPixelTriangleListStart,
			__global uint* cameraPixelTriangleListEnd,
			__global uint* cameraPixelTriangleList,

			__global uint* sampleIdP,
			__global uint* sampleCountP,

			__global float3* vertex,

			__global uint* triangleCountP,
			__global int3* triangleVertexIndex,
			__global int* triangleMaterialId,
			__global float2* triangleUv, // 3x triangleCount
			__global float3* triangleNormal, // 3x triangleCount

			__global int* axesDivCountP,
			__global float3* sceneBoxMin,
			__global uint* scenePixelTriangleListStart,
			__global uint* scenePixelTriangleList,

			__global uint2* materialImageSize, // materialCount * _MATERIAL_CHANNEL_COUNT
			__global int* materialImageStart, // materialCount * _MATERIAL_CHANNEL_COUNT
			__global uchar3* textures,

			__global uint* lightCountP,
			__global int* lightType,
			__global float3* lightPosition,
			__global float3* lightDirection,
			__global float3* lightColour,
			__global float* lightRadius,
			__global float* lightHalfAttenuationDistance,

			__global ushort* outputRed,
			__global ushort* outputGreen,
			__global ushort* outputBlue) {

	uint2 cameraImageDimension = cameraImageDimensionP[0];
	float3 cameraEye = cameraEyeP[0];
	float3 cameraEyeToTopLeftVector = cameraEyeToTopLeftVectorP[0];
	float3 cameraLeftToRightPixelSizeVector = cameraLeftToRightPixelSizeVectorP[0];
	float3 cameraTopToBottomPixelSizeVector = cameraTopToBottomPixelSizeVectorP[0];
	float cameraPixelSizeInv = cameraPixelSizeInvP[0];
	float3 commulativeOutputColour;
	int cursorBegin;
	int cursorEnd;
	int maxBounces[MAX_RAY_BUFFER_SIZE];
	uint excludedTriangleIndex[MAX_RAY_BUFFER_SIZE];
	float3 rayOrigin[MAX_RAY_BUFFER_SIZE];
	float3 rayVector[MAX_RAY_BUFFER_SIZE];
	float3 rayMultiplier[MAX_RAY_BUFFER_SIZE];
	bool fromCamera[MAX_RAY_BUFFER_SIZE];
	float minDistance[MAX_RAY_BUFFER_SIZE];
	float maxDistance[MAX_RAY_BUFFER_SIZE];
	float tmp;
	float nextPixelX = (float)get_global_id(0);
	float nextPixelY = (float)get_global_id(1);
	uint nextPixel = (uint)get_global_id(1) * cameraImageDimension.x + (uint)get_global_id(0);
	uint sampleCount = *sampleCountP;
	ulong rndSeed = (ulong)nextPixel * (ulong)sampleCount + (ulong)++(sampleIdP[0]);
	int axesDivCount = *axesDivCountP;
	uint lightCount = *lightCountP;
	if (nextPixelX < 0.f) {
		nextPixel = *nextPixelId;
		nextPixelX = (float)(nextPixel % cameraImageDimension.x);
		nextPixelY = (float)(nextPixel / cameraImageDimension.x);
		rndSeed = (ulong)nextPixel * (ulong)sampleCount + (ulong)sampleIdP[0];
	}
	if (cameraImageDimension.x <= nextPixelX || cameraImageDimension.y <= nextPixelY) {
		return;
	}
	if (353 == nextPixelX && 504 == nextPixelY) {
		nextPixelX = nextPixelX;
	}
	
	cursorBegin = 0;
	commulativeOutputColour = (float3) { 0.f, 0.f, 0.f };
	maxBounces[0] = 12;
	excludedTriangleIndex[0] = (uint)-1;
	rayOrigin[0] = cameraEye;
	rayVector[0] = cameraEyeToTopLeftVector;
	tmp = nextPixelX + randF(&rndSeed, 0.f, 1.f);
	rayVector[0].x += cameraLeftToRightPixelSizeVector.x * tmp;
	rayVector[0].y += cameraLeftToRightPixelSizeVector.y * tmp;
	rayVector[0].z += cameraLeftToRightPixelSizeVector.z * tmp;
	tmp = nextPixelY + randF(&rndSeed, 0.f, 1.f);
	rayVector[0].x += cameraTopToBottomPixelSizeVector.x * tmp;
	rayVector[0].y += cameraTopToBottomPixelSizeVector.y * tmp;
	rayVector[0].z += cameraTopToBottomPixelSizeVector.z * tmp;
	rayMultiplier[0] = (float3) { 1.f, 1.f, 1.f };
	fromCamera[0] = true;
	minDistance[0] = 0.f;
	maxDistance[0] = INFINITY;
	cursorEnd = 1;
	for (; cursorBegin != cursorEnd; cursorBegin = (cursorBegin + 1) % MAX_RAY_BUFFER_SIZE) {
		float closestTriangleMult = maxDistance[cursorBegin];
		uint closestTriangleIndex = (uint)-1;
		float closestTriangleABL = 0.f, closestTriangleACL = 0.f;
		uint i, j;
		if (fromCamera[cursorBegin]) {
			for (i = cameraPixelTriangleListStart[nextPixel]; i < cameraPixelTriangleListEnd[nextPixel]; ++i) {
				uint index = cameraPixelTriangleList[i];
				if (excludedTriangleIndex[cursorBegin] != index) {
					int3 t = triangleVertexIndex[index];
					float mult, abL, acL;
					if (RayIntersectsTriangle(rayOrigin[cursorBegin], rayVector[cursorBegin], minDistance[cursorBegin], closestTriangleMult, vertex[t.x], vertex[t.y], vertex[t.z], &mult, &abL, &acL)) {
						closestTriangleMult = mult;
						closestTriangleIndex = index;
						closestTriangleABL = abL;
						closestTriangleACL = acL;
					}
				}
			}
		}
		else {
			closestTriangleIndex = RayIntersectsTriangles(rayOrigin[cursorBegin], rayVector[cursorBegin], minDistance[cursorBegin], maxDistance[cursorBegin], excludedTriangleIndex[cursorBegin], triangleVertexIndex, triangleNormal, triangleMaterialId, triangleUv, vertex, materialImageSize, materialImageStart, textures, cameraTopToBottomPixelSizeVector, cameraLeftToRightPixelSizeVector, cameraPixelSizeInv, axesDivCount, sceneBoxMin, scenePixelTriangleListStart, scenePixelTriangleList, &closestTriangleMult, &closestTriangleABL, &closestTriangleACL);
		}
		if ((uint)-1 != closestTriangleIndex) {
			float3 localTexture = { 0.f, 0.f, 0.f };
			float3 localTransparency = { 0.f, 0.f, 0.f };
			float3 localReflectance = { 0.f, 0.f, 0.f };
			float3 localLuminance = { 0.f, 0.f, 0.f };
			float3 localLight = { 0.f, 0.f, 0.f };

			// light effects RAY TO ALL LIGHTS, IF RAY NOT INTERCEPTED, ADD LIGHT
			float3 faceLights[2] = { { 0.1f, 0.1f, 0.1f }, { 0.1f, 0.1f, 0.1f } };

			int m = triangleMaterialId[closestTriangleIndex];

			float3 closestLocation, closestTriangleNormal;
			closestLocation.x = rayOrigin[cursorBegin].x + closestTriangleMult * rayVector[cursorBegin].x;
			closestLocation.y = rayOrigin[cursorBegin].y + closestTriangleMult * rayVector[cursorBegin].y;
			closestLocation.z = rayOrigin[cursorBegin].z + closestTriangleMult * rayVector[cursorBegin].z;
			closestTriangleNormal = GetTriangleNormal(closestLocation, rayOrigin[cursorBegin], rayVector[cursorBegin], closestTriangleIndex, closestTriangleABL, closestTriangleACL, triangleVertexIndex, triangleNormal, triangleMaterialId, triangleUv, vertex, materialImageSize, materialImageStart, textures, cameraTopToBottomPixelSizeVector, cameraLeftToRightPixelSizeVector, cameraPixelSizeInv);

			if (0 <= m) {
				int mc = _MATERIAL_CHANNEL_COUNT * m;
				__global float2* uv = &triangleUv[3 * closestTriangleIndex + 0];
				if (0 < materialImageSize[mc + _MATERIAL_CHANNEL_COLOR].x)
					localTexture = Get2dTableValue3(&textures[materialImageStart[mc + _MATERIAL_CHANNEL_COLOR]], materialImageSize[mc + _MATERIAL_CHANNEL_COLOR], uv[0], uv[1], uv[2], closestTriangleABL, closestTriangleACL);
				if (0 < materialImageSize[mc + _MATERIAL_CHANNEL_TRANSPARENCY].x)
					localTransparency = Get2dTableValue3(&textures[materialImageStart[mc + _MATERIAL_CHANNEL_TRANSPARENCY]], materialImageSize[mc + _MATERIAL_CHANNEL_TRANSPARENCY], uv[0], uv[1], uv[2], closestTriangleABL, closestTriangleACL);
				if (0 < materialImageSize[mc + _MATERIAL_CHANNEL_REFLECTION].x)
					localReflectance = Get2dTableValue3(&textures[materialImageStart[mc + _MATERIAL_CHANNEL_REFLECTION]], materialImageSize[mc + _MATERIAL_CHANNEL_REFLECTION], uv[0], uv[1], uv[2], closestTriangleABL, closestTriangleACL);
				if (0 < materialImageSize[mc + _MATERIAL_CHANNEL_LUMINANCE].x)
					localLuminance = Get2dTableValue3(&textures[materialImageStart[mc + _MATERIAL_CHANNEL_LUMINANCE]], materialImageSize[mc + _MATERIAL_CHANNEL_LUMINANCE], uv[0], uv[1], uv[2], closestTriangleABL, closestTriangleACL);
			}

			for (j = 0; j < lightCount; ++j) {
				float3 toLightVectorUnit = {0.f, 0.f, 0.f}, attenuation = { 1.f, 1.f, 1.f };
				float toLightVectorMinLength = 0.f, toLightVectorMaxLength = 0.f;
				switch (lightType[j]) {
				case _LIGHT_TYPE_SPOT:
				case _LIGHT_TYPE_SPOTRECT:
				case _LIGHT_TYPE_TUBE:
				case _LIGHT_TYPE_AREA:
				case _LIGHT_TYPE_PHOTOMETRIC: {
						float3 randomLocation;
						randomLocation = GetSpherePoint(&rndSeed, lightRadius[j]);
						toLightVectorUnit.x = randomLocation.x + lightPosition[j].x - closestLocation.x;
						toLightVectorUnit.y = randomLocation.y + lightPosition[j].y - closestLocation.y;
						toLightVectorUnit.z = randomLocation.z + lightPosition[j].z - closestLocation.z;
						toLightVectorMinLength = 0.f;
						toLightVectorMaxLength = (float)sqrt(dot(toLightVectorUnit, toLightVectorUnit));
						tmp = 1.f / toLightVectorMaxLength;
						toLightVectorUnit.x *= tmp;
						toLightVectorUnit.y *= tmp;
						toLightVectorUnit.z *= tmp;
					}
					break;
				case _LIGHT_TYPE_OMNI:
					toLightVectorMinLength = 0.f;
					toLightVectorMaxLength = 0.f;
					break;
				case _LIGHT_TYPE_DISTANT:
				case _LIGHT_TYPE_PARALLEL:
				case _LIGHT_TYPE_PARSPOT:
				case _LIGHT_TYPE_PARSPOTRECT: {
						float randomLocationDistance, lightVectorLenInv;
						randomLocationDistance = (float)(sin((lightRadius[j] / 2.f) * M_PI / 180.f) * sqrt(dot(lightDirection[j], lightDirection[j])));
						toLightVectorUnit = GetSpherePoint(&rndSeed, randomLocationDistance);
						toLightVectorUnit.x -= lightDirection[j].x;
						toLightVectorUnit.y -= lightDirection[j].y;
						toLightVectorUnit.z -= lightDirection[j].z;
						lightVectorLenInv = 1.f / (float)sqrt(dot(toLightVectorUnit, toLightVectorUnit));
						toLightVectorUnit.x *= lightVectorLenInv;
						toLightVectorUnit.y *= lightVectorLenInv;
						toLightVectorUnit.z *= lightVectorLenInv;
						toLightVectorMinLength = 0.f;
						toLightVectorMaxLength = INFINITY;
					}
					break;
				}
				if (toLightVectorMinLength < toLightVectorMaxLength) {
					while(true) {
						float mult, abL, acL;
						uint triangleIndex = RayIntersectsTriangles(closestLocation, toLightVectorUnit, toLightVectorMinLength, toLightVectorMaxLength, closestTriangleIndex, triangleVertexIndex, triangleNormal, triangleMaterialId, triangleUv, vertex, materialImageSize, materialImageStart, textures, cameraTopToBottomPixelSizeVector, cameraLeftToRightPixelSizeVector, cameraPixelSizeInv, axesDivCount, sceneBoxMin, scenePixelTriangleListStart, scenePixelTriangleList, &mult, &abL, &acL);
						if (triangleIndex == (uint)-1) break;
						{
							__global float2* uv = &triangleUv[3 * triangleIndex + 0];
							int m = triangleMaterialId[triangleIndex];
							int mc = _MATERIAL_CHANNEL_COUNT * m;
							float3 transp = { 0.f, 0.f, 0.f };
							if (0 <= m && 0 < materialImageSize[mc + _MATERIAL_CHANNEL_TRANSPARENCY].x)
								transp = Get2dTableValue3(&textures[materialImageStart[mc + _MATERIAL_CHANNEL_TRANSPARENCY]], materialImageSize[mc + _MATERIAL_CHANNEL_TRANSPARENCY], uv[0], uv[1], uv[2], abL, acL);
							attenuation.x *= transp.x;
							attenuation.y *= transp.y;
							attenuation.z *= transp.z;
							if (!(0.f < attenuation.x && 0.f < attenuation.y && 0.f < attenuation.z)) break;
							toLightVectorMinLength = mult;
						}
					}
				}
				{
					float lightEffect = (float)fabs(dot(closestTriangleNormal, toLightVectorUnit));
					int index = (int)(0.f <= dot(closestTriangleNormal, toLightVectorUnit));
					float att = (float)pow(0.5f, toLightVectorMaxLength / lightHalfAttenuationDistance[j]);
					float absLightEffect = lightEffect * (att == att ? att : 1.f);
					faceLights[index].x += (1.f - faceLights[index].x) * attenuation.x * absLightEffect * (float)lightColour[j].x;
					faceLights[index].y += (1.f - faceLights[index].y) * attenuation.y * absLightEffect * (float)lightColour[j].y;
					faceLights[index].z += (1.f - faceLights[index].z) * attenuation.z * absLightEffect * (float)lightColour[j].z;
				}
			}
					
			commulativeOutputColour.x += (1.f - commulativeOutputColour.x) * localLuminance.x * rayMultiplier[cursorBegin].x;
			commulativeOutputColour.y += (1.f - commulativeOutputColour.y) * localLuminance.y * rayMultiplier[cursorBegin].y;
			commulativeOutputColour.z += (1.f - commulativeOutputColour.z) * localLuminance.z * rayMultiplier[cursorBegin].z;

			{
				float3 multiplier;
				float3 localDiffuse = { 0.f, 0.f, 0.f };
				float total;
				int isFrontFacing = (int)(dot(closestTriangleNormal, rayVector[cursorBegin]) <= 0.f);
				localLight.x = faceLights[isFrontFacing].x;
				localLight.y = faceLights[isFrontFacing].y;
				localLight.z = faceLights[isFrontFacing].z;
				commulativeOutputColour.x += (1.f - commulativeOutputColour.x) * rayMultiplier[cursorBegin].x * (1.f - localTransparency.x) * localTexture.x * localLight.x;
				commulativeOutputColour.y += (1.f - commulativeOutputColour.y) * rayMultiplier[cursorBegin].y * (1.f - localTransparency.y) * localTexture.y * localLight.y;
				commulativeOutputColour.z += (1.f - commulativeOutputColour.z) * rayMultiplier[cursorBegin].z * (1.f - localTransparency.z) * localTexture.z * localLight.z;

				// non glossy opaque light reflection (diffuse)
				if (maxBounces[cursorBegin] <= 0) continue;

				total = max(max(localReflectance.x + localTransparency.x, localReflectance.y + localTransparency.y), localReflectance.z + localTransparency.z);
				if (total < 1.f) {
					localDiffuse.x = 1.f - total;
					localDiffuse.y = 1.f - total;
					localDiffuse.z = 1.f - total;
				}
				multiplier.x = rayMultiplier[cursorBegin].x * localTexture.x * localDiffuse.x;
				multiplier.y = rayMultiplier[cursorBegin].y * localTexture.y * localDiffuse.y;
				multiplier.z = rayMultiplier[cursorBegin].z * localTexture.z * localDiffuse.z;
				if (3.f / 256.f <= multiplier.x + multiplier.y + multiplier.z) {
					maxBounces[cursorEnd] = 0;
					excludedTriangleIndex[cursorEnd] = closestTriangleIndex;
					rayOrigin[cursorEnd] = closestLocation;
					rayVector[cursorEnd] = GetSpherePoint(&rndSeed, 1.f);
					if (isFrontFacing != (0 <= dot(rayVector[cursorEnd], closestTriangleNormal))) {
						rayVector[cursorEnd].x = -rayVector[cursorEnd].x;
						rayVector[cursorEnd].y = -rayVector[cursorEnd].y;
						rayVector[cursorEnd].z = -rayVector[cursorEnd].z;
					}
					rayMultiplier[cursorEnd] = multiplier;
					fromCamera[cursorEnd] = false;
					minDistance[cursorEnd] = 0.f;
					maxDistance[cursorEnd] = INFINITY;
					cursorEnd = (cursorEnd + 1) % MAX_RAY_BUFFER_SIZE;
					if ((cursorEnd + 1) % MAX_RAY_BUFFER_SIZE == cursorBegin) continue;
				}


				// glossy opaque light reflection (mirror effect)
				multiplier.x = rayMultiplier[cursorBegin].x * localTexture.x * localReflectance.x;
				multiplier.y = rayMultiplier[cursorBegin].y * localTexture.y * localReflectance.y;
				multiplier.z = rayMultiplier[cursorBegin].z * localTexture.z * localReflectance.z;
				if (3.f / 256.f <= multiplier.x + multiplier.y + multiplier.z) {
					maxBounces[cursorEnd] = maxBounces[cursorBegin] - 1;
					excludedTriangleIndex[cursorEnd] = closestTriangleIndex;
					rayOrigin[cursorEnd] = closestLocation;
					// Compute the cosinus of the angle between the ray and the normal
					tmp = -2.f * dot(closestTriangleNormal, rayVector[cursorBegin]);
					rayVector[cursorEnd].x = rayVector[cursorBegin].x + tmp * closestTriangleNormal.x;
					rayVector[cursorEnd].y = rayVector[cursorBegin].y + tmp * closestTriangleNormal.y;
					rayVector[cursorEnd].z = rayVector[cursorBegin].z + tmp * closestTriangleNormal.z;
					rayMultiplier[cursorEnd] = multiplier;
					fromCamera[cursorEnd] = false;
					minDistance[cursorEnd] = 0.f;
					maxDistance[cursorEnd] = INFINITY;
					cursorEnd = (cursorEnd + 1) % MAX_RAY_BUFFER_SIZE;
					if ((cursorEnd + 1) % MAX_RAY_BUFFER_SIZE == cursorBegin) continue;
				}

				// glossy non opaque light reflection (glass effect)
				multiplier.x = rayMultiplier[cursorBegin].x * localTexture.x * localTransparency.x;
				multiplier.y = rayMultiplier[cursorBegin].y * localTexture.y * localTransparency.y;
				multiplier.z = rayMultiplier[cursorBegin].z * localTexture.z * localTransparency.z;
				if (3.f / 256.f <= multiplier.x + multiplier.y + multiplier.z) {
					maxBounces[cursorEnd] = maxBounces[cursorBegin] - 1;
					excludedTriangleIndex[cursorEnd] = closestTriangleIndex;
					rayOrigin[cursorEnd] = rayOrigin[cursorBegin];
					rayVector[cursorEnd] = rayVector[cursorBegin];
					rayMultiplier[cursorEnd] = multiplier;
					fromCamera[cursorEnd] = fromCamera[cursorBegin];
					minDistance[cursorEnd] = closestTriangleMult;
					maxDistance[cursorEnd] = INFINITY;
					cursorEnd = (cursorEnd + 1) % MAX_RAY_BUFFER_SIZE;
					if ((cursorEnd + 1) % MAX_RAY_BUFFER_SIZE == cursorBegin) continue;
				}
			}
		}
	}
	{
		int outputColor;
		tmp = (float)(0xFFFF) / (float)sampleCount;
		outputColor = (int)outputRed[nextPixel] + (int)(commulativeOutputColour.x * tmp);
		if (outputColor < 0) outputColor = 0;
		if (0xFFFF < outputColor) outputColor = 0xFFFF;
		outputRed[nextPixel] = (ushort)outputColor;
		outputColor = (int)outputGreen[nextPixel] + (int)(commulativeOutputColour.y * tmp);
		if (outputColor < 0) outputColor = 0;
		if (0xFFFF < outputColor) outputColor = 0xFFFF;
		outputGreen[nextPixel] = (ushort)outputColor;
		outputColor = (int)outputBlue[nextPixel] + (int)(commulativeOutputColour.z * tmp);
		if (outputColor < 0) outputColor = 0;
		if (0xFFFF < outputColor) outputColor = 0xFFFF;
		outputBlue[nextPixel] = (ushort)outputColor;
	}
}
