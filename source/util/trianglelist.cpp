#pragma warning(disable : 4756)

#include <map>
#include <set>
#include <stdlib.h>

#include "c4d.h"
#include "trianglelist.h"

char debug_symbol(cl_uint s) {
	if (0 <= s && s <= 9) {
		return (char)('0' + s);
	}
	if (10 <= s && s <= 35) {
		return (char)('a' + s - 10);
	}
	if (36 <= s && s <= 61) {
		return (char)('A' + s - 36);
	}
	return '#';
}

// From http://www.algolist.net/Algorithms/Sorting/Quicksort
template<typename SortedType>
void QuickSort(SortedType* arr, ptrdiff_t  left, ptrdiff_t  right) {
	ptrdiff_t  i = left, j = right;
	SortedType tmp;
	SortedType pivot = arr[(left + right) / 2];

	/* partition */
	while (i <= j) {
		while (arr[i] < pivot)
			i++;
		while (arr[j] > pivot)
			j--;
		if (i <= j) {
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;
			i++;
			j--;
		}
	};

	/* recursion */
	if (left < j)
		QuickSort(arr, left, j);
	if (i < right)
		QuickSort(arr, i, right);
}

// Get 2 vectors that are orthogonal to the given vector and between themselves
void GetXYUnitOrthogonal(cl_float3 vector, cl_float3* xUnit, cl_float3* yUnit) {
	cl_float3 normVector = normalize(vector);
	cl_float3 zBase = { 0.f, 0.f, 1.f };
	if (1.f <= dot(normVector, zBase)) zBase.s[1] = zBase.s[2]--;
	*xUnit = normalize(cross(normVector, zBase));
	*yUnit = normalize(cross(normVector, *xUnit));
}

// Not used anymore
cl_float2 Get2dPosition(cl_float3 eye, cl_float3 leftToRightPixelSizeVector, cl_float3 topToBottomPixelSizeVector, cl_float2 pixelSizeInv, cl_float3 v) {
	cl_float3 eyeToV;
	eyeToV.s[0] = v.s[0] - eye.s[0];
	eyeToV.s[1] = v.s[1] - eye.s[1];
	eyeToV.s[2] = v.s[2] - eye.s[2];
	cl_float2 screenPosition;
	screenPosition.s[0] = dot(leftToRightPixelSizeVector, eyeToV) * pixelSizeInv.s[0] * pixelSizeInv.s[0];
	screenPosition.s[1] = dot(topToBottomPixelSizeVector, eyeToV) * pixelSizeInv.s[1] * pixelSizeInv.s[1];
	return screenPosition;
}

// Get location of a 3D point in the frustum, projected on the 2D base of the pyramid
cl_float2 GetCameraPosition(cl_float cameraPixelSizeInv, cl_float3 cameraEye, cl_float3 cameraEyeToTopLeftVector, cl_float3 cameraLeftToRightPixelSizeVector, cl_float3 cameraTopToBottomPixelSizeVector, cl_float3 v) {
	cl_float pixelSizeInvSq = cameraPixelSizeInv * cameraPixelSizeInv;
	cl_float3 eyeToV;
	eyeToV.s[0] = v.s[0] - cameraEye.s[0];
	eyeToV.s[1] = v.s[1] - cameraEye.s[1];
	eyeToV.s[2] = v.s[2] - cameraEye.s[2];
	cl_float3 screenNormal = cross(cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector);
	cl_float eyeToVscale = dot(cameraEyeToTopLeftVector, screenNormal) / dot(eyeToV, screenNormal);
	cl_float3 topLeftToV;
	topLeftToV.s[0] = eyeToVscale * eyeToV.s[0] - cameraEyeToTopLeftVector.s[0];
	topLeftToV.s[1] = eyeToVscale * eyeToV.s[1] - cameraEyeToTopLeftVector.s[1];
	topLeftToV.s[2] = eyeToVscale * eyeToV.s[2] - cameraEyeToTopLeftVector.s[2];
	cl_float2 screenPosition;
	screenPosition.s[0] = dot(cameraLeftToRightPixelSizeVector, topLeftToV) * pixelSizeInvSq;
	screenPosition.s[1] = dot(cameraTopToBottomPixelSizeVector, topLeftToV) * pixelSizeInvSq;
	return screenPosition;
}

// Distance from a point to a 3D triangle squared
cl_float GetPointToTriangleSqLen(cl_float3 a, cl_float3 b, cl_float3 c, cl_float3 point) {
	cl_float3 ab, bc, ca, normal, apoint, bpoint, cpoint;
	cl_float dotabpoint, dotbcpoint, dotcapoint, sqLen;
	ab.s[0] = b.s[0] - a.s[0];
	ab.s[1] = b.s[1] - a.s[1];
	ab.s[2] = b.s[2] - a.s[2];
	bc.s[0] = c.s[0] - b.s[0];
	bc.s[1] = c.s[1] - b.s[1];
	bc.s[2] = c.s[2] - b.s[2];
	ca.s[0] = a.s[0] - c.s[0];
	ca.s[1] = a.s[1] - c.s[1];
	ca.s[2] = a.s[2] - c.s[2];
	normal = cross(ab, bc);
	apoint.s[0] = point.s[0] - a.s[0];
	apoint.s[1] = point.s[1] - a.s[1];
	apoint.s[2] = point.s[2] - a.s[2];
	bpoint.s[0] = point.s[0] - b.s[0];
	bpoint.s[1] = point.s[1] - b.s[1];
	bpoint.s[2] = point.s[2] - b.s[2];
	cpoint.s[0] = point.s[0] - c.s[0];
	cpoint.s[1] = point.s[1] - c.s[1];
	cpoint.s[2] = point.s[2] - c.s[2];
	dotabpoint = dot(apoint, cross(ab, normal));
	dotbcpoint = dot(bpoint, cross(bc, normal));
	dotcapoint = dot(cpoint, cross(ca, normal));
	if (0.f <= (dotabpoint * dotbcpoint) && 0.f <= (dotbcpoint * dotcapoint)) {
		sqLen = dot(apoint, normal) / (cl_float)sqrt(dot(normal, normal));
	}
	else {
		cl_float abp = GetPointToLineSqLen(a, b, point);
		cl_float bcp = GetPointToLineSqLen(b, c, point);
		cl_float cap = GetPointToLineSqLen(c, a, point);
		sqLen = min(abp, min(bcp, cap));
	}
	return sqLen;
}

// FillRectangle is used to rasterize a triangle in an array to form camera pixel triangle lists
size_t FillRectangle(cl_uint2 imageDimension, cl_ulong* trianglePairs, cl_float2 aPos, cl_float2 bPos, cl_float2 cPos, cl_uint triangleId, cl_uint triangleCount) {
	size_t trianglePairCursor = 0;

	// Precompute vectors between a, b and c
	cl_float2 ab, bc, ca;
	ab.s[0] = bPos.s[0] - aPos.s[0];
	ab.s[1] = bPos.s[1] - aPos.s[1];
	bc.s[0] = cPos.s[0] - bPos.s[0];
	bc.s[1] = cPos.s[1] - bPos.s[1];
	ca.s[0] = aPos.s[0] - cPos.s[0];
	ca.s[1] = aPos.s[1] - cPos.s[1];

	// Precompute vector slopes N.B. divs by 0 will fail the checks by design
	cl_float2 abSlope, bcSlope, caSlope;
	abSlope.s[0] = ab.s[0] / ab.s[1];
	abSlope.s[1] = 1.f / abSlope.s[0];
	bcSlope.s[0] = bc.s[0] / bc.s[1];
	bcSlope.s[1] = 1.f / bcSlope.s[0];
	caSlope.s[0] = ca.s[0] / ca.s[1];
	caSlope.s[1] = 1.f / caSlope.s[0];

	// Precompute the minimum square containing a, b and c
	cl_uint2 tmin, tmax;
	tmin.s[0] = (cl_uint)fmax(0.f, fmin(fmin(aPos.s[0], bPos.s[0]), fmin(cPos.s[0], (cl_float)(imageDimension.s[0] - 1))));
	tmin.s[1] = (cl_uint)fmax(0.f, fmin(fmin(aPos.s[1], bPos.s[1]), fmin(cPos.s[1], (cl_float)(imageDimension.s[1] - 1))));
	tmax.s[0] = (cl_uint)fmin((cl_float)(imageDimension.s[0] - 1), fmax(fmax(aPos.s[0], bPos.s[0]), fmax(cPos.s[0], 0.f)));
	tmax.s[1] = (cl_uint)fmin((cl_float)(imageDimension.s[1] - 1), fmax(fmax(aPos.s[1], bPos.s[1]), fmax(cPos.s[1], 0.f)));

	// Fill the pixels containing a if it is inside the screen - no need to check for b or c, we just need to handle the case where the 3 points are inside the pixel
	if (0.f <= aPos.s[0] && aPos.s[0] < (cl_float)imageDimension.s[0] && 0.f <= aPos.s[1] && aPos.s[1] < (cl_float)imageDimension.s[1]) {
		trianglePairs[trianglePairCursor++] = ((((cl_ulong)floor(aPos.s[0])) + ((cl_ulong)floor(aPos.s[1])) * (cl_ulong)imageDimension.s[0]) * (cl_ulong)triangleCount + (cl_ulong)triangleId);
	}

	// Scan the whole square
	for (cl_uint x = tmin.s[0]; x <= tmax.s[0]; ++x) {
		for (cl_uint y = tmin.s[1]; y <= tmax.s[1]; ++y) {
			if (x != (cl_uint)floor(aPos.s[0]) || y != (cl_uint)floor(aPos.s[1])) {
				// Check if the triangle edges intersect the pixel edges
				cl_float4 abIntersect, bcIntersect, caIntersect;
				abIntersect.s[0] = aPos.s[0] + (y - aPos.s[1]) * abSlope.s[0];
				abIntersect.s[1] = aPos.s[1] + (x - aPos.s[0]) * abSlope.s[1];
				abIntersect.s[2] = abIntersect.s[0] + abSlope.s[0];
				abIntersect.s[3] = abIntersect.s[1] + abSlope.s[1];
				bcIntersect.s[0] = bPos.s[0] + (y - bPos.s[1]) * bcSlope.s[0];
				bcIntersect.s[1] = bPos.s[1] + (x - bPos.s[0]) * bcSlope.s[1];
				bcIntersect.s[2] = bcIntersect.s[0] + bcSlope.s[0];
				bcIntersect.s[3] = bcIntersect.s[1] + bcSlope.s[1];
				caIntersect.s[0] = cPos.s[0] + (y - cPos.s[1]) * caSlope.s[0];
				caIntersect.s[1] = cPos.s[1] + (x - cPos.s[0]) * caSlope.s[1];
				caIntersect.s[2] = caIntersect.s[0] + caSlope.s[0];
				caIntersect.s[3] = caIntersect.s[1] + caSlope.s[1];
				if ((0.f <= (aPos.s[0] - abIntersect.s[0]) * (abIntersect.s[0] - bPos.s[0])) & (x == (cl_uint)(abIntersect.s[0])) |
					(0.f <= (aPos.s[0] - abIntersect.s[2]) * (abIntersect.s[2] - bPos.s[0])) & (x == (cl_uint)(abIntersect.s[2])) |
					(0.f <= (aPos.s[1] - abIntersect.s[1]) * (abIntersect.s[1] - bPos.s[1])) & (y == (cl_uint)(abIntersect.s[1])) |
					(0.f <= (aPos.s[1] - abIntersect.s[3]) * (abIntersect.s[3] - bPos.s[1])) & (y == (cl_uint)(abIntersect.s[3])) |
					(0.f <= (bPos.s[0] - bcIntersect.s[0]) * (bcIntersect.s[0] - cPos.s[0])) & (x == (cl_uint)(bcIntersect.s[0])) |
					(0.f <= (bPos.s[0] - bcIntersect.s[2]) * (bcIntersect.s[2] - cPos.s[0])) & (x == (cl_uint)(bcIntersect.s[2])) |
					(0.f <= (bPos.s[1] - bcIntersect.s[1]) * (bcIntersect.s[1] - cPos.s[1])) & (y == (cl_uint)(bcIntersect.s[1])) |
					(0.f <= (bPos.s[1] - bcIntersect.s[3]) * (bcIntersect.s[3] - cPos.s[1])) & (y == (cl_uint)(bcIntersect.s[3])) |
					(0.f <= (cPos.s[0] - caIntersect.s[0]) * (caIntersect.s[0] - aPos.s[0])) & (x == (cl_uint)(caIntersect.s[0])) |
					(0.f <= (cPos.s[0] - caIntersect.s[2]) * (caIntersect.s[2] - aPos.s[0])) & (x == (cl_uint)(caIntersect.s[2])) |
					(0.f <= (cPos.s[1] - caIntersect.s[1]) * (caIntersect.s[1] - aPos.s[1])) & (y == (cl_uint)(caIntersect.s[1])) |
					(0.f <= (cPos.s[1] - caIntersect.s[3]) * (caIntersect.s[3] - aPos.s[1])) & (y == (cl_uint)(caIntersect.s[3]))) {
					trianglePairs[trianglePairCursor++] = ((((cl_ulong)x) + ((cl_ulong)y) * (cl_ulong)imageDimension.s[0]) * (cl_ulong)triangleCount + (cl_ulong)triangleId);
				}
				// Check if the pixel is completely inside the triangle
				else {
					cl_float2 axy, bxy, cxy;
					axy.s[0] = (cl_float)x - aPos.s[0];
					axy.s[1] = (cl_float)y - aPos.s[1];
					bxy.s[0] = (cl_float)x - bPos.s[0];
					bxy.s[1] = (cl_float)y - bPos.s[1];
					cxy.s[0] = (cl_float)x - cPos.s[0];
					cxy.s[1] = (cl_float)y - cPos.s[1];
					cl_float abaxyCross, bcbxyCross, cacxyCross;
					abaxyCross = ab.s[0] * axy.s[1] - ab.s[1] * axy.s[0];
					bcbxyCross = bc.s[0] * bxy.s[1] - bc.s[1] * bxy.s[0];
					cacxyCross = ca.s[0] * cxy.s[1] - ca.s[1] * cxy.s[0];
					if ((0 <= abaxyCross * bcbxyCross) & (0 <= bcbxyCross * cacxyCross)) {
						trianglePairs[trianglePairCursor++] = ((((cl_ulong)x) + ((cl_ulong)y) * (cl_ulong)imageDimension.s[0]) * (cl_ulong)triangleCount + (cl_ulong)triangleId);
					}
				}
			}
		}
	}
	return trianglePairCursor;
}


// FillSphere used to be used to cast rays from light, but since there are rarely pin lights, the blur effect cannot be done with this technique. Not used anymore.
cl_uint FillSphere(cl_uint longitudePixel, cl_uint latitudePixel, cl_uint2 imageDimension, cl_ulong* pixelTriangles, cl_float3 eye, cl_float3 a, cl_float3 b, cl_float3 c, cl_uint triangleId, cl_uint triangleCount) {
	std::set<std::pair<cl_uint, cl_uint>> fillPixel;
	std::set<std::pair<cl_uint, cl_uint>> explorePixel;
	fillPixel.insert(std::make_pair(longitudePixel, latitudePixel));
	while (0 < explorePixel.size()) {
		std::set<std::pair<cl_uint, cl_uint>>::const_iterator it = explorePixel.cbegin();
		longitudePixel = it->first;
		latitudePixel = it->second;
		explorePixel.erase(it);

		DebugAssert(0 <= longitudePixel && longitudePixel < imageDimension.s[0], "pixel out of range, error");
		DebugAssert(0 <= latitudePixel && latitudePixel < imageDimension.s[1], "pixel out of range, error");

		cl_bool exploreEast = CL_FALSE, exploreSouth = CL_FALSE, exploreWest = CL_FALSE, exploreNorth = CL_FALSE;

		// Precompute vectors between a, b and c
		cl_float3 ab, bc, ca;
		ab.s[0] = b.s[0] - a.s[0];
		ab.s[1] = b.s[1] - a.s[1];
		ab.s[2] = b.s[2] - a.s[2];
		bc.s[0] = c.s[0] - b.s[0];
		bc.s[1] = c.s[1] - b.s[1];
		bc.s[2] = c.s[2] - b.s[2];
		ca.s[0] = a.s[0] - c.s[0];
		ca.s[1] = a.s[1] - c.s[1];
		ca.s[2] = a.s[2] - c.s[2];

		fillPixel.insert(std::make_pair(longitudePixel, latitudePixel));
		cl_float longitudeEast = -M_PI + longitudePixel * (2 * M_PI) / ((cl_float)imageDimension.s[0]);
		cl_float latitudeSouth = -M_PI / 2 + latitudePixel * M_PI / (cl_float)imageDimension.s[1];
		cl_float longitudeWest = -M_PI + (longitudePixel + 1) * (2 * M_PI) / ((cl_float)imageDimension.s[0]);
		cl_float latitudeNorth = -M_PI / 2 + (latitudePixel + 1) * M_PI / (cl_float)imageDimension.s[1];
		cl_float3 raySouthEast, rayNorthEast, raySouthWest, rayNorthWest;
		cl_float cosLatitudeSouth = cos(latitudeSouth);
		cl_float sinLatitudeSouth = sin(latitudeSouth);
		cl_float cosLatitudeNorth = cos(latitudeNorth);
		cl_float sinLatitudeNorth = sin(latitudeNorth);
		cl_float cosLongitudeEast = cos(longitudeEast);
		cl_float sinLongitudeEast = sin(longitudeEast);
		cl_float cosLongitudeWest = cos(longitudeWest);
		cl_float sinLongitudeWest = sin(longitudeWest);
		raySouthEast.s[0] = cosLatitudeSouth * cosLongitudeEast;
		raySouthEast.s[1] = cosLatitudeSouth * sinLongitudeEast;
		raySouthEast.s[2] = sinLatitudeSouth;
		rayNorthEast.s[0] = cosLatitudeNorth * cosLongitudeEast;
		rayNorthEast.s[1] = cosLatitudeNorth * sinLongitudeEast;
		rayNorthEast.s[2] = sinLatitudeNorth;
		raySouthWest.s[0] = cosLatitudeSouth * cosLongitudeWest;
		raySouthWest.s[1] = cosLatitudeSouth * sinLongitudeWest;
		raySouthWest.s[2] = sinLatitudeSouth;
		rayNorthWest.s[0] = cosLatitudeNorth * cosLongitudeWest;
		rayNorthWest.s[1] = cosLatitudeNorth * sinLongitudeWest;
		rayNorthWest.s[2] = sinLatitudeNorth;

		cl_float outRayMult, outABL, outACL;
		cl_float3 pastProjSouthEast, pastProjNorthEast, pastProjSouthWest, pastProjNorthWest;

		if (RayIntersectsTriangle(eye, raySouthEast, 0.f, INFINITY, a, b, c, &outRayMult, &outABL, &outACL)) {
			exploreEast = CL_TRUE;
			exploreSouth = CL_TRUE;
		}
		pastProjSouthEast.s[0] = eye.s[0] + 2 * raySouthEast.s[0] * outRayMult;
		pastProjSouthEast.s[1] = eye.s[1] + 2 * raySouthEast.s[1] * outRayMult;
		pastProjSouthEast.s[2] = eye.s[2] + 2 * raySouthEast.s[2] * outRayMult;

		if (RayIntersectsTriangle(eye, rayNorthEast, 0.f, INFINITY, a, b, c, &outRayMult, &outABL, &outACL)) {
			exploreEast = CL_TRUE;
			exploreNorth = CL_TRUE;
		}
		pastProjNorthEast.s[0] = eye.s[0] + 2 * rayNorthEast.s[0] * outRayMult;
		pastProjNorthEast.s[1] = eye.s[1] + 2 * rayNorthEast.s[1] * outRayMult;
		pastProjNorthEast.s[2] = eye.s[2] + 2 * rayNorthEast.s[2] * outRayMult;

		if (RayIntersectsTriangle(eye, raySouthWest, 0.f, INFINITY, a, b, c, &outRayMult, &outABL, &outACL)) {
			exploreWest = CL_TRUE;
			exploreSouth = CL_TRUE;
		}
		pastProjSouthWest.s[0] = eye.s[0] + 2 * raySouthWest.s[0] * outRayMult;
		pastProjSouthWest.s[1] = eye.s[1] + 2 * raySouthWest.s[1] * outRayMult;
		pastProjSouthWest.s[2] = eye.s[2] + 2 * raySouthWest.s[2] * outRayMult;

		if (RayIntersectsTriangle(eye, rayNorthWest, 0.f, INFINITY, a, b, c, &outRayMult, &outABL, &outACL)) {
			exploreWest = CL_TRUE;
			exploreNorth = CL_TRUE;
		}
		pastProjNorthWest.s[0] = eye.s[0] + 2 * rayNorthWest.s[0] * outRayMult;
		pastProjNorthWest.s[1] = eye.s[1] + 2 * rayNorthWest.s[1] * outRayMult;
		pastProjNorthWest.s[2] = eye.s[2] + 2 * rayNorthWest.s[2] * outRayMult;

		cl_uint indexEast = (longitudePixel - 1 + imageDimension.s[0]) % imageDimension.s[0] + latitudePixel * imageDimension.s[0];
		if (fillPixel.cend() == fillPixel.find(std::make_pair(indexEast, triangleId))) {
			if (!exploreEast) {
				exploreEast = RayIntersectsTriangle(a, ab, 0.f, 1.f, eye, pastProjSouthEast, pastProjNorthEast, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(b, bc, 0.f, 1.f, eye, pastProjSouthEast, pastProjNorthEast, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(c, ca, 0.f, 1.f, eye, pastProjSouthEast, pastProjNorthEast, &outRayMult, &outABL, &outACL);
			}
			if (exploreEast) {
				fillPixel.insert(std::make_pair((longitudePixel - 1 + imageDimension.s[0]) % imageDimension.s[0], latitudePixel));
			}
		}

		cl_uint indexSouth = longitudePixel + (latitudePixel - 1) * imageDimension.s[0];
		if (0 < latitudePixel && fillPixel.cend() == fillPixel.find(std::make_pair(indexSouth, triangleId))) {
			if (!exploreSouth) {
				exploreSouth = RayIntersectsTriangle(a, ab, 0.f, 1.f, eye, pastProjSouthEast, pastProjSouthWest, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(b, bc, 0.f, 1.f, eye, pastProjSouthEast, pastProjSouthWest, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(c, ca, 0.f, 1.f, eye, pastProjSouthEast, pastProjSouthWest, &outRayMult, &outABL, &outACL);
			}
			if (exploreSouth) {
				fillPixel.insert(std::make_pair(longitudePixel, (latitudePixel - 1 + imageDimension.s[1]) % imageDimension.s[1]));
			}
		}

		cl_uint indexWest = (longitudePixel + 1) % imageDimension.s[0] + latitudePixel * imageDimension.s[0];
		if (fillPixel.cend() == fillPixel.find(std::make_pair(indexWest, triangleId))) {
			if (!exploreWest) {
				exploreWest = RayIntersectsTriangle(a, ab, 0.f, 1.f, eye, pastProjSouthWest, pastProjNorthWest, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(b, bc, 0.f, 1.f, eye, pastProjSouthWest, pastProjNorthWest, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(c, ca, 0.f, 1.f, eye, pastProjSouthWest, pastProjNorthWest, &outRayMult, &outABL, &outACL);
			}
			if (exploreWest) {
				fillPixel.insert(std::make_pair((longitudePixel + 1) % imageDimension.s[0], latitudePixel));
			}
		}

		cl_uint indexNorth = longitudePixel + (latitudePixel + 1) * imageDimension.s[0];
		if ((latitudePixel + 1) < imageDimension.s[1] && fillPixel.cend() == fillPixel.find(std::make_pair(indexNorth, triangleId))) {
			if (!exploreNorth) {
				exploreNorth = RayIntersectsTriangle(a, ab, 0.f, 1.f, eye, pastProjNorthWest, pastProjNorthEast, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(b, bc, 0.f, 1.f, eye, pastProjNorthWest, pastProjNorthEast, &outRayMult, &outABL, &outACL) ||
					RayIntersectsTriangle(c, ca, 0.f, 1.f, eye, pastProjNorthWest, pastProjNorthEast, &outRayMult, &outABL, &outACL);
			}
			if (exploreNorth) {
				fillPixel.insert(std::make_pair(longitudePixel, latitudePixel + 1));
			}
		}
	}
	cl_uint trianglePairCursor = 0;
	for (std::set<std::pair<cl_uint, cl_uint>>::const_iterator it = fillPixel.cbegin(); it != fillPixel.cend(); ++it) {
		pixelTriangles[trianglePairCursor++] = ((cl_ulong)longitudePixel + (cl_ulong)latitudePixel * (cl_ulong)imageDimension.s[0]) * (cl_ulong)triangleCount + (cl_ulong)triangleId;
	}
	return trianglePairCursor;
}

// These functions are flattening a int3 to a ulong and reconstructing it back to int3
// a == ConvertToLong(ConvertToInt3(a))
// a == ConvertToInt3(ConvertToLong(a))
cl_ulong ConvertToLong(cl_int3 cubeId, cl_int axesDivision) {
	return ((cl_ulong)cubeId.s[0] + (cl_ulong)cubeId.s[1] * (cl_ulong)axesDivision + (cl_ulong)cubeId.s[2] * (cl_ulong)axesDivision * (cl_ulong)axesDivision);
}
cl_int3 ConvertToInt3(cl_ulong longId, cl_int axesDivision) {
	cl_int3 cubeId;
	cubeId.s[2] = (cl_int)(longId / ((cl_ulong)axesDivision * (cl_ulong)axesDivision));
	longId %= (cl_ulong)axesDivision * (cl_ulong)axesDivision;
	cubeId.s[1] = (cl_int)(longId / axesDivision);
	cubeId.s[0] = (cl_int)(longId % axesDivision);
	return cubeId;
}

// Culls a polygon of *polygonSideCount sides on dimension dim, greater than limit if isMax is set, otherwise least than limit
cl_bool Cull(cl_bool isMax, cl_float limit, cl_int dim, cl_int* polygonSideCount, cl_float3* polygon) {
	cl_int i, j, k;
	cl_bool newPoint[16] = { CL_FALSE };
	for (i = 0; i < *polygonSideCount; ++i) {
		cl_int nexti = (i + 1) % *polygonSideCount;
		cl_float diffi = limit - polygon[i].s[dim];
		cl_float diffnexti = limit - polygon[nexti].s[dim];
		if (diffi * diffnexti < 0.f) {
			cl_float3 polygoninexti;
			cl_float polygoninextiPct;
			cl_int iplusone = i + 1;
			polygoninexti.s[0] = polygon[nexti].s[0] - polygon[i].s[0];
			polygoninexti.s[1] = polygon[nexti].s[1] - polygon[i].s[1];
			polygoninexti.s[2] = polygon[nexti].s[2] - polygon[i].s[2];
			polygoninextiPct = diffi / polygoninexti.s[dim];
			for (j = (*polygonSideCount)++; iplusone < j; --j) {
				polygon[j] = polygon[j - 1];
			}
			polygon[iplusone].s[0] = polygon[i].s[0] + polygoninextiPct * polygoninexti.s[0];
			polygon[iplusone].s[1] = polygon[i].s[1] + polygoninextiPct * polygoninexti.s[1];
			polygon[iplusone].s[2] = polygon[i].s[2] + polygoninextiPct * polygoninexti.s[2];
			newPoint[iplusone] = CL_TRUE;
			i = iplusone;
		}
	}
	if (isMax) {
		for (cl_int i = 0; i < *polygonSideCount; ++i) {
			if (!newPoint[i] && limit < polygon[i].s[dim]) {
				for (j = i-- + 1, k = (*polygonSideCount)--; j < k; ++j) {
					polygon[j - 1] = polygon[j];
					newPoint[j - 1] = newPoint[j];
				}
			}
		}
	}
	else {
		for (cl_int i = 0; i < *polygonSideCount; ++i) {
			if (!newPoint[i] && polygon[i].s[dim] < limit) {
				for (j = i-- + 1, k = (*polygonSideCount)--; j < k; ++j) {
					polygon[j - 1] = polygon[j];
					newPoint[j - 1] = newPoint[j];
				}
			}
		}
	}
	if (1 == *polygonSideCount || 2 == *polygonSideCount) {
		*polygonSideCount = *polygonSideCount;
	}
	return 0 < *polygonSideCount;
}

// Check if an orthogonal box intersects a 3D triangle
cl_bool BoxIntersectsTriangle(cl_float3 min, cl_float3 max, cl_float3 a, cl_float3 b, cl_float3 c) {
	cl_bool returnValue = CL_FALSE;
	cl_int polygonSideCount = 3;
	cl_float3 polygon[16];
	polygon[0] = a;
	polygon[1] = b;
	polygon[2] = c;
	returnValue =
		((Cull(CL_FALSE, min.s[0], 0, &polygonSideCount, polygon)) &&
		 (Cull(CL_FALSE, min.s[1], 1, &polygonSideCount, polygon)) &&
		 (Cull(CL_FALSE, min.s[2], 2, &polygonSideCount, polygon)) &&
		 (Cull(CL_TRUE, max.s[0], 0, &polygonSideCount, polygon)) &&
		 (Cull(CL_TRUE, max.s[1], 1, &polygonSideCount, polygon)) &&
		 (Cull(CL_TRUE, max.s[2], 2, &polygonSideCount, polygon)));
	DebugAssert(polygonSideCount < 16, "Rounding errors could cause us to have a polygon with more than 9 sides. Margin of error exceeded: array overflow!");
	return returnValue;
}

// Rasterize a 3D triangle in a 3D grid of arbitrary sized boxes
size_t FillCube(cl_int axesDivision, cl_uchar* filledBoxBit, cl_float3* boxMin, cl_ulong* pixelTriangles, cl_float3 a, cl_float3 b, cl_float3 c, cl_uint triangleId, cl_uint triangleCount) {
	size_t trianglePairCursor = 0;
	size_t trianglePairCount = 0;
	cl_int3 cubeId = GetBoxAddress(axesDivision, boxMin, a);
	cl_ulong longId = ConvertToLong(cubeId, axesDivision);
	memset(filledBoxBit, 0, (axesDivision * axesDivision * axesDivision) / 8);
	filledBoxBit[longId / 8] |= (0x01 << (longId % 8));
	pixelTriangles[trianglePairCount++] = longId;
	while (trianglePairCursor < trianglePairCount) {
		longId = pixelTriangles[trianglePairCursor];
		pixelTriangles[trianglePairCursor++] = longId * (cl_ulong)triangleCount + (cl_ulong)triangleId;
		cubeId = ConvertToInt3(longId, axesDivision);

		DebugAssert(0 <= cubeId.s[0] && cubeId.s[0] < axesDivision, "pixel out of range, error");
		DebugAssert(0 <= cubeId.s[1] && cubeId.s[1] < axesDivision, "pixel out of range, error");
		DebugAssert(0 <= cubeId.s[2] && cubeId.s[2] < axesDivision, "pixel out of range, error");

		cl_float3 mini, maxi;
		for (cl_int i = 0; i < 3; ++i) {
			mini.s[i] = boxMin[cubeId.s[i]].s[i];
			maxi.s[i] = boxMin[cubeId.s[i] + 1].s[i];
		}
		DebugAssert(BoxIntersectsTriangle(mini, maxi, a, b, c), "ABC doesn't intersect the current box!");

		for (cl_int i = 0; i < 3; ++i) {
			for (cl_int j = -1; j <= 1; j += 2) {
				cubeId.s[i] += j;
				if (0 <= cubeId.s[i] && cubeId.s[i] < axesDivision) {
					cl_ulong longId = ConvertToLong(cubeId, axesDivision);
					cl_ulong longIdDiv8 = longId / 8;
					cl_ulong longIdMod8Char = (0x01 << (longId % 8));
					if (!(filledBoxBit[longIdDiv8] & longIdMod8Char)) {
						mini.s[i] = boxMin[cubeId.s[i]].s[i];
						maxi.s[i] = boxMin[cubeId.s[i] + 1].s[i];
						if (BoxIntersectsTriangle(mini, maxi, a, b, c)) {
							filledBoxBit[longIdDiv8] |= longIdMod8Char;
							pixelTriangles[trianglePairCount++] = longId;
						}
					}
				}
				cubeId.s[i] -= j;
			}
			mini.s[i] = boxMin[cubeId.s[i]].s[i];
			maxi.s[i] = boxMin[cubeId.s[i] + 1].s[i];
		}
	}
	DebugAssert(0 != (filledBoxBit[(ConvertToLong(GetBoxAddress(axesDivision, boxMin, a), axesDivision)) / 8] & (0x01 << ((ConvertToLong(GetBoxAddress(axesDivision, boxMin, a), axesDivision)) % 8))), "A's cube isn't filled!");
	DebugAssert(0 != (filledBoxBit[(ConvertToLong(GetBoxAddress(axesDivision, boxMin, b), axesDivision)) / 8] & (0x01 << ((ConvertToLong(GetBoxAddress(axesDivision, boxMin, b), axesDivision)) % 8))), "B's cube isn't filled!");
	DebugAssert(0 != (filledBoxBit[(ConvertToLong(GetBoxAddress(axesDivision, boxMin, c), axesDivision)) / 8] & (0x01 << ((ConvertToLong(GetBoxAddress(axesDivision, boxMin, c), axesDivision)) % 8))), "C's cube isn't filled!");

	return trianglePairCount;
}

CameraTriangleList::CameraTriangleList(cl_uint pixelCount, ptrdiff_t pixelTriangleCount) {
	triangleListStart = new cl_uint[pixelCount];
	memset(triangleListStart, 0, pixelCount * sizeof(cl_uint));
	triangleListEnd = new cl_uint[pixelCount];
	memset(triangleListEnd, 0, pixelCount * sizeof(cl_uint));
	triangleList = new cl_uint[pixelTriangleCount];
	triangleListSize = pixelTriangleCount;
}

CameraTriangleList::~CameraTriangleList() {
	if (triangleListStart) delete[] triangleListStart; triangleListStart = NULL;
	if (triangleListEnd) delete[] triangleListEnd; triangleListEnd = NULL;
	if (triangleList) delete[] triangleList; triangleList = NULL;
}

CameraTriangleList* CameraTriangleList::New(cl_uint2 cameraImageDimension, cl_float3 cameraEye, cl_float3 cameraEyeToTopLeftVector, cl_float3 cameraLeftToRightPixelSizeVector, cl_float3 cameraTopToBottomPixelSizeVector, cl_float cameraPixelSizeInv, cl_uint vertexCount, cl_uint triangleCount, cl_float3* vertex, cl_int3* triangleVertexIndex) {
	CameraTriangleList* cameraTriangleList = NULL;
	const ptrdiff_t  maxTriangles = 256 * 1024 * 1024;
	cl_ulong* cameraPixelTriangle = new cl_ulong[maxTriangles];
	ptrdiff_t  cameraPixelTriangleCursor = 0;
	cl_bool isPixelTriangleOverflowed = false;

	DebugAssert(cameraPixelTriangle);
	if (!cameraPixelTriangle) goto cleanup;

	/**************************************************************************************************************************
	* Note (idea)
	* These sets could be ordered lists instead if we assumed the triangles can have a partial order "closer" relationship.
	* But this is impossible since triangles can be symmetric when intersecting (both are closer than the other at a location).
	* Also and more importantly, the transitivity property is not repected in the case of 3 triangles t1, t2, t3 placed
	* strategically so that t1 <= t2, t2 <= t3 and t1 </= t3 (some kind of fan pattern).
	* The bug would be very infrequent and only on the pixel level, so not very important. This would speedup the ray tracing
	* a little bit (not tried).
	***************************************************************************************************************************/
	for (cl_uint triangleCursor = 0; triangleCursor < triangleCount; ++triangleCursor) {
		cl_float2 camPos0 = GetCameraPosition(cameraPixelSizeInv, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, vertex[triangleVertexIndex[triangleCursor].s[0]]);
		cl_float2 camPos1 = GetCameraPosition(cameraPixelSizeInv, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, vertex[triangleVertexIndex[triangleCursor].s[1]]);
		cl_float2 camPos2 = GetCameraPosition(cameraPixelSizeInv, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, vertex[triangleVertexIndex[triangleCursor].s[2]]);
		if (maxTriangles <= cameraPixelTriangleCursor + cameraImageDimension.s[0] * cameraImageDimension.s[1]) {
			isPixelTriangleOverflowed = true;
		}
		// TODO: Error will happen when a triangle is partially behind the camera, we need to clip the triangle into a quadrangle, pentangle, etc.
		cameraPixelTriangleCursor += FillRectangle(cameraImageDimension, &cameraPixelTriangle[isPixelTriangleOverflowed ? 0 : cameraPixelTriangleCursor], camPos0, camPos1, camPos2, triangleCursor, triangleCount);
	}

	// When the list is overflowed, we should stop recording the triangles and just count the triangles
	if (isPixelTriangleOverflowed) {
		delete[] cameraPixelTriangle; cameraPixelTriangle = NULL;
		cameraPixelTriangle = new cl_ulong[cameraPixelTriangleCursor];
		DebugAssert(cameraPixelTriangle);
		if (!cameraPixelTriangle) goto cleanup;
		cameraPixelTriangleCursor = 0;
		for (cl_uint triangleCursor = 0; triangleCursor < triangleCount; ++triangleCursor) {
			cl_float2 camPos0 = GetCameraPosition(cameraPixelSizeInv, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, vertex[triangleVertexIndex[triangleCursor].s[0]]);
			cl_float2 camPos1 = GetCameraPosition(cameraPixelSizeInv, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, vertex[triangleVertexIndex[triangleCursor].s[1]]);
			cl_float2 camPos2 = GetCameraPosition(cameraPixelSizeInv, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, vertex[triangleVertexIndex[triangleCursor].s[2]]);
			cameraPixelTriangleCursor += FillRectangle(cameraImageDimension, &cameraPixelTriangle[cameraPixelTriangleCursor], camPos0, camPos1, camPos2, triangleCursor, triangleCount);
		}
	}

	QuickSort<cl_ulong>(cameraPixelTriangle, 0, cameraPixelTriangleCursor - 1);
	cl_uint outputImageSize = cameraImageDimension.s[0] * cameraImageDimension.s[1];
	cameraTriangleList = new CameraTriangleList(outputImageSize, cameraPixelTriangleCursor);
	if (!cameraTriangleList) goto cleanup;
	for (ptrdiff_t triangleCursor = 0; triangleCursor < cameraPixelTriangleCursor; ++triangleCursor) {
		cl_uint triangleId = (cl_uint)(cameraPixelTriangle[triangleCursor] % triangleCount);
		cl_uint pixelId = (cl_uint)(cameraPixelTriangle[triangleCursor] / triangleCount);
		cameraTriangleList->triangleList[triangleCursor] = triangleId;
		++(cameraTriangleList->triangleListEnd[pixelId]);
	}
	for (cl_uint pixelCursor = 1; pixelCursor < outputImageSize; ++pixelCursor) {
		cameraTriangleList->triangleListStart[pixelCursor] = cameraTriangleList->triangleListEnd[pixelCursor - 1];
		cameraTriangleList->triangleListEnd[pixelCursor] += cameraTriangleList->triangleListStart[pixelCursor];
	}

	// Compress the triangle list by checking the previous pixels in x and y dimension. If their list is the same, copy the location.
	cl_uint compressionTriangleCount = 0;
	for (cl_uint pixelCursor = 0; pixelCursor < outputImageSize; ++pixelCursor) {
		cl_uint pixelX = pixelCursor % cameraImageDimension.s[0];
		cl_uint pixelY = pixelCursor / cameraImageDimension.s[0];
		cl_uint triangleListSize = cameraTriangleList->triangleListEnd[pixelCursor] - cameraTriangleList->triangleListStart[pixelCursor];
		memcpy(&cameraTriangleList->triangleList[cameraTriangleList->triangleListStart[pixelCursor] - compressionTriangleCount], &cameraTriangleList->triangleList[cameraTriangleList->triangleListStart[pixelCursor]], triangleListSize * sizeof(cl_uint));
		cameraTriangleList->triangleListStart[pixelCursor] -= compressionTriangleCount;
		cameraTriangleList->triangleListEnd[pixelCursor] -= compressionTriangleCount;
		cl_bool isCompressed = CL_FALSE;

		// Check if X-1 has the exact same list as this pixel
		if (0 < pixelX) {
			if (triangleListSize == cameraTriangleList->triangleListEnd[pixelCursor - 1] - cameraTriangleList->triangleListStart[pixelCursor - 1]) {
				if (0 == memcmp(&cameraTriangleList->triangleList[cameraTriangleList->triangleListStart[pixelCursor - 1]], &cameraTriangleList->triangleList[cameraTriangleList->triangleListStart[pixelCursor]], triangleListSize * sizeof(cl_uint))) {
					compressionTriangleCount += triangleListSize;
					cameraTriangleList->triangleListStart[pixelCursor] = cameraTriangleList->triangleListStart[pixelCursor - 1];
					cameraTriangleList->triangleListEnd[pixelCursor] = cameraTriangleList->triangleListEnd[pixelCursor - 1];
					isCompressed = CL_TRUE;
				}
			}
		}
		// Check if Y-1 has the exact same list as this pixel
		if (0 < pixelY && !isCompressed) {
			if (triangleListSize == cameraTriangleList->triangleListEnd[pixelCursor - cameraImageDimension.s[0]] - cameraTriangleList->triangleListStart[pixelCursor - cameraImageDimension.s[0]]) {
				if (0 == memcmp(&cameraTriangleList->triangleList[cameraTriangleList->triangleListStart[pixelCursor - cameraImageDimension.s[0]]], &cameraTriangleList->triangleList[cameraTriangleList->triangleListStart[pixelCursor]], triangleListSize * sizeof(cl_uint))) {
					compressionTriangleCount += triangleListSize;
					cameraTriangleList->triangleListStart[pixelCursor] = cameraTriangleList->triangleListStart[pixelCursor - cameraImageDimension.s[0]];
					cameraTriangleList->triangleListEnd[pixelCursor] = cameraTriangleList->triangleListEnd[pixelCursor - cameraImageDimension.s[0]];
				}
			}
		}
	}
	cameraTriangleList->triangleListSize -= compressionTriangleCount;

	//FILE* f = fopen("camera.txt", "wb+");
	//for (cl_uint i = 0; i < cameraImageDimension.s[1]; ++i) {
	//	for (cl_uint j = 0; j < cameraImageDimension.s[0]; ++j) {
	//		fprintf(f, "%c", debug_symbol(cameraTriangleList->triangleListEnd[j + i * cameraImageDimension.s[0]] - cameraTriangleList->triangleListStart[j + i * cameraImageDimension.s[0]]));
	//	}
	//	fprintf(f, "\n");
	//}
	//fclose(f);
cleanup:
	if (cameraPixelTriangle) delete[] cameraPixelTriangle; cameraPixelTriangle = NULL;
	return cameraTriangleList;
}

cl_uint* CameraTriangleList::GetTriangleListStart() {
	return triangleListStart;
}

cl_uint* CameraTriangleList::GetTriangleListEnd() {
	return triangleListEnd;
}

cl_uint* CameraTriangleList::GetTriangleList() {
	return triangleList;
}

ptrdiff_t CameraTriangleList::GetTriangleListSize() {
	return triangleListSize;
}

SceneTriangleList::SceneTriangleList(cl_uint pixelCountInc, ptrdiff_t pixelTriangleCount) {
	triangleListStart = new cl_uint[pixelCountInc];
	memset(triangleListStart, 0, pixelCountInc * sizeof(cl_uint));
	triangleList = new cl_uint[pixelTriangleCount];
}

SceneTriangleList::~SceneTriangleList() {
	if (triangleListStart) delete[] triangleListStart; triangleListStart = NULL;
	if (triangleList) delete[] triangleList; triangleList = NULL;
}

SceneTriangleList* SceneTriangleList::New(cl_uint vertexCount, cl_uint triangleCount, cl_float3* vertex, cl_int3* triangleVertexIndex) {
	SceneTriangleList* sceneTriangleList = NULL;
	cl_float3 boxMin[AXES_DIVISION + 1] = { 0.f };
	// Order every vertex in every dimension x, y, z
	// Split boxes at the location of between ordered vertex 0 and vertexCount/AXES_DIVISION, between vertex 2*vertexCount/AXES_DIVISION, ... between vertex (AXES_DIVISION-1)*vertexCount/AXES_DIVISION and vertexCount
	if (0 < vertexCount) {
		cl_float* val = new cl_float[vertexCount];
		if (val) {
			for (cl_int w = 0; w < 3; ++w) {
				for (cl_uint vertexCursor = 0; vertexCursor < vertexCount; ++vertexCursor) {
					val[vertexCursor] = vertex[vertexCursor].s[w];
				}
				QuickSort<cl_float>(val, 0, vertexCount - 1);
				for (cl_int i = 0; i < AXES_DIVISION + 1; ++i) {
					cl_uint index = (i * (vertexCount - 1)) / AXES_DIVISION;
					if (0 < index && index < vertexCount) {
						boxMin[i].s[w] = (val[index] + val[index - 1]) / 2.f;
					}
					else boxMin[i].s[w] = val[index];
				}
			}
			if (val) delete[] val; val = NULL;
		}
	}

	// Note, this constant is arbitrarily selected. We do not know yet how many triangles references will be kept. Counting them would double the time need for this.
	const ptrdiff_t  maxTriangles = (256 * 1024 * 1024);
	cl_ulong* scenePixelTriangle = new cl_ulong[maxTriangles];
	cl_uchar* filledBoxBit = new cl_uchar[(AXES_DIVISION * AXES_DIVISION * AXES_DIVISION) / 8];
	ptrdiff_t  scenePixelTriangleCursor = 0;
	cl_bool isPixelTriangleOverflowed = false;
	DebugAssert(scenePixelTriangle);
	DebugAssert(filledBoxBit);
	if (!scenePixelTriangle || !filledBoxBit) goto cleanup;
	for (cl_uint triangleCursor = 0; triangleCursor < triangleCount; ++triangleCursor) {
		cl_int3* vId = &triangleVertexIndex[triangleCursor];
		if (maxTriangles <= scenePixelTriangleCursor + AXES_DIVISION * AXES_DIVISION * AXES_DIVISION) {
			isPixelTriangleOverflowed = true;
		}
		scenePixelTriangleCursor += FillCube(AXES_DIVISION, filledBoxBit, boxMin, &scenePixelTriangle[isPixelTriangleOverflowed ? 0 : scenePixelTriangleCursor], vertex[vId->s[0]], vertex[vId->s[1]], vertex[vId->s[2]], triangleCursor, triangleCount);
	}
	if (isPixelTriangleOverflowed) {
		delete[] scenePixelTriangle; scenePixelTriangle = NULL;
		scenePixelTriangle = new cl_ulong[scenePixelTriangleCursor];
		DebugAssert(scenePixelTriangle);
		if (!scenePixelTriangle) goto cleanup;
		scenePixelTriangleCursor = 0;
		for (cl_uint triangleCursor = 0; triangleCursor < triangleCount; ++triangleCursor) {
			cl_int3* vId = &triangleVertexIndex[triangleCursor];
			scenePixelTriangleCursor += FillCube(AXES_DIVISION, filledBoxBit, boxMin, &scenePixelTriangle[scenePixelTriangleCursor], vertex[vId->s[0]], vertex[vId->s[1]], vertex[vId->s[2]], triangleCursor, triangleCount);
		}
	}
	QuickSort<cl_ulong>(scenePixelTriangle, 0, scenePixelTriangleCursor - 1);
	sceneTriangleList = new SceneTriangleList(AXES_DIVISION * AXES_DIVISION * AXES_DIVISION + 1, scenePixelTriangleCursor);
	if (sceneTriangleList) {
		memcpy(sceneTriangleList->boxMin, boxMin, sizeof(boxMin));
		for (ptrdiff_t triangleCursor = 0; triangleCursor < scenePixelTriangleCursor; ++triangleCursor) {
			cl_uint triangleId = (cl_uint)(scenePixelTriangle[triangleCursor] % triangleCount);
			cl_uint pixelId = (cl_uint)(scenePixelTriangle[triangleCursor] / triangleCount);
			sceneTriangleList->triangleList[triangleCursor] = triangleId;
			++(sceneTriangleList->triangleListStart[pixelId + 1]);
		}
		for (cl_uint pixelCursor = 1; pixelCursor <= AXES_DIVISION * AXES_DIVISION * AXES_DIVISION; ++pixelCursor) {
			sceneTriangleList->triangleListStart[pixelCursor] += sceneTriangleList->triangleListStart[pixelCursor - 1];
		}

		//FILE* f = fopen("scene.txt", "wb+");
		//for (cl_uint i = 0; i < AXES_DIVISION; ++i) {
		//	for (cl_uint j = 0; j < AXES_DIVISION; ++j) {
		//		for (cl_uint k = 0; k < AXES_DIVISION; ++k) {
		//			fprintf(f, "%c", debug_symbol(sceneTriangleList->triangleListStart[i + j * AXES_DIVISION + k * AXES_DIVISION * AXES_DIVISION + 1] - sceneTriangleList->triangleListStart[i + j * AXES_DIVISION + k * AXES_DIVISION * AXES_DIVISION]));
		//		}
		//		fprintf(f, "\n");
		//	}
		//	fprintf(f, "\n");
		//}
		//fclose(f);
	}
cleanup:
	if (scenePixelTriangle) delete[] scenePixelTriangle; scenePixelTriangle = NULL;
	if (filledBoxBit) delete[] filledBoxBit; filledBoxBit = NULL;
	return sceneTriangleList;
}

cl_float3* SceneTriangleList::GetBoxMin() {
	return boxMin;
}

cl_uint* SceneTriangleList::GetTriangleListStart() {
	return triangleListStart;
}

cl_uint* SceneTriangleList::GetTriangleList() {
	return triangleList;
}

