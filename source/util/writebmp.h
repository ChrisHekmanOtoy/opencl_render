#ifndef WRITEBMP_H
#define WRITEBMP_H

#include <cl/cl.h>

cl_bool writebmp(cl_uint2 imageSize, cl_uchar3* outputImage);
cl_bool writebmpf(cl_uint2 imageSize, cl_float3* outputImage);
cl_bool writebmp3s(cl_uint2 imageSize, cl_ushort* outputRed, cl_ushort* outputGreen, cl_ushort* outputBlue);

#endif // WRITEBMP_H
