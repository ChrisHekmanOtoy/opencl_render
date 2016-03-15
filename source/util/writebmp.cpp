#include <cl/cl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "writebmp.h"

extern "C" {
#include "opencl/raytrace.h"
}

// Functions to export a bitmap to img.bmp
// http://stackoverflow.com/questions/2654480/writing-bmp-image-in-pure-c-c-without-other-libraries
cl_bool writebmp(cl_uint2 imageSize, cl_uchar3* outputImage) {
	int w = imageSize.s[0];
	int h = imageSize.s[1];
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3 * w * h;  //w is your image width, h is image height, both int
	img = (unsigned char *)malloc(3 * w * h);
	if (img) {
		memset(img, 0, 3 * w * h);

		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				int r = outputImage[i + j * imageSize.s[0]].s[0];
				int g = outputImage[i + j * imageSize.s[0]].s[1];
				int b = outputImage[i + j * imageSize.s[0]].s[2];
				img[(i + j*w) * 3 + 2] = (unsigned char)(r);
				img[(i + j*w) * 3 + 1] = (unsigned char)(g);
				img[(i + j*w) * 3 + 0] = (unsigned char)(b);
			}
		}

		unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
		unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
		unsigned char bmppad[3] = { 0, 0, 0 };

		bmpfileheader[2] = (unsigned char)(filesize);
		bmpfileheader[3] = (unsigned char)(filesize >> 8);
		bmpfileheader[4] = (unsigned char)(filesize >> 16);
		bmpfileheader[5] = (unsigned char)(filesize >> 24);

		bmpinfoheader[4] = (unsigned char)(w);
		bmpinfoheader[5] = (unsigned char)(w >> 8);
		bmpinfoheader[6] = (unsigned char)(w >> 16);
		bmpinfoheader[7] = (unsigned char)(w >> 24);
		bmpinfoheader[8] = (unsigned char)(h);
		bmpinfoheader[9] = (unsigned char)(h >> 8);
		bmpinfoheader[10] = (unsigned char)(h >> 16);
		bmpinfoheader[11] = (unsigned char)(h >> 24);

		if (0 == fopen_s(&f, "img.bmp", "wb")) {
			fwrite(bmpfileheader, 1, 14, f);
			fwrite(bmpinfoheader, 1, 40, f);
			for (int i = 0; i < h; i++) {
				fwrite(img + (w*(h - i - 1) * 3), 3, w, f);
				fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
			}
			fclose(f);
			free(img);
			return true;
		}
		free(img);
	}
	return false;
}

cl_bool writebmpf(cl_uint2 imageSize, cl_float3* outputImage) {
	int w = imageSize.s[0];
	int h = imageSize.s[1];
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3 * w*h;  //w is your image width, h is image height, both int
	img = (unsigned char *)malloc(3 * w * h);
	if (img) {
		memset(img, 0, 3 * w * h);

		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				int r = (int)floor(0.5f + bindf(outputImage[i + j * imageSize.s[0]].s[0], 0.f, 1.f) * 255.f);
				int g = (int)floor(0.5f + bindf(outputImage[i + j * imageSize.s[0]].s[1], 0.f, 1.f) * 255.f);
				int b = (int)floor(0.5f + bindf(outputImage[i + j * imageSize.s[0]].s[2], 0.f, 1.f) * 255.f);
				img[(i + j*w) * 3 + 2] = (unsigned char)(r);
				img[(i + j*w) * 3 + 1] = (unsigned char)(g);
				img[(i + j*w) * 3 + 0] = (unsigned char)(b);
			}
		}

		unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
		unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
		unsigned char bmppad[3] = { 0, 0, 0 };

		bmpfileheader[2] = (unsigned char)(filesize);
		bmpfileheader[3] = (unsigned char)(filesize >> 8);
		bmpfileheader[4] = (unsigned char)(filesize >> 16);
		bmpfileheader[5] = (unsigned char)(filesize >> 24);

		bmpinfoheader[4] = (unsigned char)(w);
		bmpinfoheader[5] = (unsigned char)(w >> 8);
		bmpinfoheader[6] = (unsigned char)(w >> 16);
		bmpinfoheader[7] = (unsigned char)(w >> 24);
		bmpinfoheader[8] = (unsigned char)(h);
		bmpinfoheader[9] = (unsigned char)(h >> 8);
		bmpinfoheader[10] = (unsigned char)(h >> 16);
		bmpinfoheader[11] = (unsigned char)(h >> 24);

		if (0 == fopen_s(&f, "img.bmp", "wb")) {
			fwrite(bmpfileheader, 1, 14, f);
			fwrite(bmpinfoheader, 1, 40, f);
			for (int i = 0; i < h; i++) {
				fwrite(img + (w*(h - i - 1) * 3), 3, w, f);
				fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
			}
			fclose(f);
			free(img);
			return true;
		}
		free(img);
	}
	return false;
}

cl_bool writebmp3s(cl_uint2 imageSize, cl_ushort* outputRed, cl_ushort* outputGreen, cl_ushort* outputBlue) {
	int w = imageSize.s[0];
	int h = imageSize.s[1];
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3 * w*h;  //w is your image width, h is image height, both int
	img = (unsigned char *)malloc(3 * w * h);
	if (img) {
		memset(img, 0, 3 * w * h);

		for (int i = 0; i < w; i++) {
			for (int j = 0; j < h; j++) {
				int r = outputRed[i + j * imageSize.s[0]];
				int g = outputGreen[i + j * imageSize.s[0]];
				int b = outputBlue[i + j * imageSize.s[0]];
				img[(i + j*w) * 3 + 2] = (unsigned char)(r);
				img[(i + j*w) * 3 + 1] = (unsigned char)(g);
				img[(i + j*w) * 3 + 0] = (unsigned char)(b);
			}
		}

		unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
		unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
		unsigned char bmppad[3] = { 0, 0, 0 };

		bmpfileheader[2] = (unsigned char)(filesize);
		bmpfileheader[3] = (unsigned char)(filesize >> 8);
		bmpfileheader[4] = (unsigned char)(filesize >> 16);
		bmpfileheader[5] = (unsigned char)(filesize >> 24);

		bmpinfoheader[4] = (unsigned char)(w);
		bmpinfoheader[5] = (unsigned char)(w >> 8);
		bmpinfoheader[6] = (unsigned char)(w >> 16);
		bmpinfoheader[7] = (unsigned char)(w >> 24);
		bmpinfoheader[8] = (unsigned char)(h);
		bmpinfoheader[9] = (unsigned char)(h >> 8);
		bmpinfoheader[10] = (unsigned char)(h >> 16);
		bmpinfoheader[11] = (unsigned char)(h >> 24);

		if (0 == fopen_s(&f, "img.bmp", "wb")) {
			fwrite(bmpfileheader, 1, 14, f);
			fwrite(bmpinfoheader, 1, 40, f);
			for (int i = 0; i < h; i++) {
				fwrite(img + (w*(h - i - 1) * 3), 3, w, f);
				fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
			}
			fclose(f);
			free(img);
			return true;
		}
		free(img);
	}
	return false;
}
