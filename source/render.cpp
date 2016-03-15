#pragma warning(disable : 4756)

#include <map>
#include <time.h>

#include "c4d.h"
#include "c4d_reflection.h"
#include "c4d_symbols.h"
#include "c4d_memory_mp.h"
#include "main.h"
#include "util/trianglelist.h"
#include "util/writebmp.h"

extern "C" {
#include "opencl/raytrace.h"
}

// custom gui headers
#include "customgui_htmlviewer.h"
#include "customgui_priority.h"
#include "customgui_descproperty.h"
#include "customgui_hyperlink.h"
#include "customgui_inexclude.h"
#include "c4d_graphview.h"
#include "customgui_datetime.h"

enum {
	ID_TEXT_IMAGE_SIZE_X = 1100,
	ID_EDITNUMBER_IMAGE_SIZE_X,
	ID_TEXT_IMAGE_SIZE_Y,
	ID_EDITNUMBER_IMAGE_SIZE_Y,
	ID_TEXT_SAMPLING,
	ID_EDITNUMBER_SAMPLING,
	ID_COMBOBOX_PROCESSOR,
	ID_BUTTON_REFRESH_OPENCL,
	ID_TEXT_EMPTY,
	ID_BUTTON_RENDER,
	ID_TEXT_PREPARETIME_TITLE,
	ID_TEXT_PREPARETIME,
	ID_TEXT_PROGRESS_TITLE,
	ID_TEXT_PROGRESS,
	ID_TEXT_TIMELEFT_TITLE,
	ID_TEXT_TIMELEFT,
	ID_TEXT_TIMEPAST_TITLE,
	ID_TEXT_TIMEPAST,
};

Bool parseAndRender(Matrix cameraMatrix, cl_float fov, BaseDocument* doc, cl_uint platformDevice, cl_uint2 outputImageDimension, cl_uint sampleCount, BaseBitmap** bitmap);

class CInitOpenCL;

// This class is used as a thread to initialize OpenCL
// These functions will keep track of the progress:
//   void ResetComputationType();
//   cl_bool GetIsComputationTypeUpdated();
//   size_t GetComputationTypeCount();
//   cl_bool GetComputationTypeName(size_t id, size_t strLen, cl_char* str);

class CInitOpenCL : public C4DThread {
private:
	Char name[32];
public:
	CInitOpenCL() : C4DThread() {
		strcpy(name, "CInitOpenCL");
	}
	virtual ~CInitOpenCL() {}
	virtual void Main(void) {
		InitOpenCL();
	}
	virtual const Char* GetThreadName(void) { return name; }
};

// This class is used as a thread to run the ray tracing and return the image
// input:
//  cameraMatrix: a 4x3 matrix representing the position and rotation of the camera
//  fov: the angle in radian of the camera (horizontally)
//  doc: a polygonized version of the scene. Since this thread is going to read it, it must not change for the duration of the thread
//  platformDevice: id of the platform/device used to do the computation. 0 means locally on this thread with this CPU, x=[1, n] means run it on OpenCL platform/device x-1
//  outputImageDimension: size of the resulting (ALLOCATED INSIDE THE THREAD) bitmap
//  sampleCount: number of rays cast per pixel
// output:
//  bitmap: location to return the allocated bitmap. THIS BITMAP MUST BE DEALLOCATED AFTER USAGE!
//
// These functions will keep track of the progress:
//   cl_float GetProgress();
//   void SetProgress(cl_float p);
//   clock_t GetStartTime();
//   clock_t GetEndTime();
//   void ResetTime();
//
class ComputeOpenCL : public C4DThread {
private:
	Char name[32];
	Matrix cameraMatrix;
	cl_float fov;
	BaseDocument* doc;
	cl_uint platformDevice;
	cl_uint2 outputImageDimension;
	cl_uint sampleCount;
	BaseBitmap** bitmap;
public:
	ComputeOpenCL(Matrix cameraMatrix, cl_float fov, BaseDocument* doc, cl_uint platformDevice, cl_uint2 outputImageDimension, cl_uint sampleCount, BaseBitmap** bitmap) : C4DThread() {
		strcpy(name, "ComputeOpenCL");
		this->cameraMatrix = cameraMatrix;
		this->fov = fov;
		this->doc = doc;
		this->platformDevice = platformDevice;
		this->outputImageDimension = outputImageDimension;
		this->sampleCount = sampleCount;
		this->bitmap = bitmap;
	}
	virtual ~ComputeOpenCL() {}
	virtual void Main(void) {
		parseAndRender(cameraMatrix, fov, doc, platformDevice, outputImageDimension, sampleCount, bitmap);
	}
	virtual const Char* GetThreadName(void) { return name; }
};

// Dialog element that pops up when the plugin is activated
class OpenCLDialog : public GeDialog {
private:
	BaseDocument* doc;
	CInitOpenCL* initOpenCL;
	ComputeOpenCL* computeOpenCL;
	bool isOpenCLExpanded;
	clock_t prepareStart, computeStart, lastUpdate;
	cl_float lastProgress;
	BaseBitmap* bitmap;

public:
	OpenCLDialog();
	virtual ~OpenCLDialog(void);
	virtual Bool CreateLayout();
	virtual Bool Command(Int32 id, const BaseContainer & msg);
	virtual void Timer(const BaseContainer& msg);
	void SetDocument(BaseDocument* doc);

private:
	void UpdateDialog();
};


OpenCLDialog::OpenCLDialog(void) : doc(NULL) {
	isOpenCLExpanded = false;
	initOpenCL = new CInitOpenCL();
	initOpenCL->Start();
	computeOpenCL = NULL;
	bitmap = NULL;
};

OpenCLDialog::~OpenCLDialog(void) {
	if (initOpenCL) {
		initOpenCL->End(false);
		delete initOpenCL;
		initOpenCL = NULL;
	}
	if (computeOpenCL) {
		computeOpenCL->End(false);
		delete computeOpenCL;
		computeOpenCL = NULL;
	}
	if (bitmap) {
		BaseBitmap::Free(bitmap);
		bitmap = NULL;
	}
}

// Dialog window layout
Bool OpenCLDialog::CreateLayout() {
	SetTitle("OpenCL Render");

	GroupBegin(1000, BFH_SCALEFIT | BFV_FIT, 2, 0, "", BFV_BORDERGROUP_FOLD_OPEN, 850, 100);

	AddStaticText(ID_TEXT_IMAGE_SIZE_X, BFH_LEFT | BFV_TOP, 400, 10, "Image Width", BORDER_NONE);
	AddEditNumber(ID_EDITNUMBER_IMAGE_SIZE_X, BFH_LEFT | BFV_TOP, 100, 10);
	SetInt32(ID_EDITNUMBER_IMAGE_SIZE_X, 1024, 1, 8192);
	AddStaticText(ID_TEXT_IMAGE_SIZE_Y, BFH_LEFT | BFV_TOP, 400, 10, "Image Height", BORDER_NONE);
	AddEditNumber(ID_EDITNUMBER_IMAGE_SIZE_Y, BFH_LEFT | BFV_TOP, 100, 10);
	SetInt32(ID_EDITNUMBER_IMAGE_SIZE_Y, 768, 1, 8192);
	AddStaticText(ID_TEXT_SAMPLING, BFH_LEFT | BFV_TOP, 400, 10, "Sample Count", BORDER_NONE);
	AddEditNumber(ID_EDITNUMBER_SAMPLING, BFH_LEFT | BFV_TOP, 100, 10);
	SetInt32(ID_EDITNUMBER_SAMPLING, 100, 1, 16384);

	GroupBorderSpace(10, 10, 10, 0);

	AddComboBox(ID_COMBOBOX_PROCESSOR, BFH_SCALEFIT, 100, 10, false);

	cl_char deviceName[256];
	if (GetComputationTypeName(0, sizeof(deviceName), deviceName)) {
		AddChild(ID_COMBOBOX_PROCESSOR, 0, (char*)deviceName);
	}

	AddButton(ID_BUTTON_REFRESH_OPENCL, BFH_LEFT | BFV_TOP, 100, 10, "Refresh");
	AddStaticText(ID_TEXT_EMPTY, BFH_LEFT | BFV_TOP, 100, 10, "", BORDER_NONE);
	AddButton(ID_BUTTON_RENDER, BFH_LEFT | BFV_TOP, 100, 10, "Render");

	GroupBorderSpace(10, 10, 10, 0);

	AddStaticText(ID_TEXT_PREPARETIME_TITLE, BFH_LEFT | BFV_TOP, 500, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_PREPARETIME, BFH_LEFT | BFV_TOP, 100, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_PROGRESS_TITLE, BFH_LEFT | BFV_TOP, 500, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_PROGRESS, BFH_LEFT | BFV_TOP, 100, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_TIMEPAST_TITLE, BFH_LEFT | BFV_TOP, 500, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_TIMEPAST, BFH_LEFT | BFV_TOP, 100, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_TIMELEFT_TITLE, BFH_LEFT | BFV_TOP, 500, 10, "", BORDER_NONE);
	AddStaticText(ID_TEXT_TIMELEFT, BFH_LEFT | BFV_TOP, 100, 10, "", BORDER_NONE);

	GroupEnd();

	SetTimer(500);
	isOpenCLExpanded = false;

	return true;
}

// Button action
Bool OpenCLDialog::Command(Int32 id, const BaseContainer & msg) {
	Int32 tmp;
	cl_uint2 imageSize;
	cl_uint samplingCount;
	cl_uint platformDevice;
	switch (id) {
	case ID_BUTTON_REFRESH_OPENCL:
		if (GetIsComputationTypeUpdated() && initOpenCL) {
			ResetComputationType();
			initOpenCL->End(false);
			delete initOpenCL;
			initOpenCL = new CInitOpenCL();
			initOpenCL->Start();
			FreeChildren(ID_COMBOBOX_PROCESSOR);
			cl_char deviceName[256];
			if (GetComputationTypeName(0, sizeof(deviceName), deviceName)) {
				AddChild(ID_COMBOBOX_PROCESSOR, 0, (char*)deviceName);
			}
			SetInt32(ID_COMBOBOX_PROCESSOR, 0);
			SetTimer(500);
			isOpenCLExpanded = false;
		}
		break;
	case ID_BUTTON_RENDER:
		if (!computeOpenCL || 1.f <= GetProgress()) {
			GetInt32(ID_EDITNUMBER_IMAGE_SIZE_X, tmp);
			imageSize.s[0] = (cl_uint)tmp;
			GetInt32(ID_EDITNUMBER_IMAGE_SIZE_Y, tmp);
			imageSize.s[1] = (cl_uint)tmp;
			GetInt32(ID_EDITNUMBER_SAMPLING, tmp);
			samplingCount = (cl_uint)tmp;
			GetInt32(ID_COMBOBOX_PROCESSOR, tmp);
			platformDevice = (cl_uint)tmp;
			if (computeOpenCL) {
				computeOpenCL->End(false);
				delete computeOpenCL;
				computeOpenCL = NULL;
			}
			if (bitmap) {
				BaseBitmap::Free(bitmap);
				bitmap = NULL;
			}
			BaseDocument* polygonizedDoc = doc->Polygonize();
			BaseDraw* bd = doc->GetActiveBaseDraw();
			CameraObject* cam = (CameraObject*)bd->GetSceneCamera(doc);
			if (!cam) cam = (CameraObject*)bd->GetEditorCamera();
			GeData d;
			cam->GetParameter(CAMERAOBJECT_FOV, d, DESCFLAGS_GET_0);
			cl_float fov = (cl_float)d.GetFloat();
			Matrix mg = cam->GetMg();
			computeOpenCL = new ComputeOpenCL(mg, fov, polygonizedDoc, platformDevice, imageSize, samplingCount, &bitmap);
			computeOpenCL->Start();
			ResetTime();
			SetProgress(0.f);
			lastProgress = 0.f;
			computeStart = clock();
			lastUpdate = computeStart;
			prepareStart = computeStart;
		}
		break;
	}

	return GeDialog::Command(id, msg);
}

// Write the time in ?h?m?s format
int GetTimeStr(clock_t timeMs, size_t len, char* str) {
	int hours = (int)floor(timeMs / 3600000.f);
	int minutes = (int)floor((timeMs - hours * 3600000.f) / 60000.f);
	int seconds = (int)floor((timeMs - hours * 3600000.f - minutes * 60000.f) / 1000.f);
	if (100 <= hours) {
		return sprintf_s(str, len, "%dh", hours);
	}
	if (0 < hours) {
		return sprintf_s(str, len, "%dh%dm", hours, minutes);
	}
	if (0 < minutes) {
		return sprintf_s(str, len, "%dm%ds", minutes, seconds);
	}
	return sprintf_s(str, len, "%ds", seconds);
}

// Update the text and OpenCL options every 0.5s
void OpenCLDialog::UpdateDialog() {
	if (GetIsComputationTypeUpdated() && !isOpenCLExpanded) {
		for (size_t i = 1, tc = GetComputationTypeCount(); i < tc; ++i) {
			cl_char deviceName[256];
			if (GetComputationTypeName(i, sizeof(deviceName), deviceName)) {
				AddChild(ID_COMBOBOX_PROCESSOR, (Int32)i, (char*)deviceName);
			}
			else {
				AddChild(ID_COMBOBOX_PROCESSOR, -1, "?");
			}
		}
		isOpenCLExpanded = true;
	}
	if (computeOpenCL || 1.f <= GetProgress()) {
		cl_float t, currentProgress = GetProgress();
		clock_t currentTime = clock();
		clock_t computeStart = GetStartTime();
		char str[32];
		SetString(ID_TEXT_PREPARETIME_TITLE, String("Scene & OpenCL Prepare Time"));
		t = (cl_float)(currentTime - prepareStart);
		if (0.f < computeStart) t = (cl_float)(computeStart - prepareStart);
		GetTimeStr((clock_t)t, sizeof(str), str);
		SetString(ID_TEXT_PREPARETIME, String(str));
		SetString(ID_TEXT_PROGRESS_TITLE, String("Rendering Progress"));
		sprintf(str, "%.1f %%", 100.f * currentProgress);
		SetString(ID_TEXT_PROGRESS, String(str));
		SetString(ID_TEXT_TIMELEFT_TITLE, String("Estimated Time Left"));
		currentTime = clock();
		if (lastProgress != currentProgress) {
			if (1.f <= currentProgress) {
				ShowBitmap(bitmap);
			}
			lastUpdate = currentTime;
			lastProgress = currentProgress;
		}
		t = (computeStart - currentTime) + (lastUpdate - computeStart) / currentProgress;
		if (t != t || t * 2.f == t) {
			SetString(ID_TEXT_TIMELEFT, String(""));
		}
		else {
			if (t < 0.f) t = 0.f;
			GetTimeStr((clock_t)t, sizeof(str), str);
			SetString(ID_TEXT_TIMELEFT, String(str));
		}
		SetString(ID_TEXT_TIMEPAST_TITLE, String("Rendering Time"));
		if (0.f < GetStartTime()) {
			t = (cl_float)(GetEndTime() - GetStartTime());
			if (t <= 0.f) {
				t = (cl_float)(currentTime - GetStartTime());
			}
			GetTimeStr((clock_t)t, sizeof(str), str);
			SetString(ID_TEXT_TIMEPAST, String(str));
		}
		else SetString(ID_TEXT_TIMEPAST, String(""));
	}
	else {
		SetString(ID_TEXT_PREPARETIME_TITLE, String(""));
		SetString(ID_TEXT_PREPARETIME, String(""));
		SetString(ID_TEXT_PROGRESS_TITLE, String(""));
		SetString(ID_TEXT_PROGRESS, String(""));
		SetString(ID_TEXT_TIMELEFT_TITLE, String(""));
		SetString(ID_TEXT_TIMELEFT, String(""));
		SetString(ID_TEXT_TIMEPAST_TITLE, String(""));
		SetString(ID_TEXT_TIMEPAST, String(""));
	}
}

void OpenCLDialog::SetDocument(BaseDocument* doc) {
	this->doc = doc;
}

void OpenCLDialog::Timer(const BaseContainer& msg) {
	UpdateDialog();
}




#define ID_SDK_DIALOG_COMMAND 1035350 ///< Command ID

// Class that controls the opening of the dialog
class OpenCLRendering : public CommandData {
	INSTANCEOF(OpenCLRendering, CommandData)

public:
	virtual Bool Execute(BaseDocument* doc);
	virtual Bool RestoreLayout(void* secret);

	static OpenCLRendering* Alloc() { return NewObjClear(OpenCLRendering); }

private:
	OpenCLDialog _dialog;
};

Bool OpenCLRendering::Execute(BaseDocument* doc) {
	if (_dialog.IsOpen() == false) {
		_dialog.Open(DLG_TYPE_ASYNC, ID_SDK_DIALOG_COMMAND, -1, -1, 400, 208);
		_dialog.SetDocument(doc);
	}
	return true;
}

Bool OpenCLRendering::RestoreLayout(void* secret) {
	return _dialog.RestoreLayout(ID_SDK_DIALOG_COMMAND, 0, secret);
}

Bool RegisterOpenCLRender(void) {
	// be sure to use a unique ID obtained from www.plugincafe.com
	return RegisterCommandPlugin(1000956, String("OpenCL Rendering"), 0, AutoBitmap("icon.tif"), String("OpenCL Rendering"), OpenCLRendering::Alloc());
}



// Get the 3D normal of a 3D triangle
cl_float3 getNormal(cl_float3 a, cl_float3 b, cl_float3 c) {
	cl_float3 ab, ac;
	ab.s[0] = b.s[0] - a.s[0];
	ab.s[1] = b.s[1] - a.s[1];
	ab.s[2] = b.s[2] - a.s[2];
	ac.s[0] = c.s[0] - a.s[0];
	ac.s[1] = c.s[1] - a.s[1];
	ac.s[2] = c.s[2] - a.s[2];
	return cross(ab, ac);
}

cl_float3 VectorToFloat3(Vector& v) {
	cl_float3 f3;
	f3.s[0] = (cl_float)v.x;
	f3.s[1] = (cl_float)v.y;
	f3.s[2] = (cl_float)v.z;
	return f3;
}

// Quaternion multiply and rotation. They are not used anymore, but I keep them just in case I need some rotational decomposition of the scene
// from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/code/
cl_float4 quatMul(cl_float4 q1, cl_float4 q2) {
	cl_float4 r;
	r.s[0] = q1.s[0] * q2.s[3] + q1.s[1] * q2.s[2] - q1.s[2] * q2.s[1] + q1.s[3] * q2.s[0];
	r.s[1] = -q1.s[0] * q2.s[2] + q1.s[1] * q2.s[3] + q1.s[2] * q2.s[0] + q1.s[3] * q2.s[1];
	r.s[2] = q1.s[0] * q2.s[1] - q1.s[1] * q2.s[0] + q1.s[2] * q2.s[3] + q1.s[3] * q2.s[2];
	r.s[3] = -q1.s[0] * q2.s[0] - q1.s[1] * q2.s[1] - q1.s[2] * q2.s[2] + q1.s[3] * q2.s[3];
	return r;
}
cl_float3 quatRot(cl_float4 q, cl_float3 v) {
	cl_float3 r;
	cl_float4 v0;
	v0.s[0] = v.s[0];
	v0.s[1] = v.s[1];
	v0.s[2] = v.s[2];
	v0.s[3] = 0;
	cl_float4 qv0 = quatMul(q, v0);
	cl_float4 qconj;
	qconj.s[0] = -q.s[0];
	qconj.s[1] = -q.s[1];
	qconj.s[2] = -q.s[2];
	qconj.s[3] = q.s[3];
	r = quatMul(qv0, qconj);
	return r;
}

// Finds the camera vectors and pixel sizes from the position, lookat, up, field of view and image size
void SetCamera(cl_float3* outCameraEyeToTopLeftVector, cl_float3* outCameraLeftToRightPixelSizeVector, cl_float3* outCameraTopToBottomPixelSizeVector, cl_float* outCameraPixelSizeInv, cl_float3 position, cl_float3 object, cl_float3 up, cl_float fov, cl_uint2 imageDimension) {
	float3 cameraVector, rightToLeftVector, bottomToTopVector, rightToLeftUnitVector, bottomToTopUnitVector;
	float midToLeftLength, midToTopLength, rightToLeftVectorLen, bottomToTopVectorLen;
	cameraVector.s[0] = object.s[0] - position.s[0];
	cameraVector.s[1] = object.s[1] - position.s[1];
	cameraVector.s[2] = object.s[2] - position.s[2];
	rightToLeftVector = cross(up, cameraVector);
	bottomToTopVector = up;
	midToLeftLength = (float)sqrt(dot(cameraVector, cameraVector)) * (float)tan(fov / 2.f);
	midToTopLength = midToLeftLength * (float)imageDimension.s[1] / (float)imageDimension.s[0];

	rightToLeftVectorLen = (float)sqrt(dot(rightToLeftVector, rightToLeftVector));
	bottomToTopVectorLen = (float)sqrt(dot(bottomToTopVector, bottomToTopVector));
	rightToLeftUnitVector.s[0] = rightToLeftVector.s[0] / rightToLeftVectorLen;
	rightToLeftUnitVector.s[1] = rightToLeftVector.s[1] / rightToLeftVectorLen;
	rightToLeftUnitVector.s[2] = rightToLeftVector.s[2] / rightToLeftVectorLen;
	bottomToTopUnitVector.s[0] = bottomToTopVector.s[0] / bottomToTopVectorLen;
	bottomToTopUnitVector.s[1] = bottomToTopVector.s[1] / bottomToTopVectorLen;
	bottomToTopUnitVector.s[2] = bottomToTopVector.s[2] / bottomToTopVectorLen;

	*outCameraPixelSizeInv = ((float)imageDimension.s[0]) / (2.f * midToLeftLength);
	outCameraEyeToTopLeftVector->s[0] = cameraVector.s[0] - midToLeftLength * rightToLeftUnitVector.s[0] + midToTopLength * bottomToTopUnitVector.s[0];
	outCameraEyeToTopLeftVector->s[1] = cameraVector.s[1] - midToLeftLength * rightToLeftUnitVector.s[1] + midToTopLength * bottomToTopUnitVector.s[1];
	outCameraEyeToTopLeftVector->s[2] = cameraVector.s[2] - midToLeftLength * rightToLeftUnitVector.s[2] + midToTopLength * bottomToTopUnitVector.s[2];
	outCameraLeftToRightPixelSizeVector->s[0] = rightToLeftUnitVector.s[0] / *outCameraPixelSizeInv;
	outCameraLeftToRightPixelSizeVector->s[1] = rightToLeftUnitVector.s[1] / *outCameraPixelSizeInv;
	outCameraLeftToRightPixelSizeVector->s[2] = rightToLeftUnitVector.s[2] / *outCameraPixelSizeInv;
	outCameraTopToBottomPixelSizeVector->s[0] = -bottomToTopUnitVector.s[0] / *outCameraPixelSizeInv;
	outCameraTopToBottomPixelSizeVector->s[1] = -bottomToTopUnitVector.s[1] / *outCameraPixelSizeInv;
	outCameraTopToBottomPixelSizeVector->s[2] = -bottomToTopUnitVector.s[2] / *outCameraPixelSizeInv;
}

// This function is straight from the Maxon forums, it maps a texture dependent on the projection type
// https://developers.maxon.net/docs/Cinema4DCPPSDK/html/struct_tex_data.html
Bool ShdProjectPoint(TextureTag* textureTag, const Vector& p, const Vector& n, Vector* uv) {
	GeData d;
	textureTag->GetParameter(TEXTURETAG_OFFSETX, d, DESCFLAGS_GET_0);
	Float ox = d.GetFloat();
	textureTag->GetParameter(TEXTURETAG_OFFSETY, d, DESCFLAGS_GET_0);
	Float oy = d.GetFloat();
	textureTag->GetParameter(TEXTURETAG_LENGTHX, d, DESCFLAGS_GET_0);
	Float lenx = d.GetFloat();
	textureTag->GetParameter(TEXTURETAG_LENGTHY, d, DESCFLAGS_GET_0);
	Float leny = d.GetFloat();
	textureTag->GetParameter(TEXTURETAG_PROJECTION, d, DESCFLAGS_GET_0);
	Int32 proj = d.GetInt32();
	textureTag->GetParameter(TEXTURETAG_TILE, d, DESCFLAGS_GET_0);
	Bool textile = d.GetBool();

	Float lenxinv = 0.0, lenyinv = 0.0;
	if (lenx != 0.0) lenxinv = 1.0 / lenx;
	if (leny != 0.0) lenyinv = 1.0 / leny;
	switch (proj) {
	case P_VOLUMESHADER:
	{
		*uv = p;// *tdp->im;
		return true;
	}
	case P_SPHERICAL: default:
	{
		Vector d = p;// *tdp->im;
		Float sq = Sqrt(d.x*d.x + d.z*d.z);
		if (sq == 0.0)
		{
			uv->x = 0.0;
			if (d.y>0.0)
				uv->y = +0.5;
			else
				uv->y = -0.5;
		}
		else
		{
			uv->x = ACos(d.x / sq) / PI2;
			if (d.z<0.0) uv->x = 1.0 - uv->x;
			uv->x -= ox;
			if (lenx>0.0 && uv->x<0.0)
				uv->x += 1.0;
			else if (lenx<0.0 && uv->x>0.0)
				uv->x -= 1.0;
			uv->x *= lenxinv;
			uv->y = 0.5 + ATan(d.y / sq) / PI;
		}
		uv->y = -(uv->y - oy)*lenyinv;
		break;
	}
	case P_SHRINKWRAP:
	{
		Vector d = p;// *tdp->im;
		Float   sn, cs, sq = Sqrt(d.x*d.x + d.z*d.z);
		if (sq == 0.0)
		{
			uv->x = 0.0;
			if (d.y>0.0)
				uv->y = 0.0;
			else
				uv->y = 1.0;
		}
		else
		{
			uv->x = ACos(d.x / sq) / PI2;
			if (d.z<0.0) uv->x = 1.0 - uv->x;
			uv->y = 0.5 - ATan(d.y / sq) / PI;
		}
		SinCos(uv->x*PI2, sn, cs);
		uv->x = (0.5 + 0.5*cs*uv->y - ox)*lenxinv;
		uv->y = (0.5 + 0.5*sn*uv->y - oy)*lenyinv;
		break;
	}
	case P_CYLINDRICAL:
	{
		Vector d = p;// *tdp->im;
		Float sq = Sqrt(d.x*d.x + d.z*d.z);
		if (sq == 0.0)
			uv->x = 0.0;
		else
		{
			uv->x = ACos(d.x / sq) / PI2;
			if (d.z<0.0) uv->x = 1.0 - uv->x;
			uv->x -= ox;
			if (lenx>0.0 && uv->x<0.0)
				uv->x += 1.0;
			else if (lenx<0.0 && uv->x>0.0)
				uv->x -= 1.0;
			uv->x *= lenxinv;
		}
		uv->y = -(d.y*0.5 + oy)*lenyinv;
		break;
	}
	case P_FLAT: case P_SPATIAL:
	{
		Vector d = p;// *tdp->im;
		uv->x = (d.x*0.5 - ox)*lenxinv;
		uv->y = -(d.y*0.5 + oy)*lenyinv;
		break;
	}
	case P_CUBIC:
	{
		Vector d = p;// *tdp->im;
		Vector v = n;// ^ tdp->im;
		Int32   dir;
		if (Abs(v.x)>Abs(v.y))
		{
			if (Abs(v.x)>Abs(v.z))
				dir = 0;
			else
				dir = 2;
		}
		else
		{
			if (Abs(v.y)>Abs(v.z))
				dir = 1;
			else
				dir = 2;
		}
		switch (dir)
		{
		case 0: // x axis
		{
					if (v.x<0.0)
						uv->x = (-d.z*0.5 - ox)*lenxinv;
					else
						uv->x = (d.z*0.5 - ox)*lenxinv;
					uv->y = -(d.y*0.5 + oy)*lenyinv;
					break;
		}
		case 1:  // y axis
		{
						if (v.y<0.0)
							uv->y = (d.z*0.5 - oy)*lenyinv;
						else
							uv->y = (-d.z*0.5 - oy)*lenyinv;
						uv->x = (d.x*0.5 - ox)*lenxinv;
						break;
		}
		case 2: // z axis
		{
					if (v.z<0.0)
						uv->x = (d.x*0.5 - ox)*lenxinv;
					else
						uv->x = (-d.x*0.5 - ox)*lenxinv;
					uv->y = -(d.y*0.5 + oy)*lenyinv;
					break;
		}
		}
		break;
	}
	case P_FRONTAL:
	{
		DebugAssert(false, "not handled yet :(");
		/*RayParameter *param = sd->GetRayParameter();
		Float ox = 0.0, oy = 0.0, ax = param->xres, ay = param->yres;
		Int32 curr_x, curr_y, scl;
		sd->GetXY(&curr_x, &curr_y, &scl);
		uv->x = ((Float(curr_x) / Float(scl) - ox) / ax - tdp->ox)*lenxinv;
		uv->y = ((Float(curr_y) / Float(scl) - ox) / ay - tdp->oy)*lenyinv;*/
		break;
	}
	case P_UVW:
	{
		DebugAssert(false, "not handled yet :(");
		/*RayObject *op = sd->ID_to_Obj(lhit, nullptr);
		if (op && tdp->uvwind<op->uvwcnt && op->uvwadr[tdp->uvwind])
			*uv = sd->GetPointUVW(tdp, lhit, p);
		else
			uv->x = uv->y = 0.0;*/
		break;
	}
	}
	if (textile)
		return true;
	else
		return uv->x >= 0.0 && uv->x <= 1.0 && uv->y >= 0.0 && uv->y <= 1.0;
}

// This function goes recursively through the whole scene to count triangles, vertices and lights (for allocation size)
void CountPolygonsRecursive(BaseObject* op, cl_uint* vertexCount, cl_uint* triangleCount, cl_uint* lightCount) {
	while (op) {
		Int32 objectType = op->GetType();
		switch (objectType) {
		case Opolygon: {
				PolygonObject* po = static_cast<PolygonObject*>(op);
				*vertexCount += po->GetPointCount();
				*triangleCount += 2 * po->GetPolygonCount();
			}
			break;
		case Olight:
			++(*lightCount);
			break;
		case Ocamera:
			break;
		case Onull:
			break;
		default:
			DebugAssert(false, "Unrecognized object");
		}
		BaseObject* child = op->GetDown();
		if (child) {
			CountPolygonsRecursive(child, vertexCount, triangleCount, lightCount);
		}
		op = op->GetNext();
	}
}

// This function recursively goes through the scene to grab vertices, triangles, normals, uvs, lights, materialIds
// Should be worked on and refactored (TODO)
// When I can't find normals, I will make the triangle face the camera. It's a cheap tactic that will cause bugs, but until I find how to do it... (TODO)
void AddPolygonsRecursive(BaseObject* op, cl_uint* vertexCursor, cl_uint* triangleCursor, cl_uint* lightCursor, cl_float3* vertex, cl_int3* triangleVertexIndex, cl_int* triangleMaterialId, cl_float2* triangleUv, cl_float3* triangleNormal, cl_int* lightType, cl_float3* lightPosition, cl_float3* lightDirection, cl_float3* lightColour, cl_float* lightRadius, cl_float* lightHalfAttenuationDistance, std::map<BaseMaterial*, Int32>& materialDatabase, cl_float3 cameraEye, cl_float* materialTextureRatio) {
	while (op) {
		Int32 objectType = op->GetType();
		const Matrix& globalMatrix = op->GetMg();
		Matrix globalRotation = globalMatrix;
		globalRotation.off.SetZero();
		switch (objectType) {
		case Opolygon: {
				PolygonObject* po = (PolygonObject*)(op);
				cl_uint firstVertex = *vertexCursor;
				cl_uint firstTriangle = *triangleCursor;

				const Vector* vertexLocal = po->GetPointR();
				for (Int32 i = 0, im = po->GetPointCount(); i < im; ++i) {
					Vector tv = globalMatrix * vertexLocal[i];
					vertex[*vertexCursor].s[0] = (cl_float)tv.x;
					vertex[*vertexCursor].s[1] = (cl_float)tv.y;
					vertex[(*vertexCursor)++].s[2] = (cl_float)tv.z;
				}

				Vector32* phongNormals = po->CreatePhongNormals();
				const Int32 polygonCount = po->GetPolygonCount();
				const CPolygon* p = po->GetPolygonR();
				*triangleCursor = firstTriangle;
				for (Int32 i = 0; i < polygonCount; ++i) {
					Bool isSquare = (p[i].c != p[i].d);
					triangleVertexIndex[*triangleCursor].s[0] = firstVertex + p[i].a;
					triangleVertexIndex[*triangleCursor].s[1] = firstVertex + p[i].b;
					triangleVertexIndex[*triangleCursor].s[2] = firstVertex + p[i].c;
					if (phongNormals) {
						Vector an(phongNormals[4 * i + 0].x, phongNormals[4 * i + 0].y, phongNormals[4 * i + 0].z);
						Vector bn(phongNormals[4 * i + 1].x, phongNormals[4 * i + 1].y, phongNormals[4 * i + 1].z);
						Vector cn(phongNormals[4 * i + 2].x, phongNormals[4 * i + 2].y, phongNormals[4 * i + 2].z);
						an = globalRotation * an;
						bn = globalRotation * bn;
						cn = globalRotation * cn;
						an.Normalize();
						bn.Normalize();
						cn.Normalize();
						triangleNormal[3 * *triangleCursor + 0].s[0] = (cl_float)an.x;
						triangleNormal[3 * *triangleCursor + 0].s[1] = (cl_float)an.y;
						triangleNormal[3 * *triangleCursor + 0].s[2] = (cl_float)an.z;
						triangleNormal[3 * *triangleCursor + 1].s[0] = (cl_float)bn.x;
						triangleNormal[3 * *triangleCursor + 1].s[1] = (cl_float)bn.y;
						triangleNormal[3 * *triangleCursor + 1].s[2] = (cl_float)bn.z;
						triangleNormal[3 * *triangleCursor + 2].s[0] = (cl_float)cn.x;
						triangleNormal[3 * *triangleCursor + 2].s[1] = (cl_float)cn.y;
						triangleNormal[3 * *triangleCursor + 2].s[2] = (cl_float)cn.z;
					}
					else {
						cl_float3 tn = getNormal(vertex[triangleVertexIndex[*triangleCursor].s[0]], vertex[triangleVertexIndex[*triangleCursor].s[1]], vertex[triangleVertexIndex[*triangleCursor].s[2]]);
						cl_float lenInv = 1.f / (cl_float)sqrt(dot(tn, tn));
						// If we don't have normals, make the triangles face the camera...
						if (0 <= dot(vector(cameraEye, vertex[triangleVertexIndex[*triangleCursor].s[0]]), tn)) lenInv = -lenInv;
						tn.s[0] *= lenInv;
						tn.s[1] *= lenInv;
						tn.s[2] *= lenInv;
						triangleNormal[3 * *triangleCursor + 0].s[0] = tn.s[0];
						triangleNormal[3 * *triangleCursor + 0].s[1] = tn.s[1];
						triangleNormal[3 * *triangleCursor + 0].s[2] = tn.s[2];
						triangleNormal[3 * *triangleCursor + 1].s[0] = tn.s[0];
						triangleNormal[3 * *triangleCursor + 1].s[1] = tn.s[1];
						triangleNormal[3 * *triangleCursor + 1].s[2] = tn.s[2];
						triangleNormal[3 * *triangleCursor + 2].s[0] = tn.s[0];
						triangleNormal[3 * *triangleCursor + 2].s[1] = tn.s[1];
						triangleNormal[3 * *triangleCursor + 2].s[2] = tn.s[2];
					}
					++(*triangleCursor);
					if (isSquare) {
						triangleVertexIndex[*triangleCursor].s[0] = firstVertex + p[i].a;
						triangleVertexIndex[*triangleCursor].s[1] = firstVertex + p[i].c;
						triangleVertexIndex[*triangleCursor].s[2] = firstVertex + p[i].d;
						if (phongNormals) {
							Vector an(phongNormals[4 * i + 0].x, phongNormals[4 * i + 0].y, phongNormals[4 * i + 0].z);
							Vector cn(phongNormals[4 * i + 2].x, phongNormals[4 * i + 2].y, phongNormals[4 * i + 2].z);
							Vector dn(phongNormals[4 * i + 3].x, phongNormals[4 * i + 3].y, phongNormals[4 * i + 3].z);
							an = globalRotation * an;
							cn = globalRotation * cn;
							dn = globalRotation * dn;
							an.Normalize();
							cn.Normalize();
							dn.Normalize();
							triangleNormal[3 * *triangleCursor + 0].s[0] = (cl_float)an.x;
							triangleNormal[3 * *triangleCursor + 0].s[1] = (cl_float)an.y;
							triangleNormal[3 * *triangleCursor + 0].s[2] = (cl_float)an.z;
							triangleNormal[3 * *triangleCursor + 1].s[0] = (cl_float)cn.x;
							triangleNormal[3 * *triangleCursor + 1].s[1] = (cl_float)cn.y;
							triangleNormal[3 * *triangleCursor + 1].s[2] = (cl_float)cn.z;
							triangleNormal[3 * *triangleCursor + 2].s[0] = (cl_float)dn.x;
							triangleNormal[3 * *triangleCursor + 2].s[1] = (cl_float)dn.y;
							triangleNormal[3 * *triangleCursor + 2].s[2] = (cl_float)dn.z;
						}
						else {
							cl_float3 tn = getNormal(vertex[triangleVertexIndex[*triangleCursor].s[0]], vertex[triangleVertexIndex[*triangleCursor].s[1]], vertex[triangleVertexIndex[*triangleCursor].s[2]]);
							cl_float lenInv = 1.f / (cl_float)sqrt(dot(tn, tn));
							// If we don't have normals, make the triangles face the camera...
							if (0 <= dot(vector(cameraEye, vertex[triangleVertexIndex[*triangleCursor].s[0]]), tn)) lenInv = -lenInv;
							tn.s[0] *= lenInv;
							tn.s[1] *= lenInv;
							tn.s[2] *= lenInv;
							triangleNormal[3 * *triangleCursor + 0].s[0] = tn.s[0];
							triangleNormal[3 * *triangleCursor + 0].s[1] = tn.s[1];
							triangleNormal[3 * *triangleCursor + 0].s[2] = tn.s[2];
							triangleNormal[3 * *triangleCursor + 1].s[0] = tn.s[0];
							triangleNormal[3 * *triangleCursor + 1].s[1] = tn.s[1];
							triangleNormal[3 * *triangleCursor + 1].s[2] = tn.s[2];
							triangleNormal[3 * *triangleCursor + 2].s[0] = tn.s[0];
							triangleNormal[3 * *triangleCursor + 2].s[1] = tn.s[1];
							triangleNormal[3 * *triangleCursor + 2].s[2] = tn.s[2];
						}
						++(*triangleCursor);
					}
				}
				if (phongNormals) DeleteMem(phongNormals); phongNormals = NULL;

				std::map<std::string, BaseSelect*> selectionMap;
				BaseTag* tags = op->GetTag(Tpolygonselection);
				while (tags) {
					if (tags->GetType() == Tpolygonselection) {
						SelectionTag* selectionTag = (SelectionTag*)tags;
						Char buff[2048];
						selectionTag->GetName().GetCString(buff, 2048);
						selectionMap.insert(std::make_pair(std::string(buff), selectionTag->GetBaseSelect()));
					}
					tags = tags->GetNext();
				}

				std::map<cl_int, TextureTag*> ttmap;
				tags = op->GetTag(Ttexture);
				while (tags) {
					if (tags->GetType() == Ttexture) {
						TextureTag* textureTag = (TextureTag*)tags;
						BaseMaterial* m = textureTag->GetMaterial();
						BaseChannel* channel = m->GetChannel(CHANNEL_COLOR);
						if (channel) {
							Char* str = channel->GetData().GetString(BASECHANNEL_TEXTURE).GetCStringCopy();
							str = str;
						}
						GeData d;
						textureTag->GetParameter(TEXTURETAG_RESTRICTION, d, DESCFLAGS_GET_0);
						const String& selection = d.GetString();
						Char buff[2048];
						selection.GetCString(buff, 2048);
						std::string selectionStr(buff);

						std::map<BaseMaterial*, Int32>::const_iterator materialDatabaseIt = materialDatabase.find(m);
						if (materialDatabase.cend() != materialDatabaseIt) {
							ttmap.insert(std::make_pair(materialDatabaseIt->second, textureTag));
							std::map<std::string, BaseSelect*>::const_iterator selectionMapIt = selectionMap.find(selectionStr);
							BaseSelect* bs = selectionMap.cend() == selectionMapIt ? NULL : selectionMapIt->second;
							*triangleCursor = firstTriangle;
							for (Int32 i = 0; i < polygonCount; ++i) {
								Bool isSquare = (p[i].c != p[i].d);
								if (!bs || bs->IsSelected(i)) {
									for (int t = 0; t < ((isSquare) ? 2 : 1); ++t) {
										cl_float3 a = vertex[triangleVertexIndex[*triangleCursor].s[0]];
										cl_float3 b = vertex[triangleVertexIndex[*triangleCursor].s[1]];
										cl_float3 c = vertex[triangleVertexIndex[*triangleCursor].s[2]];
										cl_float3 ab = vector(a, b); cl_float3 bc = vector(b, c); cl_float3 ca = vector(c, a);
										cl_float maxVertexDistanceSq = fmax(dot(ab, ab), fmax(dot(bc, bc), dot(ca, ca)));
										cl_float maxTextureX = fmax(triangleUv[3 * *triangleCursor + 0].s[0], fmax(triangleUv[3 * *triangleCursor + 1].s[0], triangleUv[3 * *triangleCursor + 2].s[0]));
										cl_float maxTextureY = fmax(triangleUv[3 * *triangleCursor + 0].s[1], fmax(triangleUv[3 * *triangleCursor + 1].s[1], triangleUv[3 * *triangleCursor + 2].s[1]));
										cl_float maxTextureSq = maxTextureX * maxTextureX + maxTextureY * maxTextureY;
										cl_float spaceToTextureRatioSq = maxVertexDistanceSq / (0.f < maxTextureSq ? maxTextureSq : 1.f);
										Int32 materialId = materialDatabaseIt->second;
										materialTextureRatio[materialId] =
											fmax(materialTextureRatio[materialId],
												sqrt(spaceToTextureRatioSq / GetPointToTriangleSqLen(a, b, c, cameraEye)));
										triangleMaterialId[(*triangleCursor)++] = materialId;
									}
								}
								else {
									*triangleCursor += (isSquare) ? 2 : 1;
								}
							}
						}
					}
					tags = tags->GetNext();
				}

				Int32 uvwTagDataCount = 0;
				tags = op->GetTag(Tuvw);
				while (tags) {
					if (tags->GetType() == Tuvw) {
						UVWTag* uvwTag = (UVWTag*)tags;
						uvwTagDataCount = uvwTag->GetDataCount();
						*triangleCursor = firstTriangle;
						for (Int32 i = 0; i < uvwTagDataCount; i++)
						{
							Bool isSquare = (p[i].c != p[i].d);
							UVWStruct s = uvwTag->GetSlow(i);
							triangleUv[3 * *triangleCursor + 0].s[0] = (cl_float)s.a.x;
							triangleUv[3 * *triangleCursor + 0].s[1] = (cl_float)s.a.y;
							triangleUv[3 * *triangleCursor + 1].s[0] = (cl_float)s.b.x;
							triangleUv[3 * *triangleCursor + 1].s[1] = (cl_float)s.b.y;
							triangleUv[3 * *triangleCursor + 2].s[0] = (cl_float)s.c.x;
							triangleUv[3 * *triangleCursor + 2].s[1] = (cl_float)s.c.y;
							++(*triangleCursor);
							if (isSquare) {
								triangleUv[3 * *triangleCursor + 0].s[0] = (cl_float)s.a.x;
								triangleUv[3 * *triangleCursor + 0].s[1] = (cl_float)s.a.y;
								triangleUv[3 * *triangleCursor + 1].s[0] = (cl_float)s.c.x;
								triangleUv[3 * *triangleCursor + 1].s[1] = (cl_float)s.c.y;
								triangleUv[3 * *triangleCursor + 2].s[0] = (cl_float)s.d.x;
								triangleUv[3 * *triangleCursor + 2].s[1] = (cl_float)s.d.y;
								++(*triangleCursor);
							}
						}
					}
					tags = tags->GetNext();
				}
				if (uvwTagDataCount != polygonCount) {
					*triangleCursor = firstTriangle;
					for (Int32 i = 0; i < polygonCount; i++)
					{
						Bool isSquare = (p[i].c != p[i].d);
						Int32 triangleCount = isSquare ? 2 : 1;
						for (Int32 t = 0; t < triangleCount; ++t) {
							std::map<cl_int, TextureTag*>::const_iterator ttIt = ttmap.find(triangleMaterialId[*triangleCursor]);
							if (ttmap.cend() != ttIt) {
								cl_int v[3];
								v[0] = triangleVertexIndex[*triangleCursor].s[0];
								v[1] = triangleVertexIndex[*triangleCursor].s[1];
								v[2] = triangleVertexIndex[*triangleCursor].s[2];
								Vector p[3] = {
									vertexLocal[v[0] - firstVertex],
									vertexLocal[v[1] - firstVertex],
									vertexLocal[v[2] - firstVertex]
								};
								Vector n = Cross(p[1] - p[0], p[2] - p[0]);
								for (Int32 i = 0; i < 3; ++i) {
									Vector uv;
									ShdProjectPoint(ttIt->second, p[i], n, &uv);
									triangleUv[3 * *triangleCursor + i].s[0] = (cl_float)uv.x;
									triangleUv[3 * *triangleCursor + i].s[1] = (cl_float)uv.y;
								}
								++(*triangleCursor);
							}
							else {
								triangleUv[3 * *triangleCursor + 0].s[0] = 0.f;
								triangleUv[3 * *triangleCursor + 0].s[1] = 0.f;
								triangleUv[3 * *triangleCursor + 1].s[0] = 0.f;
								triangleUv[3 * *triangleCursor + 1].s[1] = 1.f;
								triangleUv[3 * *triangleCursor + 2].s[0] = 1.f;
								triangleUv[3 * *triangleCursor + 2].s[1] = 1.f;
								++(*triangleCursor);
							}
						}
					}
				}
			}
			break;
		case Olight: {
				// The Sun's angular size is 0.52 degrees
				cl_float sunAngleDegrees = 0.52f;
				Vector absPosLocal = op->GetAbsPos();
				Vector absPosGlobal = globalMatrix * absPosLocal;
				Vector absLightDirection = globalMatrix * Vector(0, 0, 1.);
				GeData d;
				op->GetParameter(LIGHT_COLOR, d, DESCFLAGS_GET_0);
				Vector lightColor = d.GetVector();
				op->GetParameter(LIGHT_BRIGHTNESS, d, DESCFLAGS_GET_0);
				Float lightBrightness = d.GetFloat();
				op->GetParameter(LIGHT_TYPE, d, DESCFLAGS_GET_0);
				Int32 type = d.GetInt32();
				lightType[(*lightCursor)] = type;
				lightPosition[(*lightCursor)].s[0] = (cl_float)absPosGlobal.x;
				lightPosition[(*lightCursor)].s[1] = (cl_float)absPosGlobal.y;
				lightPosition[(*lightCursor)].s[2] = (cl_float)absPosGlobal.z;
				lightHalfAttenuationDistance[(*lightCursor)] = INFINITY;
				cl_float absLightDirectionLen = (cl_float)sqrt(Dot(absLightDirection, absLightDirection));
				lightDirection[(*lightCursor)].s[0] = (cl_float)absLightDirection.x / absLightDirectionLen;
				lightDirection[(*lightCursor)].s[1] = (cl_float)absLightDirection.y / absLightDirectionLen;
				lightDirection[(*lightCursor)].s[2] = (cl_float)absLightDirection.z / absLightDirectionLen;
				lightColour[(*lightCursor)].s[0] = (cl_float)lightColor.x * (cl_float)lightBrightness;
				lightColour[(*lightCursor)].s[1] = (cl_float)lightColor.y * (cl_float)lightBrightness;
				lightColour[(*lightCursor)].s[2] = (cl_float)lightColor.z * (cl_float)lightBrightness;
				//light[(*lightCursor)].rayDimension.s[0] = 128;
				//light[(*lightCursor)].rayDimension.s[1] = 128;
				lightRadius[(*lightCursor)] = sunAngleDegrees;
				++(*lightCursor);
			}
			break;
		case Ocamera:
			break;
		case Onull:
			break;
		default:
			DebugAssert(false, "Unrecognized object");
		}
		BaseObject* child = op->GetDown();
		if (child) {
			AddPolygonsRecursive(child, vertexCursor, triangleCursor, lightCursor, vertex, triangleVertexIndex, triangleMaterialId, triangleUv, triangleNormal, lightType, lightPosition, lightDirection, lightColour, lightRadius, lightHalfAttenuationDistance, materialDatabase, cameraEye, materialTextureRatio);
		}
		op = op->GetNext();
	}
}

// This function sets up the camera, gets the materials and textures (not well) and the lauches the ray tracing. It will return when the rendering is done.
// This can be improved a lot, was only tested on one scene and material import is plain bad (TODO)
// Should be worked on and refactored (TODO)
// In other words, TODO TODO TODO!
Bool parseAndRender(Matrix camMatrix, cl_float fov, BaseDocument* doc, cl_uint platformDevice, cl_uint2 outputImageDimension, cl_uint sampleCount, BaseBitmap** bitmap) {
//	GeShowMouse(MOUSE_BUSY);

	cl_bool computationDone = false;
	cl_uint outputImageSize = outputImageDimension.s[0] * outputImageDimension.s[1];
	cl_ushort* outputRed = new cl_ushort[outputImageSize];
	memset(outputRed, 0, outputImageSize * sizeof(cl_ushort));
	cl_ushort* outputGreen = new cl_ushort[outputImageSize];
	memset(outputGreen, 0, outputImageSize * sizeof(cl_ushort));
	cl_ushort* outputBlue = new cl_ushort[outputImageSize];
	memset(outputBlue, 0, outputImageSize * sizeof(cl_ushort));

	cl_uint vertexCount = 0;
	cl_float3* vertex = NULL;

	cl_uint triangleCount = 0;
	cl_int3* triangleVertexIndex = NULL;
	cl_int* triangleMaterialId = NULL;
	cl_float2* triangleUv = NULL; // 3x triangleCount
	cl_float3* triangleNormal = NULL; // 3x triangleCount

	CameraTriangleList* cameraTriangleList = NULL;
	SceneTriangleList* sceneTriangleList = NULL;

	cl_uint materialCount = 0;
	cl_uint2* materialImageSize = NULL; // materialCount * _MATERIAL_CHANNEL_COUNT
	cl_int* materialImageStart = NULL; // materialCount * _MATERIAL_CHANNEL_COUNT
	cl_float* materialTextureRatio = NULL; // materialCount

	cl_uint lightCount = 0;
	cl_int* lightType = NULL;
	cl_float3* lightPosition = NULL;
	cl_float3* lightDirection = NULL;
	cl_float3* lightColour = NULL;
	cl_float* lightRadius = NULL;
	cl_float* lightHalfAttenuationDistance = NULL;

	cl_float3 cameraEye;
	cl_float3 cameraEyeToTopLeftVector;
	cl_float3 cameraLeftToRightPixelSizeVector;
	cl_float3 cameraTopToBottomPixelSizeVector;
	cl_float cameraPixelSizeInv;
	Vector eye = camMatrix * Vector(0., 0., 0.);
	Vector obj = camMatrix * Vector(0., 0., 1.);
	Vector up = camMatrix * Vector(0., 1., 0.);
	Vector localUp = up - eye;
	SetCamera(&cameraEyeToTopLeftVector, &cameraLeftToRightPixelSizeVector, &cameraTopToBottomPixelSizeVector, &cameraPixelSizeInv, VectorToFloat3(eye), VectorToFloat3(obj), VectorToFloat3(localUp), fov, outputImageDimension);
	cameraEye = VectorToFloat3(eye);

	BaseDocument* polygonizedDoc = doc;
	BaseObject* op = polygonizedDoc->GetFirstObject();
	BaseMaterial* bm = polygonizedDoc->GetFirstMaterial();
	Material* material = (Material*)bm;
	std::map<BaseMaterial*, Int32> materialDatabase;
	while (material) {
		materialDatabase.insert(std::make_pair(material, materialCount++));
		material = (Material*)material->GetNext();
	}
	materialImageSize = new cl_uint2[materialCount * _MATERIAL_CHANNEL_COUNT]; // materialCount * _MATERIAL_CHANNEL_COUNT
	if (!materialImageSize) goto cleanup;
	memset(materialImageSize, 0, materialCount * _MATERIAL_CHANNEL_COUNT * sizeof(cl_uint2));
	materialImageStart = new cl_int[materialCount * _MATERIAL_CHANNEL_COUNT + 1]; // materialCount * _MATERIAL_CHANNEL_COUNT
	if (!materialImageStart) goto cleanup;
	memset(materialImageStart, 0, (materialCount * _MATERIAL_CHANNEL_COUNT + 1) * sizeof(cl_int));
	materialTextureRatio = new cl_float[materialCount]; // materialCount
	if (!materialTextureRatio) goto cleanup;
	memset(materialTextureRatio, 0, materialCount * sizeof(cl_float));

	cl_uint texturesCursor = 0;
	cl_uchar3* texturesTmp = new cl_uchar3[256 * 1024 * 1024];
	if (!texturesTmp) goto cleanup;

	if (op) {
		Vector absPos = op->GetAbsPos();
		Vector froPos = op->GetFrozenPos();
		Vector relPos = op->GetRelPos();
		BaseContainer bc = op->GetData();
		CountPolygonsRecursive(op, &vertexCount, &triangleCount, &lightCount);
		if (0 == vertexCount || 0 == triangleCount) goto cleanup;

		vertex = new cl_float3[vertexCount];
		if (!vertex) goto cleanup;
		memset(vertex, 0, vertexCount * sizeof(cl_float3));

		triangleVertexIndex = new cl_int3[triangleCount];
		if (!triangleVertexIndex) goto cleanup;
		memset(triangleVertexIndex, 0, triangleCount * sizeof(cl_int3));
		triangleMaterialId = new cl_int[triangleCount];
		if (!triangleMaterialId) goto cleanup;
		memset(triangleMaterialId, -1, triangleCount * sizeof(cl_int));
		triangleUv = new cl_float2[3 * triangleCount]; // 3x triangleCount
		if (!triangleUv) goto cleanup;
		memset(triangleUv, 0, 3 * triangleCount * sizeof(cl_float2));
		triangleNormal = new cl_float3[3 * triangleCount]; // 3x triangleCount
		if (!triangleNormal) goto cleanup;
		memset(triangleNormal, 0, 3 * triangleCount * sizeof(cl_float3));

		lightType = new cl_int[lightCount];
		if (!lightType) goto cleanup;
		memset(lightType, 0, lightCount * sizeof(cl_int));
		lightPosition = new cl_float3[lightCount];
		if (!lightPosition) goto cleanup;
		memset(lightPosition, 0, lightCount * sizeof(cl_float3));
		lightDirection = new cl_float3[lightCount];
		if (!lightDirection) goto cleanup;
		memset(lightDirection, 0, lightCount * sizeof(cl_float3));
		lightColour = new cl_float3[lightCount];
		if (!lightColour) goto cleanup;
		memset(lightColour, 0, lightCount * sizeof(cl_float3));
		lightRadius = new cl_float[lightCount];
		if (!lightRadius) goto cleanup;
		memset(lightRadius, 0, lightCount * sizeof(cl_float));
		lightHalfAttenuationDistance = new cl_float[lightCount];
		if (!lightHalfAttenuationDistance) goto cleanup;
		memset(lightHalfAttenuationDistance, 0, lightCount * sizeof(cl_float));

		{
			cl_uint vertexCursor = 0, triangleCursor = 0, lightCursor = 0;
			AddPolygonsRecursive(op, &vertexCursor, &triangleCursor, &lightCursor, vertex, triangleVertexIndex, triangleMaterialId, triangleUv, triangleNormal, lightType, lightPosition, lightDirection, lightColour, lightRadius, lightHalfAttenuationDistance, materialDatabase, cameraEye, materialTextureRatio);
			DebugAssert(vertexCursor == vertexCount, "Didn't add the right amount of vertices");
			vertexCount = vertexCursor;
			DebugAssert(triangleCursor <= triangleCount, "Didn't add the right amount of triangles");
			triangleCount = triangleCursor;
			DebugAssert(lightCursor == lightCount, "Didn't add the right amount of lights");
			lightCount = lightCursor;
		}

		Int32 channels[] = { CHANNEL_COLOR, CHANNEL_REFLECTION, CHANNEL_TRANSPARENCY, CHANNEL_BUMP, CHANNEL_LUMINANCE };

		material = (Material*)bm;
		while (material) {
			std::map<BaseMaterial*, Int32>::const_iterator matIt = materialDatabase.find(material);
			if (materialDatabase.cend() != matIt) {
				for (Int32 ci = 0; ci < sizeof(channels) / sizeof(Int32); ++ci) {
					BaseChannel* channel = material->GetChannel(channels[ci]);
					Bool enabled = material->GetChannelState(channels[ci]);
					if (!channel || !enabled) {
						materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] = 0;
						materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1] = 0;
						materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + ci] = texturesCursor;
					}
					else {
						Filename path;
						if (GenerateTexturePath(doc->GetDocumentPath(), Filename(channel->GetData().GetString(BASECHANNEL_TEXTURE)), Filename(), &path)) {
							Char* texName = path.GetString().GetCStringCopy();
							texName = texName;
						}
						InitRenderStruct irs;
						irs.version = GetC4DVersion();
						irs.fps = doc->GetFps();
						irs.time = doc->GetTime();
						irs.docpath = doc->GetDocumentPath();
						INITRENDERRESULT res = channel->InitTexture(irs);
						if (INITRENDERRESULT_OK == res) {
							BaseBitmap* bitmap = channel->GetBitmap();
							BaseShader* shader = channel->GetShader();
							if (bitmap) {
								Int32 bw = bitmap->GetBw();
								Int32 bh = bitmap->GetBh();
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] = bw;
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1] = bh;
								materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + ci] = texturesCursor;
								cl_uchar3* imageCursor = &(texturesTmp[texturesCursor]);
								for (Int32 y = 0; y < bh; ++y) {
									for (Int32 x = 0; x < bw; ++x) {
										UInt16 r, g, b;
										bitmap->GetPixel(x, y, &r, &g, &b);
										imageCursor->s[0] = (cl_uchar)r;
										imageCursor->s[1] = (cl_uchar)g;
										imageCursor->s[2] = (cl_uchar)b;
										imageCursor++;
									}
								}
								texturesCursor += bw * bh;
							}
							else if (shader) {
								ChannelData cd;
								cd.p = Vector(0., 0., 0.);
								cd.n = Vector(0., 0., 1.);
								cd.d = Vector(0., 0., 0.);
								cd.t = 0.;
								cd.texflag = 0;
								cd.vd = nullptr;
								cd.off = 0.;
								cd.scale = 0.;

								// Make the texture pixel size proportional to the closest point of the closest triangle with this mat
								cl_float pixelSizeAtDist = materialTextureRatio[matIt->second] / (cl_float)sin(fov / fmax(outputImageDimension.s[0], outputImageDimension.s[1]));
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] = (cl_uint)bindf(pixelSizeAtDist, 1.f, 4096.f);
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1] = materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0];
								materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + ci] = texturesCursor;
								cl_uchar3* imageCursor = &(texturesTmp[texturesCursor]);
								Float ds0 = 1. / (Float)(materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0]);
								Float ds1 = 1. / (Float)(materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1]);
								Float s1 = 0.5 * ds1;
								for (cl_uint y = 0; y < materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1]; ++y) {
									Float s0 = 0.5 * ds0;
									for (cl_uint x = 0; x < materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0]; ++x) {
										cd.p.x = s0;
										s0 += ds0;
										cd.p.y = s1;
										Vector color = shader->Sample(&cd);
										imageCursor->s[0] = (cl_uchar)floor(0.5f + (255.f * bindf((cl_float)color.x, 0.f, 1.f)));
										imageCursor->s[1] = (cl_uchar)floor(0.5f + (255.f * bindf((cl_float)color.y, 0.f, 1.f)));
										imageCursor->s[2] = (cl_uchar)floor(0.5f + (255.f * bindf((cl_float)color.z, 0.f, 1.f)));
										imageCursor++;
									}
									s1 += ds1;
								}
								texturesCursor += materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] * materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1];
							}
							else if (CHANNEL_REFLECTION == channels[ci]) {
								cl_float reflectance = 0.2f;
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] = 1;
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1] = 1;
								materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + ci] = texturesCursor;
								texturesTmp[texturesCursor].s[0] = (cl_uchar)floor(0.5f + reflectance * 255.f);
								texturesTmp[texturesCursor].s[1] = (cl_uchar)floor(0.5f + reflectance * 255.f);
								texturesTmp[texturesCursor].s[2] = (cl_uchar)floor(0.5f + reflectance * 255.f);
								texturesCursor += 1;
							}
							else if (CHANNEL_TRANSPARENCY == channels[ci]) {
								cl_float transparency = 1.f;
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] = 1;
								materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1] = 1;
								materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + ci] = texturesCursor;
								texturesTmp[texturesCursor].s[0] = (cl_uchar)floor(0.5f + transparency * 255.f);
								texturesTmp[texturesCursor].s[1] = (cl_uchar)floor(0.5f + transparency * 255.f);
								texturesTmp[texturesCursor].s[2] = (cl_uchar)floor(0.5f + transparency * 255.f);
								texturesCursor += 1;
							}
							channel->FreeTexture();
						}
					}
					if (CHANNEL_COLOR != channels[ci] && 0 == materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0]) {
						materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[0] = 1;
						materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + ci].s[1] = 1;
						materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + ci] = texturesCursor;
						texturesTmp[texturesCursor].s[0] = 0;
						texturesTmp[texturesCursor].s[1] = 0;
						texturesTmp[texturesCursor].s[2] = 0;
						texturesCursor += 1;
					}
				}

				if (0 == materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + 0].s[0]) {
					cl_float3 color = { 1.f, 1.f, 1.f };
					GeData d;
					if (material->GetParameter(DescID(MATERIAL_COLOR_COLOR), d, DESCFLAGS_GET_0)) {
						Vector c = d.GetVector();
						color.s[0] = (cl_float)c.x;
						color.s[1] = (cl_float)c.y;
						color.s[2] = (cl_float)c.z;
						if (material->GetParameter(DescID(MATERIAL_COLOR_BRIGHTNESS), d, DESCFLAGS_GET_0)) {
							Float b = d.GetFloat();
							color.s[0] *= (cl_float)b;
							color.s[1] *= (cl_float)b;
							color.s[2] *= (cl_float)b;
						}
					}
					materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + 0].s[0] = 1;
					materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + 0].s[1] = 1;
					materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + 0] = texturesCursor;
					texturesTmp[texturesCursor].s[0] = (cl_uchar)floor(0.5f + (Float)color.s[0] * 255.f);
					texturesTmp[texturesCursor].s[1] = (cl_uchar)floor(0.5f + (Float)color.s[1] * 255.f);
					texturesTmp[texturesCursor].s[2] = (cl_uchar)floor(0.5f + (Float)color.s[2] * 255.f);
					texturesCursor += 1;
				}
				if (0 == materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + 2].s[0]) {
					cl_float3 color = { 1.f, 1.f, 1.f };
					GeData d;
					if (material->GetParameter(DescID(MATERIAL_TRANSPARENCY_COLOR), d, DESCFLAGS_GET_0)) {
						Vector c = d.GetVector();
						color.s[0] = (cl_float)c.x;
						color.s[1] = (cl_float)c.y;
						color.s[2] = (cl_float)c.z;
						if (material->GetParameter(DescID(MATERIAL_TRANSPARENCY_BRIGHTNESS), d, DESCFLAGS_GET_0)) {
							Float b = d.GetFloat();
							color.s[0] *= (cl_float)b;
							color.s[1] *= (cl_float)b;
							color.s[2] *= (cl_float)b;
						}
					}
					materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + 2].s[0] = 1;
					materialImageSize[_MATERIAL_CHANNEL_COUNT * matIt->second + 2].s[1] = 1;
					materialImageStart[_MATERIAL_CHANNEL_COUNT * matIt->second + 2] = texturesCursor;
					texturesTmp[texturesCursor].s[0] = (cl_uchar)floor(0.5f + (Float)color.s[0] * 255.f);
					texturesTmp[texturesCursor].s[1] = (cl_uchar)floor(0.5f + (Float)color.s[1] * 255.f);
					texturesTmp[texturesCursor].s[2] = (cl_uchar)floor(0.5f + (Float)color.s[2] * 255.f);
					texturesCursor += 1;
				}
			}
			material = (Material*)material->GetNext();
		}
	}
	polygonizedDoc->Free(polygonizedDoc);
	
	materialImageStart[materialCount * _MATERIAL_CHANNEL_COUNT] = texturesCursor;
	cl_uchar3* textures = new cl_uchar3[texturesCursor];
	memcpy(textures, texturesTmp, texturesCursor * sizeof(cl_uchar3));
	if (texturesTmp) delete[] texturesTmp; texturesTmp = NULL;

	cameraTriangleList = CameraTriangleList::New (outputImageDimension, cameraEye, cameraEyeToTopLeftVector, cameraLeftToRightPixelSizeVector, cameraTopToBottomPixelSizeVector, cameraPixelSizeInv, vertexCount, triangleCount, vertex, triangleVertexIndex);
	sceneTriangleList = SceneTriangleList::New (vertexCount, triangleCount, vertex, triangleVertexIndex);

	computationDone =
	RaytraceAll(platformDevice,
				outputImageDimension,
				cameraEye,
				cameraEyeToTopLeftVector,
				cameraLeftToRightPixelSizeVector,
				cameraTopToBottomPixelSizeVector,
				cameraPixelSizeInv,
				cameraTriangleList->GetTriangleListStart(),
				cameraTriangleList->GetTriangleListEnd(),
				cameraTriangleList->GetTriangleList(),
				cameraTriangleList->GetTriangleListSize(),
				sampleCount,
				vertexCount,
				vertex,
				triangleCount,
				triangleVertexIndex,
				triangleMaterialId,
				triangleUv,
				triangleNormal,
				SceneTriangleList::AXES_DIVISION,
				sceneTriangleList->GetBoxMin(),
				sceneTriangleList->GetTriangleListStart(),
				sceneTriangleList->GetTriangleList(),
				materialCount,
				materialImageSize,
				materialImageStart,
				texturesCursor,
				textures,
				lightCount,
				lightType,
				lightPosition,
				lightDirection,
				lightColour,
				lightRadius,
				lightHalfAttenuationDistance,
				outputRed,
				outputGreen,
				outputBlue);

cleanup:
	if (vertex) delete[] vertex; vertex = NULL;
	if (triangleVertexIndex) delete[] triangleVertexIndex; triangleVertexIndex = NULL;
	if (triangleMaterialId) delete[] triangleMaterialId; triangleMaterialId = NULL;
	if (triangleUv) delete[] triangleUv; triangleUv = NULL;
	if (triangleNormal) delete[] triangleNormal; triangleNormal = NULL;
	if (materialImageSize) delete[] materialImageSize; materialImageSize = NULL;
	if (materialImageStart) delete[] materialImageStart; materialImageStart = NULL;
	if (materialTextureRatio) delete[] materialTextureRatio; materialTextureRatio = NULL;
	if (lightType) delete[] lightType; lightType = NULL;
	if (lightPosition) delete[] lightPosition; lightPosition = NULL;
	if (lightDirection) delete[] lightDirection; lightDirection = NULL;
	if (lightColour) delete[] lightColour; lightColour = NULL;
	if (lightRadius) delete[] lightRadius; lightRadius = NULL;
	if (lightHalfAttenuationDistance) delete[] lightHalfAttenuationDistance; lightHalfAttenuationDistance = NULL;
	if (cameraTriangleList) delete cameraTriangleList; cameraTriangleList = NULL;
	if (sceneTriangleList) delete sceneTriangleList; sceneTriangleList = NULL;

	if (computationDone) {
		writebmp3s(outputImageDimension, outputRed, outputGreen, outputBlue);
		*bitmap = BaseBitmap::Alloc();
		(*bitmap)->Init(outputImageDimension.s[0], outputImageDimension.s[1], 24, INITBITMAPFLAGS_0);
		if (bitmap) {
			for (cl_uint x = 0; x < outputImageDimension.s[0]; x++) {
				for (cl_uint y = 0; y < outputImageDimension.s[1]; y++) {
					(*bitmap)->SetPixel(x, y,
						(int)(outputRed[x + y * outputImageDimension.s[0]] / 256),
						(int)(outputGreen[x + y * outputImageDimension.s[0]] / 256),
						(int)(outputBlue[x + y * outputImageDimension.s[0]] / 256));
				}
			}
		}
	}
	else computationDone = false;

	if (outputRed) delete[] outputRed; outputRed = NULL;
	if (outputGreen) delete[] outputGreen; outputGreen = NULL;
	if (outputBlue) delete[] outputBlue; outputBlue = NULL;

//	GeShowMouse(MOUSE_NORMAL);

//	MessageDialog("Done");

	SetProgress(1.f);

	return computationDone;
}
