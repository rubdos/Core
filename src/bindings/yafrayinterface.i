%module yafrayinterface

%include "cpointer.i"
%pointer_functions(float, floatp);
%pointer_functions(int, intp);
%pointer_functions(unsigned int, uintp);

%include "carrays.i"
%include "std_string.i"
%include "std_vector.i"

%array_functions(float, floatArray);

namespace std
{
	%template(StrVector) vector<string>;
}

#ifdef SWIGPYTHON  // Begining of python specific code

%{
#include <yafraycore/monitor.h>
#include <core_api/output.h>
#include <interface/yafrayinterface.h>

struct yafTilePixel4_t
{
	float r;
	float g;
	float b;
	float a;
};

struct yafTilePixel3_t
{
	float x;
	float y;
	float z;
};

struct yafTilePixel1_t
{
	float v;
};

struct YafTile4Object_t
{
	PyObject_HEAD
	int resx, resy;
	int x0, x1, y0, y1;
	int w, h;
	yafTilePixel4_t *mem;
};

struct YafTile3Object_t
{
	PyObject_HEAD
	int resx, resy;
	int x0, x1, y0, y1;
	int w, h;
	yafTilePixel3_t *mem;
};

struct YafTile1Object_t
{
	PyObject_HEAD
	int resx, resy;
	int x0, x1, y0, y1;
	int w, h;
	yafTilePixel1_t *mem;
};


static Py_ssize_t yaf_tile_4_length(YafTile4Object_t *self)
{
	self->w = (self->x1 - self->x0);
	self->h = (self->y1 - self->y0);

	return self->w * self->h;
}

static PyObject *yaf_tile_4_subscript_int(YafTile4Object_t *self, int keynum)
{
	// Check boundaries and fill w and h
	if (keynum >= yaf_tile_4_length(self) || keynum < 0)
		return Py_BuildValue("[f,f,f,f]", 1, 0, 0, 1);

	// Calc position of the tile in the image region
	int vy = keynum / self->w;
	int vx = keynum - vy * self->w;

	// Map tile position to image buffer
	vx = self->x0 + vx;
	vy = (self->y0 + self->h - 1) - vy;

	// Get pixel
	yafTilePixel4_t &pix = self->mem[ self->resx * vy + vx ];

	return Py_BuildValue("[f,f,f,f]", pix.r, pix.g, pix.b, pix.a);
}

static void yaf_tile_4_dealloc(YafTile4Object_t *self)
{
	PyObject_Del(self);
}

PySequenceMethods sequence_methods_4 =
{
	( lenfunc ) yaf_tile_4_length,
	NULL,
	NULL,
	( ssizeargfunc ) yaf_tile_4_subscript_int
};

PyTypeObject yafTile4_Type =
{
	PyVarObject_HEAD_INIT(NULL, 0)
	"yaf_tile_4",						/* tp_name */
	sizeof(YafTile4Object_t),			/* tp_basicsize */
	0,									/* tp_itemsize */
	( destructor ) yaf_tile_4_dealloc,	/* tp_dealloc */
	NULL,                       		/* printfunc tp_print; */
	NULL,								/* getattrfunc tp_getattr; */
	NULL,								/* setattrfunc tp_setattr; */
	NULL,								/* tp_compare */ /* DEPRECATED in python 3.0! */
	NULL,								/* tp_repr */
	NULL,                       		/* PyNumberMethods *tp_as_number; */
	&sequence_methods_4,				/* PySequenceMethods *tp_as_sequence; */
	NULL,								/* PyMappingMethods *tp_as_mapping; */
	NULL,								/* hashfunc tp_hash; */
	NULL,								/* ternaryfunc tp_call; */
	NULL,                       		/* reprfunc tp_str; */
	NULL,								/* getattrofunc tp_getattro; */
	NULL,								/* setattrofunc tp_setattro; */
	NULL,                       		/* PyBufferProcs *tp_as_buffer; */
	Py_TPFLAGS_DEFAULT,         		/* long tp_flags; */
};


static Py_ssize_t yaf_tile_1_length(YafTile1Object_t *self)
{
	self->w = (self->x1 - self->x0);
	self->h = (self->y1 - self->y0);

	return self->w * self->h;
}

static PyObject *yaf_tile_1_subscript_int(YafTile1Object_t *self, int keynum)
{
	// Check boundaries and fill w and h
	if (keynum >= yaf_tile_1_length(self) || keynum < 0)
		return Py_BuildValue("[f]", 0);

	// Calc position of the tile in the image region
	int vy = keynum / self->w;
	int vx = keynum - vy * self->w;

	// Map tile position to image buffer
	vx = self->x0 + vx;
	vy = (self->y0 + self->h - 1) - vy;

	// Get pixel
	yafTilePixel1_t &pix = self->mem[ self->resx * vy + vx ];

	return Py_BuildValue("[f]", pix.v);
}

static void yaf_tile_1_dealloc(YafTile1Object_t *self)
{
	PyObject_Del(self);
}

PySequenceMethods sequence_methods_1 =
{
	( lenfunc ) yaf_tile_1_length,
	NULL,
	NULL,
	( ssizeargfunc ) yaf_tile_1_subscript_int
};

PyTypeObject yafTile1_Type =
{
	PyVarObject_HEAD_INIT(NULL, 0)
	"yaf_tile_1",						/* tp_name */
	sizeof(YafTile1Object_t),			/* tp_basicsize */
	0,									/* tp_itemsize */
	( destructor ) yaf_tile_1_dealloc,	/* tp_dealloc */
	NULL,                       		/* printfunc tp_print; */
	NULL,								/* getattrfunc tp_getattr; */
	NULL,								/* setattrfunc tp_setattr; */
	NULL,								/* tp_compare */ /* DEPRECATED in python 3.0! */
	NULL,								/* tp_repr */
	NULL,                       		/* PyNumberMethods *tp_as_number; */
	&sequence_methods_1,				/* PySequenceMethods *tp_as_sequence; */
	NULL,								/* PyMappingMethods *tp_as_mapping; */
	NULL,								/* hashfunc tp_hash; */
	NULL,								/* ternaryfunc tp_call; */
	NULL,                       		/* reprfunc tp_str; */
	NULL,								/* getattrofunc tp_getattro; */
	NULL,								/* setattrofunc tp_setattro; */
	NULL,                       		/* PyBufferProcs *tp_as_buffer; */
	Py_TPFLAGS_DEFAULT,         		/* long tp_flags; */
};


static Py_ssize_t yaf_tile_3_length(YafTile3Object_t *self)
{
	self->w = (self->x1 - self->x0);
	self->h = (self->y1 - self->y0);

	return self->w * self->h;
}

static PyObject *yaf_tile_3_subscript_int(YafTile3Object_t *self, int keynum)
{
	// Check boundaries and fill w and h
	if (keynum >= yaf_tile_3_length(self) || keynum < 0)
		return Py_BuildValue("[f,f,f]", 0, 0, 0);

	// Calc position of the tile in the image region
	int vy = keynum / self->w;
	int vx = keynum - vy * self->w;

	// Map tile position to image buffer
	vx = self->x0 + vx;
	vy = (self->y0 + self->h - 1) - vy;

	// Get pixel
	yafTilePixel3_t &pix = self->mem[ self->resx * vy + vx ];

	return Py_BuildValue("[f,f,f]", pix.x, pix.y, pix.z);
}

static void yaf_tile_3_dealloc(YafTile3Object_t *self)
{
	PyObject_Del(self);
}

PySequenceMethods sequence_methods_3 =
{
	( lenfunc ) yaf_tile_3_length,
	NULL,
	NULL,
	( ssizeargfunc ) yaf_tile_3_subscript_int
};

PyTypeObject yafTile3_Type =
{
	PyVarObject_HEAD_INIT(NULL, 0)
	"yaf_tile_3",						/* tp_name */
	sizeof(YafTile3Object_t),			/* tp_basicsize */
	0,									/* tp_itemsize */
	( destructor ) yaf_tile_3_dealloc,	/* tp_dealloc */
	NULL,                       		/* printfunc tp_print; */
	NULL,								/* getattrfunc tp_getattr; */
	NULL,								/* setattrfunc tp_setattr; */
	NULL,								/* tp_compare */ /* DEPRECATED in python 3.0! */
	NULL,								/* tp_repr */
	NULL,                       		/* PyNumberMethods *tp_as_number; */
	&sequence_methods_3,				/* PySequenceMethods *tp_as_sequence; */
	NULL,								/* PyMappingMethods *tp_as_mapping; */
	NULL,								/* hashfunc tp_hash; */
	NULL,								/* ternaryfunc tp_call; */
	NULL,                       		/* reprfunc tp_str; */
	NULL,								/* getattrofunc tp_getattro; */
	NULL,								/* setattrofunc tp_setattro; */
	NULL,                       		/* PyBufferProcs *tp_as_buffer; */
	Py_TPFLAGS_DEFAULT,         		/* long tp_flags; */
};


class pyOutput_t : public yafaray::colorOutput_t
{

public:

	pyOutput_t(int x, int y, int borderStartX, int borderStartY, bool prev, PyObject *drawAreaCallback, PyObject *flushCallback) :
	resx(x),
	resy(y),
	bsX(borderStartX),
	bsY(borderStartY),
	preview(prev),
	mDrawArea(drawAreaCallback),
	mFlush(flushCallback)
	{
		tile = PyObject_New(YafTile4Object_t, &yafTile4_Type);
		dTile = PyObject_New(YafTile1Object_t, &yafTile1_Type);

		tile->mem = new yafTilePixel4_t[x*y];
		tile->resx = x;
		tile->resy = y;

		dTile->mem = new yafTilePixel1_t[x*y];
		dTile->resx = x;
		dTile->resy = y;
	}

	virtual ~pyOutput_t()
	{
		delete [] tile->mem;
		delete [] dTile->mem;
		Py_XDECREF(tile);
		Py_XDECREF(dTile);
	}

	virtual bool putPixel(int x, int y, const float *c, bool alpha = true, bool depth = false, float z = 0.f)
	{
		yafTilePixel4_t &pix= tile->mem[resx * y + x];
		pix.r = c[0];
		pix.g = c[1];
		pix.b = c[2];
		pix.a = alpha ? c[3] : 1.0f;

        yafTilePixel1_t &dpix = dTile->mem[resx * y + x];
        dpix.v = depth ? z : 1.f;

		return true;
	}

	virtual void flush()
	{
		tile->x0 = 0;
		tile->x1 = resx;
		tile->y0 = 0;
		tile->y1 = resy;

		dTile->x0 = 0;
		dTile->x1 = resx;
		dTile->y0 = 0;
		dTile->y1 = resy;

		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyEval_CallObject(mFlush, Py_BuildValue("ii(OO)", resx, resy, tile, dTile));
		PyGILState_Release(gstate);
	}

	virtual void flushArea(int x0, int y0, int x1, int y1)
	{
		// Do nothing if we are rendering preview renders
		if(preview) return;

		tile->x0 = x0 - bsX;
		tile->x1 = x1 - bsX;
		tile->y0 = y0 - bsY;
		tile->y1 = y1 - bsY;

		dTile->x0 = x0 - bsX;
		dTile->x1 = x1 - bsX;
		dTile->y0 = y0 - bsY;
		dTile->y1 = y1 - bsY;

		int w = x1 - x0;
		int h = y1 - y0;

		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyEval_CallObject(mDrawArea, Py_BuildValue("iiii(OO)", tile->x0, resy - tile->y1, w, h, tile, dTile));
		PyGILState_Release(gstate);
	}

	virtual void highliteArea(int x0, int y0, int x1, int y1)
	{
		// Do nothing if we are rendering preview renders
		if(preview) return;

		tile->x0 = x0 - bsX;
		tile->x1 = x1 - bsX;
		tile->y0 = y0 - bsY;
		tile->y1 = y1 - bsY;

		dTile->x0 = x0 - bsX;
		dTile->x1 = x1 - bsX;
		dTile->y0 = y0 - bsY;
		dTile->y1 = y1 - bsY;

		int w = x1 - x0;
		int h = y1 - y0;
		int lineL = std::min( 4, std::min( h - 1, w - 1 ) );

		drawCorner(tile->x0, tile->y0, lineL, TL_CORNER);
		drawCorner(tile->x1, tile->y0, lineL, TR_CORNER);
		drawCorner(tile->x0, tile->y1, lineL, BL_CORNER);
		drawCorner(tile->x1, tile->y1, lineL, BR_CORNER);

		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyEval_CallObject(mDrawArea, Py_BuildValue("iiii(OO)", tile->x0, resy - tile->y1, w, h, tile, dTile));
		PyGILState_Release(gstate);
	}

private:

	enum corner
	{
		TL_CORNER,
		TR_CORNER,
		BL_CORNER,
		BR_CORNER
	};

	void drawCorner(int x, int y, int len, corner pos)
	{
		int minX = 0;
		int minY = 0;
		int maxX = 0;
		int maxY = 0;

		switch(pos)
		{
			case TL_CORNER:
				minX = x;
				minY = y;
				maxX = x + len;
				maxY = y + len;
				break;

			case TR_CORNER:
				minX = x - len - 1;
				minY = y;
				maxX = x - 1;
				maxY = y + len;
				--x;
				break;

			case BL_CORNER:
				minX = x;
				minY = y - len - 1;
				maxX = x + len;
				maxY = y - 1;
				--y;
				break;

			case BR_CORNER:
				minX = x - len - 1;
				minY = y - len - 1;
				maxX = x;
				maxY = y - 1;
				--x;
				--y;
				break;
		}

		for(int i = minX; i < maxX; ++i)
		{
			yafTilePixel4_t &pix = tile->mem[resx * y + i];
			pix.r = 0.625f;
			pix.g = 0.f;
			pix.b = 0.f;
			pix.a = 1.f;
		}

		for(int j = minY; j < maxY; ++j)
		{
			yafTilePixel4_t &pix = tile->mem[resx * j + x];
			pix.r = 0.625f;
			pix.g = 0.f;
			pix.b = 0.f;
			pix.a = 1.f;
		}
	}

	int resx, resy;
	int bsX, bsY;
	bool preview;
	PyObject *mDrawArea;
	PyObject *mFlush;
	YafTile4Object_t *tile;
	YafTile1Object_t *dTile;
};

class pyProgress : public yafaray::progressBar_t
{

public:

	pyProgress(PyObject *callback) : callb(callback) {}

	void report_progress(float percent)
	{
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyEval_CallObject(callb, Py_BuildValue("sf", "progress", percent));
		PyGILState_Release(gstate);
	}

	virtual void init(int totalSteps)
	{
		steps_to_percent = 1.f / (float) totalSteps;
		doneSteps = 0;
		report_progress(0.f);
	}

	virtual void update(int steps = 1)
	{
		doneSteps += steps;
		report_progress(doneSteps * steps_to_percent);
	}

	virtual void done()
	{
		report_progress(1.f);
	}

	virtual void setTag(const char* text)
	{
		PyGILState_STATE gstate;
		gstate = PyGILState_Ensure();
		PyEval_CallObject(callb, Py_BuildValue("ss", "tag", text));
		PyGILState_Release(gstate);
	}

private:

	PyObject *callb;
	float steps_to_percent;
	int doneSteps;
};

%}

%init %{
	PyType_Ready(&yafTile4_Type);
	PyType_Ready(&yafTile1_Type);
	PyType_Ready(&yafTile3_Type);
%}

%typemap(in) PyObject *pyfunc
{
	if (!PyCallable_Check($input))
	{
		PyErr_SetString(PyExc_TypeError, "Need a callback method.");
		return NULL;
	}

	$1 = $input;
}

%extend yafaray::yafrayInterface_t
{
	void render(int x, int y, int borderStartX, int borderStartY, bool prev, PyObject *drawAreaCallBack, PyObject *flushCallBack, PyObject *progressCallback)
	{
		pyOutput_t output_wrap(x, y, borderStartX, borderStartY, prev, drawAreaCallBack, flushCallBack);
		pyProgress *pbar_wrap = new pyProgress(progressCallback);

		Py_BEGIN_ALLOW_THREADS;
		self->render(output_wrap, pbar_wrap);
		Py_END_ALLOW_THREADS;
	}
    void addTriangles(PyObject *pyArray)
    {
        if(!PyList_Check(pyArray)) {
            PyErr_SetString(PyExc_TypeError,
                "yafrayInterface_t: "
                "addTriangles expects a list.");
            return;
        }
        const size_t count = PyList_Size(pyArray);
        PyObject *obj;
        for(size_t i = 0; i < count; ++i) {
            obj = PyList_GetItem(pyArray, i);
            if(!PyTuple_Check(obj)) {
                PyErr_SetString(PyExc_TypeError,
                    "yafrayInterface_t: "
                    "addVertices expects a list of tuples.");
                std::cout << "Element " << i << " is not a list" << std::endl;
                return;
            }
            const size_t innersize = PyTuple_Size(obj); 
            if(innersize == 2) {
                PyObject *verts = PyTuple_GetItem(obj, 0);
                PyObject *ymaterial = PyTuple_GetItem(obj, 1);
                const material_t *mat = reinterpret_cast<material_t *>(ymaterial);
                self->addTriangle(
                        PyInt_AsLong(PyList_GetItem(verts, 0)),
                        PyInt_AsLong(PyList_GetItem(verts, 1)),
                        PyInt_AsLong(PyList_GetItem(verts, 2)),
                        mat
                        );
                if (PyList_Size(verts) == 4) { // We got a quad instead of a triangle.
                                               // Triangulate
                    self->addTriangle(
                            PyInt_AsLong(PyList_GetItem(verts, 0)),
                            PyInt_AsLong(PyList_GetItem(verts, 2)),
                            PyInt_AsLong(PyList_GetItem(verts, 3)),
                            mat
                            );
                }
            } else if(innersize == 3) {
                PyObject *verts = PyTuple_GetItem(obj, 0);
                PyObject *uv = PyTuple_GetItem(obj, 1);
                PyObject *ymaterial = PyTuple_GetItem(obj, 2);
                const material_t *mat = reinterpret_cast<material_t *>(ymaterial);
                self->addTriangle(
                        PyInt_AsLong(PyList_GetItem(verts, 0)),
                        PyInt_AsLong(PyList_GetItem(verts, 1)),
                        PyInt_AsLong(PyList_GetItem(verts, 2)),
                        PyInt_AsLong(PyTuple_GetItem(uv, 0)),
                        PyInt_AsLong(PyTuple_GetItem(uv, 1)),
                        PyInt_AsLong(PyTuple_GetItem(uv, 2)),
                        mat
                        );
                if (PyList_Size(verts) == 4) { // We got a quad instead of a triangle.
                                               // Triangulate
                    self->addTriangle(
                            PyInt_AsLong(PyList_GetItem(verts, 0)),
                            PyInt_AsLong(PyList_GetItem(verts, 2)),
                            PyInt_AsLong(PyList_GetItem(verts, 3)),
                            PyInt_AsLong(PyTuple_GetItem(uv, 0)),
                            PyInt_AsLong(PyTuple_GetItem(uv, 2)),
                            PyInt_AsLong(PyTuple_GetItem(uv, 3)),
                            mat
                            );
                }
            } else {
                PyErr_SetString(PyExc_TypeError,
                    "yafrayInterface_t: "
                    "addTriangles expects inner tuples with sizes 2 or 3.");
                return;
            }
        }
    }
    void addVertices(PyObject *pyArray)
    {
        // Structure should be:
        // list of tuples.
        // The inner tuples should have 3 or 6 items:
        // three coordinates and optionally three Orco coordinates.
        if(!PyList_Check(pyArray)) {
            PyErr_SetString(PyExc_TypeError,
                "yafrayInterface_t: "
                "addVertices expects a list.");
            return;
        }
        const size_t count = PyList_Size(pyArray);
        PyObject *obj;
        for(size_t i = 0; i < count; ++i) {
            obj = PyList_GetItem(pyArray, i);
            if(!PyTuple_Check(obj)) {
                PyErr_SetString(PyExc_TypeError,
                    "yafrayInterface_t: "
                    "addVertices expects a list of tuples.");
                std::cout << "Element " << i << " is not a tuple" << std::endl;
                return;
            }
            const size_t innersize = PyTuple_Size(obj); 
            if (innersize < 3) {
                PyErr_SetString(PyExc_TypeError,
                    "yafrayInterface_t: "
                    "addVertices expects inner tuples with sizes 3 or 6.");
                return;
            }
            for(size_t i = 0; i < innersize; ++i) {
                if(!PyFloat_Check(PyTuple_GetItem(obj, i))) {
                    PyErr_SetString(PyExc_TypeError,
                        "yafrayInterface_t: "
                        "addVertices expects tuples of Floats.");
                    return;
                }
            }
            if (innersize >= 6) {
                self->addVertex(
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 0)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 1)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 2)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 3)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 4)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 5))
                        );
            } else {
                self->addVertex(
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 0)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 1)),
                        PyFloat_AsDouble(PyTuple_GetItem(obj, 2))
                        );
            }
        }
    }
}

%exception yafaray::yafrayInterface_t::loadPlugins
{
	Py_BEGIN_ALLOW_THREADS
	$action
	Py_END_ALLOW_THREADS
}

#endif // End of python specific code

%{
#include <yafray_constants.h>
#include <interface/yafrayinterface.h>
#include <interface/xmlinterface.h>
#include <yafraycore/imageOutput.h>
#include <yafraycore/memoryIO.h>
#include <core_api/matrix4.h>
using namespace yafaray;
%}

#ifdef SWIGRUBY
%feature("director") colorOutput_t::putPixel;
%feature("director") colorOutput_t::flush;
%feature("director") colorOutput_t::flushArea;
#endif

namespace yafaray
{
	// Required abstracts

	class colorOutput_t
	{
		public:
			virtual ~colorOutput_t() {};
			virtual bool putPixel(int x, int y, const float *c, bool alpha = true, bool depth = false, float z = 0.f)=0;
			virtual void flush()=0;
			virtual void flushArea(int x0, int y0, int x1, int y1)=0;
	};

	class imageHandler_t
	{
		public:
			virtual void initForOutput(int width, int height, bool withAlpha = false, bool withDepth = true) = 0;
			virtual ~imageHandler_t() {};
			virtual bool loadFromFile(const std::string &name) = 0;
			virtual bool loadFromMemory(const yByte *data, size_t size) {return false; };
			virtual bool saveToFile(const std::string &name) = 0;
			virtual void putPixel(int x, int y, const colorA_t &rgba, float depth = 0.f) = 0;
			virtual colorA_t getPixel(int x, int y) = 0;

		protected:
			std::string handlerName;
			int m_width;
			int m_height;
			bool m_hasAlpha;
			bool m_hasDepth;
			rgba2DImage_nw_t *m_rgba;
			gray2DImage_nw_t *m_depth;
	};

	// Outputs

	class imageOutput_t : public colorOutput_t
	{
		public:
			imageOutput_t(imageHandler_t *handle, const std::string &name, int bx, int by);
			imageOutput_t(); //!< Dummy initializer
			virtual ~imageOutput_t();
			virtual bool putPixel(int x, int y, const float *c, bool alpha = true, bool depth = false, float z = 0.f);
			virtual void flush();
			virtual void flushArea(int x0, int y0, int x1, int y1) {}; // not used by images... yet
		private:
			imageHandler_t *image;
			std::string fname;
	};

	class memoryIO_t : public colorOutput_t
	{
		public:
			memoryIO_t(int resx, int resy, float* iMem);
			virtual bool putPixel(int x, int y, const float *c, bool alpha = true, bool depth = false, float z = 0.f);
			void flush();
			virtual void flushArea(int x0, int y0, int x1, int y1) {}; // no tiled file format used...yet
			virtual ~memoryIO_t();
		protected:
			int sizex, sizey;
			float* imageMem;
	};

	// Utility classes

	class matrix4x4_t
	{
		public:
			matrix4x4_t() {};
			matrix4x4_t(const float init);
			matrix4x4_t(const matrix4x4_t & source);
			matrix4x4_t(const float source[4][4]);
			matrix4x4_t(const double source[4][4]);
			~matrix4x4_t() {};
			/*! attention, this function can cause the matrix to become invalid!
				unless you are sure the matrix is invertible, check invalid() afterwards! */
			matrix4x4_t & inverse();
			matrix4x4_t & transpose();
			void identity();
			void translate(float dx,float dy,float dz);
			void rotateX(float degrees);
			void rotateY(float degrees);
			void rotateZ(float degrees);
			void scale(float sx, float sy, float sz);
			int invalid() const { return _invalid; }
			// ignored by swig
			//const PFLOAT * operator [] (int i) const { return matrix[i]; }
			//PFLOAT * operator [] (int i) { return matrix[i]; }
			void setVal(int row, int col, float val)
			{
				matrix[row][col] = val;
			};

			float getVal(int row, int col)
			{
				return matrix[row][col];
			};

		protected:

			float  matrix[4][4];
			int _invalid;
		};

		// Interfaces

	class yafrayInterface_t
	{
		public:
			yafrayInterface_t();
			virtual ~yafrayInterface_t();
			// directly related to scene_t:
			virtual void loadPlugins(const char *path); //!< load plugins from path, if NULL load from default path, if available.
			virtual bool startGeometry(); //!< call before creating geometry; only meshes and vmaps can be created in this state
			virtual bool endGeometry(); //!< call after creating geometry;
			/*! start a triangle mesh
				in this state only vertices, UVs and triangles can be created
				\param id returns the ID of the created mesh
			*/
			virtual unsigned int getNextFreeID();
			virtual bool startTriMesh(unsigned int id, int vertices, int triangles, bool hasOrco, bool hasUV=false, int type=0);
			virtual bool startCurveMesh(unsigned int id, int vertices);
			virtual bool startTriMeshPtr(unsigned int *id, int vertices, int triangles, bool hasOrco, bool hasUV=false, int type=0);
			virtual bool endTriMesh(); //!< end current mesh and return to geometry state
			virtual bool endCurveMesh(const material_t *mat, float strandStart, float strandEnd, float strandShape); //!< end current mesh and return to geometry state
			virtual int  addVertex(double x, double y, double z); //!< add vertex to mesh; returns index to be used for addTriangle
			virtual int  addVertex(double x, double y, double z, double ox, double oy, double oz); //!< add vertex with Orco to mesh; returns index to be used for addTriangle
			virtual void addNormal(double nx, double ny, double nz); //!< add vertex normal to mesh; the vertex that will be attached to is the last one inserted by addVertex method
			virtual bool addTriangle(int a, int b, int c, const material_t *mat); //!< add a triangle given vertex indices and material pointer
			virtual bool addTriangle(int a, int b, int c, int uv_a, int uv_b, int uv_c, const material_t *mat); //!< add a triangle given vertex and uv indices and material pointer
			virtual int  addUV(float u, float v); //!< add a UV coordinate pair; returns index to be used for addTriangle
			virtual bool smoothMesh(unsigned int id, double angle); //!< smooth vertex normals of mesh with given ID and angle (in degrees)
			virtual bool addInstance(unsigned int baseObjectId, matrix4x4_t objToWorld);
			// functions to build paramMaps instead of passing them from Blender
			// (decouling implementation details of STL containers, paraMap_t etc. as much as possible)
			virtual void paramsSetPoint(const char* name, double x, double y, double z);
			virtual void paramsSetString(const char* name, const char* s);
			virtual void paramsSetBool(const char* name, bool b);
			virtual void paramsSetInt(const char* name, int i);
			virtual void paramsSetFloat(const char* name, double f);
			virtual void paramsSetColor(const char* name, float r, float g, float b, float a=1.f);
			virtual void paramsSetColor(const char* name, float *rgb, bool with_alpha=false);
			virtual void paramsSetMatrix(const char* name, float m[4][4], bool transpose=false);
			virtual void paramsSetMatrix(const char* name, double m[4][4], bool transpose=false);
			virtual void paramsSetMemMatrix(const char* name, float* matrix, bool transpose=false);
			virtual void paramsSetMemMatrix(const char* name, double* matrix, bool transpose=false);
			virtual void paramsClearAll(); 	//!< clear the paramMap and paramList
			virtual void paramsStartList(); //!< start writing parameters to the extended paramList (used by materials)
			virtual void paramsPushList(); 	//!< push new list item in paramList (e.g. new shader node description)
			virtual void paramsEndList(); 	//!< revert to writing to normal paramMap
			// functions directly related to renderEnvironment_t
			virtual light_t* 		createLight			(const char* name);
			virtual texture_t* 		createTexture		(const char* name);
			virtual material_t* 	createMaterial		(const char* name);
			virtual camera_t* 		createCamera		(const char* name);
			virtual background_t* 	createBackground	(const char* name);
			virtual integrator_t* 	createIntegrator	(const char* name);
			virtual VolumeRegion* 	createVolumeRegion	(const char* name);
			virtual imageHandler_t*	createImageHandler	(const char* name, bool addToTable = true); //!< The addToTable parameter, if true, allows to avoid the interface from taking ownership of the image handler
			virtual unsigned int 	createObject		(const char* name);
			virtual void clearAll(); //!< clear the whole environment + scene, i.e. free (hopefully) all memory.
			virtual void render(colorOutput_t &output, progressBar_t *pb = 0); //!< render the scene...
			virtual bool startScene(int type=0); //!< start a new scene; Must be called before any of the scene_t related callbacks!
			virtual void setInputGamma(float gammaVal, bool enable);
			virtual void abort();
			virtual paraMap_t* getRenderParameters() { return params; }
			virtual bool getRenderedImage(colorOutput_t &output); //!< put the rendered image to output
			virtual std::vector<std::string> listImageHandlers();
			virtual std::vector<std::string> listImageHandlersFullName();
			virtual std::string getImageFormatFromFullName(const std::string &fullname);
			virtual std::string getImageFullNameFromFormat(const std::string &format);

			virtual void setVerbosityLevel(int vlevel);
			virtual void setVerbosityInfo();
			virtual void setVerbosityWarning();
			virtual void setVerbosityError();
			virtual void setVerbosityMute();

			virtual void setDrawParams(bool on = true);
			virtual bool getDrawParams();

			virtual char* getVersion() const; //!< Get version to check aginst the exporters

			/*! Console Printing wrappers to report in color with yafaray's own console coloring */
			void printInfo(const std::string &msg);
			void printWarning(const std::string &msg);
			void printError(const std::string &msg);
			void printLog(const std::string &msg);

		protected:
			paraMap_t *params;
			std::list<paraMap_t> *eparams; //! for materials that need to define a whole shader tree etc.
			paraMap_t *cparams; //! just a pointer to the current paramMap, either params or a eparams element
			renderEnvironment_t *env;
			scene_t *scene;
			imageFilm_t *film;
			float inputGamma;
			bool gcInput;
	};


	class xmlInterface_t: public yafrayInterface_t
	{
		public:
			xmlInterface_t();
			// directly related to scene_t:
			virtual void loadPlugins(const char *path);
			virtual bool startGeometry();
			virtual bool endGeometry();
			virtual unsigned int getNextFreeID();
			virtual bool startTriMesh(unsigned int id, int vertices, int triangles, bool hasOrco, bool hasUV=false, int type=0);
			virtual bool startCurveMesh(unsigned int id, int vertices);
			virtual bool startTriMeshPtr(unsigned int *id, int vertices, int triangles, bool hasOrco, bool hasUV=false, int type=0);
			virtual bool endTriMesh();
			virtual bool endCurveMesh(const material_t *mat, float strandStart, float strandEnd, float strandShape); //!< end current mesh and return to geometry state
			virtual int  addVertex(double x, double y, double z); //!< add vertex to mesh; returns index to be used for addTriangle
			virtual int  addVertex(double x, double y, double z, double ox, double oy, double oz); //!< add vertex with Orco to mesh; returns index to be used for addTriangle
			virtual void addNormal(double nx, double ny, double nz); //!< add vertex normal to mesh; the vertex that will be attached to is the last one inserted by addVertex method
			virtual bool addTriangle(int a, int b, int c, const material_t *mat);
			virtual bool addTriangle(int a, int b, int c, int uv_a, int uv_b, int uv_c, const material_t *mat);
			virtual int  addUV(float u, float v);
			virtual bool smoothMesh(unsigned int id, double angle);

			// functions directly related to renderEnvironment_t
			virtual light_t* 		createLight		(const char* name);
			virtual texture_t* 		createTexture	(const char* name);
			virtual material_t* 	createMaterial	(const char* name);
			virtual camera_t* 		createCamera	(const char* name);
			virtual background_t* 	createBackground(const char* name);
			virtual integrator_t* 	createIntegrator(const char* name);
			virtual unsigned int 	createObject	(const char* name);
			virtual void clearAll(); //!< clear the whole environment + scene, i.e. free (hopefully) all memory.
			virtual void render(colorOutput_t &output); //!< render the scene...
			virtual bool startScene(int type=0); //!< start a new scene; Must be called before any of the scene_t related callbacks!

			virtual void setOutfile(const char *fname);
		protected:
			void writeParamMap(const paramMap_t &pmap, int indent=1);
			void writeParamList(int indent);

			std::map<const material_t *, std::string> materials;
			std::ofstream xmlFile;
			std::string xmlName;
			const material_t *last_mat;
			size_t nmat;
			int n_uvs;
			unsigned int nextObj;
	};

}
