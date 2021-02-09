/* File: read_dcdmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Tue Feb  9 03:03:26 2021
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <string.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *read_dcd_error;
static PyObject *read_dcd_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
typedef char * string;

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
\
#define FAILNULL(p) do {                                            \
    if ((p) == NULL) {                                              \
        PyErr_SetString(PyExc_MemoryError, "NULL pointer found");   \
        goto capi_fail;                                             \
    }                                                               \
} while (0)

#define STRINGMALLOC(str,len)\
    if ((str = (string)malloc(sizeof(char)*(len+1))) == NULL) {\
        PyErr_SetString(PyExc_MemoryError, "out of memory");\
        goto capi_fail;\
    } else {\
        (str)[len] = '\0';\
    }

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#define STRINGFREE(str) do {if (!(str == NULL)) free(str);} while (0)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#define STRINGCOPYN(to,from,buf_size)                           \
    do {                                                        \
        int _m = (buf_size);                                    \
        char *_to = (to);                                       \
        char *_from = (from);                                   \
        FAILNULL(_to); FAILNULL(_from);                         \
        (void)strncpy(_to, _from, sizeof(char)*_m);             \
        _to[_m-1] = '\0';                                      \
        /* Padding with spaces instead of nulls */              \
        for (_m -= 2; _m >= 0 && _to[_m] == '\0'; _m--) {      \
            _to[_m] = ' ';                                      \
        }                                                       \
    } while (0)


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int
string_from_pyobj(string *str,int *len,const string inistr,PyObject *obj,const char *errmess)
{
    PyArrayObject *arr = NULL;
    PyObject *tmp = NULL;
#ifdef DEBUGCFUNCS
fprintf(stderr,"string_from_pyobj(str='%s',len=%d,inistr='%s',obj=%p)\n",(char*)str,*len,(char *)inistr,obj);
#endif
    if (obj == Py_None) {
        if (*len == -1)
            *len = strlen(inistr); /* Will this cause problems? */
        STRINGMALLOC(*str,*len);
        STRINGCOPYN(*str,inistr,*len+1);
        return 1;
    }
    if (PyArray_Check(obj)) {
        if ((arr = (PyArrayObject *)obj) == NULL)
            goto capi_fail;
        if (!ISCONTIGUOUS(arr)) {
            PyErr_SetString(PyExc_ValueError,"array object is non-contiguous.");
            goto capi_fail;
        }
        if (*len == -1)
            *len = (PyArray_ITEMSIZE(arr))*PyArray_SIZE(arr);
        STRINGMALLOC(*str,*len);
        STRINGCOPYN(*str,PyArray_DATA(arr),*len+1);
        return 1;
    }
    if (PyBytes_Check(obj)) {
        tmp = obj;
        Py_INCREF(tmp);
    }
    else if (PyUnicode_Check(obj)) {
        tmp = PyUnicode_AsASCIIString(obj);
    }
    else {
        PyObject *tmp2;
        tmp2 = PyObject_Str(obj);
        if (tmp2) {
            tmp = PyUnicode_AsASCIIString(tmp2);
            Py_DECREF(tmp2);
        }
        else {
            tmp = NULL;
        }
    }
    if (tmp == NULL) goto capi_fail;
    if (*len == -1)
        *len = PyBytes_GET_SIZE(tmp);
    STRINGMALLOC(*str,*len);
    STRINGCOPYN(*str,PyBytes_AS_STRING(tmp),*len+1);
    Py_DECREF(tmp);
    return 1;
capi_fail:
    Py_XDECREF(tmp);
    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = read_dcd_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}

static int
int_from_pyobj(int* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;

    if (PyLong_Check(obj)) {
        *v = Npy__PyLong_AsInt(obj);
        return !(*v == -1 && PyErr_Occurred());
    }

    tmp = PyNumber_Long(obj);
    if (tmp) {
        *v = Npy__PyLong_AsInt(tmp);
        Py_DECREF(tmp);
        return !(*v == -1 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj, 0);
    if (tmp) {
        PyErr_Clear();
        if (int_from_pyobj(v, tmp, errmess)) {
            Py_DECREF(tmp);
            return 1;
        }
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = read_dcd_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC_US(read_dcd,READ_DCD)(string,int*,int*,double*,float*,float*,float*,size_t);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************** read_dcd **********************************/
static char doc_f2py_rout_read_dcd_read_dcd[] = "\
box,x,y,z = read_dcd(fname,nstr,natom)\n\nWrapper for ``read_dcd``.\
\n\nParameters\n----------\n"
"fname : input string(len=256)\n"
"nstr : input int\n"
"natom : input int\n"
"\nReturns\n-------\n"
"box : rank-2 array('d') with bounds (nstr,3)\n"
"x : rank-2 array('f') with bounds (nstr,natom)\n"
"y : rank-2 array('f') with bounds (nstr,natom)\n"
"z : rank-2 array('f') with bounds (nstr,natom)";
/* extern void F_FUNC_US(read_dcd,READ_DCD)(string,int*,int*,double*,float*,float*,float*,size_t); */
static PyObject *f2py_rout_read_dcd_read_dcd(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(string,int*,int*,double*,float*,float*,float*,size_t)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

  string fname = NULL;
  int slen(fname);
  PyObject *fname_capi = Py_None;
  int nstr = 0;
  PyObject *nstr_capi = Py_None;
  int natom = 0;
  PyObject *natom_capi = Py_None;
  double *box = NULL;
  npy_intp box_Dims[2] = {-1, -1};
  const int box_Rank = 2;
  PyArrayObject *capi_box_tmp = NULL;
  int capi_box_intent = 0;
  float *x = NULL;
  npy_intp x_Dims[2] = {-1, -1};
  const int x_Rank = 2;
  PyArrayObject *capi_x_tmp = NULL;
  int capi_x_intent = 0;
  float *y = NULL;
  npy_intp y_Dims[2] = {-1, -1};
  const int y_Rank = 2;
  PyArrayObject *capi_y_tmp = NULL;
  int capi_y_intent = 0;
  float *z = NULL;
  npy_intp z_Dims[2] = {-1, -1};
  const int z_Rank = 2;
  PyArrayObject *capi_z_tmp = NULL;
  int capi_z_intent = 0;
    static char *capi_kwlist[] = {"fname","nstr","natom",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOO|:read_dcd.read_dcd",\
        capi_kwlist,&fname_capi,&nstr_capi,&natom_capi))
        return NULL;
/*frompyobj*/
  /* Processing variable fname */
  slen(fname) = 256;
  f2py_success = string_from_pyobj(&fname,&slen(fname),"",fname_capi,"string_from_pyobj failed in converting 1st argument `fname' of read_dcd.read_dcd to C string");
  if (f2py_success) {
  /* Processing variable nstr */
    f2py_success = int_from_pyobj(&nstr,nstr_capi,"read_dcd.read_dcd() 2nd argument (nstr) can't be converted to int");
  if (f2py_success) {
  /* Processing variable natom */
    f2py_success = int_from_pyobj(&natom,natom_capi,"read_dcd.read_dcd() 3rd argument (natom) can't be converted to int");
  if (f2py_success) {
  /* Processing variable box */
  box_Dims[0]=nstr,box_Dims[1]=3;
  capi_box_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_box_tmp = array_from_pyobj(NPY_DOUBLE,box_Dims,box_Rank,capi_box_intent,Py_None);
  if (capi_box_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : read_dcd_error,"failed in converting hidden `box' of read_dcd.read_dcd to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    box = (double *)(PyArray_DATA(capi_box_tmp));

  /* Processing variable x */
  x_Dims[0]=nstr,x_Dims[1]=natom;
  capi_x_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_x_tmp = array_from_pyobj(NPY_FLOAT,x_Dims,x_Rank,capi_x_intent,Py_None);
  if (capi_x_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : read_dcd_error,"failed in converting hidden `x' of read_dcd.read_dcd to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    x = (float *)(PyArray_DATA(capi_x_tmp));

  /* Processing variable y */
  y_Dims[0]=nstr,y_Dims[1]=natom;
  capi_y_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_y_tmp = array_from_pyobj(NPY_FLOAT,y_Dims,y_Rank,capi_y_intent,Py_None);
  if (capi_y_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : read_dcd_error,"failed in converting hidden `y' of read_dcd.read_dcd to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    y = (float *)(PyArray_DATA(capi_y_tmp));

  /* Processing variable z */
  z_Dims[0]=nstr,z_Dims[1]=natom;
  capi_z_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_z_tmp = array_from_pyobj(NPY_FLOAT,z_Dims,z_Rank,capi_z_intent,Py_None);
  if (capi_z_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : read_dcd_error,"failed in converting hidden `z' of read_dcd.read_dcd to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    z = (float *)(PyArray_DATA(capi_z_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(fname,&nstr,&natom,box,x,y,z,slen(fname));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("NNNN",capi_box_tmp,capi_x_tmp,capi_y_tmp,capi_z_tmp);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  }  /*if (capi_z_tmp == NULL) ... else of z*/
  /* End of cleaning variable z */
  }  /*if (capi_y_tmp == NULL) ... else of y*/
  /* End of cleaning variable y */
  }  /*if (capi_x_tmp == NULL) ... else of x*/
  /* End of cleaning variable x */
  }  /*if (capi_box_tmp == NULL) ... else of box*/
  /* End of cleaning variable box */
  } /*if (f2py_success) of natom*/
  /* End of cleaning variable natom */
  } /*if (f2py_success) of nstr*/
  /* End of cleaning variable nstr */
    STRINGFREE(fname);
  }  /*if (f2py_success) of fname*/
  /* End of cleaning variable fname */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/****************************** end of read_dcd ******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"read_dcd",-1,{{-1}},0,(char *)F_FUNC_US(read_dcd,READ_DCD),(f2py_init_func)f2py_rout_read_dcd_read_dcd,doc_f2py_rout_read_dcd_read_dcd},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "read_dcd",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_read_dcd(void) {
  int i;
  PyObject *m,*d, *s, *tmp;
  m = read_dcd_module = PyModule_Create(&moduledef);
  Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module read_dcd (failed to import numpy)"); return m;}
  d = PyModule_GetDict(m);
  s = PyUnicode_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);
  s = PyUnicode_FromString(
    "This module 'read_dcd' is auto-generated with f2py (version:2).\nFunctions:\n"
"  box,x,y,z = read_dcd(fname,nstr,natom)\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  Py_DECREF(s);
  s = PyUnicode_FromString("1.20.0");
  PyDict_SetItemString(d, "__f2py_numpy_version__", s);
  Py_DECREF(s);
  read_dcd_error = PyErr_NewException ("read_dcd.error", NULL, NULL);
  /*
   * Store the error object inside the dict, so that it could get deallocated.
   * (in practice, this is a module, so it likely will not and cannot.)
   */
  PyDict_SetItemString(d, "_read_dcd_error", read_dcd_error);
  Py_DECREF(read_dcd_error);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
    tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
    PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
    Py_DECREF(tmp);
  }

/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"read_dcd");
#endif
  return m;
}
#ifdef __cplusplus
}
#endif