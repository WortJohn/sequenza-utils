#include "common.h"
#include "parsers.h"

static char module_doc[] = "This is a simple, module "
    "wrapping a C pileup parser";
static char pileup_doc[] = "This function parse the quality and "
    "pileup string to return a python dictionary";

static PyObject *main_parser(PyObject *module, PyObject *args,
    PyObject *kwargs) {

    static char *kwlist[] = {
        "pileup",
        "quality",
        "depth",
        "reference",
        "qlimit",
        "noend",
        "nostart",
    };

    int nucleot_list[4] = {0, 0, 0, 0};
    int strand_list[4]  = {0, 0, 0, 0};

    PyObject *ret = PyDict_New(); // new reference
    assert(PyDict_Check(ret));

    PyObject *z_list = PyList_New(4);
    assert(PyList_Check(z_list));

    char* quality;
    char* pileup;


    char reference;
    int depth, qlimit, noend, nostart, i;
    qlimit    = 53;
    nostart   = 0;
    noend     = 0;
    #if PY_MAJOR_VERSION >= 3
        if (! PyArg_ParseTupleAndKeywords(args, kwargs, "zziC|iii", kwlist,
            &pileup, &quality, &depth, &reference, &qlimit,
            &noend, &nostart)) {
            goto except;
        }
    #else
        if (! PyArg_ParseTupleAndKeywords(args, kwargs, "zzic|iii", kwlist,
            &pileup, &quality, &depth, &reference, &qlimit,
            &noend, &nostart)) {
            goto except;
        }
    #endif


    pileup_parse(nucleot_list, strand_list, pileup, quality, depth,
        reference, qlimit, noend, nostart);


    for(i = 0; i < 4; ++i){
        PyList_SetItem(z_list, i, Py_BuildValue("i", strand_list[i]));
    }

    PyDict_SetItemString(ret, "A",
        Py_BuildValue("i", nucleot_list[0]));
    PyDict_SetItemString(ret, "C",
        Py_BuildValue("i", nucleot_list[1]));
    PyDict_SetItemString(ret, "G",
        Py_BuildValue("i", nucleot_list[2]));
    PyDict_SetItemString(ret, "T",
        Py_BuildValue("i", nucleot_list[3]));
    PyDict_SetItemString(ret, "Z", z_list);

    Py_DECREF(z_list);


    assert(! PyErr_Occurred());
    goto finally;
except:
    Py_XDECREF(ret);
    ret = NULL;
finally:
    return ret;
}

static PyMethodDef module_methods[] = {
    {
        "acgt",
        (PyCFunction)main_parser,
        METH_VARARGS | METH_KEYWORDS,
        pileup_doc
    },
    {
        NULL,
        NULL,
        0,
        NULL
    }
};


MOD_INIT(c_pileup)
{
    PyObject *m;

    MOD_DEF(m, "c_pileup", module_doc,
            module_methods)

    if (m == NULL)
        return MOD_ERROR_VAL;

    return MOD_SUCCESS_VAL(m);
}
