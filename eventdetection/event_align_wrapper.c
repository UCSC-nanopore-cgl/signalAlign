//
// Created by Andrew Bailey on 6/27/18.
//

#include <Python.h>
#include <hdf5.h>

#include <eventAligner.h>



static PyObject *eventAlignError;

static PyObject *wrap_load_from_raw(PyObject *self, PyObject *args, PyObject *keywds) {
    /* Parse the input tuple */
    char* fast5_path;
    char* template_model_file;
    char* nuc_sequence;
    char* path_in_fast5;
    static char *kwlist[] = {"fast5_path", "template_model_file", "nuc_sequence", "path_in_fast5", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "ssss", kwlist,
                                     &fast5_path, &template_model_file, &nuc_sequence, &path_in_fast5))
        return NULL;

    int status = load_from_raw(fast5_path, template_model_file, nuc_sequence, path_in_fast5);

    return Py_BuildValue("i", status);
}



fast5_raw_scaling get_channel_params(const char* hdf5_path){

    hid_t hdf_file = fast5_open(hdf5_path);
    fast5_raw_scaling scaling_params = fast5_get_channel_params(hdf_file);
    fast5_close(hdf_file);
    return scaling_params;
}


static PyObject *wrap_get_channel_params(PyObject *self, PyObject *args) {
    /* Parse the input tuple */
    const char *command;
//    PyArrayObject *scalings;
    fast5_raw_scaling scalings;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;

    scalings = get_channel_params(command);

    return Py_BuildValue("{s:f, s:f, s:f, s:f}",
                         "offset", scalings.offset,
                         "sample_rate", scalings.sample_rate,
                         "range", scalings.range,
                         "digitisation", scalings.digitisation);
}


char* get_raw_read_name(const char* hdf5_path){

    hid_t hdf_file = fast5_open(hdf5_path);
    char* name = fast5_get_raw_read_name(hdf_file);
    fast5_close(hdf_file);
    return name;
}


static PyObject *wrap_get_raw_read_name(PyObject *self, PyObject *args)
{
    /* Parse the input tuple */
    const char* command;
    char* name;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;

    name = get_raw_read_name(command);

    return Py_BuildValue("s", name);
}

int* str_to_int(const char* name) {
    return (int*) name;
}

static PyObject *wrap_str_to_int(PyObject *self, PyObject *args)
{
    /* Parse the input tuple */
    const char* command;
    int* char_to_int;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;

    if (strcmp(command ,"error") == 0){
        PyErr_SetString(eventAlignError, "No number for you! Try again!");
        return NULL;
    } else{
    char_to_int = str_to_int(command);
    }
// convert the c int to python object
    return Py_BuildValue("i", char_to_int);
//    return PyBytes_FromString(name);
}


static PyObject *version(PyObject *self){
    return Py_BuildValue("s", "KmerAlign Version 0.1");
}


/* Module specification */
static PyMethodDef eventalign_methods[] = {
        {"str_to_int", wrap_str_to_int, METH_VARARGS, "Cast a string to an integer"},
        {"get_id", wrap_get_raw_read_name, METH_VARARGS, "Get the Read ID from a fast5 file."},
        {"version", (PyCFunction)version, METH_NOARGS, "Returns the version kmerAligner"},
        {"get_channel_params", wrap_get_channel_params, METH_VARARGS, "Returns the channel parameters of a read"},
        {"load_from_raw", wrap_load_from_raw, METH_VARARGS|METH_KEYWORDS, "Performs banded alignment between kmers of a sequence "
                                                                         "and the raw current readings given a "
                                                                         "model file"},
        {NULL, NULL, 0, NULL} /* Sentinel */
};


static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "kmeralign",
        "Create kmer to events alignment",
        -1,
//        sizeof(struct module_state),
        eventalign_methods,
};


PyMODINIT_FUNC PyInit_kmeralign(void) {
    PyObject *m;

    m = PyModule_Create(&moduledef);
    if (m == NULL)
        return NULL;

    eventAlignError = PyErr_NewException("eventAlign.Error", NULL, NULL);
    Py_INCREF(eventAlignError);
    PyModule_AddObject(m, "error", eventAlignError);
    return m;
}

