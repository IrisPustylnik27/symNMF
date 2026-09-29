#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "symnmf.h"
#include "matrix.h"

static double **python_list_to_matrix_array(PyObject *pyList, int n, int d);
static PyObject *matrix_array_to_python_list(double **matrix, int n, int d);
static double **parsing_PyObject_to_matrix_array(PyObject *py_X, int *n, int *d);

static PyObject *symnmf_wrapper(PyObject *self, PyObject *args)
{
    (void)self;

    PyObject *py_W = NULL, *py_H_input = NULL, *py_H_output = NULL;
    double **W = NULL, **H = NULL, **result = NULL;
    int n = 0, k = 0, W_rows = 0, W_col = 0;

    if (!PyArg_ParseTuple(args, "OO", &py_H_input, &py_W))
        return NULL;

    H = parsing_PyObject_to_matrix_array(py_H_input, &n, &k);
    if (H == NULL)
        return NULL;

    if (!(W = parsing_PyObject_to_matrix_array(py_W, &W_rows, &W_col)))
        goto error;

    if (W_rows != n || W_col != n)
    {
        PyErr_SetString(PyExc_ValueError, "An Error Has Occurred");
        goto error;
    }

    if (!(result = optimize_H_matrix(H, W, n, k)))
        goto error;

    py_H_output = matrix_array_to_python_list(result, n, k);
    goto cleanup;

error:
    if (!PyErr_Occurred())
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
    py_H_output = NULL;

cleanup:
    free_matrix(H, n);
    free_matrix(W, W_rows);
    free_matrix(result, n);

    return py_H_output;
}

static PyObject *sym_wrapper(PyObject *self, PyObject *args)
{
    (void)self; // Unused parameter

    PyObject *py_X = NULL, *py_A = NULL;
    double **X = NULL, **result = NULL;
    int n = 0, d = 0;

    if (!PyArg_ParseTuple(args, "O", &py_X))
        return NULL;

    if (!(X = parsing_PyObject_to_matrix_array(py_X, &n, &d)))
        return NULL;

    if (!(result = similarity_matrix(X, n, d)))
        goto error;

    py_A = matrix_array_to_python_list(result, n, n);
    goto cleanup;

error:
    if (!PyErr_Occurred())
    {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
    }
    py_A = NULL;

cleanup:
    free_matrix(X, n);
    free_matrix(result, n);

    return py_A;
}

static PyObject *ddg_wrapper(PyObject *self, PyObject *args)
{
    (void)self; // Unused parameter

    PyObject *py_X = NULL, *py_D = NULL;
    double **A = NULL, **X = NULL, **result = NULL;
    int n = 0, d = 0;

    if (!PyArg_ParseTuple(args, "O", &py_X))
        return NULL;

    if (!(X = parsing_PyObject_to_matrix_array(py_X, &n, &d)))
        return NULL;

    if (!(A = similarity_matrix(X, n, d)))
        goto error;

    if (!(result = diagonal_matrix(A, n)))
        goto error;

    py_D = matrix_array_to_python_list(result, n, n);
    goto cleanup;

error:
    if (!PyErr_Occurred())
    {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
    }
    py_D = NULL;

cleanup:
    free_matrix(X, n);
    free_matrix(A, n);
    free_matrix(result, n);

    return py_D;
}

static PyObject *norm_wrapper(PyObject *self, PyObject *args)
{
    (void)self;

    PyObject *py_X = NULL, *py_W = NULL;
    double **X = NULL, **result = NULL;
    int n = 0, d = 0;

    if (!PyArg_ParseTuple(args, "O", &py_X))
        return NULL;

    if (!(X = parsing_PyObject_to_matrix_array(py_X, &n, &d)))
        return NULL;

    if (!(result = make_W_matrix(X, n, d)))
        goto error;

    py_W = matrix_array_to_python_list(result, n, n);
    goto cleanup;

error:
    if (!PyErr_Occurred())
    {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
    }
    py_W = NULL;

cleanup:
    free_matrix(X, n);
    free_matrix(result, n);

    return py_W;
}

static double **make_W_matrix(double **X, int n, int d)
{
    double **A = NULL, **D = NULL, **result = NULL;

    if (!(A = similarity_matrix(X, n, d)))
        goto error;

    if (!(D = diagonal_matrix(A, n)))
        goto error;

    if (!(result = normalized_similarity_matrix(D, A, n)))
        goto error;
    goto no_error;

error:
    free_matrix(result, n);
    result = NULL;

no_error:
    free_matrix(A, n);
    free_matrix(D, n);
    return result;
}

static PyMethodDef symnmfMethods[] = {
    {"symnmf", symnmf_wrapper, METH_VARARGS, "SymNMF wrapper"},
    {"sym", sym_wrapper, METH_VARARGS, "Sym wrapper"},
    {"ddg", ddg_wrapper, METH_VARARGS, "DDG wrapper"},
    {"norm", norm_wrapper, METH_VARARGS, "Norm wrapper"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef symnmfModule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",
    "Python wrapper for SymNMF",
    -1,
    symnmfMethods};

PyMODINIT_FUNC PyInit_symnmfmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&symnmfModule);
    if (!m)
    {
        return NULL;
    }
    return m;
}

static double **parsing_PyObject_to_matrix_array(PyObject *py_X, int *n, int *d)
{
    double **X;
    PyObject *first_row;

    if (!PyList_Check(py_X))
    {
        PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
        return NULL;
    }

    *n = (int)PyList_Size(py_X);
    if (*n <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "An Error Has Occurred");
        return NULL;
    }

    first_row = PyList_GetItem(py_X, 0);

    if (first_row == NULL || !PyList_Check(first_row))
    {
        PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
        return NULL;
    }
    *d = (int)PyList_Size(first_row);
    if (*d <= 0)
    {
        PyErr_SetString(PyExc_ValueError, "An Error Has Occurred");
        return NULL;
    }

    X = python_list_to_matrix_array(py_X, *n, *d);

    return X;
}

static double **python_list_to_matrix_array(PyObject *py_list, int n, int d)
{
    double **list = NULL;
    PyObject *vector, *item;
    int i, j;

    if (!(list = allocate_matrix(n, d, 0)))
    {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL;
    }
    for (i = 0; i < n; i++)
    {
        vector = PyList_GetItem(py_list, i);
        if (vector == NULL || !PyList_Check(vector))
        {
            PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
            goto error;
        }
        if ((int)PyList_Size(vector) != d)
        {
            PyErr_SetString(PyExc_ValueError, "An Error Has Occurred");
            goto error;
        }
        for (j = 0; j < d; j++)
        {
            if (!(item = PyList_GetItem(vector, j)))
            {
                PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
                goto error;
            }
            list[i][j] = PyFloat_AsDouble(item);
            if (PyErr_Occurred())
                goto error;
        }
    }
    goto no_error;

error:
    free_matrix(list, n);
    list = NULL;

no_error:
    return list;
}

static PyObject *matrix_array_to_python_list(double **matrix, int n, int d)
{
    int i, j;
    PyObject *output_list, *item, *vector_list;

    if (!(output_list = PyList_New(n)))
        return NULL;

    for (i = 0; i < n; i++)
    {
        if (!(vector_list = PyList_New(d)))
        {
            Py_DECREF(output_list);
            return NULL;
        }
        for (j = 0; j < d; j++)
        {
            if (!(item = PyFloat_FromDouble(matrix[i][j])))
            {
                PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
                Py_DECREF(output_list);
                Py_DECREF(vector_list);
                return NULL;
            }
            PyList_SET_ITEM(vector_list, j, item);
        }
        PyList_SET_ITEM(output_list, i, vector_list);
    }
    return output_list;
}
