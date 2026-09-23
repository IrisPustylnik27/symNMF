
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include "symnmf.h"

double **python_list_to_matrix_array(PyObject *pyList, int N, int d);
PyObject *matrix_array_to_python_list(double **matrix, int N, int d);
void free_matrix(double **matrix, int n);

static PyObject *symnmf_wrapper (PyObject *self, PyObject *args){

    PyObject *py_W, *py_H;
    double **W;
    double **H;
    double **result;
    int n, k;

    if (!PyArg_ParseTuple(args, "OO", &py_H, &py_W))
    {
        return NULL;
    }

    n = (int)PyList_Size(py_W);
    k = (int)PyList_Size(PyList_GetItem(py_H, 0));

    W = python_list_to_matrix_array(py_W, n, n);
    if (W == NULL)
    {
        return NULL;
    }

    H = python_list_to_matrix_array(py_H, n, k);
    if (H == NULL)
    {
        free_matrix(W, n);
        return NULL;
    }

    result = optimize_H_matrix(H, W, n, k);

    py_H = matrix_array_to_python_list(H, n, k);

    free_matrix(W, n);
    free_matrix(H, n);
    free_matrix(result, n);

    return py_H;

}

static PyObject *sym_wrapper (PyObject *self, PyObject *args){
    
    PyObject *py_X, *py_A;
    double **X;
    double **result;
    int n, d;

    if (!PyArg_ParseTuple(args, "O", &py_X))
    {
        return NULL;
    }

    n = (int)PyList_Size(py_X);
    d = (int)PyList_Size(PyList_GetItem(py_X, 0));

    X = python_list_to_matrix_array(py_X, n, d);
    if (X == NULL)
    {
        return NULL;
    }
    result = similarity_matrix(X, n, d);
    if (result == NULL)
    {
        free_matrix(X, n);
        return NULL;
    }

    py_A = matrix_array_to_python_list(result, n, n);

    free_matrix(X, n);
    free_matrix(result, n);

    return py_A;
    
}

static PyObject *ddg_wrapper (PyObject *self, PyObject *args){
    
    PyObject *py_X, *py_D;
    double **A, **X;
    double **result;
    int n, d;

    if (!PyArg_ParseTuple(args, "O", &py_X))
    {
        return NULL;
    }

    n = (int)PyList_Size(py_X);
    d = (int)PyList_Size(PyList_GetItem(py_X, 0));

    X = python_list_to_matrix_array(py_X, n, d);
    if (X == NULL)
    {
        return NULL;
    }

    A = similarity_matrix(X, n, d);
    if (A == NULL)
    {
        free_matrix(X, n);
        return NULL;
    }

    result = diagonal_matrix(A, n, n);
    if (result == NULL)
    {
        free_matrix(X, n);
        free_matrix(A, n);
        return NULL;
    }

    py_D = matrix_array_to_python_list(result, n, n);

    free_matrix(X, n);
    free_matrix(A, n);
    free_matrix(result, n);


    return py_D;
    
}

static PyObject *norm_wrapper (PyObject *self, PyObject *args){
    
    PyObject *py_X, *py_W;
    double **A, **D, **X;
    double **result;
    int n, d;

    if (!PyArg_ParseTuple(args, "O", &py_X))
    {
        return NULL;
    }

    n = (int)PyList_Size(py_X);
    d = (int)PyList_Size(PyList_GetItem(py_X, 0));

    X = python_list_to_matrix_array(py_X, n, d);
    if (X == NULL)
    {
        return NULL;
    }

    A = similarity_matrix(X, n, d);
    if (A == NULL)
    {
        free_matrix(X, n);
        return NULL;
    }

    D = diagonal_matrix(A, n, n);
    if (D == NULL)
    {
        free_matrix(X, n);
        free_matrix(A, n);
        return NULL;
    }

    result = norm_matrix(D, A, n, n);
    if (result == NULL)
    {
        free_matrix(X, n);
        free_matrix(A, n);
        free_matrix(D, n);
        return NULL;
    }
    py_W = matrix_array_to_python_list(result, n, n);

    free_matrix(X, n);
    free_matrix(A, n);
    free_matrix(D, n);
    free_matrix(result, n);

    return py_W;
}


static PyMethodDef symnmfMethods[] = {
    {"symnmf", symnmf_wrapper, METH_VARARGS, "SymNMF wrapper"},
    {"sym", sym_wrapper, METH_VARARGS, "Sym wrapper"},
    {"ddg", ddg_wrapper, METH_VARARGS, "DDG wrapper"},
    {"norm", norm_wrapper, METH_VARARGS, "Norm wrapper"},
    {NULL, NULL, 0, NULL}
};


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

double **python_list_to_matrix_array(PyObject *pyList, int N, int d)
{
    double **pList = NULL;
    int i, j;

    pList = malloc(N * sizeof(double *));
    if (pList == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; // Return NULL if memory allocation fails
    }

    for (i = 0; i < N; i++)
    {
        PyObject *vector = PyList_GetItem(pyList, i);
        if (vector == NULL || !PyList_Check(vector))
            {
                PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
                free_matrix(pList, i);
                return NULL; 
            }

            if ((int)PyList_Size(vector) != d)
            {
                PyErr_SetString(PyExc_ValueError, "An Error Has Occurred");
                free_matrix(pList, i);
                return NULL; 
            }
        
        pList[i] = malloc(d * sizeof(double));
        if (pList[i] == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
            free_matrix(pList, i);

            return NULL; // Return NULL if memory allocation fails
        }
        for (j = 0; j < d; j++)
        {
            PyObject *item = PyList_GetItem(vector, j);

            if (item == NULL)
            {
                PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
                free_matrix(pList, i+1);
                return NULL; 
            }

            pList[i][j] = PyFloat_AsDouble(item);

            if (PyErr_Occurred())
            {
                PyErr_SetString(PyExc_TypeError, "An Error Has Occurred");
                free_matrix(pList, i+1);
                return NULL; // Return NULL if conversion fails
            }
        }
    }
    return pList;
}

PyObject *matrix_array_to_python_list(double **matrix, int N, int d)
{
    int i, j;

    PyObject *outputList = PyList_New(N);
    if (outputList == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
        return NULL; // Return NULL if memory allocation fails
    }

    PyObject *item;
    PyObject *vectorList;

    for (i = 0; i < N; i++)
    {
        vectorList = PyList_New(d);
        if (vectorList == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
            Py_DECREF(outputList);
            return NULL;
        }
        for (j = 0; j < d; j++)
        {
            item = PyFloat_FromDouble(matrix[i][j]);
            if (item == NULL)
            {
                PyErr_SetString(PyExc_MemoryError, "An Error Has Occurred");
                Py_DECREF(outputList);
                Py_DECREF(vectorList);
                return NULL;
            }
            PyList_SET_ITEM(vectorList, j, item);
        }
        PyList_SET_ITEM(outputList, i, vectorList);
    }
    return outputList;
}

void free_matrix(double **matrix, int n){
    int i;
    if (matrix == NULL){
        return;
    }
    for (i = 0; i < n; i++)
        {
            free(matrix[i]);
        }
        free(matrix);
}

