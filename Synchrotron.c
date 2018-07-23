#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "Synchrotron.h"
#include <stdio.h>
/*
 * Follow the paper A&A 367, 809-825 (2001)
 * The multifrequency emission of Mrk 501 from radio to TeV gamma-rays
 * The emission and absorption coefficient calculation using the approximation (Appendix A)
 * The equation (3) has some difficult in high energy, donot use it.
 */

double U_B(double B){
    return B * B / (8 * M_PI);
}

double Nu_B(double B){
    return e * B / (2 * M_PI * m_e * c);
}

double t(double nu, double gamma, double B){
    double nu_B = Nu_B(B);
    return nu / (3 * pow(gamma, 2) * nu_B);
}


double c1 = 0.78;
double c2 = 0.25;
double c3 = 2.175;

double g(double nu, double t, double n, double K, double nu_B){
    double a = c2 + (n-1)/2;
    double x = c3*t;
    double fx = gsl_sf_gamma_inc(a, x);

    return K * pow((nu / (3 * nu_B)), -n/2) * pow(t, c2 + (n-1)/2) * pow(c3*t, (1-n)/2 - c2) * fx;
}


//Power-Law
double j_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tmax = t(nu, gamma_max, B);
    
    return 9 * sigma_T * c * u_B * c1 / (24 * pow(M_PI, 2) * nu_B) * sqrt(nu/nu_B)\
             * (g(nu, tmax, n, K, nu_B) - g(nu, tmin, n, K, nu_B));
}

double J_PowerLaw(double nu, double B, double* pars){
    return j_PowerLaw(nu, B, pars[0], pars[1], pars[2], pars[3]);
}

/*
void J_PowerLaw(double* nu, double B, double K, double n, double gamma_min, double gamma_max, double* j, int size){
    for(int i=0; i<size; i++){
        j[i] = j_PowerLaw(nu[i], B, K, n, gamma_min, gamma_max);
    }
}
*/



//Break Power-Law
double j_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tb = t(nu, gamma_b, B);
    double tmax = t(nu, gamma_max, B);
    
    return 9 * sigma_T * c * u_B * c1 / (24 * pow(M_PI, 2) * nu_B) * sqrt(nu/nu_B)\
             * (g(nu, tb, n1, K1, nu_B) - g(nu, tmin, n1, K1, nu_B) + g(nu, tmax, n2, K2, nu_B) - g(nu, tb, n2, K2, nu_B));
}

double J_BrokenPowerLaw(double nu, double B, double* pars){
    return j_BrokenPowerLaw(nu, B, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6]);
}

/*
void J_BrokenPowerLaw(double* nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max, double* j, int size){
    for(int i=0; i<size; i++){
        j[i] = j_BrokenPowerLaw(nu[i], B, K1, K2, n1, n2, gamma_min, gamma_b, gamma_max);
    }
}
*/

//SpectrumType: 1-PowerLaw; 2-BrokenPowerLaw
static PyObject* J(PyObject* self, PyObject* args){
    PyArrayObject* nu;
    double B;
    int SpectrumType;
    PyTupleObject* SpectrumPars;
    PyArg_ParseTuple(args, "O!diO!", &PyArray_Type, &nu, &B, &SpectrumType, &PyTuple_Type, &SpectrumPars);
    
    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ++ipar)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    PyArrayObject* in_array = nu;
    PyObject* out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        if(SpectrumType == 1)
            *out_dataptr = J_PowerLaw(*in_dataptr, B, pars);
        if(SpectrumType == 2)
            *out_dataptr = J_BrokenPowerLaw(*in_dataptr, B, pars);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}


double h(double nu, double t, double n, double K, double nu_B){
    double a = c2 + n/2;
    double x = c3*t;
    double fx = gsl_sf_gamma_inc(a, x);

    return K * (n+2) * pow(3, n/2) * pow(nu/nu_B, -n/2) * pow(t, c2+n/2) * pow(c3*t, -c2-n/2) * fx;
}


//Power-Law
double k_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tmax = t(nu, gamma_max, B);
    
    return 3 * sqrt(3) * sigma_T * c * u_B * c1 / (16 * pow(M_PI, 2) * m_e * pow(nu, 2) * nu_B)\
             * (h(nu, tmax, n, K, nu_B) - h(nu, tmin, n, K, nu_B));
}

double K_PowerLaw(double nu, double B, double* pars){
    return k_PowerLaw(nu, B, pars[0], pars[1], pars[2], pars[3]);
}

/*
void K_PowerLaw(double* nu, double B, double K, double n, double gamma_min, double gamma_max, double* k, int size){
    for(int i=0; i<size; i++){
        k[i] = k_PowerLaw(nu[i], B, K, n, gamma_min, gamma_max);
    }
}
*/

//Broken Power-Law
double k_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max){
    double u_B = U_B(B);
    double nu_B = Nu_B(B);
    double tmin = t(nu, gamma_min, B);
    double tb = t(nu, gamma_b, B);
    double tmax = t(nu, gamma_max, B);
    
    return 3 * sqrt(3) * sigma_T * c * u_B * c1 / (16 * pow(M_PI, 2) * m_e * pow(nu, 2) * nu_B)\
             * (h(nu, tb, n1, K1, nu_B) - h(nu, tmin, n1, K1, nu_B) + h(nu, tmax, n2, K2, nu_B) - h(nu, tb, n2, K2, nu_B));
}

double K_BrokenPowerLaw(double nu, double B, double* pars){
    return k_BrokenPowerLaw(nu, B, pars[0], pars[1], pars[2], pars[3], pars[4], pars[5], pars[6]);
}

/*
void K_BrokenPowerLaw(double* nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max, double* k, int size){
    for(int i=0; i<size; i++){
        k[i] = k_BrokenPowerLaw(nu[i], B, K1, K2, n1, n2, gamma_min, gamma_b, gamma_max);
    }
}
*/

//SpectrumType: 1-PowerLaw; 2-BrokenPowerLaw
static PyObject* K(PyObject* self, PyObject* args){
    PyArrayObject* nu;
    double B;
    int SpectrumType;
    PyTupleObject* SpectrumPars;
    PyArg_ParseTuple(args, "O!diO!", &PyArray_Type, &nu, &B, &SpectrumType, &PyTuple_Type, &SpectrumPars);

    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ipar++)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    PyArrayObject* in_array = nu;
    PyObject* out_array;
    out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        if(SpectrumType == 1)
            *out_dataptr = K_PowerLaw(*in_dataptr, B, pars);
        if(SpectrumType == 2)
            *out_dataptr = K_BrokenPowerLaw(*in_dataptr, B, pars);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    Py_DECREF(in_iter);
    Py_DECREF(out_iter);
    Py_INCREF(out_array);
    return out_array;
}


static PyMethodDef Synchrotron[] = 
{
    {"J", (PyCFunction)J, METH_VARARGS, "calculate the emission coefficient of Synchrotron"},
    {"K", (PyCFunction)K, METH_VARARGS, "calculate the absorption coefficient of Synchrotron"},
    {NULL, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "Synchrotron",
        NULL,
        -1,
        Synchrotron,
        NULL,
        NULL,
        NULL,
        NULL
};

PyMODINIT_FUNC
PyInit_Synchrotron(void)
{
    PyObject *module = PyModule_Create(&moduledef);
    import_array();
    return module;
}
