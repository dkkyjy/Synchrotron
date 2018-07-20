#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "Synchrotron.h"
#include "Python.h"
#include <numpy/arrayobject.h>

double c = 3.00e10;
double e = 4.803e-10;
double m_e = 9.11e-28;
double sigma_T = 6.65e-25;

double U_B(double B){
    return pow(B, 2) / (8 * M_PI);
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

//SpectrumType: 1-PowerLaw; 2-BrokenPowerLaw
PyObject* J(PyObject* nu, double B, int SpectrumType, PyObject* SpectrumPars){
    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ipar++)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    PyArrayObject* in_array;
    PyObject* out_array;
    PyArg_ParseTuple(nu, "O!", &PyArray_Type, &in_array);
    out_array = PyArray_NewLikeArray(in_array, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(in_array);
    PyArrayIterObject* out_iter = PyArray_IterNew(out_array);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        if (SpectrumType == 1):
            double* out_dataptr = J_PowerLaw(*in_dataptr);
        if (SpectrumType == 2):
            double* out_dataptr = J_BrokenPowerLaw(*in_dataptr);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    PyArray_DECREF(in_iter);
    PyArray_DECREF(out_iter);
    PyArray_INCREF(out_array);
    return out_array;
}

    /*
    Py_ssize_t num = PyList_Size(nuList);
    PyObject* jList = PyList_New(num);
    for(int i=0; i<num; i++){
        double J;
        double nu = PyFloat_AsDouble(PyList_GetItem(nuList, i));
        if (SpectrumType == 1)
            J = J_PowerLaw(nu, B, pars);
        if (SpectrumType == 2)
            J = J_BrokenPowerLaw(nu, B, pars);
        PyList_SetItem(jList, i, PyFloat_FromDouble(J));
    }
    return jList;
    */


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


//SpectrumType: 1-PowerLaw; 2-BrokenPowerLaw
PyArrayObject* K(PyArrayObject* nu, double B, int SpectrumType, PyObject* SpectrumPars){
    Py_ssize_t parNum = PyTuple_Size(SpectrumPars);
    double pars[parNum];
    for(int ipar=0; ipar<parNum; ipar++)
        pars[ipar] = PyFloat_AsDouble(PyTuple_GetItem(SpectrumPars, ipar));

    PyArrayObject* k = PyArray_NewLikeArray(nu, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(nu);
    PyArrayIterObject* out_iter = PyArray_IterNew(k);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        if (SpectrumType == 1):
            double* out_dataptr = K_PowerLaw(*in_dataptr);
        if (SpectrumType == 2):
            double* out_dataptr = K_BrokenPowerLaw(*in_dataptr);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    PyArray_DECREF(in_iter);
    PyArray_DECREF(out_iter);
    PyArray_INCREF(k);
    return k;
}

    /*
    Py_ssize_t num = PyList_Size(nuList);
    PyObject* kList = PyList_New(num);
    for(int i=0; i<num; i++){
        double K;
        double nu = PyFloat_AsDouble(PyList_GetItem(nuList, i));
        if (SpectrumType == 1)
            K = K_PowerLaw(nu, B, pars);
        if (SpectrumType == 2)
            K = K_BrokenPowerLaw(nu, B, pars);
        PyList_SetItem(kList, i, PyFloat_FromDouble(K));
    }
    return kList;
    */

PyArrayObject* Tau(PyArrayObject* k, double R){
    PyArrayObject* tau = PyArray_NewLikeArray(k, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(nu);
    PyArrayIterObject* out_iter = PyArray_IterNew(k);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = 2 * R * *in_dataptr;
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    PyArray_DECREF(in_iter);
    PyArray_DECREF(out_iter);
    PyArray_INCREF(tau);
    return tau;
}

PyArrayObject* I(PyArrayObject* j, PyArrayObject* k, double R){
    PyArrayObject* tau = Tau(k, R)
    PyArrayObject* I = PyArray_NewLikeArray(j, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* j_iter = PyArray_IterNew(j);
    PyArrayIterObject* k_iter = PyArray_IterNew(k);
    PyArrayIterObject* tau_iter = PyArray_IterNew(tau);
    PyArrayIterObject* I_iter = PyArray_IterNew(I);
    while(j_iter->index < j_iter->size && I_iter->index < I_iter->size){
        double* j_dataptr = j_iter->dataptr;
        double* k_dataptr = k_iter->dataptr;
        double* tau_dataptr = tau_iter->dataptr;
        double* I_dataptr = I_iter->dataptr;
        *I_dataptr = *j_dataptr / *k_dataptr * (1 - 2/pow(*tau_dataptr, 2) * (1 - exp(-(*tau_dataptr)*(*tau_dataptr + 1))));
        PyArray_ITER_NEXT(j_iter);
        PyArray_ITER_NEXT(k_iter);
        PyArray_ITER_NEXT(tau_iter);
        PyArray_ITER_NEXT(I_iter);
    }
    PyArray_DECREF(j_iter);
    PyArray_DECREF(k_iter);
    PyArray_DECREF(tau_iter);
    PyArray_DECREF(I_iter);
    PyArray_INCREF(I);
    return I;
}

    /*
    Py_ssize_t num = PyList_Size(nuList);
    PyObject* IList = PyList_New(num);
    for(int i=0; i<num; i++){
        double j = PyFloat_AsDouble(PyList_GetItem(jList, i));
        double k = PyFloat_AsDouble(PyList_GetItem(kList, i));
        double nu = PyFloat_AsDouble(PyList_GetItem(nuList, i));
        double tau = 2 * R * k;
        double I = j / k * (1 - 2 / pow(tau, 2) * (1 - exp(-tau) * (tau + 1)));
        PyList_SetItem(IList, i, PyFloat_FromDouble(I));
    }
    return IList;
    */

/*
PyObject* L(PyObject* nuList, double B, double R, int SpectrumType, PyObject* SpectrumPars){
    PyObject* IList = I(nuList, B, R, SpectrumType, SpectrumPars);

    Py_ssize_t num = PyList_Size(nuList);
    PyObject* LList = PyList_New(num);
    for(int i=0; i<num; i++){
        double I = PyFloat_AsDouble(PyList_GetItem(IList, i));
        double L = 4 * pow(M_PI, 2) * R * I;
        PyList_SetItem(LList, i, PyFloat_FromDouble(L));
    }
    return LList;
}

PyObject* F(PyObject* nuList, double B, double R, double delta, double dL, double z, int SpectrumType, PyObject* SpectrumPars){
    PyObject* LList = L(nuList, B, R, SpectrumType, SpectrumPars);

    Py_ssize_t num = PyList_Size(nuList);
    PyObject* FList = PyList_New(num);
    for(int i=0; i<num; i++){
        double L = PyFloat_AsDouble(PyList_GetItem(LList, i));
        double F = L / (4 * M_PI * pow(dL, 2)) * pow(delta, 3) * (1+z);
        PyList_SetItem(FList, i, PyFloat_FromDouble(F));
    }
    return FList;
}
*/

PyArrayObject* Nu_obs(PyArrayObject* nu, double delta, double z){
    PyArrayObject* nu_obs = PyArray_NewLikeArray(nu, NPY_ANYORDER, NULL, 0);
    PyArrayIterObject* in_iter = PyArray_IterNew(nu);
    PyArrayIterObject* out_iter = PyArray_IterNew(nu_obs);
    while(in_iter->index < in_iter->size && out_iter->index < out_iter->size){
        double* in_dataptr = in_iter->dataptr;
        double* out_dataptr = out_iter->dataptr;
        *out_dataptr = *in_dataptr * delta / (1+z);
        PyArray_ITER_NEXT(in_iter);
        PyArray_ITER_NEXT(out_iter);
    }
    PyArray_DECREF(in_iter);
    PyArray_DECREF(out_iter);
    PyArray_INCREF(nu_obs);
    return nu_obs;
}