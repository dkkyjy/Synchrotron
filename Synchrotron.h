#include <Python.h>
#include <numpy/arrayobject.h>

const double c = 3.00e10;
const double e = 4.803e-10;
const double m_e = 9.11e-28;
const double sigma_T = 6.65e-25;

double U_B(double B);

double j_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max);
double J_PowerLaw(double nu, double B, double* pars);
double j_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max);
double J_BrokenPowerLaw(double nu, double B, double* pars);
static PyObject* J(PyObject* self, PyObject* args);

double k_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max);
double K_PowerLaw(double nu, double B, double* pars);
double k_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max);
double K_BrokenPowerLaw(double nu, double B, double* pars);
static PyObject* K(PyObject* self, PyObject* args);

//PyObject* Tau(PyObject* k, double R);
//PyObject* I(PyObject* j, PyObject* k, double R);
//PyObject* Nu_obs(PyObject* nu, double delta, double z);
