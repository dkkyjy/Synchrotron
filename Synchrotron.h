double U_B(double B);
double Nu_B(double B);
double t(double nu, double gamma, double B);

double g(double nu, double t, double n, double K, double nu_B);
double j_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max);
double J_PowerLaw(double nu, double B, double* pars);
double j_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max);
double J_BrokenPowerLaw(double nu, double B, double* pars);
PyObject* J(PyObject* nu, double B, int SpectrumType, PyObject* SpectrumPars){

double h(double nu, double t, double n, double K, double nu_B);
double k_PowerLaw(double nu, double B, double K, double n, double gamma_min, double gamma_max);
double K_PowerLaw(double nu, double B, double* pars);
double k_BrokenPowerLaw(double nu, double B, double K1, double K2, double n1, double n2, double gamma_min, double gamma_b, double gamma_max);
double K_BrokenPowerLaw(double nu, double B, double* pars);
PyArrayObject* K(PyArrayObject* nu, double B, int SpectrumType, PyObject* SpectrumPars);

PyArrayObject* Tau(PyArrayObject* k, double R);
PyArrayObject* I(PyArrayObject* j, PyArrayObject* k, double R);
PyArrayObject* Nu_obs(PyArrayObject* nu, double delta, double z);