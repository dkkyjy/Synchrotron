# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.8
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Synchrotron', [dirname(__file__)])
        except ImportError:
            import _Synchrotron
            return _Synchrotron
        if fp is not None:
            try:
                _mod = imp.load_module('_Synchrotron', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Synchrotron = swig_import_helper()
    del swig_import_helper
else:
    import _Synchrotron
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0



def U_B(B):
    return _Synchrotron.U_B(B)
U_B = _Synchrotron.U_B

def Nu_B(B):
    return _Synchrotron.Nu_B(B)
Nu_B = _Synchrotron.Nu_B

def t(nu, gamma, B):
    return _Synchrotron.t(nu, gamma, B)
t = _Synchrotron.t

def g(nu, t, n, K, nu_B):
    return _Synchrotron.g(nu, t, n, K, nu_B)
g = _Synchrotron.g

def j_PowerLaw(nu, B, K, n, gamma_min, gamma_max):
    return _Synchrotron.j_PowerLaw(nu, B, K, n, gamma_min, gamma_max)
j_PowerLaw = _Synchrotron.j_PowerLaw

def J_PowerLaw(nu, B, pars):
    return _Synchrotron.J_PowerLaw(nu, B, pars)
J_PowerLaw = _Synchrotron.J_PowerLaw

def j_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_b, gamma_max):
    return _Synchrotron.j_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_b, gamma_max)
j_BrokenPowerLaw = _Synchrotron.j_BrokenPowerLaw

def J_BrokenPowerLaw(nu, B, pars):
    return _Synchrotron.J_BrokenPowerLaw(nu, B, pars)
J_BrokenPowerLaw = _Synchrotron.J_BrokenPowerLaw

def J(nu, B, SpectrumType, SpectrumPars):
    return _Synchrotron.J(nu, B, SpectrumType, SpectrumPars)
J = _Synchrotron.J

def h(nu, t, n, K, nu_B):
    return _Synchrotron.h(nu, t, n, K, nu_B)
h = _Synchrotron.h

def k_PowerLaw(nu, B, K, n, gamma_min, gamma_max):
    return _Synchrotron.k_PowerLaw(nu, B, K, n, gamma_min, gamma_max)
k_PowerLaw = _Synchrotron.k_PowerLaw

def K_PowerLaw(nu, B, pars):
    return _Synchrotron.K_PowerLaw(nu, B, pars)
K_PowerLaw = _Synchrotron.K_PowerLaw

def k_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_b, gamma_max):
    return _Synchrotron.k_BrokenPowerLaw(nu, B, K1, K2, n1, n2, gamma_min, gamma_b, gamma_max)
k_BrokenPowerLaw = _Synchrotron.k_BrokenPowerLaw

def K_BrokenPowerLaw(nu, B, pars):
    return _Synchrotron.K_BrokenPowerLaw(nu, B, pars)
K_BrokenPowerLaw = _Synchrotron.K_BrokenPowerLaw

def K(nu, B, SpectrumType, SpectrumPars):
    return _Synchrotron.K(nu, B, SpectrumType, SpectrumPars)
K = _Synchrotron.K

def Tau(k, R):
    return _Synchrotron.Tau(k, R)
Tau = _Synchrotron.Tau

def Nu_obs(nu, delta, z):
    return _Synchrotron.Nu_obs(nu, delta, z)
Nu_obs = _Synchrotron.Nu_obs
# This file is compatible with both classic and new-style classes.


