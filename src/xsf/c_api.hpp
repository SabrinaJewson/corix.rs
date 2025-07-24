#include "vendored/include/xsf/airy.h"
#include <complex.h>

#ifndef PREFIX
# define PREFIX
#endif
// WTF C????
#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)
#define PREFIXED(ident) CONCAT(PREFIX, ident)

extern "C" void PREFIXED(airy)(double x, double *ai, double *aid, double *bi, double *bid) { xsf::airy(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airye)(double x, double *ai, double *aid, double *bi, double *bid) { xsf::airye(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airy_complex)(std::complex<double> x, std::complex<double> *ai, std::complex<double> *aid, std::complex<double> *bi, std::complex<double> *bid) { xsf::airy(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airye_complex)(std::complex<double> x, std::complex<double> *ai, std::complex<double> *aid, std::complex<double> *bi, std::complex<double> *bid) { xsf::airye(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airyzo)(int nt, int kf, double *xa, double *xb, double *xc, double *xd) { xsf::airyzo(nt, kf, xa, xb, xc, xd); }
extern "C" void PREFIXED(itairy)(double x, double *apt, double *bpt, double *ant, double *bnt) { xsf::itairy(x, *apt, *bpt, *ant, *bnt); };
