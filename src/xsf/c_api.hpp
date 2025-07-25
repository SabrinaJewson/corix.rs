#include "vendored/include/xsf/airy.h"
#include "vendored/include/xsf/ellip.h"

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

extern "C" void PREFIXED(ellipj)(double u, double m, double *sn, double *cn, double *dn, double *ph) { xsf::ellipj(u, m, *sn, *cn, *dn, *ph); }
extern "C" double PREFIXED(ellipk)(double m) { return xsf::ellipk(m); }
extern "C" double PREFIXED(ellipkm1)(double m) { return xsf::ellipkm1(m); }
extern "C" double PREFIXED(ellipkinc)(double φ, double m) { return xsf::ellipkinc(φ, m); }
extern "C" double PREFIXED(ellipe)(double m) { return xsf::ellipe(m); }
extern "C" double PREFIXED(ellipeinc)(double φ, double m) { return xsf::ellipeinc(φ, m); }
