#include "vendored/include/xsf/airy.h"
#include "vendored/include/xsf/bessel.h"
#include "vendored/include/xsf/cephes/kn.h"
#include "vendored/include/xsf/ellip.h"
#include "vendored/include/xsf/specfun/specfun.h"
#include "vendored/include/xsf/sph_bessel.h"
#include "vendored/include/xsf/wright_bessel.h"

using std::complex;

#ifndef PREFIX
# define PREFIX
#endif
// WTF C????
#define CONCAT_(a, b) a ## b
#define CONCAT(a, b) CONCAT_(a, b)
#define PREFIXED(ident) CONCAT(PREFIX, ident)

extern "C" void PREFIXED(airy)(double x, double *ai, double *aid, double *bi, double *bid) { xsf::airy(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airye)(double x, double *ai, double *aid, double *bi, double *bid) { xsf::airye(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airy_complex)(complex<double> x, complex<double> *ai, complex<double> *aid, complex<double> *bi, complex<double> *bid) { xsf::airy(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airye_complex)(complex<double> x, complex<double> *ai, complex<double> *aid, complex<double> *bi, complex<double> *bid) { xsf::airye(x, *ai, *aid, *bi, *bid); }
extern "C" void PREFIXED(airyzo)(int nt, int kf, double *xa, double *xb, double *xc, double *xd) { xsf::airyzo(nt, kf, xa, xb, xc, xd); }
extern "C" void PREFIXED(itairy)(double x, double *apt, double *bpt, double *ant, double *bnt) { xsf::itairy(x, *apt, *bpt, *ant, *bnt); };

extern "C" void PREFIXED(ellipj)(double u, double m, double *sn, double *cn, double *dn, double *ph) { xsf::ellipj(u, m, *sn, *cn, *dn, *ph); }
extern "C" double PREFIXED(ellipk)(double m) { return xsf::ellipk(m); }
extern "C" double PREFIXED(ellipkm1)(double m) { return xsf::ellipkm1(m); }
extern "C" double PREFIXED(ellipkinc)(double φ, double m) { return xsf::ellipkinc(φ, m); }
extern "C" double PREFIXED(ellipe)(double m) { return xsf::ellipe(m); }
extern "C" double PREFIXED(ellipeinc)(double φ, double m) { return xsf::ellipeinc(φ, m); }

extern "C" double PREFIXED(cyl_bessel_j)(double v, double x) { return xsf::cyl_bessel_j(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_j_complex)(double v, complex<double> x) { return xsf::cyl_bessel_j(v, x); }
extern "C" double PREFIXED(cyl_bessel_je)(double v, double x) { return xsf::cyl_bessel_je(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_je_complex)(double v, complex<double> x) { return xsf::cyl_bessel_je(v, x); }
extern "C" double PREFIXED(cyl_bessel_j0)(double x) { return xsf::cyl_bessel_j0(x); }
extern "C" double PREFIXED(cyl_bessel_j1)(double x) { return xsf::cyl_bessel_j1(x); }

extern "C" double PREFIXED(cyl_bessel_y)(double v, double x) { return xsf::cyl_bessel_y(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_y_complex)(double v, complex<double> x) { return xsf::cyl_bessel_y(v, x); }
extern "C" double PREFIXED(cyl_bessel_ye)(double v, double x) { return xsf::cyl_bessel_ye(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_ye_complex)(double v, complex<double> x) { return xsf::cyl_bessel_ye(v, x); }
extern "C" double PREFIXED(cyl_bessel_y0)(double x) { return xsf::cyl_bessel_y0(x); }
extern "C" double PREFIXED(cyl_bessel_y1)(double x) { return xsf::cyl_bessel_y1(x); }
extern "C" double PREFIXED(cyl_bessel_yn)(int n, double x) { return xsf::cephes::yn(n, x); }

extern "C" double PREFIXED(cyl_bessel_i)(double v, double x) { return xsf::cyl_bessel_i(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_i_complex)(double v, complex<double> x) { return xsf::cyl_bessel_i(v, x); }
extern "C" double PREFIXED(cyl_bessel_ie)(double v, double x) { return xsf::cyl_bessel_ie(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_ie_complex)(double v, complex<double> x) { return xsf::cyl_bessel_ie(v, x); }
extern "C" double PREFIXED(cyl_bessel_i0)(double x) { return xsf::cyl_bessel_i0(x); }
extern "C" double PREFIXED(cyl_bessel_i0e)(double x) { return xsf::cyl_bessel_i0e(x); }
extern "C" double PREFIXED(cyl_bessel_i1)(double x) { return xsf::cyl_bessel_i1(x); }
extern "C" double PREFIXED(cyl_bessel_i1e)(double x) { return xsf::cyl_bessel_i1e(x); }

extern "C" double PREFIXED(cyl_bessel_k)(double v, double x) { return xsf::cyl_bessel_k(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_k_complex)(double v, complex<double> x) { return xsf::cyl_bessel_k(v, x); }
extern "C" double PREFIXED(cyl_bessel_ke)(double v, double x) { return xsf::cyl_bessel_ke(v, x); }
extern "C" complex<double> PREFIXED(cyl_bessel_ke_complex)(double v, complex<double> x) { return xsf::cyl_bessel_ke(v, x); }
extern "C" double PREFIXED(cyl_bessel_k0)(double x) { return xsf::cyl_bessel_k0(x); }
extern "C" double PREFIXED(cyl_bessel_k0e)(double x) { return xsf::cyl_bessel_k0e(x); }
extern "C" double PREFIXED(cyl_bessel_k1)(double x) { return xsf::cyl_bessel_k1(x); }
extern "C" double PREFIXED(cyl_bessel_k1e)(double x) { return xsf::cyl_bessel_k1e(x); }
extern "C" double PREFIXED(cyl_bessel_kn)(int n, double x) { return xsf::cephes::kn(n, x); }

extern "C" complex<double> PREFIXED(cyl_hankel_1)(double v, complex<double> z) { return xsf::cyl_hankel_1(v, z); }
extern "C" complex<double> PREFIXED(cyl_hankel_1e)(double v, complex<double> z) { return xsf::cyl_hankel_1e(v, z); }
extern "C" complex<double> PREFIXED(cyl_hankel_2)(double v, complex<double> z) { return xsf::cyl_hankel_2(v, z); }
extern "C" complex<double> PREFIXED(cyl_hankel_2e)(double v, complex<double> z) { return xsf::cyl_hankel_2e(v, z); }

extern "C" double PREFIXED(wright_bessel)(double a, double b, double x) { return xsf::wright_bessel(a, b, x); }
extern "C" double PREFIXED(log_wright_bessel)(double a, double b, double x) { return xsf::log_wright_bessel(a, b, x); }

extern "C" void PREFIXED(lamv)(double v, double x, double *vm, double *vl, double *dl) { xsf::specfun::lamv(v, x, vm, vl, dl); }
extern "C" void PREFIXED(lamn)(int n, double x, int *nm, double *vl, double *dl) { xsf::specfun::lamn(n, x, nm, vl, dl); }

extern "C" void PREFIXED(it1j0y0)(double x, double *j0int, double *y0int) { xsf::it1j0y0(x, *j0int, *y0int); }
extern "C" void PREFIXED(it2j0y0)(double x, double *j0int, double *y0int) { xsf::it2j0y0(x, *j0int, *y0int); }
extern "C" void PREFIXED(it1i0k0)(double x, double *i0int, double *k0int) { xsf::it1i0k0(x, *i0int, *k0int); }
extern "C" void PREFIXED(it2i0k0)(double x, double *i0int, double *k0int) { xsf::it2i0k0(x, *i0int, *k0int); }
extern "C" double PREFIXED(besselpoly)(double a, double λ, double v) { return xsf::besselpoly(a, λ, v); }

extern "C" double PREFIXED(sph_bessel_j)(uint64_t n, double x) { return xsf::sph_bessel_j(n, x); }
extern "C" double PREFIXED(sph_bessel_y)(uint64_t n, double x) { return xsf::sph_bessel_y(n, x); }
extern "C" double PREFIXED(sph_bessel_i)(uint64_t n, double x) { return xsf::sph_bessel_i(n, x); }
extern "C" double PREFIXED(sph_bessel_k)(uint64_t n, double x) { return xsf::sph_bessel_k(n, x); }
extern "C" double PREFIXED(sph_bessel_j_jac)(uint64_t n, double x) { return xsf::sph_bessel_j_jac(n, x); }
extern "C" double PREFIXED(sph_bessel_y_jac)(uint64_t n, double x) { return xsf::sph_bessel_y_jac(n, x); }
extern "C" double PREFIXED(sph_bessel_i_jac)(uint64_t n, double x) { return xsf::sph_bessel_i_jac(n, x); }
extern "C" double PREFIXED(sph_bessel_k_jac)(uint64_t n, double x) { return xsf::sph_bessel_k_jac(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_j_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_j(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_y_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_y(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_i_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_i(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_k_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_k(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_j_jac_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_j_jac(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_y_jac_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_y_jac(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_i_jac_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_i_jac(n, x); }
extern "C" complex<double> PREFIXED(sph_bessel_k_jac_complex)(uint64_t n, complex<double> x) { return xsf::sph_bessel_k_jac(n, x); }

// Currently unused… (defined in specfun.h)
// extern "C" xsf::specfun::Status PREFIXED(jdzo)(int nt, double *zo, int *n, int *m, int *p) { return xsf::specfun::jdzo(nt, zo, n, m, p); }
// extern "C" void PREFIXED(jyzo)(int n, int nt, double *rj0, double *rj1, double *ry0, double *ry1) { xsf::specfun::jyzo(n, nt, rj0, rj1, ry0, ry1); }
// extern "C" void PREFIXED(cyzo)(int nt, int kf, int kc, complex<double> *zo, complex<double> *zv) { xsf::specfun::cyzo(nt, kf, kc, zo, zv); }

template <typename T>
struct Slice {
	T *ptr;
	int len;
	int extent(int n) { return this->len; }
	T& operator[](size_t i) { return this->ptr[i]; }
};
extern "C" void PREFIXED(rctj)(double x, int *nm, double *rj, double *dj, int len) { return xsf::rctj(x, nm, (Slice<double>){ rj, len }, (Slice<double>){ dj, len }); }
extern "C" void PREFIXED(rcty)(double x, int *nm, double *rj, double *dj, int len) { return xsf::rcty(x, nm, (Slice<double>){ rj, len }, (Slice<double>){ dj, len }); }
