#if !defined(__calerf_h)
#define __calerf_h

// -------- netlib package by J. W. Cody ---------

namespace calerfpack
{
long double calerf_LD(const long double arg, const int jint);
long double erf_LD(const long double x);
long double erfc_LD(const long double x);
long double erex_LD(const long double x);
};
// -----------------------------------------------

#endif
