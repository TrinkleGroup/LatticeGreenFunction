#ifndef _GSL_COMPLEX_OP_H_
#define _GSL_COMPLEX_OP_H_
#include <gsl/gsl_complex.h>

gsl_complex operator-(const gsl_complex &s1, const gsl_complex &s2);
gsl_complex operator-(const gsl_complex &s1, const double &s2);
gsl_complex operator-(const double &s1, const gsl_complex &s2);
gsl_complex operator-(const gsl_complex &s);
gsl_complex operator+(const gsl_complex &s1, const gsl_complex &s2);
gsl_complex operator+(const gsl_complex &s1, const double &s2);
gsl_complex operator+(const double &s1, const gsl_complex &s2);
gsl_complex operator*(const gsl_complex &s1, const gsl_complex &s2);
gsl_complex operator/(const gsl_complex &s1, const gsl_complex &s2);
gsl_complex operator*(const gsl_complex &s1, const double &s2);
gsl_complex operator*(const double &s1, const gsl_complex &s2);
gsl_complex operator/(const gsl_complex &s1, const double &s2);

#endif
