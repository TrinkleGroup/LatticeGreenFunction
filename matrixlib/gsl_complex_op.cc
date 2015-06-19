#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "gsl_complex_op.hh"

gsl_complex operator-(const gsl_complex &s1, const gsl_complex &s2) {
	return gsl_complex_sub(s1,s2);
}

gsl_complex operator-(const gsl_complex &s1, const double &s2) {
	return gsl_complex_sub_real(s1,s2);
}

gsl_complex operator-(const double &s1, const gsl_complex &s2) {
	return gsl_complex_negative(gsl_complex_sub_real(s2,s1));
}

gsl_complex operator-(const gsl_complex &s) {
	return gsl_complex_negative(s);
}

gsl_complex operator+(const gsl_complex &s1, const gsl_complex &s2) {
	return gsl_complex_add(s1,s2);
}

gsl_complex operator+(const gsl_complex &s1, const double &s2) {
	return gsl_complex_add_real(s1,s2);
}

gsl_complex operator+(const double &s1, const gsl_complex &s2) {
	return gsl_complex_add_real(s2,s1);
}

gsl_complex operator*(const gsl_complex &s1, const gsl_complex &s2) {
	return gsl_complex_mul(s1,s2);
}

gsl_complex operator/(const gsl_complex &s1, const gsl_complex &s2) {
	return gsl_complex_div(s1,s2);
}

gsl_complex operator*(const gsl_complex &s1, const double &s2) {
	return gsl_complex_mul_real(s1,s2);
}

gsl_complex operator*(const double &s1, const gsl_complex &s2) {
	return gsl_complex_mul_real(s2,s1);
}

gsl_complex operator/(const gsl_complex &s1, const double &s2) {
	return gsl_complex_div_real(s1,s2);
}
