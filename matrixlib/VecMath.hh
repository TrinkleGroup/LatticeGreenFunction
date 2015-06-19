#ifndef _VEC_MATH_H_
#define _VEC_MATH_H_

inline double vecmag(double x[3]) {
	return std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
}

inline double vecmag2(double x[3]) {
	return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}

/* n = t x m */
inline void crossprod(double t[3], double m[3], double n[3]) {
	n[0] = t[1]*m[2] - t[2]*m[1];
	n[1] = t[2]*m[0] - t[0]*m[2];
	n[2] = t[0]*m[1] - t[1]*m[0];
}

inline double vecdot(double a[3], double b[3]) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

#endif
