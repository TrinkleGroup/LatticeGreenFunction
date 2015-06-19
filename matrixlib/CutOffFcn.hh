/* file: CutOffFcn.hh
 * Defines the cut-off function for the divergent pieces
 * of the lattice Green function.
 * Where doing analytical integration, cut-off the 1/k^2, i/k, and k^0
 * pieces with this function.  Semi-continuum gets 1-fcutspline
 * wrapped into it.  Cubic spline from ALPHA_CUT to the minimum Wigner-Seitz
 * cell distance.
 */
#ifndef _CUT_OFF_FCN_H_
#define _CUT_OFF_FCN_H_

#define ALPHA_CUT 0.5

inline double fcutspline(double x) {
	return (1.0-x)*(1.0-x)/(1.0-ALPHA_CUT)/(1.0-ALPHA_CUT)*(3.0 - 2.0*(1.0 - x)/(1.0-ALPHA_CUT));
}

#endif
