#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
		int imax,
		int jmax,
		double **U,
		double **V,
		const int wl,
		const int wr,
		const int wt,
		const int wb
);
/**
 * Function to set special boundary values.
 */
void spec_boundary_val(
		char *problem,
		int imax,
		int jmax,
		double **U,
		double **V
);

#endif
