#include "boundary_val.h"
#include "helper.h"

/* ----------------------------------------------------------------------- */
/*                             Set Boundary Conditions                     */
/* ----------------------------------------------------------------------- */

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
) {

	int i,j;
	/*Set values for the left boundary*/
	switch(wl){
	case NO_SLIP:
		for (j = 1; j < jmax + 1; j++){
			/*U velocity on left boundary */
			U[0][j] = 0;
			/*V velocity left boundary */
			V[0][j]=-1*V[1][j];
		}
		break;
	case FREE_SLIP:
		for (j = 1; j < jmax + 1; j++){
			/*U velocity on left boundary */
			U[0][j] = 0;
			/*V velocity left boundary */
			V[0][j] = V[1][j];
		}
		break;
	case OUTFLOW:
		for (j = 1; j < jmax + 1; j++){
			/*U velocity on left boundary */
			U[0][j] = U[1][j];
			/*V velocity left boundary */
			V[0][j]= V[1][j];
		}
		break;
	default:
		break;
	}

	/*Set values for the right boundary*/
	switch(wr){
	case NO_SLIP:
		for (j = 1; j < jmax + 1; j++){
			U[imax][j] = 0;
			V[imax+1][j]=-1*V[imax][j];
		}
		break;
	case FREE_SLIP:
		for (j = 1; j < jmax + 1; j++){
			U[imax][j] = 0;
			V[imax+1][j] = V[imax][j];
		}
		break;
	case OUTFLOW:
		for (j = 1; j < jmax + 1; j++){
			U[imax][j] = U[imax-1][j];
			V[imax+1][j]= V[imax][j];
		}
		break;
	default:
		break;
	}

	/*Set values for the top boundary*/
	switch(wt){
	case NO_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = 0;
			U[i][jmax+1]= -1*U[i][jmax];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = 0;
			U[i][jmax+1] = U[i][jmax];
		}
		break;
	case OUTFLOW:
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = V[i][jmax-1];
			U[i][jmax+1] = U[i][jmax];
		}
		break;
	default:
		break;
	}

	/*Set values for the bottom boundary*/
	switch(wb){
	case NO_SLIP:
		for (i = 1; i <= imax; i++){
			V[i][0] = 0;
			U[i][0]=-1*U[i][1];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i <= imax; i++){
			V[i][0] = 0;
			U[i][0] = U[i][1];
		}
		break;
	case OUTFLOW:
		for (i = 1; i <= imax; i++){
			V[i][0] = V[i][1];
			U[i][0]= U[i][1];
		}
		break;
	default:
		break;
	}
}
