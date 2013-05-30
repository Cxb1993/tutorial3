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
		const int wb,
		int **Flag
) {

	int i,j;
	/*Initialize corners*/
	U[0][0]=0.0;
	U[0][jmax+1]=0.0;
	U[imax+1][0]=0.0;
	U[imax+1][jmax+1]=0.0;
	V[0][0]=0.0;
	V[0][jmax+1]=0.0;
	V[imax+1][0]=0.0;
	V[imax+1][jmax+1]=0.0;

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
		for (i = 1; i < imax + 1; i++){
			V[i][0] = 0;
			U[i][0]=-1*U[i][1];
		}
		break;
	case FREE_SLIP:
		for (i = 1; i < imax + 1; i++){
			V[i][0] = 0;
			U[i][0] = U[i][1];
		}
		break;
	case OUTFLOW:
		for (i = 1; i < imax + 1; i++){
			V[i][0] = V[i][1];
			U[i][0]= U[i][1];
		}
		break;
	default:
		break;
	}

	/**
	 * Loop to check for boundary cells in the inner domain
	 */
	for(i = 1; i < imax+1; i++){
		for(j = 1; j < jmax+1; j++){
			if((Flag[i][j]&31)==B_N){
				V[i][j]=0;
				U[i-1][j]=-1*U[i-1][j+1];
				U[i][j]=-1*U[i][j+1];
			}
			else if((Flag[i][j]&31)==B_S){
				V[i][j-1]=0;
				U[i-1][j]=-1*U[i-1][j-1];
				U[i][j]=-1*U[i][j-1];
			}
			else if((Flag[i][j]&31)==B_W){
				U[i-1][j]=0;
				V[i][j-1]=-1*V[i-1][j-1];
				V[i][j]=-1*V[i-1][j];
			}
			else if((Flag[i][j]&31)==B_O){
				U[i][j]=0;
				V[i][j-1]=-1*V[i+1][j-1];
				V[i][j]=-1*V[i+1][j];
			}
			else if((Flag[i][j]&31)==B_NO){
				U[i][j]=0;
				U[i-1][j]=-1*U[i-1][j+1];
				V[i][j]=0;
				V[i][j-1]=-1*V[i+1][j-1];
			}
			else if((Flag[i][j]&31)==B_NW){
				U[i-1][j]=0;
				U[i][j]=-1*U[i][j+1];
				V[i][j]=0;
				V[i][j-1]=-1*V[i-1][j-1];
			}
			else if((Flag[i][j]&31)==B_SO){
				U[i][j]=0;
				U[i-1][j]=-1*U[i-1][j-1];
				V[i][j-1]=0;
				V[i][j]=-1*V[i+1][j];
			}
			else if((Flag[i][j]&31)==B_SW){
				U[i-1][j]=0;
				U[i][j]=-1*U[i][j-1];
				V[i][j]=-1*V[i-1][j];
				V[i][j-1]=0;
			}
		}
	}

}

void spec_boundary_val(
		char *problem,
		int imax,
		int jmax,
		double **U,
		double **V
){
	int i,j;
	if(strcmp(problem,"cavity")==0){
		for (i = 1; i < imax + 1; i++){
			V[i][jmax] = 0;
			U[i][jmax+1]= 2.0-1*U[i][jmax];
		}
	}
	if(strcmp(problem,"KarmanVortexStreet")==0){
		for (j = 1; j < jmax + 1; j++){
			U[0][j] = 1.0;
			V[0][j]= 0.0;
		}
	}
	/*if(strcmp(problem,"PlaneShearFlow")==0){
		for (j = 1; j < jmax + 1; j++){
			U[0][j] = 0.0;
			V[0][j]= 0.0;
		}
	}*/
	if(strcmp(problem,"FlowOverStep")==0){
		for (j = 1; j <= (jmax)/2; j++){
			U[0][j] = 0.0;
			V[0][j]= 0.0;
		}
		for (j = ((jmax)/2) + 1; j < jmax + 1; j++){
			U[0][j] = 1.0;
			V[0][j]= 0.0;
		}
	}
}
