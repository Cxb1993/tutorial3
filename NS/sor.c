#include "sor.h"
#include <math.h>
#include "helper.h"

void sor(
		double omg,
		double dx,
		double dy,
		int    imax,
		int    jmax,
		double **P,
		double **RS,
		double *res,
		double lp,
		double rp,
		double dp,
		int **Flag
) {
	int i,j;
	double rloc;
	double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
	int count;

	count = 0;

	/* SOR iteration */
	for(i = 1; i <= imax; i++) {
		for(j = 1; j<=jmax; j++) {
			if((Flag[i][j]&B_C)==B_C){
				P[i][j] = (1.0-omg)*P[i][j] + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) +
						( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
				count++;
			}
			else if((Flag[i][j]&31)==B_N){
				P[i][j]=P[i][j+1];
			}
			else if((Flag[i][j]&31)==B_S){
				P[i][j]=P[i][j-1];
			}
			else if((Flag[i][j]&31)==B_W){
				P[i][j]=P[i-1][j];
			}
			else if((Flag[i][j]&31)==B_O){
				P[i][j]=P[i+1][j];
			}
			else if((Flag[i][j]&31)==B_NO){
				P[i][j]=(P[i+1][j]+P[i][j+1])/2.0;
			}
			else if((Flag[i][j]&31)==B_NW){
				P[i][j]=(P[i][j+1]+P[i-1][j])/2.0;
			}
			else if((Flag[i][j]&31)==B_SO){
				P[i][j]=(P[i][j-1]+P[i+1][j])/2.0;
			}
			else if((Flag[i][j]&31)==B_SW){
				P[i][j]=(P[i][j-1]+P[i-1][j])/2.0;
			}
		}
	}


	/* set boundary values */
			for(i = 1; i <= imax; i++) {

				P[i][0] = P[i][1];
				P[i][jmax+1] = P[i][jmax];
			}
			for(j = 1; j <= jmax; j++) {
				/*left*/

				if((Flag[0][j]&P_L)==P_L){
					if(lp>=0){
						P[0][j] = 2*lp-P[1][j];
					}
					else if(dp!=0){
						P[0][j] = 2*dp-P[1][j];
					}
				}
				else{
					P[0][j] = P[1][j];
				}
				/*right*/

				if((Flag[imax+1][j]&P_R)==P_R){
					if(rp>=0){
						P[imax+1][j] = 2*rp-P[imax][j];
					}
					else if(dp!=0){
						P[imax+1][j] = -P[imax][j];
					}
				}
				else{
					P[imax+1][j] = P[imax][j];
				}
			}

	/* compute the residual */
	rloc = 0.0;
	for(i = 1; i <= imax; i++) {
		for(j = 1; j <= jmax; j++) {
			if((Flag[i][j]&B_C)==B_C){
				rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
						( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
			}
		}
	}
	rloc = rloc/((double)count);
	rloc = sqrt(rloc);
	/* set residual */
	*res = rloc;
}
