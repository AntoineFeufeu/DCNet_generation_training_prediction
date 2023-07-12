/* changes
 17 MAR 2022 - in computing the envelope , change "float x, y" to "double x, y"
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;

void gsac_set_param_env(int ncmd, char **cmdstr)
{
		
}

void gsac_exec_env(void)
{
	int i, j, k, ntrc, npts;
	int n2;
	float depmax, depmin, depmen;
	int indmax, indmin;
	float dt, df;
	float env;
	double x, y;

	/* only apply the filtering if permitted */

	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;

	/* process the traces */
	for ( k=0 ; k < ntrc ; k ++){
		npts = sacdata[k].sachdr.ihdr[H_NPTS];
		dt = sacdata[k].sachdr.rhdr[H_DELTA];
		if(npts > 0){
			n2 = npow2(npts);
			/* allocate space for the spectra array */
			sacdata[k].sac_spectra = 
				(float *)realloc(sacdata[k].sac_spectra,
						 (n2+n2)*sizeof(float));
			/* fill the array with the time series */
			/* The even elements are the time series, the odd are 
			 * the imaginary part of the time series which is
			 * zero. The trace is effectively zero filled at 
			 * the end */
			for(i=0, j=0; i < n2 ; i++){
				if(i<npts)
					sacdata[k].sac_spectra[j] = 
						sacdata[k].sac_data[i];
				else
					sacdata[k].sac_spectra[j] = 0.0;
				j++;
				sacdata[k].sac_spectra[j++] = 0.0;
			}
			/* get the DFT */
			df = 0.0;
			four(sacdata[k].sac_spectra, n2, -1, &dt, &df);
			sacdata[k].df = df;
			sacdata[k].npow2 = n2;
			/* to get the Hilbert transform we zero out the
			 * negative frequencies and inverse transform */
			for(j= 0; j < n2 + n2 ; j++){
				if(j == 0)
					sacdata[k].sac_spectra[j] *= 1.0;
				else if(j == 1)
					sacdata[k].sac_spectra[j] *= 0.0;
				else if(j == n2)
					sacdata[k].sac_spectra[j] *= 1.0;
				else if(j == n2 +1)
					sacdata[k].sac_spectra[j] *= 0.0;
				else if(j < n2 && j > 0)
					sacdata[k].sac_spectra[j] *= 2.0;
				else
					sacdata[k].sac_spectra[j] = 0.0;
			}
			/* get the inverse FFT */
			four(sacdata[k].sac_spectra, n2, +1, &dt, &df);
			/* refill the original array */
			for(i=0, j=0; i < npts; i++){
				x = sacdata[k].sac_spectra[j++];
				y = sacdata[k].sac_spectra[j++];
				env = (float)sqrt(x*x + y*y);
					sacdata[k].sac_data[i] = env;
			}
			/* redetermine the depmax depmin depmen */
			getmxmn(sacdata[k].sac_data, npts,&depmax, &depmin, &depmen,&indmax,&indmin);
			sacdata[k].sachdr.rhdr[H_TIMMAX] = sacdata[k].sachdr.rhdr[H_B]  + ( indmax)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_TIMMIN] = sacdata[k].sachdr.rhdr[H_B]  + ( indmin)*sacdata[k].sachdr.rhdr[H_DELTA] ;
			sacdata[k].sachdr.rhdr[H_DEPMIN] = depmin;
			sacdata[k].sachdr.rhdr[H_DEPMAX] = depmax;
			sacdata[k].sachdr.rhdr[H_DEPMEN] = depmen;
		}
	}
}
