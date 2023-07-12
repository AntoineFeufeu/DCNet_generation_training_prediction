#include        <stdio.h>
#include        <string.h>
#include "gsac.h"
#include "gsac_docommand.h"
#include "gsac_sac.h"
#include "gsac_sachdr.h"
#include "gsac_arg.h"
#include "csstim.h"

/* CHANGES:
        05 MAR 2018 - in December 2017, gsac_ch was changed so that
		the values ids DIST GCARC are not overwritten when
		lcalda = false.  This was to permit the sac file to be
		used for exploration ithout the subterfuge of making 
		up lat/lon.  A side effect is that lcalda must be set before the 
		test for proper stla/stlo/evla/evlo
	
		For safety so that the sequence must not be exactly
			ch lcalda true
			ch evla 12 evlo 34 stla 45 stlo 56
		but can be
			ch evla 12 evlo 34 stla 45 stlo 56
			ch lcalda true
		The wh checks this again
	30 OCT 2022 Support NVHDR = 7
		
                
*/

extern struct sacfile_ *sacdata;
extern int *sortptr;


void gsac_set_param_wh(int ncmd, char **cmdstr)
{
	int i;
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
}

void gsac_exec_wh(void)
{
	int k, ntrchdr;
	double stla, stlo, evla, evlo;
	double az, baz, dist, gcarc;
	/* safety - do not write header if within a cut - force user to
	 * overwrite */
	if(gsac_control.docut){
		printf("Cannot write headers while CUT is on\n");
		printf("Note that a WRITE will replace original file\n");
	} else {
		ntrchdr = gsac_control.number_iheaders;
		if(ntrchdr < 1)
			return;
		for ( k=0 ; k < ntrchdr ; k ++){
			/* safety introduced 05 MAR 2018 */
			stla = sacdata[k].sachdr.rhdr[H_STLA];
			stlo = sacdata[k].sachdr.rhdr[H_STLO];
			evla = sacdata[k].sachdr.rhdr[H_EVLA];
			evlo = sacdata[k].sachdr.rhdr[H_EVLO];
			if(stla != -12345. && stlo != -12345. 
				&& evla != -12345. && evlo != -12345.
				&& sacdata[k].sachdr.ihdr[H_LCALDA] == 1 ){
				delaz( evla,  evlo,  stla,  stlo,  &gcarc,  &az,  &baz,  &dist);
				sacdata[k].sachdr.rhdr[H_DIST] = dist;
				sacdata[k].sachdr.rhdr[H_AZ] = az;
				sacdata[k].sachdr.rhdr[H_BAZ] = baz;
				sacdata[k].sachdr.rhdr[H_GCARC] = gcarc;
			}
			bwsach(sacdata[k].sac_ofile_name,sacdata[k].sachdr);
		}
	}
}
