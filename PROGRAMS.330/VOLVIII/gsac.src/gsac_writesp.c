/* CHANGES:
	17 MAR 2016 - added ability to output phase, real and imaginary parts of spectra
		      to previous amplitude spectrum
	20 NOV 2019 - DEPMAX DEPMIN DEPMEN were not correctly set. The three lines were moved after
		      the call getmxmn
        17 MAR 2022 - in computing the amplitude spectra, change "float tr, ti" to "double tr, ti"
                      
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

#define    NOUT   4

static  int write_sp_ctl[NOUT];

#define	WRITESP_DFLT	0
#define	WRITESP_AM	1
#define	WRITESP_PH	2
#define	WRITESP_RL	3
#define	WRITESP_IM	4

#define W_AM 0
#define W_PH 1
#define W_RL 2
#define W_IM 3


struct arghdr writesparg[] = {
	{WRITESP_DFLT, "DEFAULT", IHDR, NO, 0, NO, "", -1},
	{WRITESP_AM  , "AM"  , RHDR, NO, 0, NO, "", 1},
	{WRITESP_PH  , "PH"  , RHDR, NO, 0, NO, "", 1},
	{WRITESP_RL  , "RL"  , RHDR, NO, 0, NO, "", 1},
	{WRITESP_IM  , "IM"  , RHDR, NO, 0, NO, "", 1},
	{0           ,	""   , IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float writesp_real[10];
int   writesp_int [10];
int   writesp_yn;
int   writesp_num;

static float *temp = (float *)NULL;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_writesp(int ncmd, char **cmdstr)
{
	int i, j;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* initialize */
        for ( i= 0 ; i < NOUT ; i++)
		write_sp_ctl[i] = NO;
 	/* if no option on the command line, continue but force amplitude spectra output */
	if(ncmd == 1){
		write_sp_ctl[W_AM] = YES;
		return;
	}
	if(testarg(ncmd, cmdstr, writesparg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; writesparg[i].key[0] != '\0' ; i++){
		if(writesparg[i].used > 0){
			if(writesparg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, writesp_real);
			} else if(writesparg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, writesp_int );
			} else if(writesparg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, &writesp_yn );
			} else if(writesparg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, writesparg[i].key, 
					writesparg[i].mfit,writesparg[i].narg, &writesp_num );
			}
			switch(writesparg[i].id){
				case WRITESP_AM:
					write_sp_ctl[W_AM] = YES;
					break;
				case WRITESP_PH:
					write_sp_ctl[W_PH] = YES;
					break;
				case WRITESP_RL:
					write_sp_ctl[W_RL] = YES;
					break;
				case WRITESP_IM:
					write_sp_ctl[W_IM] = YES;
					break;
			}
		}
	}
	/* if nothing has been set force output of amplitude spectrum */
        for(i=0, j=0 ; i < NOUT ; i++){
		j+= write_sp_ctl[i];
	}
	if(j == 0 )
		write_sp_ctl[W_AM] = YES;
}

void gsac_exec_writesp(void)
{
	int i, j, k, ntrc, n2, n21, ncmd;
	float depmax, depmin, depmen;
        double tr, ti;
	int indmax, indmin;
	char t[1000];
	struct sacfile_ tsac;
	/* if there are no traces return */
	ntrc = gsac_control.number_itraces;
	if(ntrc < 1)
		return;
        if(gsac_control.fft == NO){
                printf("Execute FFT first before trying to write spectra\n");
                return;
        }
	/* first get the spectra even if it has already been - this is not
		too time consuming and also ensures that the spectra is of the
		current trace */
	gsac_exec_fft();
	/* now create the create the output and write the spectra - 
		allocate a structure for sac
		fill trace with spectra
		bwsac
	*/
	for(ncmd = 0 ; ncmd < NOUT ;ncmd ++){
		if(write_sp_ctl[ncmd] == YES){
			/* ensure that the sac_data file is nulled - we will free it later */
			/* get a temporary sacfile_ structure */
			for ( k=0 ; k < ntrc ; k ++){
				strcpy(t, sacdata[k].sac_ifile_name);
				switch (ncmd){
					case W_AM:  
						strcat(t,".am");
						break;
					case W_PH:  ;;
						strcat(t,".ph");
						break;
					case W_RL:  ;;
						strcat(t,".rl");
						break;
					case W_IM:  ;;
						strcat(t,".im");
						break;
				}
				n2 = sacdata[k].npow2/2 ;
				n21 = n2 + 1 ;
				if(temp == NULL){
					if((temp= (float *)calloc(n21,sizeof(float)) ) == NULL){
						fprintf(stderr,"could not allocate memory for spectra file\n");
						return;
					}
				} else {
					if((temp=(float *)realloc(temp,n21*sizeof(float)))==NULL){
						fprintf(stderr,"could not allocate memory for spectra file\n");
						return;
					}
				}
				for(i=0, j=0; i<=n2 ; i++){
					tr = sacdata[k].sac_spectra[j++];
					ti = sacdata[k].sac_spectra[j++];
					switch (ncmd){
						case W_AM:  
							temp[i] = (float)sqrt(tr*tr + ti*ti);
							break;
						case W_PH:  ;;
							temp[i] = atan2(ti,tr);;
							break;
						case W_RL:  ;;
							temp[i] = tr;
							break;
						case W_IM:  ;;
							temp[i] = ti;
							break;
					}
				}
				/* copy the header from sacdata[k] */
				for(i=0 ; i < 70 ; i++)
					tsac.sachdr.rhdr[i] = sacdata[k].sachdr.rhdr[i];
				for(i=0 ; i < 40 ; i++)
					tsac.sachdr.ihdr[i] = sacdata[k].sachdr.ihdr[i];
				for(i=0 ; i < 24 ; i++)
					strncpy(tsac.sachdr.chdr[i],sacdata[k].sachdr.chdr[i],9);
				tsac.sachdr.rhdr[H_B] =  0.0 ;
				tsac.sachdr.ihdr[H_NPTS] =  n21 ;
				tsac.sachdr.rhdr[H_E] =  0.5/sacdata[k].sachdr.rhdr[H_DELTA] ;
				tsac.sachdr.rhdr[H_DELTA] = sacdata[k].df ;
				/* unset time values */
				tsac.sachdr.rhdr[H_O] = -12345. ;
				tsac.sachdr.rhdr[H_A] = -12345. ;
				tsac.sachdr.rhdr[H_T0] = -12345. ;
				tsac.sachdr.rhdr[H_T1] = -12345. ;
				tsac.sachdr.rhdr[H_T2] = -12345. ;
				tsac.sachdr.rhdr[H_T3] = -12345. ;
				tsac.sachdr.rhdr[H_T4] = -12345. ;
				tsac.sachdr.rhdr[H_T5] = -12345. ;
				tsac.sachdr.rhdr[H_T6] = -12345. ;
				tsac.sachdr.rhdr[H_T7] = -12345. ;
				tsac.sachdr.rhdr[H_T8] = -12345. ;
				tsac.sachdr.rhdr[H_T9] = -12345. ;
				/* define as IFTYPE (15) = IXY = 4 */
				tsac.sachdr.ihdr[H_IFTYPE] =  ENUM_IXY ;
				/* define LEVEN as true */
				tsac.sachdr.ihdr[35] =  1 ;
				tsac.sac_data = temp;
				getmxmn(tsac.sac_data, n21,&depmax, &depmin, &depmen,&indmax,&indmin);
				tsac.sachdr.rhdr[H_DEPMAX] = depmax;
				tsac.sachdr.rhdr[H_DEPMIN] = depmin;
				tsac.sachdr.rhdr[H_DEPMEN] = depmen;
				tsac.sachdr.rhdr[H_TIMMAX] = tsac.sachdr.rhdr[H_B]  + ( indmax)*tsac.sachdr.rhdr[H_DELTA] ;
				tsac.sachdr.rhdr[H_TIMMIN] = tsac.sachdr.rhdr[H_B]  + ( indmin)*tsac.sachdr.rhdr[H_DELTA] ;
				bwsac(t,tsac.sachdr,tsac.sac_data);
		
			}
		}
	}
/*
	free(temp);
*/

}
