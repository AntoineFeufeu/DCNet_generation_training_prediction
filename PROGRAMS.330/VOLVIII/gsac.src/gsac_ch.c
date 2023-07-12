/* DEBUG later recompute dist az etc */
/* PUT IN SPECIAL FOR ch o gmt */
/* put in dummy header values */
/* CHANGES:
	26 JAN 2005 - problem of setting GCARC traced to the fact that
		there was no check for charg[i].id >= 0 in gsac_set_param
		This meant that I inadvertantly wrote over the rhdr 
	30 DEC 2017 - changed the lines from
			if(stla != -12345. && stlo != -12345. 
				&& evla != -12345. && evlo != -12345.){
                  to
			if(stla != -12345. && stlo != -12345. 
				&& evla != -12345. && evlo != -12345.
				&& sacdata[k].sachdr.ihdr[H_LCALDA] == 1 ){
		  to prevent DIST AZ BAZ being computed. This is required,
		  for example if the user wants to enter km for coordinates
		  instead of degrees
	30 OCT 2022 - permit ch NVHDR 6 and ch NHVDR 7 else no change
		
*/
#include	<stdio.h>
#include	<string.h>
#include	"gsac.h"
#include	"gsac_docommand.h"
#include	"gsac_sac.h"
#include	"gsac_sachdr.h"
#include	"gsac_arg.h"
#include	"csstim.h"

extern struct sacfile_ *sacdata;

/* variables */

#define NZYEAR 0
#define NZJDAY 1
#define NZHOUR 2
#define NZMIN  3
#define NZSEC  4
#define NZMSEC 5

/* define options for setting O A Tn - we can do both GMT and CAL 
 * notation. To make code cleaner in appearance, we note that
 * the array positions of O A Tn are 7 8 10...19 - you will see
 * a pattern here */
#define  OGMT -17
#define  AGMT -18
#define T0GMT -20
#define T1GMT -21
#define T2GMT -22
#define T3GMT -23
#define T4GMT -24
#define T5GMT -25
#define T6GMT -26
#define T7GMT -27
#define T8GMT -28
#define T9GMT -29
#define  OCAL -37
#define  ACAL -38
#define T0CAL -40
#define T1CAL -41
#define T2CAL -42
#define T3CAL -43
#define T4CAL -44
#define T5CAL -45
#define T6CAL -46
#define T7CAL -47
#define T8CAL -48
#define T9CAL -49
#define KEVNM 102

/* permit change to NVHDR need unique ID for this */
#define NVHDR -200




struct arghdr charg[] = {
	{ OGMT, "OGMT"  , IHDR, 0, 6, YES, "OGMT YY DOY HH MM SS MSEC",-1},
	{ AGMT, "AGMT"  , IHDR, 0, 6, YES, "AGMT YY DOY HH MM SS MSEC",-1},
	{T0GMT, "T0GMT" , IHDR, 0, 6, YES, "T0GMT YY DOY HH MM SS MSEC",-1},
	{T1GMT, "T1GMT" , IHDR, 0, 6, YES, "T1GMT YY DOY HH MM SS MSEC",-1},
	{T2GMT, "T2GMT" , IHDR, 0, 6, YES, "T2GMT YY DOY HH MM SS MSEC",-1},
	{T3GMT, "T3GMT" , IHDR, 0, 6, YES, "T3GMT YY DOY HH MM SS MSEC",-1},
	{T4GMT, "T4GMT" , IHDR, 0, 6, YES, "T4GMT YY DOY HH MM SS MSEC",-1},
	{T5GMT, "T5GMT" , IHDR, 0, 6, YES, "T5GMT YY DOY HH MM SS MSEC",-1},
	{T6GMT, "T6GMT" , IHDR, 0, 6, YES, "T6GMT YY DOY HH MM SS MSEC",-1},
	{T7GMT, "T7GMT" , IHDR, 0, 6, YES, "T7GMT YY DOY HH MM SS MSEC",-1},
	{T8GMT, "T8GMT" , IHDR, 0, 6, YES, "T8GMT YY DOY HH MM SS MSEC",-1},
	{T9GMT, "T9GMT" , IHDR, 0, 6, YES, "T9GMT YY DOY HH MM SS MSEC",-1},
	{ OCAL, "OCAL"  , IHDR, 0, 7, YES, "OCAL YY MM DD HH MM SS MSEC",-1},
	{ ACAL, "ACAL"  , IHDR, 0, 7, YES, "ACAL YY MM DD HH MM SS MSEC",-1},
	{T0CAL, "T0CAL" , IHDR, 0, 7, YES, "T0CAL YY MM DD HH MM SS MSEC",-1},
	{T1CAL, "T1CAL" , IHDR, 0, 7, YES, "T1CAL YY MM DD HH MM SS MSEC",-1},
	{T2CAL, "T2CAL" , IHDR, 0, 7, YES, "T2CAL YY MM DD HH MM SS MSEC",-1},
	{T3CAL, "T3CAL" , IHDR, 0, 7, YES, "T3CAL YY MM DD HH MM SS MSEC",-1},
	{T4CAL, "T4CAL" , IHDR, 0, 7, YES, "T4CAL YY MM DD HH MM SS MSEC",-1},
	{T5CAL, "T5CAL" , IHDR, 0, 7, YES, "T5CAL YY MM DD HH MM SS MSEC",-1},
	{T6CAL, "T6CAL" , IHDR, 0, 7, YES, "T6CAL YY MM DD HH MM SS MSEC",-1},
	{T7CAL, "T7CAL" , IHDR, 0, 7, YES, "T7CAL YY MM DD HH MM SS MSEC",-1},
	{T8CAL, "T8CAL" , IHDR, 0, 7, YES, "T8CAL YY MM DD HH MM SS MSEC",-1},
	{T9CAL, "T9CAL" , IHDR, 0, 7, YES, "T9CAL YY MM DD HH MM SS MSEC",-1},
	{  0, "DELTA"   , DHDR, 0, 1, YES, "DELTA dt" ,-1},
	{  1, "DEPMIN"  , RHDR, 0, 1, YES, "DEPMAX depmax" ,-1},
	{  2, "DEPMAX"  , RHDR, 0, 1, YES, "DEPMIN depmin" ,-1},
	{  3, "SCALE"   , RHDR, 0, 1, YES, "" ,-1},
	{  4, "ODELTA"  , RHDR, 0, 1, YES, "" ,-1},
	{  5, "B"       , DHDR, 0, 1, YES, "B btime" ,-1},
	{  6, "E"       , DHDR, 0, 1, YES, "E etime" ,-1},
	{  7, "O"       , DHDR, 0, 1, YES, "" ,-1},
	{  8, "A"       , DHDR, 0, 1, YES, "" ,-1},
	{  9, "FMT"     , RHDR, 0, 1, YES, "" ,-1},
	{ 10, "T0"      , DHDR, 0, 1, YES, "" ,-1},
	{ 11, "T1"      , DHDR, 0, 1, YES, "" ,-1},
	{ 12, "T2"      , DHDR, 0, 1, YES, "" ,-1},
	{ 13, "T3"      , DHDR, 0, 1, YES, "" ,-1},
	{ 14, "T4"      , DHDR, 0, 1, YES, "" ,-1},
	{ 15, "T5"      , DHDR, 0, 1, YES, "" ,-1},
	{ 16, "T6"      , DHDR, 0, 1, YES, "" ,-1},
	{ 17, "T7"      , DHDR, 0, 1, YES, "" ,-1},
	{ 18, "T8"      , DHDR, 0, 1, YES, "" ,-1},
	{ 19, "T9"      , DHDR, 0, 1, YES, "" ,-1},
	{ 20, "F"       , DHDR, 0, 1, YES, "" ,-1},
	{ 21, "RESP0"   , RHDR, 0, 1, YES, "" ,-1},
	{ 22, "RESP1"   , RHDR, 0, 1, YES, "" ,-1},
	{ 23, "RESP2"   , RHDR, 0, 1, YES, "" ,-1},
	{ 24, "RESP3"   , RHDR, 0, 1, YES, "" ,-1},
	{ 25, "RESP4"   , RHDR, 0, 1, YES, "" ,-1},
	{ 26, "RESP5"   , RHDR, 0, 1, YES, "" ,-1},
	{ 27, "RESP6"   , RHDR, 0, 1, YES, "" ,-1},
	{ 28, "RESP7"   , RHDR, 0, 1, YES, "" ,-1},
	{ 29, "RESP8"   , RHDR, 0, 1, YES, "" ,-1},
	{ 30, "RESP9"   , RHDR, 0, 1, YES, "" ,-1},
	{ 31, "STLA"    , DHDR, 0, 1, YES, "" ,-1},
	{ 32, "STLO"    , DHDR, 0, 1, YES, "" ,-1},
	{ 33, "STEL"    , RHDR, 0, 1, YES, "" ,-1},
	{ 34, "STDP"    , RHDR, 0, 1, YES, "" ,-1},
	{ 35, "EVLA"    , DHDR, 0, 1, YES, "" ,-1},
	{ 36, "EVLO"    , DHDR, 0, 1, YES, "" ,-1},
	{ 37, "EVEL"    , RHDR, 0, 1, YES, "" ,-1},
	{ 38, "EVDP"    , RHDR, 0, 1, YES, "" ,-1},
	{ 39, "MAG"     , RHDR, 0, 1, YES, "" ,-1},
	{ 40, "USER0"   , RHDR, 0, 1, YES, "" ,-1},
	{ 41, "USER1"   , RHDR, 0, 1, YES, "USER1 user1" ,-1},
	{ 42, "USER2"   , RHDR, 0, 1, YES, "" ,-1},
	{ 43, "USER3"   , RHDR, 0, 1, YES, "" ,-1},
	{ 44, "USER4"   , RHDR, 0, 1, YES, "" ,-1},
	{ 45, "USER5"   , RHDR, 0, 1, YES, "" ,-1},
	{ 46, "USER6"   , RHDR, 0, 1, YES, "" ,-1},
	{ 47, "USER7"   , RHDR, 0, 1, YES, "" ,-1},
	{ 48, "USER8"   , RHDR, 0, 1, YES, "" ,-1},
	{ 49, "USER9"   , RHDR, 0, 1, YES, "" ,-1},
	{ 50, "DIST"    , RHDR, 0, 1, YES, "" ,-1},
	{ 51, "AZ"      , RHDR, 0, 1, YES, "" ,-1},
	{ 52, "BAZ"     , RHDR, 0, 1, YES, "" ,-1},
	{ 53, "GCARC"   , RHDR, 0, 1, YES, "" ,-1},
	{ 54, "SB"      , DHDR, 0, 1, YES, "" ,-1},
	{ 55, "SDELTA"  , DHDR, 0, 1, YES, "" ,-1},
	{ 56, "DEPMEN"  , RHDR, 0, 1, YES, "DEPMEN depmen" ,-1},
	{ 57, "CMPAZ"   , RHDR, 0, 1, YES, "" ,-1},
	{ 58, "CMPINC"  , RHDR, 0, 1, YES, "" ,-1},
	{ 59, "XMINIMUM", RHDR, 0, 1, YES, "" ,-1},
	{ 60, "XMAXIMUM", RHDR, 0, 1, YES, "" ,-1},
	{ 61, "YMINIMUM", RHDR, 0, 1, YES, "" ,-1},
	{ 62, "YMAXIMUM", RHDR, 0, 1, YES, "" ,-1},
	{ 63, "ADJTM"   , RHDR, 0, 1, YES, "" ,-1},
	{ 64, "TIMMAX"  , RHDR, 0, 1, YES, "" ,-1},
	{ 65, "TIMMIN"  , RHDR, 0, 1, YES, "" ,-1},
	{ 66, "FHDR67"  , RHDR, 0, 1, YES, "" ,-1},
	{ 67, "FHDR68"  , RHDR, 0, 1, YES, "" ,-1},
	{ 68, "FHDR69"  , RHDR, 0, 1, YES, "" ,-1},
	{ 69, "FHDR70"  , RHDR, 0, 1, YES, "" ,-1},
	{  0, "NZYEAR"  , IHDR, 0, 1, YES, "NZYEAR nzyear" ,-1},
	{  1, "NZJDAY"  , IHDR, 0, 1, YES, "" ,-1},
	{  2, "NZHOUR"  , IHDR, 0, 1, YES, "" ,-1},
	{  3, "NZMIN"   , IHDR, 0, 1, YES, "" ,-1},
	{  4, "NZSEC"   , IHDR, 0, 1, YES, "" ,-1},
	{  5, "NZMSEC"  , IHDR, 0, 1, YES, "" ,-1},
	{  6, "NVHDR"   , IHDR, 0, 1, YES, "" ,-1}, /* must be 6 or 7 */
	{  7, "NINF"    , IHDR, 0, 1, YES, "" ,-1},
	{  8, "NHST"    , IHDR, 0, 1, YES, "" ,-1},
	{  9, "NPTS"    , IHDR, 0, 0, YES, "NPTS npts" ,-1}, /* never change*/
	{ 10, "NSNPTS"  , IHDR, 0, 1, YES, "" ,-1},
	{ 11, "NSN"     , IHDR, 0, 1, YES, "" ,-1},
	{ 12, "NXSIZE"  , IHDR, 0, 1, YES, "" ,-1},
	{ 13, "NYSIZE"  , IHDR, 0, 1, YES, "" ,-1},
	{ 14, "NHDR15"  , IHDR, 0, 1, YES, "" ,-1},
	{ 15, "IFTYPE"  , EHDR, 0, 1, YES, "" ,-1},
	{ 16, "IDEP"    , EHDR, 0, 1, YES, "" ,-1},
	{ 17, "IZTYPE"  , EHDR, 0, 1, YES, "" ,-1},
	{ 18, "IHDR4"   , EHDR, 0, 1, YES, "" ,-1},
	{ 19, "IINST"   , EHDR, 0, 1, YES, "" ,-1},
	{ 20, "ISTREG"  , IHDR, 0, 1, YES, "" ,-1},
	{ 21, "IEVREG"  , IHDR, 0, 1, YES, "" ,-1},
	{ 22, "IEVTYP"  , EHDR, 0, 1, YES, "" ,-1},
	{ 23, "IQUAL"   , EHDR, 0, 1, YES, "" ,-1},
	{ 24, "ISYNTH"  , EHDR, 0, 1, YES, "" ,-1},
	{ 25, "IHDR11"  , IHDR, 0, 1, YES, "" ,-1},
	{ 26, "IHDR12"  , IHDR, 0, 1, YES, "" ,-1},
	{ 27, "IHDR13"  , IHDR, 0, 1, YES, "" ,-1},
	{ 28, "IHDR14"  , IHDR, 0, 1, YES, "" ,-1},
	{ 29, "IHDR15"  , IHDR, 0, 1, YES, "" ,-1},
	{ 30, "IHDR16"  , IHDR, 0, 1, YES, "" ,-1},
	{ 31, "IHDR17"  , IHDR, 0, 1, YES, "" ,-1},
	{ 32, "IHDR18"  , IHDR, 0, 1, YES, "" ,-1},
	{ 33, "IHDR19"  , IHDR, 0, 1, YES, "" ,-1},
	{ 34, "IHDR20"  , IHDR, 0, 1, YES, "" ,-1},
	{ 36, "LPSPOL"  , LHDR, 0, 1, YES, "" ,-1},
	{ 37, "LOVROK"  , LHDR, 0, 1, YES, "LOVROK [T|F]" ,-1},
	{ 38, "LCALDA"  , LHDR, 0, 1, YES, "LCALDA [T|F]" ,-1},
	{ 39, "LHDR5"   , LHDR, 0, 1, YES, "" ,-1},
	{  0, "KSTNM"   , CHDR, 0, 1, YES, "" ,-1},
	{  2, "KEVNM"  , CHDL, 0, 1, YES, "" ,-1},
	{  3, "KHOLE"   , CHDR, 0, 1, YES, "" ,-1},
	{  4, "KO"      , CHDR, 0, 1, YES, "" ,-1},
	{  5, "KA"      , CHDR, 0, 1, YES, "" ,-1},
	{  6, "KT0"     , CHDR, 0, 1, YES, "" ,-1},
	{  7, "KT1"     , CHDR, 0, 1, YES, "" ,-1},
	{  8, "KT2"     , CHDR, 0, 1, YES, "" ,-1},
	{  9, "KT3"     , CHDR, 0, 1, YES, "" ,-1},
	{ 10, "KT4"     , CHDR, 0, 1, YES, "" ,-1},
	{ 11, "KT5"     , CHDR, 0, 1, YES, "" ,-1},
	{ 12, "KT6"     , CHDR, 0, 1, YES, "" ,-1},
	{ 13, "KT7"     , CHDR, 0, 1, YES, "" ,-1},
	{ 14, "KT8"     , CHDR, 0, 1, YES, "" ,-1},
	{ 15, "KT9"     , CHDR, 0, 1, YES, "" ,-1},
	{ 16, "KF"      , CHDR, 0, 1, YES, "" ,-1},
	{ 17, "KUSER0"  , CHDR, 0, 1, YES, "" ,-1},
	{ 18, "KUSER1"  , CHDR, 0, 1, YES, "" ,-1},
	{ 19, "KUSER2"  , CHDR, 0, 1, YES, "" ,-1},
	{ 20, "KCMPNM"  , CHDR, 0, 1, YES, "" ,-1},
	{ 21, "KNETWK"  , CHDR, 0, 1, YES, "" ,-1},
	{ 22, "KDATRD"  , CHDR, 0, 1, YES, "" ,-1},
	{ 23, "KINST"   , CHDR, 0, 1, YES, "" ,-1},
	{-100, ""        , CHDR, 0, 1, YES, "" ,-1}
};	

static float ch_real[10];
static double ch_double[10];
static int   ch_int [10];
static int   ch_logic [10];
static int chcoord;
static int choat;
static int chdatetime;

/* array for time changes */
/* temporary array for epoch times for new O A Tn */
#define NOAT 20
static double oatepoch[NOAT];
static int    oateind[NOAT];
static char chinstr[100];

void gsac_set_param_ch(int ncmd, char **cmdstr)
{
	/* code to change header values */
	int i, k, m, ntrchdr;
	double tztmp;
	char str1[9], str2[9];
	int numgmt, numcal;
	chcoord = NO;
	chdatetime = NO;
	choat = NO;
	for(i=0 ; i < NOAT; i++)
		oateind[i] = NO;
	if(ncmd == 1)
		return;
	numgmt = gsac_countgmt(ncmd, cmdstr,"GMT") ;
	if(numgmt > 0){
		for(i=0;i<numgmt;i++){
			gsac_mergegmt(&ncmd,cmdstr,"GMT");
		}
	}
	numcal = gsac_countgmt(ncmd, cmdstr,"CAL") ;
	if(numcal > 0){
		for(i=0;i<numcal;i++){
			gsac_mergegmt(&ncmd,cmdstr,"CAL");
		}
	}
	if(numgmt > 0 || numcal > 0){
		printf("Executing: ");
		for(i=0;i<ncmd;i++)
			printf("%s ",cmdstr[i]);
		printf("\n");
	}
	/* is the command syntax correct ? */
	if(testarg(ncmd, cmdstr, charg, NO, YES))
		return;

	/* NOTE THIS IS ONE TIME THAT THE EXEC IS DONE IN THE SET PARAM
	 * THIS IS BECAUSE I NEED THE ncmd and cmdstr
	 * */
	/* loop through traces and then look through header values */
	/* for speed loop through options and then set traces */
	ntrchdr = gsac_control.number_iheaders;
	if ( ntrchdr < 1)
		return;
	/* check for presence of a sole GMT - if present give the
	 * number of occurrences */

	for(i=0; charg[i].key[0] != '\0'  ;i++){
		if(charg[i].used >= 1 ){
			if(charg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, charg[i].key, charg[i].mfit, charg[i].narg, ch_real);
				/* special logic to handle coordinate change 
				 * which then requires a recompute of distances
				 */
				if(charg[i].id == 31 || charg[i].id == 32 ||
					charg[i].id == 35 || charg[i].id == 36){
					chcoord = YES;
				}
			} else if(charg[i].ricell == DHDR){
				getargd(ncmd, cmdstr, charg[i].key, charg[i].mfit, charg[i].narg, ch_double);
				/* special logic to handle coordinate change 
				 * which then requires a recompute of distances
				 */
				if(charg[i].id == 31 || charg[i].id == 32 ||
					charg[i].id == 35 || charg[i].id == 36){
					chcoord = YES;
				}
			} else if(charg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, charg[i].key, charg[i].mfit, charg[i].narg, ch_int );
				if(charg[i].id <= -37 && charg[i].id >= -49){
					/* CAL conversion */
					htoe2(ch_int[0], ch_int[1],ch_int[2],ch_int[3],
						ch_int[4],ch_int[5],ch_int[6],&tztmp);
					m = -30 - charg[i].id ;
					oateind[m] = YES;
					oatepoch[m] = tztmp;
					choat = YES;
				}
				if(charg[i].id <= -17 && charg[i].id >= -29){
					/* GMT conversion */
					htoe1(ch_int[0], ch_int[1],ch_int[2],ch_int[3],
						ch_int[4],ch_int[5],&tztmp);
					m = -10 - charg[i].id ;
					oateind[m] = YES;
					oatepoch[m] = tztmp;
					choat = YES;
				}
				if(charg[i].id >= 0 && charg[i].id <=5){
					chdatetime =YES;
				}
				if(charg[i].id == NVHDR ){
					/* check value for nvhdr - only 6 or 7 permitted */
					if(ch_int[0] >= 6 && ch_int[0] <= 7 ){
						fprintf(stderr,"Changing NVHDR to %d\n",ch_int[0]);
					} else {
						fprintf(stderr,"ch NVHDR %d has incorrect value. Must be either 6 (original) or 7 (extended)\n",ch_int[0]);
					}
				}

			} else if(charg[i].ricell == LHDR){
				getargl(ncmd, cmdstr, charg[i].key, charg[i].mfit, charg[i].narg, ch_logic );
			} else if(charg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, charg[i].key, charg[i].mfit, charg[i].narg, chinstr );
				fillstr1(chinstr,str1);
			} else if(charg[i].ricell == CHDL){
				getargs(ncmd, cmdstr, charg[i].key, charg[i].mfit, charg[i].narg, chinstr );
				fillstr12(chinstr,str1,str2);
			}
			if(charg[i].id >= 0) {
				for ( k=0 ; k < ntrchdr ; k ++){
					if(charg[i].ricell == RHDR){
						sacdata[k].sachdr.rhdr[charg[i].id] = ch_real[0];
						/* note interannly sahdr.rhdr is doublesothis was just  cst */
					} if(charg[i].ricell == DHDR){
						sacdata[k].sachdr.rhdr[charg[i].id] = ch_double[0];
					} else if(charg[i].ricell == IHDR){
						sacdata[k].sachdr.ihdr[charg[i].id] = ch_int[0];
					} else if(charg[i].ricell == LHDR){
						sacdata[k].sachdr.ihdr[charg[i].id] = ch_logic[0];
					} else if(charg[i].ricell == CHDR){
						/* update in header for WH and also display */
						strncpy(sacdata[k].sachdr.chdr[charg[i].id],str1,8); 
						strcpy(sacdata[k].schdr[charg[i].id],str1); 
					} else if(charg[i].ricell == CHDL){
						/* update in header for WH and also display */
						strncpy(sacdata[k].sachdr.chdr[1],str1,8); 
						strncpy(sacdata[k].sachdr.chdr[2],str2,8); 
						strcpy(sacdata[k].schdr[1],str1); 
						strcpy(sacdata[k].schdr[2],str2); 
					}
				}
			}
		}
	}

}

void gsac_exec_ch(void)
{
	double evla, evlo, stla, stlo;
	
	double gcarc, az, baz, dist;
	int i, k;
	int ntrchdr;
	int month, day;
	/* if there are no traces teturn */
	ntrchdr = gsac_control.number_iheaders;
	if(ntrchdr < 1)
		return;

	/* special updating 
	 * if station or source coordinates are changed recompute distance 
	 * 	and azimuth
	 * if any of the arrival times are set, such as the 
	 * OGMT, AGMT, TnGMT, OCAL, ACAL, TnCAL then recompute the time values
	 * */

	/* check for coordinate change */
	if(chcoord == YES ){
		for ( k=0 ; k < ntrchdr ; k ++){
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
		}
	}
	if(chdatetime == YES ){
		/* put in the date -time change here 
		 * we first computer the correct epoch and then convert
		 * to human time - this will take care fo thinks like
		 * ch nzhour 25 */
		for ( k=0 ; k < ntrchdr ; k ++){
			htoe1(sacdata[k].sachdr.ihdr[NZYEAR], 
				sacdata[k].sachdr.ihdr[NZJDAY], 
				sacdata[k].sachdr.ihdr[NZHOUR],
				sacdata[k].sachdr.ihdr[NZMIN],
				sacdata[k].sachdr.ihdr[NZSEC],
				sacdata[k].sachdr.ihdr[NZMSEC],
				&sacdata[k].tzref);
		/* NOW RESET THE OTHER TIME FIELDS */
		etoh(sacdata[k].tzref, &sacdata[k].sachdr.ihdr[NZYEAR], 
				&sacdata[k].sachdr.ihdr[NZJDAY], &month, &day,
				&sacdata[k].sachdr.ihdr[NZHOUR], 
				&sacdata[k].sachdr.ihdr[NZMIN],
				&sacdata[k].sachdr.ihdr[NZSEC], 
				&sacdata[k].sachdr.ihdr[NZMSEC]);
		}
	}
	if(choat ==YES ){
		/* change the time of the O A Tn fields which
		 * are relative to the reference time */
		for ( k=0 ; k < ntrchdr ; k ++){
			for(i=0; i < NOAT; i++){
				if(oateind[i] == YES)
					sacdata[k].sachdr.rhdr[i] 
						= oatepoch[i]-sacdata[k].tzref;
			}
		}
	}

}
