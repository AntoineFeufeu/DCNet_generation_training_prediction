/* changes
	19 NOV 2019 added TIMMIN TIMMAX  
	24 FEB 2020 corrected logic so that
		GSAC> fileid fname   or
		GSAC> fileid dist
	    will not work. the "list" is required, e.g.,
		GSAC> fileid list
	07 AUG 2020  added DEPMAX DEPMIN to list
	16 SEP 2020  Fixed CONCAT ON|OFF and FORMAT COLON to require 3 character match
                     for CONcat  to avoid confusion with COlon
        04 JAN 2021  replaced strncat at line 650 by strcat
	20 SEP 2021  implemented a deblack to eliminate a lot of blanks in the annotation
	25 OCT 2021  correct to recognize depmin depmax - problem was the DEfault required
                     a two character match which conflicted with DEPMAX DEPMIN -
		     now DEFault requires a 3 character match
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"
#include	"csstim.h"
#include	<libgen.h>

extern struct sacfile_ *sacdata;
extern int *sortptr;
extern void gnomesortint(int n, int ar[], int key[] ) ;


#define	FILEID_DEFAULT	-1
#define	FILEID_ON	-2
#define	FILEID_OFF	-3
#define	FILEID_LOCATION	-4
#define	FILEID_FORMAT	-5
#define	FILEID_TYPE	-6
#define	FILEID_NAME	-7
#define	FILEID_LIST	-8
#define FILEID_DUMMY	-9
#define FILEID_UR	-10
#define FILEID_UL	-11
#define FILEID_LR	-12
#define FILEID_LL	-13
#define FILEID_UC	-14
#define FILEID_LC	-15
#define FILEID_EQUALS	-16
#define FILEID_COLONS	-17
#define FILEID_NONAMES	-18
#define FILEID_CONCAT   -19
#define KZDATE  100
#define KZTIME  101
#define KEVNM   102
#define FILEID_FNAME   103
#define FILEID_BNAME   104

void gsac_plot_fileid(float x0,float y0,float xlen,float ylen, int k, int kmax);
void gsac_set_param_fileid(int ncmd, char **cmdstr);
void gsac_exec_fileid(void);
void gsac_deblank(char *source, char *result);


struct arghdr fileidarg[] = {
	{FILEID_DEFAULT ,"DEFAULT"  , IHDR, NO , 0, NO , "", 3},
	{FILEID_ON      , "ON"      , IHDR, NO , 0, NO , "",2},
	{FILEID_OFF     , "OFF"     , IHDR, NO , 0, NO , "",2},
	{FILEID_LIST    , "LIST"    , IHDR, NO , 0, NO , "",2},
	{FILEID_FORMAT  , "FORMAT"  , CHDR, NO , 1, NO , "FOrmat EQUAls|COlons|NOnames",2},
	{FILEID_TYPE    , "TYPE"    , IHDR, NO , 0, NO , "",2},
	{FILEID_NAME    , "NAME"    , IHDR, NO , 0, NO , "", 1},
	{FILEID_LOCATION,"LOCATION" , CHDR, NO , 1, NO , "LOcation UL|UC|UR|LL|LC|LR", 2},
	{FILEID_CONCAT  , "CONCAT"  , YHDR, NO , 1, NO , "CONcat [ON|OFF] ", 3},
	{ KZDATE        , "KZDATE"  , IHDR, YES, 0, YES, "" ,-1},
	{ KZTIME        , "KZTIME"  , IHDR, YES, 0, YES, "" ,-1},
	{ KEVNM         , "KEVNM"   , IHDR, YES, 0, YES, "" ,-1},
	{ FILEID_FNAME  , "FNAME"   , IHDR, YES, 0, YES, "" ,-1},
	{ FILEID_BNAME  , "BNAME"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_STLA        , "STLA"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_STLO        , "STLO"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_STEL        , "STEL"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_EVLA        , "EVLA"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_EVLO        , "EVLO"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_EVDP        , "EVDP"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER0       , "USER0"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER1       , "USER1"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER2       , "USER2"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER3       , "USER3"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER4       , "USER4"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER5       , "USER5"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER6       , "USER6"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER7       , "USER7"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER8       , "USER8"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_USER9       , "USER9"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_TIMMIN      , "TIMMIN"  , IHDR, YES, 0, YES, "" ,-1},
	{ H_TIMMAX      , "TIMMAX"  , IHDR, YES, 0, YES, "" ,-1},
	{ H_DIST        , "DIST"    , IHDR, YES, 0, YES, "" ,-1},
	{ H_AZ          , "AZ"      , IHDR, YES, 0, YES, "" ,-1},
	{ H_BAZ         , "BAZ"     , IHDR, YES, 0, YES, "" ,-1},
	{ H_GCARC       , "GCARC"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_KSTNM       , "KSTNM"   , IHDR, YES, 0, YES, "" ,-1},
	{ H_KCMPNM      , "KCMPNM"  , IHDR, YES, 0, YES, "" ,-1},
	{ H_DEPMAX      , "DEPMAX"  , IHDR, YES, 0, YES, "" ,-1},
	{ H_DEPMIN      , "DEPMIN"  , IHDR, YES, 0, YES, "" ,-1},
	{0, ""  , IHDR, NO, 0, NO, "", -1}
};

/* these are temporary variables only used here */
float fileid_real[10];
int   fileid_int [10];
int   fileid_yn;
int   fileid_num;

/* these are prototypes for global variables to be used by the routine */
#define FILEID_MAXLIST 10
int fileid_use = YES;
int fileid_location = FILEID_UR ;
int fileid_format = FILEID_NONAMES;
int fileid_type = FILEID_LIST;
int fileid_n = 2;
int fileid_concat = NO;
int fileid_list[FILEID_MAXLIST] =
	{ 0, 20, -1, -1, -1, -1, -1, -1, -1, -1
	};	/* this is the actual list to be displayed */
char fileid_clist[FILEID_MAXLIST][8] ={
	"KSTNM", "KCMPNM", "", "", "", "", "", "", "", ""
	};	/* this is the list of title strings to be displayed */
int fileid_ar[FILEID_MAXLIST] = 
	{ 0, 1, -1, -1, -1, -1, -1, -1, -1, -1
	};	/* this is used for the sort to ensure the list is in order given */
int fileid_key[FILEID_MAXLIST] =
	{ 0, 1, -1, -1, -1, -1, -1, -1, -1, -1
	}; /* this created by gnomesortint */
#define MAXCONCATSTR 301
char fileid_concatstr[MAXCONCATSTR];
char fileid_tconcatstr[MAXCONCATSTR];
#define MAXSTR 40
char fileid_ostr[FILEID_MAXLIST][MAXSTR];

/* temporary space */
	

void gsac_set_param_fileid(int ncmd, char **cmdstr)
{
	int i;
	char instr[100];
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, fileidarg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; fileidarg[i].key[0] != '\0' ; i++){
		if(fileidarg[i].used > 0){
			if(fileidarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, fileid_real);
			} else if(fileidarg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, fileid_int );
			} else if(fileidarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, &fileid_yn );
			} else if(fileidarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit,fileidarg[i].narg, &fileid_num );
			} else if(fileidarg[i].ricell == CHDR){
				getargs(ncmd, cmdstr, fileidarg[i].key, 
					fileidarg[i].mfit, fileidarg[i].narg, instr );
			}
			switch(fileidarg[i].id){
				case FILEID_LIST:
					fileid_n = 0;
					fileid_type = FILEID_LIST;
				case FILEID_ON:
					fileid_use = YES;
					break;
				case FILEID_OFF:
					fileid_use = NO;
					break;
				case FILEID_DEFAULT:
					fileid_use = YES;
					fileid_location = FILEID_UR ;
					fileid_format = FILEID_NONAMES;
					fileid_type = FILEID_DEFAULT;
					fileid_n = 4;
					fileid_list[0]=H_KSTNM;
					strcpy(fileid_clist[0],"KSTNM");
					fileid_list[1]=H_KCMPNM;
					strcpy(fileid_clist[1],"KCMPNM");
					fileid_list[2]=KZDATE;
					strcpy(fileid_clist[2],"KZDATE");
					fileid_list[3]=KZTIME;
					strcpy(fileid_clist[3],"KZTIME");
					break;
				case FILEID_FORMAT:
					gsac_strupr(instr);
					if(strncmp(instr,"COL",3)==0)
						fileid_format = FILEID_COLONS;
					else if(strncmp(instr,"EQ",2)==0)
						fileid_format = FILEID_EQUALS;
					else if(strncmp(instr,"NO",2)==0)
						fileid_format = FILEID_NONAMES;
					break;
				case FILEID_LOCATION:
					gsac_strupr(instr);
					printf("fileid_location %s\n",instr);
					if(strncmp(instr,"UL",2)==0)
						fileid_location = FILEID_UL;
					else if(strncmp(instr,"UR",2)==0)
						fileid_location = FILEID_UR;
					else if(strncmp(instr,"UC",2)==0)
						fileid_location = FILEID_UC;
					else if(strncmp(instr,"LL",2)==0)
						fileid_location = FILEID_LL;
					else if(strncmp(instr,"LR",2)==0)
						fileid_location = FILEID_LR;
					else if(strncmp(instr,"LC",2)==0)
						fileid_location = FILEID_LC;
					break;
                                case FILEID_CONCAT:
                                        if(fileid_yn == NO)
                                                fileid_concat = NO;
                                        else if(fileid_yn == YES)
                                                fileid_concat = YES;
                                        break;
				case FILEID_NAME:
					fileid_type = FILEID_NAME;
					break;
				case KZDATE:
				case KZTIME:
				case H_KCMPNM:
				case H_KSTNM:
				case KEVNM:
				case FILEID_FNAME:
				case FILEID_BNAME:
				case H_GCARC:
				case H_DIST:
				case H_AZ:
				case H_BAZ:
				case H_STLA:
				case H_STLO:
				case H_STEL:
				case H_EVLA:
				case H_EVLO:
				case H_EVDP:
				case H_USER0:
				case H_USER1:
				case H_USER2:
				case H_USER3:
				case H_USER4:
				case H_USER5:
				case H_USER6:
				case H_USER7:
				case H_USER8:
				case H_USER9:
				case H_TIMMIN:
				case H_TIMMAX:
				case H_DEPMAX:
				case H_DEPMIN:
					if(fileid_type == FILEID_LIST) {
					/* safety so never exceed array dimension  */
					if(fileid_n < FILEID_MAXLIST  ){
					fileid_list[fileid_n]=fileidarg[i].id;;
					/* now annotate the labeling string */
					switch(fileidarg[i].id){
						case KZDATE:
						strcpy(fileid_clist[fileid_n],"KZDATE");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KZDATE");
						break;
						case KZTIME:
						strcpy(fileid_clist[fileid_n],"KZTIME");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KZTIME");
						break;
						case H_KCMPNM:
						strcpy(fileid_clist[fileid_n],"KCMPNM");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KCMPNM");
						break;
						case H_KSTNM:
						strcpy(fileid_clist[fileid_n],"KSTNM");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KSTNM");
						break;
						case KEVNM:
						strcpy(fileid_clist[fileid_n],"KEVNM");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"KEVNM");
						break;
						case FILEID_FNAME:
						strcpy(fileid_clist[fileid_n],"FNAME");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"FNAME");
						break;
						case FILEID_BNAME:
						strcpy(fileid_clist[fileid_n],"BNAME");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"BNAME");
						break;
						case H_GCARC:
						strcpy(fileid_clist[fileid_n],"GCARC");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"GCARC");
						break;
						case H_DIST:
						strcpy(fileid_clist[fileid_n],"DIST");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"DIST");
						break;
						case H_AZ:
						strcpy(fileid_clist[fileid_n],"AZ");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"AZ");
						break;
						case H_BAZ:
						strcpy(fileid_clist[fileid_n],"BAZ");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"BAZ");
						break;
						case H_STLA:
						strcpy(fileid_clist[fileid_n],"STLA");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"STLA");
						break;
						case H_STLO:
						strcpy(fileid_clist[fileid_n],"STLO");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"STLO");
						break;
						case H_STEL:
						strcpy(fileid_clist[fileid_n],"STEL");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"STEL");
						break;
						case H_EVLA:
						strcpy(fileid_clist[fileid_n],"EVLA");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"EVLA");
						break;
						case H_EVLO:
						strcpy(fileid_clist[fileid_n],"EVLO");
						fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"EVLO");
						break;
						case H_EVDP:
							strcpy(fileid_clist[fileid_n],"EVDP");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"EVDP");
						break;
						case H_USER0:
							strcpy(fileid_clist[fileid_n],"USER0");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER0");
						break;
						case H_USER1:
							strcpy(fileid_clist[fileid_n],"USER1");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER1");
						break;
						case H_USER2:
							strcpy(fileid_clist[fileid_n],"USER2");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER2");
						break;
						case H_USER3:
							strcpy(fileid_clist[fileid_n],"USER3");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER3");
						break;
						case H_USER4:
							strcpy(fileid_clist[fileid_n],"USER4");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER4");
						break;
						case H_USER5:
							strcpy(fileid_clist[fileid_n],"USER5");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER5");
						break;
						case H_USER6:
							strcpy(fileid_clist[fileid_n],"USER6");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER6");
						break;
						case H_USER7:
							strcpy(fileid_clist[fileid_n],"USER7");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER7");
						break;
						case H_USER8:
							strcpy(fileid_clist[fileid_n],"USER8");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER8");
						break;
						case H_USER9:
							strcpy(fileid_clist[fileid_n],"USER9");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"USER9");
						break;
						case H_DEPMIN:
							strcpy(fileid_clist[fileid_n],"DEPMIN");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"DEPMIN");
						break;
						case H_DEPMAX:
							strcpy(fileid_clist[fileid_n],"DEPMAX");
							fileid_ar[fileid_n]=
							whicharg(ncmd,cmdstr,"DEPMAX");
						break;
					}
					fileid_key[fileid_n] = fileid_n;
					fileid_n++;
					}
					break;
				}
			}
		}
	}
}

void gsac_exec_fileid(void)
{
}

void gsac_plot_fileid(float x0,float y0,float xlen,float ylen, int k, int kmax)
{
	/* implement the fileid labeling 
	the bounding box is given by (x0,y0) for the lower left 
		and (x0+xlen,y0+ylen) for the upper right corner
	k is the index of the trace file - note we need a pointer to the
		sacdata structure at the top of this routine
	*/
	int i,j;
	float ht;
	float xp, yp;
	float dy;
	int LCR;	/* Left < 0 , Center = 0 , Right > 0 */
	int clength;	/* maximum length of character string for placement */
	float cw;	/* width of character string*/
	char outstr[MAXSTR];
	int ls;
	if(fileid_use == YES){
		/* change order of output of used LIST */
		if(fileid_type == FILEID_LIST){
			if(fileid_concat == NO){
				/* sort the list */
				gnomesortint(fileid_n, fileid_ar, fileid_key);
				ht = MIN(0.4*ylen/fileid_n,0.1);
			} else {
			ht = MIN(0.10*ylen,0.1);
			}
		} else {
			ht = MIN(0.10*ylen,0.1);
		}
		switch(fileid_location){
			case FILEID_UR:
				xp = x0 + xlen ;
				dy = -1.25*ht;
				yp = y0 + ylen + 2.5*dy;
				LCR = 1;
				break;
			case FILEID_UC:
				xp = x0 + 0.5*xlen;
				dy = -1.25*ht;
				yp = y0 + ylen + 2.5*dy;
				LCR = 0;
				break;
			case FILEID_UL:
				xp = x0 + 2.*ht;
				dy = -1.25*ht;
				yp = y0 + ylen + 2.5*dy;
				LCR = -1;
				break;
			case FILEID_LR:
				xp = x0 + xlen ;
				dy = 1.25*ht;
				yp =  y0 + dy + kmax*dy;
				LCR = 1;
				break;
			case FILEID_LC:
				xp = x0 + 0.5*xlen;
				dy = 1.25*ht;
				yp =  y0 + dy + kmax*dy;
				LCR = 0;
				break;
			case FILEID_LL:
				xp = x0 + 2.*ht;
				dy = 1.25*ht;
				yp =  y0 + dy + kmax*dy;
				LCR = -1;
				break;
		}
		if(fileid_type == FILEID_NAME){
			clength = strlen(sacdata[k].sac_ofile_name);
			cw = clength*ht;
			if(LCR < 0)
				gleft(xp, yp,ht,sacdata[k].sac_ofile_name,0.0);
			else if(LCR == 0)
				gleft(xp-0.5*cw,yp,ht,sacdata[k].sac_ofile_name,0.0);
			else if(LCR > 0)
				gleft(xp -cw-ht-ht,yp,ht,sacdata[k].sac_ofile_name,0.0);
		} else if(fileid_type == FILEID_LIST){
			/* this is complicated by the size of the fields, the
				fact that the fields are numeric or string
				and the formatting - so be patient 
				for SAC compatability the field width depends on the
				numbers - yuck */
			clength = 0;
			for(j = 0 ; j < fileid_n; j++){
				i = fileid_key[j];
				switch(fileid_list[i]){
				case KZDATE:
					printkdatestr(sacdata[k].sachdr.ihdr[0],
                                                        sacdata[k].sachdr.ihdr[1],outstr);
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-18s",outstr);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%8s:%-18s",fileid_clist[i],outstr);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-18s",fileid_clist[i],outstr);
							break;
					}
					break;
				case KZTIME:
					printktimestr(sacdata[k].sachdr.ihdr[2],
                                                        sacdata[k].sachdr.ihdr[3],
                                                        sacdata[k].sachdr.ihdr[4],
                                                        sacdata[k].sachdr.ihdr[5],outstr);
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-12s",outstr);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%8s:%-12s",fileid_clist[i],outstr);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-12s",fileid_clist[i],outstr);
							break;
					}
					break;
				case H_KCMPNM:
				case H_KSTNM:
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-8s",sacdata[k].schdr[fileid_list[i]]);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%8s:%-8s",fileid_clist[i],sacdata[k].schdr[fileid_list[i]]);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-8s",fileid_clist[i],sacdata[k].schdr[fileid_list[i]]);
							break;
					}
					break;
				case KEVNM:
					sprintf(fileid_ostr[i],"%-8s%-8s",
						sacdata[k].schdr[1],
						sacdata[k].schdr[2]);
					break;
				case FILEID_FNAME:
					/* be careful here since the fileid_ostr is
						only 40 characters long */
					ls = strlen(sacdata[k].sac_ofile_name);
					if(ls >= MAXSTR)
						ls = MAXSTR-1 ;
					outstr[0] = '\0';
					strncat(outstr,sacdata[k].sac_ofile_name,ls);
					sprintf(fileid_ostr[i],"%-s",outstr);
					break;
				case FILEID_BNAME:
					/* be careful here since the fileid_ostr is
						only 40 characters long */
					ls = strlen(basename(sacdata[k].sac_ofile_name));
					if(ls >= MAXSTR)
						ls = MAXSTR-1 ;
					outstr[0] = '\0';
					strncat(outstr,basename(sacdata[k].sac_ofile_name),ls);
					sprintf(fileid_ostr[i],"%-s",outstr);
					break;
				case H_GCARC:
				case H_DIST:
				case H_AZ:
				case H_BAZ:
				case H_STLA:
				case H_STLO:
				case H_STEL:
				case H_EVLA:
				case H_EVLO:
				case H_EVDP:
				case H_USER0:
				case H_USER1:
				case H_USER2:
				case H_USER3:
				case H_USER4:
				case H_USER5:
				case H_USER6:
				case H_USER7:
				case H_USER8:
				case H_USER9:
				case H_DEPMAX:
				case H_DEPMIN:
					switch (fileid_format){
						case FILEID_NONAMES:
							sprintf(fileid_ostr[i],"%-9.3g",sacdata[k].sachdr.rhdr[fileid_list[i]]);
							break;
						case FILEID_COLONS:
							sprintf(fileid_ostr[i],"%8s:%-9.3g",fileid_clist[i],sacdata[k].sachdr.rhdr[fileid_list[i]]);
							break;
						case FILEID_EQUALS:
							sprintf(fileid_ostr[i],"%8s=%-9.3g",fileid_clist[i],sacdata[k].sachdr.rhdr[fileid_list[i]]);
							break;
					}
					break;
				}
					/* safety so never exceed array dimension  */
				clength=MAX(clength,strlen(fileid_ostr[i]));
			}
			if(fileid_concat == NO){
				/* the output strings are placed vertically */
				/* now that we have the output strings and the length we
					can properly place them on the plot */
				cw = clength*ht;
				for(j = 0 ; j < fileid_n; j++){
					i = fileid_key[j];
					if(LCR < 0)
						gleft(xp, yp,ht,fileid_ostr[i],0.0);
					else if(LCR == 0)
						gleft(xp-0.5*cw,yp,ht,fileid_ostr[i],0.0);
					else if(LCR > 0)
						gleft(xp -cw-ht-ht,yp,ht,fileid_ostr[i],0.0);
					yp += dy;
				}
			} else {
				/* the output strings are plotted horizontally */
				fileid_concatstr[0] = '\0';
				for(j = 0 ; j < fileid_n; j++){
					i = fileid_key[j];
					ls = strlen(fileid_ostr[i] );
					/* put in test for the 300 */
					if(strlen(fileid_concatstr) + ls + 1 < MAXCONCATSTR  ){
						strncat(fileid_concatstr, fileid_ostr[i],ls);
						/* add a blank */
						strcat(fileid_concatstr, " " );
					}

				}
				ls = strlen(fileid_concatstr);
				gsac_deblank(fileid_concatstr,fileid_tconcatstr);
				if(LCR < 0)
					gleft (xp + ht,   yp,ht,fileid_tconcatstr,0.0);
				else if(LCR == 0)
					gcent (xp        ,yp,ht,fileid_tconcatstr,0.0);
				else if(LCR > 0)
					gright(xp -ht -ht,yp,ht,fileid_tconcatstr,0.0);
		
			}
		}

	}

}

void gsac_deblank(char *source, char *result)
{
	int haveblank;
	haveblank = NO;
	while(*source){
		if(isblank(*source)){
			*result = *source; 
			source++;
			result++;
			while(isblank(*source))
				source++ ;
		} else {
			*result = *source; 
			source++;
			result++;
		}
	}
        *result = '\0';
}

