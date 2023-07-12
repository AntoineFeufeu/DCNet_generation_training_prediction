/* This is the C version of the older shwsac.f 
   It differs in that it supports the extended SAC header filed NVHDR=7
   in addition to the older NVHRD = 6
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>
#include "sacsubc.h"
#include "calplot.h"

/* parameters */
#define TRUE 1
#define FALSE 0
#define LLN 10000

/* prototypes */
void plotit(float *x, int npts, float rhdr1, float rhdr2, float rhdr3);
void usage(void);
void  gcmdln(int argc, char**argv,int *sacbin, int *onlyhead, char *fname);
void parsesac(int LN, float *x, int onlyhead);

/* parameter from SACIO */
extern struct sachdr_ sachdr;
extern struct dsachdr_ dsachdr;

/* external values defines in sacsubc.h */
extern char *rstr[];
extern char *dstr[] ;
extern char *istr[] ;
extern char *cstr[] ;
extern char *estr[] ;
extern char *Istr[] ;
extern int NRSTR, NISTR, NCSTR, NDSTR;


int main(int argc, char **argv)
{
	int sacbin, nerr;
	int onlyhead;
	char fname[81];
float *x;
	/*
	parse command line arguments
	*/
        gcmdln(argc,argv,&sacbin,&onlyhead,fname);
	if(strlen(fname) == 0){
		fprintf(stderr,"Input file not defined\n");
        	return(EXIT_SUCCESS);
	}
	/*
	Read input waveform data in SAC binary format.
	*/
	if(onlyhead == TRUE){
		if(sacbin == TRUE){
                	brsach(fname,&nerr);
		} else {
                	arsach(fname,&nerr);
		}
	} else {
		if(sacbin == TRUE){
                	brsac(LLN,fname,&x,&nerr);
		} else {
                	arsac(LLN,fname,&x,&nerr);
		}
	}
        if(sachdr.ihdr[6] == 6){
		printf("NVHDR=%d - original Sac file %s\n",sachdr.ihdr[6],fname);
        } else if(sachdr.ihdr[6] == 7){
		printf("File version: NVHDR=%d - extended Sac file %s\n",sachdr.ihdr[6],fname);
	} else {
		printf("NVHDR=%d - not a Sac file %s\n",sachdr.ihdr[6],fname);
        	return(EXIT_SUCCESS);
	}
	if(nerr != -1){
printf("Tabulation of header values showing index(C language), name VALUE from header \n");
printf("and value obtained using getfhv, getnhv, getlhv, getkhv and getdhv(new). \n");
printf("This is followed by a tabulation of a few data values. The plot SHWSAC.PLT displays\n");
printf("the entire trace\n");
                parsesac(LLN,x,onlyhead);
		if(onlyhead == FALSE)free(x);
	}
        return(EXIT_SUCCESS);
}

void  gcmdln(int argc, char**argv,int *sacbin, int *onlyhead, char *fname){
/* parse command line argiments, e.g., 
	shwsac -B | -A  sacfile
*/
	char *cp;

	*sacbin = TRUE;
	*onlyhead = FALSE;

	while(--argc > 0 ){
			cp = argv[1];
		if(*argv[1] == '-'){
			if(strcmp("-A",cp) == 0) {
				*sacbin = FALSE;
			} else if (strcmp("-B",cp) == 0){
				*sacbin = TRUE;
			} else if (strcmp("-H",cp) == 0){
				*onlyhead = TRUE;
			} else if (strcmp("-h",cp) == 0){
				usage();
			} else if (strcmp("-?",cp) == 0){
				usage();
			}
		} else {
			strcpy(fname,argv[1]);
		}
		argv++ ;
	}
}

void usage(void){
	fprintf(stderr,"Usage: shwsac [-B] [-A] [-?] file_name\n");
	fprintf(stderr," -B    file is a binary SAC file (default)\n");
	fprintf(stderr," -A    file is an ASCII SAC file (default)\n");
	fprintf(stderr," -?    Command help - do not run\n");
	fprintf(stderr," -h    Command help - do not run\n");
	fprintf(stderr," file_name Name of SAC file (required)\n");
	exit(EXIT_SUCCESS);
}


/* string patterns */

/*
static int  NRSTR=70;
char *rstr[] = {
"DELTA   ", "DEPMIN  ", "DEPMAX  ", "SCALE   ", "ODELTA  ", 
"B       ", "E       ", "O       ", "A       ", "FMT     ", 
"T0      ", "T1      ", "T2      ", "T3      ", "T4      ", 
"T5      ", "T6      ", "T7      ", "T8      ", "T9      ", 
"F       ", "RESP0   ", "RESP1   ", "RESP2   ", "RESP3   ", 
"RESP4   ", "RESP5   ", "RESP6   ", "RESP7   ", "RESP8   ", 
"RESP9   ", "STLA    ", "STLO    ", "STEL    ", "STDP    ", 
"EVLA    ", "EVLO    ", "EVEL    ", "EVDP    ", "MAG     ", 
"USER0   ", "USER1   ", "USER2   ", "USER3   ", "USER4   ",
"USER5   ", "USER6   ", "USER7   ", "USER8   ", "USER9   ", 
"DIST    ", "AZ      ", "BAZ     ", "GCARC   ", "SB      ", 
"SDELTA  ", "DEPMEN  ", "CMPAZ   ", "CMPINC  ", "XMINIMUM", 
"XMAXIMUM", "YMINIMUM", "YMAXIMUM", "ADJTM   ", "TIMMAX  ", 
"TIMMIN  ", "FHDR67  ", "FHDR68  ", "FHDR69  ", "FHDR70  " 
	};
static int  NDSTR=22;
char *dstr[] = {
"DELTA   ", "B       ", "E       ", "O       ", "A       ", 
"T0      ", "T1      ", "T2      ", "T3      ", "T4      ", 
"T5      ", "T6      ", "T7      ", "T8      ", "T9      ", 
"F       ", "EVLO    ", "EVLA    ", "STLO    ", "STLA    ",
"SDELTA  ", "SB      "
	};
static int NISTR=40;
char *istr[] = {
"NZYEAR  ", "NZJDAY  ", "NZHOUR  ", "NZMIN   ", "NZSEC   ", 
"NZMSEC  ", "NVHDR   ", "NINF    ", "NHST    ", "NPTS    ", 
"NSNPTS  ", "NSN     ", "NXSIZE  ", "NYSIZE  ", "NHDR15  ", 
"IFTYPE  ", "IDEP    ", "IZTYPE  ", "IHDR4   ", "IINST   ", 
"ISTREG  ", "IEVREG  ", "IEVTYP  ", "IQUAL   ", "ISYNTH  ", 
"IHDR11  ", "IHDR12  ", "IHDR13  ", "IHDR14  ", "IHDR15  ", 
"IHDR16  ", "IHDR17  ", "IHDR18  ", "IHDR19  ", "IHDR20  ", 
"LEVEN   ", "LPSPOL  ", "LOVROK  ", "LCALDA  ", "LHDR5   " 
	};
static int NCSTR=24;
char *cstr[] = {
"KSTNM   ", "KEVNM   ", "KEVNMC  ", "KHOLE   ", 
"KO      ", "KA      ", "KT0     ", "KT1     ", 
"KT2     ", "KT3     ", "KT4     ", "KT5     ", 
"KT6     ", "KT7     ", "KT8     ", "KT9     ", 
"KF      ", "KUSER0  ", "KUSER1  ", "KUSER2  ", 
"KCMPNM  ", "KNETWK  ", "KDATRD  ", "KINST   "
	};

static int NESTR=50;
char *estr[] = {
	"ITIME   ", "IRLIM   ", "IAMPH   ", "IXY     ", "IUNKN   ", 
	"IDISP   ", "IVEL    ", "IACC    ", "IB      ", "IDAY    ", 
	"IO      ", "IA      ", "IT0     ", "IT1     ", "IT2     ", 
	"IT3     ", "IT4     ", "IT5     ", "IT6     ", "IT7     ", 
	"IT8     ", "IT9     ", "IRADNV  ", "ITANNV  ", "IRADEV  ", 
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ", 
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ", 
	"ISRO    ", "INUCL   ", "IPREN   ", "IPOSTN  ", "IQUAKE  ", 
	"IPREQ   ", "IPOSTQ  ", "ICHEM   ", "IOTHER  ", "IGOOD   ", 
	"IGLCH   ", "IDROP   ", "ILOWSN  ", "IRLDTA  ", "IVOLTS  "
	} ;
static int NIVAL=50;
char *Istr[] = {
	"IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IFTYPE  ", "IDEP    ",
	if(onlyhead == FALSE){
	"IDEP    ", "IDEP    ", "IDEP    ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ", "IZTYPE  ",
	"IZTYPE  ", "IZTYPE  ", "IRADNV  ", "ITANNV  ", "IRADEV  ",
	"ITANEV  ", "INORTH  ", "IEAST   ", "IHORZA  ", "IDOWN   ",
	"IUP     ", "ILLLBB  ", "IWWSN1  ", "IWWSN2  ", "IHGLP   ",
	"ISRO    ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IEVTYP  ",
	"IEVTYP  ", "IEVTYP  ", "IEVTYP  ", "IQUAL   ", "IQUAL   ",
	"IQUAL   ", "IQUAL   ", "IQUAL   ", "ISYNTH  ", "IDEP    "    
};
*/

#define NPIXEL 1000

void parsesac(int LN, float *x, int onlyhead){
	int i, nerr, npts;
	float rval;
	double dval;
	int ival, lval;
	char cval[9];
	int inc;
        char ts1[9], ts2[9];



	/* output float header */
	printf("\nREAL         INDEX  NAME            VALUE  getfhv VALUE\n");
        for(i=0 ; i < NRSTR ; i++){
		getfhv(rstr[i],&rval,&nerr);
		printf("             %5d  %8s%13.6e %13.6e\n",i,rstr[i],sachdr.rhdr[i], rval);
	}
	/* output integer header */
	printf("\nINTEGER      INDEX  NAME            VALUE  getnkv VALUE\n");
	for (i=0 ; i < NISTR ; i++){
		getnhv(istr[i],&ival,&nerr);
		printf("             %5d  %8s%13d %13d\n",i,istr[i],sachdr.ihdr[i], ival);
	}
	/* output enumerated header */
	printf("\nENNUMERATED  INDEX  NAME            VALUE     ENU VALUE\n");
	for(i=15; i< 25 ; i++){
		getihv(istr[i],cval,&nerr);
		printf("             %5d  %8s%13d     %8s\n",i,istr[i],sachdr.ihdr[i], cval);
	}
	/* output logical header */
	printf("\nLOGICAL      INDEX  NAME            VALUE  getlhv VALUE (1=TRUE 0=FALSE)\n");
	for(i=35; i < 40 ; i++){
		getlhv(istr[i],&lval,&nerr);
		printf("             %5d  %8s%13d     %8d\n",i,istr[i],sachdr.ihdr[i], lval);
	}
	/* output character header */
	printf("\nCHARACTER    INDEX  NAME            VALUE  getkhv VALUE\n");
        for(i=0; i < NCSTR ; i++){
		getkhv(cstr[i],cval,&nerr);
		strncpy(ts1,sachdr.chdr[i],8);
		ts1[8] = '\0';
		strncpy(ts2,cval,8);
		ts2[8] = '\0';
		printf("             %5d  %8s     %8s     %8s\n",i,cstr[i],ts1,ts2);
	}
	/* if NVHDR = 7 output extended header */
	if(sachdr.ihdr[6] == 7 ){
		printf("\nEXTENDED REAL INDEX  NAME            VALUE         getdhv VALUE\n");
        	for(i=0 ; i < NDSTR ; i++){
			getdhv(dstr[i],&dval,&nerr);
			printf("             %5d  %8s%21.14e %21.14e\n",i,dstr[i],dsachdr.dhdr[i], dval);
		}
	}

	if(onlyhead == FALSE){
		/* plot it */
        	getnhv("NPTS    ",&npts,&nerr);
		/* output some representative values */
		/* for plotting decimate to have 1000 pixels in length */
        	inc = npts/10 ;
		if(inc == 0){
			inc = 1 ;
		}
		printf("Selected trace values\n");
		printf("    Sample     Trace value \n");
			for(i=0; i < npts; i+= inc){
				printf("%10d %15.7f\n",i,x[i]);
			}
                pinitf("SHWSAC.PLT");
        	plotit(x,npts,sachdr.rhdr[0],sachdr.rhdr[1],sachdr.rhdr[2]);
                pend();
	}
}

void plotit(float *x, int npts, float rhdr1, float rhdr2, float rhdr3){
float ymx, ymn, ymxv, ymnv, dx, xx, yy;
int i, ipen;
	/* get extreme values of time series */
	ymx = x[1] ;
	ymn = x[1] ;
	for (i=0;i < npts; i++){
		if(x[i] > ymx){
			ymx = x[i] ;
		}
		if(x[i] < ymn){
			ymn = x[i] ;
		}
	}
	/* safety for plot */
        if(ymx == ymn){
            ymn = -ymx ;
        }
        if( ymn == 0.0 && ymx == 0 ){
            ymnv = -1.0 ;
            ymxv =  1.0 ;
        } else {
            ymnv = ymn ;
            ymxv = ymx ;
	}
        if(npts > 0){
        	dx = 8.0/npts ;
        	ipen = 3 ;
        	xx = 1.0 ;
		for(i=0;i< npts;i++){
            	xx = xx + dx ;
            	yy = 3.5 + 2.0*(x[i] - ymnv)/(ymxv - ymnv) ;
            	plot(xx,yy,ipen);
            	ipen= 2 ;
		}
        	symbol(  1.0,0.6,0.10,"NPTS=",0.0,5) ;
        	number(999.0,0.6,0.10,(float)npts,0.0,-1) ;
        	symbol(  3.0,0.6,0.10," MAX=",0.0,5) ;
        	number(999.0,0.6,0.10,ymx,0.0,1003) ;
        	symbol(  1.0,0.4,0.10,"DT=",0.0,3) ;
        	number(999.0,0.4,0.10,rhdr1,0.0,1003) ;
        	symbol(  3.0,0.4,0.10," MIN=",0.0,5) ;
        	number(999.0,0.4,0.10,ymn,0.0,1003) ;
	}
}
