/* convert SAC binary SUN to SAC binary Intel */
/* Robert B. Herrmann
   May 8, 2001
   Seoul National University

	Changes:

	27 APRIL 2003 - build in intelligence to determine if byte order
		should be reversed. the criteria is that 
		o at least one integer header = -12345
		o and NO integer header values are < -99999
		  note we could also try to look at the year doy field
		  for time series (but may not work for spectra other stuff
	06 MAY 2014 - return value for main (Larry Baker, USGS Menlo Park)
	17 OCT 2022 - Updated to accommodate upgrade of SAC which 
		adds support for NHVDR=7 which appends 22 8-dyne 
		double precision after the
		waveform to support  DELTA, B, E, O, A, T0 ... T9, F, 
		EVLO, EVLA, STLO, STLA, SB, and SDELTA

		This means implementing ieee4intel and the redoing 
		the data waveform read

		The -I means that if the sac file is in the correct local
		machine byte order, then the input from stdin is
		copied directy to stdout, otherwise it is byte swapped.
 
		If -I is not specified then the byte swap always occurs

		Finally the program will stop if NVHDR != 6 or != 7 
		and if there is no -12345 in the integer header or if there is 
		some integer < -99999  which would be a strange numbevr given the
		usage of the integer field numbers

		If the byte order for 4 byte integer or floats n IEEE is A B C D then the
		byte order for INTEL is D C B A
		For doubles *8 byte) the comparison is A B C D E F G H <=> H G F E D C B A
		Note the character entries in the character portion of the header is NOT swapped

		Interesting info is at
			https://docs.oracle.com/cd/E19205-01/819-5265/bjbds/index.html
		and
			https://en.wikipedia.org/wiki/SPARC

*/

#include	<stdio.h>
#include	<stdlib.h>
#ifdef MSDOS
#include	<fcntl.h>
#endif
#include	<strings.h>

#define BUFFLOAT	280
#define NUMFLOAT	70
#define BUFINTEGER	160
#define NUMINTEGER	40
#define BUFCHAR		192
#define NUMCHARACTER	192
#define BUFDOUBLE	176
#define NUMDOUBLE	22
#define BUF		8192

#define DO_FLOAT        1
#define DO_INTEGER      2
#define DO_DOUBLE       3

#define YES             1
#define NO		0


/* selected SAC header values */
#define H_DELTA  0
#define H_NVHDR  6
#define H_NPTS 9
#define INT int
/* headers */
/* function prototypes */
void ieee2intel( char *arr, int nbytes, int type);
void ieee4intel( char *arr, int nbytes, int type);
void gcmdln(int argc, char **argv, int *do_intelligence);
void usage();
int issacfile(int *npts, int *nvhdr, int *swapin);
void doconvert(void);
int  issacheader(void);

/* global variable */

struct sachdr_  {
        float rhdr[70];
        INT ihdr[40];
        char chdr[24][8];
        } sacdata  ;
char arr[BUF] ;


int aargc;
char **aargv;



void main(int argc, char **argv)
{
int do_intelligence;	/* 1 force conversion */
int npts;	/* number of dta samples */
int nvhdr;
int nread;	/* nvhdr = 6 old sac file
		   nvhdr = 7 new sac with doubles at end */
int nbufread;
int doconvert;
int swapin;
int ntoread;


	/* parse the command line */
	gcmdln(argc, argv, & do_intelligence);

	/* is the input stream a sac file? */
	if( issacfile(&npts,&nvhdr,&swapin) == NO){
		fprintf(stderr,"Input does not have a valid sac header \n");
		exit(1);
	}

fprintf(stderr,"npts %d nvhdr %d swapin %d \n",npts, nvhdr,swapin);
	if(do_intelligence){
		doconvert = swapin;		/* output in current machine byteorder */
	} else {
		doconvert = YES;		/* force a byteswap */
	}
	/* perform the conversion by by starting at the beginning of the file  */
	fseek(stdin, 0L, SEEK_SET);
	/* always get input in machine byte order */
	fread(&sacdata,sizeof(struct sachdr_), 1,stdin);
	/* output the header */
	if(doconvert == YES){
		ieee2intel((char *)sacdata.rhdr, NUMFLOAT   , DO_FLOAT);
		ieee2intel((char *)sacdata.ihdr, NUMINTEGER , DO_INTEGER);
	}
	fwrite(&sacdata,sizeof(struct sachdr_), 1,stdout);
	/* now handle the waveform data */
	nbufread = npts*4/BUF;
	ntoread= (npts%BUF)*4;
	while(nbufread > 0 ){
		fread(arr, sizeof(char), BUF, stdin);
		if(doconvert == YES)
			ieee2intel( arr, BUF/4, DO_FLOAT );
		fwrite(arr, sizeof(char), BUF, stdout) ;
		nbufread--;
	}
	if(ntoread > 0){
		fread(arr, sizeof(char), ntoread, stdin);
		if(doconvert == YES)
			ieee2intel( arr,ntoread/4, DO_FLOAT );
		fwrite(arr, sizeof(char), ntoread, stdout) ;
	}
	/* IF nvhdr == 7 read and output the double header at the end */
	if(nvhdr == 7) {
		fread(arr, sizeof(char), BUFDOUBLE, stdin);
			ieee4intel( arr, BUFDOUBLE, DO_DOUBLE );
		if(doconvert == YES)
			ieee4intel( arr, BUFDOUBLE, DO_DOUBLE );
		fwrite(arr, sizeof(char), BUFDOUBLE, stdout) ;
	}
	exit(0);

}


int issacfile(int *npts, int *nvhdr, int *swapin)
{
	int i;
	if(fread(&sacdata,sizeof(struct sachdr_), 1,stdin) != 1 ){
		fprintf(stderr, "Input is too short for sac header \n");
		exit(1);
	}
	/* tests*/
	/* NVHDR must be 6 or 7, but
           since this could be a random occurrence, reject the file if
	   there is NO -12345 and at least one ihrd < -99999
           so it is sacheader?
	   if not byte swap
           so if it is not sacheader
             abort
           else
             copy to stdout
	*/
	*swapin = NO;
	if( issacheader() == NO){
		/* try a swap */
		ieee2intel((char *)sacdata.rhdr, NUMFLOAT   , DO_FLOAT);
		ieee2intel((char *)sacdata.ihdr, NUMINTEGER , DO_INTEGER);
		*swapin = YES;
		if( issacheader() == NO){
			return (NO);
		} else {
		*npts = sacdata.ihdr[H_NPTS ];
		*nvhdr = sacdata.ihdr[H_NVHDR];
			return (YES);
		}
	} else {
		*npts = sacdata.ihdr[H_NPTS ];
		*nvhdr = sacdata.ihdr[H_NVHDR];
		return (YES);
	}
}

int issacheader(void)
{
	/* A valid sac integer header muhst have NVHDR = 6 or 7
		AND
	   no integer value < -99999
		AND
	   at least one -12345
	This last one may be too strict
	*/

	int i;
	int neg12345 = 0 ;
	int havebigneg = 0 ;;
	fprintf(stderr,"NVHDR    %d\n",sacdata.ihdr[H_NVHDR]);
	if(sacdata.ihdr[H_NVHDR] == 6 || sacdata.ihdr[H_NVHDR] == 7){
		for(i=0 ; i < NUMINTEGER; i++){
			if(sacdata.ihdr[i] == -12345) neg12345++;
			if(sacdata.ihdr[i] <  -99999) havebigneg++;
		}
		if(neg12345 ==0 && havebigneg > 0 ){
			return(NO);
		} else {
			return(YES);
		}

	} else {
		fprintf(stderr," NVHDR != 6 && != 7\n");
		return (NO);
	}
}


void ieee2intel( char *arr, int nword, int type)
/*	arr	array of bytes to be converted in place 
	nbytes	number of  bytes to be converted
	type	0 float
*/
{
	char a, b, c, d;
	int i;
	for (i=0; i < nword*4 ; i+=4 ){
		a= arr[i]  ;
		b= arr[i+1];
		c= arr[i+2];
		d= arr[i+3];
		arr[i  ] = d;
		arr[i+1] = c;
		arr[i+2] = b;
		arr[i+3] = a;
	}
}

void ieee4intel( char *arr, int nword, int type)
/*	arr	array of bytes to be converted in place 
	nbytes	number of  bytes to be converted
	type	0 float
*/
{
	char a, b, c, d, e, f, g, h;
	char *cp;
	int i;
	for (i=0; i < nword*8 ; i+=8 ){
		a= arr[i]  ;
		b= arr[i+1];
		c= arr[i+2];
		d= arr[i+3];
		e= arr[i+4];
		f= arr[i+5];
		g= arr[i+6];
		h= arr[i+7];
		arr[i  ] = h;
		arr[i+1] = g;
		arr[i+2] = f;
		arr[i+3] = e;
		arr[i+4] = d;
		arr[i+5] = c;
		arr[i+6] = b;
		arr[i+7] = a;
	}
}

void gcmdln(int argc, char **argv, int *do_intelligence){

	/* parse the command line flags */
	*do_intelligence = NO;
	while(argc-- > 1){
		if(*argv[1] == '-'){
			switch(argv[1][1]){ 
			case 'I':
				*do_intelligence = YES; /* ensure output compatible with actual hardware */
				break;
			case '?':
                                usage();
                                break;
			case 'h':
                                usage();
                                break;
                        default:
                                break;                                                                  }
                }
                argv++;
        }                                 
}

void usage()
{
fprintf(stderr,	"Usage: saccvt [-I] [-h] [-?]\n");
fprintf(stderr, " Example: saccvt < SAC_BINARY > tmp ; mv tmp SAC_BINARY\n");
fprintf(stderr,	"    Convert SAC binary  IEEE to INTEL format \n");
fprintf(stderr,	"    and  SAC binary  INTEL to IEEE format \n");
fprintf(stderr,	"    All 4 byte integers and floats (a,b,c,d) are\n");
fprintf(stderr,	"    transposed to (d,c,b,a)\n");
fprintf(stderr, "    and doubles (a,b,c,d,e,f,g,h) to (h,g,f,e,d,c,b,a) for NVHDR=7\n"); 
fprintf(stderr,	" -I   (default none)  intelligently guess whether to convert to match use architecture\n");
fprintf(stderr, "    else always convert\n");
fprintf(stderr,	" -h   (default none)  this help message\n");
fprintf(stderr,	" -?   (default none)  this help message\n");
exit(0);
}
