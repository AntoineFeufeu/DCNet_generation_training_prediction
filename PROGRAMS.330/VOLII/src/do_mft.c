/* do_mft - main driver
 * 	Changes:
 * 	04 SEP 2002 - put azimuth on first screen for reference
 *	30 OCT 2002 - extended selections to very short periods for exploration
 *	23 JUL 2003 Following Email from Meijian An, 
 *	Department of Geophysics, Institute of Astronomics, 
 *	Geophysical and Atmospheric Sciences, *	University of Sao Paulo
 *	Sao Paulo, Brazil, carefully free memory
 *
 * 	13 SEP 2003 - put in controls for sacmft96 contour plot presentation
 *		in do_mft3.c sclstr[] and absstr[]
 * 	12 FEB 2004 - put in error return and termination if sacmft96 cannot
 * 	be found
 * 	23 OCT 2004 - minor clean up in clearscreen in do_mft4
 *	11 AUG 2005 - added a distance sort and started command line 
 *              arguments so that
 *		we can reset the number per page e.g., -N10 is the default
 *	10 JUN 2006 - added -G flag to that the edited dispersion is called
 *		sacfile.dsp
 *	01 OCT 2006 - fixed allocation error in pathfind in do_mft3.c
 *		I did not have the +2 in the calloc and realloc lines with
 *		strlen(name)+2
 *      13 NOV 2012 - added -T flag to fork off a special script to
 *      	overlay tomography group velocity on top of the MFT96.PLT output
 *      27 MAY 2014 - in insertnode in fmenunod.c use strncpy and null terminate
 *              big error in do_mft4.c line 1620 - got rid of test OVERLAYTOMO and change 
 *                  for (i = 0 ; i < sizeof(m)/sizeof(struct menu) ; i++ ) to
 *                  for (i = 0 ; i < sizeof(mp)/sizeof(struct menu) ; i++ )
 *                  This fixed a problem in insertnode
 *      12 JUN 2014 - fixed so that if distance is negative, e.g., not defined 
 *                  if -12345, the program will not crach, but will rather 
 *                  flag the entry with a RED instead of BLUE color in
 *                  the routine Pick_file
 *      07 JUN 2016 - refined the Automatic mode to automatically advance by one trace
 *                  and skip the Units/Display parameter pages. this speed up processing
 *                  To change the parameters, get out of automatic mode and then Exit
 *                  and then go back to the trace desired.
 *                  Also when in Automatic mode the previous mode selection is used. Again
 *                  this saves some menu time.
 *      10 FEB 2020 - when initializing use the XORArrow instead of the
 *                    hardware Arrow since the hardware version was sluggish
 */ 	
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

/* Turn this on if using mtrace */
#ifdef MTRACE
#include <mcheck.h> /* get the mtrace definitions */
#endif


#include <sys/stat.h>
#include <unistd.h>

#define 	PNAME 		"do_mft"
#include	"nfmenu.h"
#include	"nmenu.h"

/* DEFINES */

#define		FILE_SAC_BINARY	1
#define		FILE_SAC_ASCII 	2
#define		FILE_UNKNOWN	0

#define		MENU_FILE_QUIT	-1
#define		MENU_FILE_NEXT	-2
#define		MENU_FILE_PREV	-3
#define		MENU_FILE_NUMB	10
 
extern struct menu menu_p1[] ;

#define ON      1
#define OFF     0
/* plot window for MFT96 graphs 
	   The plot will be placed between (XL,YL) and (XH,YH) 
	   The MFT96.PLT will be shifted by adding (XOFF,YOFF) */
#define XL      0.0
#define XH      10.00
#define YL      1.5
#define YH      8.0
#define XOFF	0.0
#define YOFF	-0.5;
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define LIN 0
#define LOG 1
#include "sacsubc.h"
#include "csstime.h"
#define MAXSACARR 132000

char *ftype[] = { "UNK", "BIN", "ASC"  }; 
char Strout[100], str1out[100];
#include "calplot.h"
#include "grphsubc.h"

/* GLOBAL VARIABLE DECLARATIONS */
	/* display information */
int HasMouse; 
float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
int Color;

int black = 1;	/* 1 white background, 0 black background */
int kolor = 1;	/* 1 gray scale, 2 color */

fmenu **file_menu;
int ndfiles;

struct date_time dt_begin, dt_origin, dt_ptime, dt_stime;
struct date_time dt_refer;
float evla, evlo, evdp, stla, stlo, delta, 
	dist, baz, az, gcarc, b, e, o, a, t0;
int npts, nzyear, nzjday, nzhour, nzmin, nzsec, nzmsec;
int ihdr11;
char kstnm[9], kcmpnm[9], kevnm[9], kevnmc[9];

/* prototypes */
int main(int argc, char **argv , char **arge);
void Inquire_file(int argc, char **cp, int *npage);
void DInquire_file(int *argc, char **cp);
int type_file(char *cstr,int *nsamp, int *fsize, char *datetime, char *kstnm, char *kcmpnm);
void mgwrtxt(float x0, float y0, char *str, int cmd, int color);
int Pick_file(int npage, int argc, int *curpage);
	/* ret =  */
int do_page1(int npage, int *curpage);
	/*  ret = */
int do_page2(char *fname);
	/*  ret = */
extern int do_sacmft(char *fname, int Units);
extern int do_units(void);
int do_page3(char *fname, char *unitstr);
	/*  ret = */
int do_page4(char *fname, char *type);
	/*  ret = */
void do_check(float xl, float yl, float xh, float yh);
void do_reject(float xl, float yl, float xh, float yh);
void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);
void usage(void);

/* structures to implement distance sort */
void gnomesort(int n, float ar[], int key[] ) ;
	/* http://www.cs.vu.nl/~dick/gnomesort.html */
int *sortptr;
float *sortfloat;

/* integers to handle command flags */
#define MAXNUMPERPAGE 10
int numperpage = MAXNUMPERPAGE;
int interstationgreen = 0;
char *ccp;
int overlaytomography = 0;

float distmin = 0.0 ;
float distmax = 1.0e+37 ;
int ihdr11min = -123456.;
int ihdr11max =  123456.;

int do_auto_feed = OFF;
int do_auto_i_next = -1 ;
int do_auto_curfile;
int do_auto_nfiles;
int do_auto_curpage;
int do_auto_npages;
int do_auto_i_max ;
int do_auto_ihdr11 ;

int do_dist = OFF;
int do_ihdr11 = OFF;
int do_phvel = OFF;

int AutoMode = -1 ;
int Mode  = -1;
int iMode  = -1;





int Units = 0;
int dotype = 0;
extern int doshade;
extern int doabs;
extern int doscl;
extern int dovrb;
extern int XaxisPeriod;
extern int XaxisLog;

/* **cp is pointer to waveform names - this is created parsing the command flags
   from what is left over
   **tcp is pointer to original char **argv
*/
char **cp = NULL;
char **tcp = NULL;

fmenu *p;
extern fmenu *Start, *End;

int main(int argc, char **argv, char **arge)
{
	int i;
	int npage;
	int aargc, taargc;
	int curpage = 0;
	int ret, line, page;
	float xl, yl, xh, yh;

	/* Call to mtrace routines in glibc library */
#ifdef MTRACE
	    mtrace();  /* Turn on mtrace function */
#endif

	taargc = argc;
	/* make private copy of file names to be examined */
	cp = (char **) calloc(taargc, sizeof (char *));
	tcp = (char **) calloc(taargc, sizeof (char *));
	/* argc[0] is the name of the executable */

#ifdef DEBUG
		fprintf(stderr,"argc: %d aargc: %d\n",argc,aargc);
#endif

	for (i=1 ; i < taargc; i++){
		tcp[i-1] = (char *)calloc(strlen(argv[i])+1, sizeof (char));
		strcpy(tcp[i-1], argv[i]);
		strcat(tcp[i-1],"\0");
	}
	aargc = 0;
	/* parse command line */
	for(i=0 ; i < taargc -1 ; i++){
		/* note need error check for atoi and atof */
		if(strncmp(tcp[i],"-11MIN",6) ==0){
			i++;
			ihdr11min = atoi(tcp[i]);
			if(ihdr11min < 0)
				ihdr11min = 0;
			do_ihdr11 = ON;
		} else if(strncmp(tcp[i],"-11MAX",6) ==0){
			i++;
			ihdr11max = atoi(tcp[i]);
			if(ihdr11max < 0)
				ihdr11max = 0;
			do_ihdr11 = ON;
		} else if(strncmp(tcp[i],"-DMIN",5) ==0){
			i++;
			distmin = atof(tcp[i]);
			if(distmin < 0 ) distmin = 0.0;
			do_dist = ON;
		} else if(strncmp(tcp[i],"-DMAX",5) ==0){
			i++;
			distmax = atof(tcp[i]);
			do_dist = ON;
		} else if(strncmp(tcp[i],"-IG",3) ==0){
			do_phvel = 1;
		} else if(strncmp(tcp[i],"-P0",3) ==0){
			do_phvel = 2;
		} else if(*tcp[i]  == '-'){
			switch(tcp[i][1]){
				case 'N':
					ccp = tcp[i];
					ccp++;
					ccp++;
					numperpage = atoi(ccp);
					break;
				case 'G':
					interstationgreen = 1;
					break;
				case 'T':
					overlaytomography = 1;
					break;
				case 'h':
				case '?':
					usage();
					return (0);;
				default:
					break;
			}
		} else  {
			cp[aargc] = (char *)calloc(strlen(tcp[i])+1, sizeof (char));
			strcpy(cp[aargc++],tcp[i]);
		}
	}
	/* open graphics - note the calxvig window size can be
		fixed by adding simething like this to
		.Xdefaults:
		plotxvig.calxvig.nm*geometry: 800x600+00+00 
		where nm is temporarily the PNAME
	*/
	ginitf("INTER", PNAME);
	/* get information about the display */
	ginfo(&HasMouse, &XminDev, &YminDev, 
		&XmaxDev, &YmaxDev, &XminClip, 
		&YminClip, &XmaxClip, &YmaxClip,&Color);
		if(Color >= 4)
			black = 0;
		else
			black = 1;
		kolor = Color%4;
	gcursor("XORArrow");

	gmesg("Examining Files ");
	/* initialize the nodes for the file menu */
		Start = End = getnode();
	/* file_menu relates to those displayed by page */
		file_menu = (fmenu **)calloc(14,sizeof(fmenu *));
	/* initialize distance sort */
		sortfloat = (float *)calloc(aargc, sizeof(float));
		sortptr = (int *)calloc(aargc, sizeof(int));
		DInquire_file(&aargc,cp);
	/* grab information about all files on the command line */
		Inquire_file(aargc-1,cp,&npage);
	/* OK choose the file for use */
		ret = 0;
		while(ret >= 0){
			/* recompute the page and line entries */
			p = Start;
			ndfiles = 0;
			npage = 0;
			while(p != End){
				ndfiles++;
				line = (ndfiles-1)%numperpage; 
				page = (ndfiles-1)/numperpage; 
				if(page > npage)
					npage = page;
				xl = 0.1;
				yl = 6.0 -0.5*line;
				xh = 9.9;
				yh = yl + 0.5;
				p->xl = xl;
				p->yl = yl;
				p->xh = xh;
				p->yh = yh;
				p->page = page;
				#ifdef DEBUG
				printf("%d %d %d %d %s\n",p->lstrmx, p->page, p->line, p->type,p->str);
				#endif
				p = p->next;
			}
			/* select an ooption */
			do_auto_nfiles =ndfiles;
			ret = Pick_file(npage, aargc, &curpage );
		}
	/* terminate the program */
	/* eventually invoke the signal mechanism to cause a
		clean up of all memory allocated

		also this may be a place to save the processing
		summary
		so that the command will be nnn Dir
		or nnn will look for the .xxxx file int he current directory
	*/
	p = Start;
	while(p != End)
		deletenode(p);
	/*
	 * Remove free the cp memory
	 */
	for (i=1 ; i < aargc; i++){
		free(cp[i-1]);
	}
	free (cp);
	free(file_menu);

	gend(1);
	return (0);
}

int Pick_file(int npage, int aargc, int *curpage)
{
	int i, ls;
        int dist_is_neg;
	int page_entry;
	i = 0;
	gframe(2);
	ginfo(&HasMouse, &XminDev, &YminDev, 
		&XmaxDev, &YmaxDev, &XminClip, 
		&YminClip, &XmaxClip, &YmaxClip,&Color);
	newpen(1);
	if( do_auto_feed == OFF){
		gwidth(0.02);
		gcent(5.0,7.2,0.2,"MFT96 File Selection",0.0);
		gwidth(0.00);
		if(do_phvel != OFF){
			mgwrtxt(0.2,6.8,
			"Type File                                                             Ihdr11", 1,   1);
			mgwrtxt(0.2,6.6,
			"Stnm    Cmpnm         Npts    Bytes   First Sample Time          Dist     Az ", 0,   4);
		} else {
			mgwrtxt(0.2,6.8,
			"Type File                                                             ", 1,   1);
			mgwrtxt(0.2,6.6,
			"Stnm    Cmpnm         Npts    Bytes   First Sample Time          Dist     Az Proc", 0,   4);
		}
	}
	page_entry = 2;
	p = Start;
	while(p != End){
#ifdef DEBUG
		printf("%d %s %s %d %d %d %s %s %s %d %f\n",
			p->type, ftype[p->type], 
			p->str, p->line, 
			p->fsize, p->nsamp, 
			p->kstnm, p->kcmpnm, 
			p->datetime, p->page,p->dist);
#endif
		if(p->xl > 0.0 && *curpage == p->page){
			if(do_phvel != OFF){
			sprintf(Strout ,"%3s %-65s %6d", ftype[p->type], p->str,p->ihdr11);
			} else {
			sprintf(Strout ,"%3s %-85s", ftype[p->type], p->str);
			}
			ls = strlen(Strout);
			if(ls > 79)ls = 79;
			Strout[ls] = '\0';
			if(p->dist < 1.0 && p->dist >= 0.0){
			/* special fix for exploration processing */
				if(do_phvel != OFF)
 				sprintf(str1out,
				"%-8s %-8s %8d %8d %23s %10.5f %5.0f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist, p->az);
				else
 				sprintf(str1out,
				"%-8s %-8s %8d %8d %23s %10.5f %5.0f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist, p->az);
				dist_is_neg = 4 ;
			} else if(p->dist < 0.0 ){
				if(do_phvel != OFF)
 				sprintf(str1out,
				"%-8s %-8s %8d %8d %23s %10.3f %5.0f ",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist, p->az);
				else
 				sprintf(str1out,
				"%-8s %-8s %8d %8d %23s %10.3f %5.0f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist, p->az);
				dist_is_neg = 2 ;
			} else {
				if(do_phvel != OFF)
 				sprintf(str1out,
				"%-8s %-8s %8d %8d %23s %10.3f %5.0f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist, p->az);
				else
 				sprintf(str1out,
				"%-8s %-8s %8d %8d %23s %10.3f %5.0f",
					p->kstnm, p->kcmpnm, p->nsamp,
					p->fsize, p->datetime, p->dist, p->az);
				dist_is_neg = 4 ;
			}
				
			ls = strlen(str1out);
			if(ls > 79)ls = 79;
			str1out[ls] = '\0';
			if( do_auto_feed == OFF){
				mgwrtxt(p->xl+0.1,p->yl+0.25,Strout, 0,  1);
				mgwrtxt(p->xl+0.1,p->yl+0.05,str1out,0,  dist_is_neg);
				gbox(0.1,p->yl,9.9,p->yh);
				if(p->used > 0){
					do_check(p->xl,p->yl, p->xh,p->yh);
				} else if(p->used < 0){
					do_reject(p->xl,p->yl, p->xh,p->yh);
				}
			}
			page_entry++;
			menu_p1[page_entry].fileptr = i;
			menu_p1[page_entry].action = page_entry -2;
			menu_p1[page_entry].xl = p->xl;
			menu_p1[page_entry].yl = p->yl;
			menu_p1[page_entry].xh = p->xh;
			menu_p1[page_entry].yh = p->yh;
			file_menu[page_entry] = p;
		
		}
		p = p->next;
	}
	/* set menu entries with no data */
	for( i= page_entry +1; i < MENU_FILE_NUMB+3; i++){
		
		menu_p1[i].xl = -1;
		menu_p1[i].yl = -1;
		menu_p1[i].xh = -1;
		menu_p1[i].yh = -1;
		menu_p1[i].action = 0;
	}

	/* now display the menu */
	
	return(do_page1(npage, curpage));
}



/* put in a check mark */
void do_check(float xl, float yl, float xh, float yh)
{
	float x1,y1,x2,y2,x3,y3,yc,xc,ht,hht;
	newpen(2);
	hht = yh - yl;
	ht = 0.6*hht;
	yc = yl + 0.5*hht - 0.5*ht;
	xc = xh - hht;
	x1 = xc - 0.2*ht; y1 = yc + 0.2*ht;
	x2 = xc; y2 = yc;
	x3 = xc + 0.2*ht; y3 = yc + 0.8*ht;
	plot(x1,y1,3);plot(x2,y2,2);plot(x3,y3,2);plot(x3,y3,3);
	newpen(1);
}

/* put in a rejection mark */
void do_reject(float xl, float yl, float xh, float yh)
{
	float x1,y1,x2,y2,x3,y3,x4,y4,yc,xc,ht,hht;
	newpen(1);
	hht = yh - yl;
	ht = 0.4*hht;
	yc = yl + 0.5*hht ;
	xc = xh - hht;
	x1 = xc - 0.5*ht; y1 = yc - 0.5*ht;
	x2 = xc + 0.5*ht; y2 = yc + 0.5*ht;
	x3 = xc - 0.5*ht; y3 = yc + 0.5*ht;
	x4 = xc + 0.5*ht; y4 = yc - 0.5*ht;
	plot(x1,y1,3);plot(x2,y2,2);plot(x3,y3,3);plot(x4,y4,2);plot(x4,y4,3);
	newpen(1);
}


/* select file to be processed from list valid files 
	this takes time since we read the SAC files */
void Inquire_file(int aargc, char **cp, int *npage)
{
	int i,k,retval;
	int l, lstr;
	struct stat statbuf;
	float xl, yl, xh, yh;
	int action, lstrmx, type, line, fsize, nsamp, page, used;
	char kstnm[9], kcmpnm[9], datetime[24];
	/* get the maximum string length in this menu category */
	action = -1;
	lstr = 0;
	for(i=0 ; i <= aargc ; i++){
		l = strlen(cp[i]);
		if(l > lstr)lstr = l;
	}
	ndfiles = 0;
	*npage = 0;
	for(k=0;k<= aargc; k++){
		{
		i = sortptr[k];
		retval = type_file(cp[i],&nsamp,&fsize,datetime,kstnm,kcmpnm);
		if(retval == 1 || retval == 2){
			if( stat(cp[i],&statbuf)==0)
				fsize = statbuf.st_size;
			else
				fsize = -1 ;
			/* get rid of dregs */
			used = 0;
			lstrmx = lstr;
			type = retval;
			ndfiles++;
			line = (ndfiles-1)%10; 
			page = (ndfiles-1)/10; 
			if(page > *npage)
				*npage = page;
			xl = 0.1;
			yl = 6.0 -0.5*line;
			xh = 9.9;
			yh = yl + 0.5;
		
			appendnode(xl, yl, xh, yh,
				cp[i], action, lstrmx, type, line, fsize,
				nsamp, kstnm, kcmpnm, datetime,
				page, used, dist, az, baz, ihdr11);
		}
	}
	}
}
void DInquire_file(int *aargc, char **cp)
{
	int i,retval;
	int l, lstr;
	/* struct stat statbuf; */
	/* float xl, yl, xh, yh; */
	int fsize, nsamp;
	char kstnm[9], kcmpnm[9], datetime[24];
	/* float dval; */
	int do_use_trace ;
	/* get the maximum string length in this menu category */
	/*action = -1; */
	lstr = 0;
	for(i=0 ; i < *aargc ; i++){
		l = strlen(cp[i]);
		if(l > lstr)lstr = l;
	}
	ndfiles = 0;
	for(i=0;i< *aargc; i++){
		retval = type_file(cp[i],&nsamp,&fsize,datetime,kstnm,kcmpnm);
		/* if(retval == 1 || retval == 2){
			dval = dist ;
		} else {
			dval = -12345.;
		} */
		do_use_trace = ON;
		/* apply selection criteria */
		if(do_ihdr11 == ON)
			if( ihdr11 < ihdr11min || ihdr11 > ihdr11max)
				do_use_trace = OFF;
		if(do_dist == ON)
			if( dist < distmin || dist > distmax)
				do_use_trace = OFF;
		if(do_use_trace == ON){
			sortptr[ndfiles] = i;
			sortfloat[ndfiles] = dist;
			ndfiles++;
		};
	}
	gnomesort(ndfiles, sortfloat, sortptr) ;
	*aargc = ndfiles;
}


int type_file(char *cp,int *nsamp, int *fsize, char *datetime, char *ksnm, char *kcmnm)
{
	int nberr, naerr, nerr;
	float *data;
	char ostr[80];
	
/*
	struct date_time dt_begin, dt_origin, dt_ptime, dt_stime;
	struct date_time dt_refer;
	float evla, evlo, evdp, stla, stlo, delta, 
		dist, baz, az, gcarc, b, e, o, a, t0;
	int npts, nzyear, nzjday, nzhour, nzmin, 
		nzsec, nzmsec;
	char kstnm[9], kcmpnm[9], kevnm[9], kevnmc[9];
*/


	/* if the first character is a - then this is a command flag
		and not a file name */
	if(cp[0] == '-')
		return(FILE_UNKNOWN);
	/* attempt to open as a SAC file */
	nberr = -100;
	naerr = -100;
	brsac(MAXSACARR,cp,  &data, &nberr);
	if(nberr >= 0){
		free(data);
		
	}
	if(nberr < 0){
		arsac(MAXSACARR,cp,  &data, &naerr);
		if(naerr >= 0){
		free(data);
		}
	}
	if(naerr >= 0 || nberr >=0){
		getfhv("EVLA",&evla,&nerr);
		getfhv("EVLO",&evlo,&nerr);
		getfhv("EVDP",&evdp,&nerr);
		getfhv("STLA",&stla,&nerr);
		getfhv("STLO",&stlo,&nerr);
		getfhv("BAZ",&baz,&nerr);
		getfhv("AZ",&az,&nerr);
		getfhv("DELTA",&delta,&nerr);
		getfhv("DIST",&dist,&nerr);
		getfhv("GCARC",&gcarc,&nerr);
		getfhv("B",&b,&nerr);
		getfhv("O",&o,&nerr);
		getfhv("A",&a,&nerr);
		getfhv("T0",&t0,&nerr);
		getnhv("NPTS",&npts,&nerr);
		getnhv("NZYEAR",&nzyear,&nerr);
		getnhv("NZJDAY",&nzjday,&nerr);
		getnhv("NZHOUR",&nzhour,&nerr);
		getnhv("NZMIN",&nzmin,&nerr);
		getnhv("NZSEC",&nzsec,&nerr);
		getnhv("NZMSEC",&nzmsec,&nerr);
		getnhv("IHDR11",&ihdr11,&nerr);
		getkhv("KSTNM",kstnm,&nerr);

		getkhv("KCMPNM",kcmpnm,&nerr);
		getkhv("KEVNM",kevnm,&nerr);
		getkhv("KEVNMC",kevnmc,&nerr);


		/* convert the reference to epoch time */
		dt_refer.date = 1000L*nzyear + nzjday;
		dt_refer.hour = nzhour;
		dt_refer.minute = nzmin;
		dt_refer.second = (float)nzsec + (float)nzmsec/1000.0;
		/* convert to epoch */
		htoe(&dt_refer);
		etoh(&dt_refer);
		timeprintstr(&dt_refer,ostr);  
		/* create an entry for origin time */
		dt_origin.epoch = dt_refer.epoch + o;
		etoh(&dt_origin);
		/* now create one for first sample time */
		dt_begin.epoch = dt_refer.epoch + b;
		etoh(&dt_begin);
		/* now create one for first P time */
		if(a != -12345.){
			dt_ptime.epoch = dt_refer.epoch + a;
			etoh(&dt_ptime);
		}
		/* now create one for first S time */
		if(t0 != -12345.){
			dt_stime.epoch = dt_refer.epoch + t0;
			etoh(&dt_stime);
		}
		timestr(&dt_begin,ostr);  
		*nsamp = npts;
		strcpy(datetime, ostr);
		strcpy(ksnm,kstnm);
		strcpy(kcmnm,kcmpnm);
	}
	/* safety check for bad file */
	if(npts <= 0)
		return(FILE_UNKNOWN);
	if(nberr >= 0){
		return(FILE_SAC_BINARY);
	} else if(naerr >= 0 ){
		return(FILE_SAC_ASCII);
	} else {
		return(FILE_UNKNOWN);
	}
}

/* Gnome Sort - The Simplest Sort Algorithm
 * http://www.cs.vu.nl/~dick/gnomesort.html
 * */
void gnomesort(int n, float ar[], int key[]) {
	int i;
	float tmp;
	int   kmp;
	i = 1;
	while (i < n) {
		if( i == 0){
		       i++;
		} else if( ar[i-1] <= ar[i]) {
			i++;
		}
		else {
			tmp = ar[i]; ar[i] = ar[i-1]; ar[i-1] = tmp;
			kmp = key[i]; key[i] = key[i-1]; key[i-1] = kmp;
			i--;
		}
	}
}

void usage(void)
{
fprintf(stderr,"Usage: do_mft [-Nnumberperpage] [-G] [-T] [-11MIN min11] [-11MAX max11] [-DMIN dmin] [-DMAX dmax] sacfiles\n");
fprintf(stderr,"  e.g.:  do_mft -N2 -G -T *Z\n");
fprintf(stderr,"  -Nnumberperpage  (default 10) number of file menu items per page\n");
fprintf(stderr,"        this option is useful when using a slow connection since writing \n");
fprintf(stderr,"        a complete menu takes time\n");
fprintf(stderr,"  -G              (default off) The default dispersion file name is of the form \n");
fprintf(stderr,"        StationComponent.dsp, e.g., SLMBHZ.dsp \n");
fprintf(stderr,"        When working with ground-noise crosscorrelation for interstation Green functions, the\n");
fprintf(stderr,"        naming is Station1Component1Station2Component2.dsp , e.g., SLMBHZFVMBHZ.dsp\n");
fprintf(stderr,"  -T              (default off)  run script MFTDOOVERLAY\n");
fprintf(stderr,"  -11MIN min11 \n");
fprintf(stderr,"  -11MAX max11 \n");
fprintf(stderr,"  -DMIN  dmin \n");
fprintf(stderr,"  -DMAX  dmax \n");
fprintf(stderr,"   These options control the selection of the files for MFT\n");
fprintf(stderr,"   analysis. The first two use the number in the IHDR11 field,\n   which is the number of waveforms stacked for cross-correlation \n   of ground noise. The last two select the distances \n");
fprintf(stderr,"  -IG             (default false)  Inter-station phase velocity from noise cross-correlation (has pi/4 correction) \n");
fprintf(stderr,"  -P0             (default false)  Inter-station phase velocity from cross-correlation \n");
fprintf(stderr,"  -h              (default false)  Usage \n");
}
