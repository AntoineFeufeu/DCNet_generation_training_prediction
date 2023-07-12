/* Changes
  29 JUL 2015 - created this to permit interactive picks of amplitude spectra
		only the L command is implemented
  06 SEP 2015 - fixed tics to be white when the backgound color is not white
  10 SEP 2015 - implemented picking of smoothed trace in psppk_getfreq
		implemented a linear regression to get t*  using the Z Z sequence to
		get the frequency range
  17 MAR 2022 - in computing the amplitude spectra, change "float tr, ti" to "double tr, ti"
*/
#include	<stdio.h>
#include        "gsac.h"

#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include	"csstim.h"
#include	<libgen.h>

extern struct sacfile_ *sacdata;
extern int  sfgetline(FILE *fp, char s[], int lim);

static float *y = (float *)NULL;
static float *x = (float *)NULL;
extern int *sortptr;
/* save the depmin depmax of the y axes */
static float *save_depmin = (float *)NULL ;
static float *save_depmax = (float *)NULL ;


#define PSPPK_AM 0
#define PSPPK_PH 1
#define PSPPK_OVERLAY 2
#define PSPPK_PERPLOT 3
#define PSPPK_XLOG 4
#define PSPPK_XLIN 5
#define PSPPK_YLOG 6
#define PSPPK_YLIN 7
#define PSPPK_FMIN 8
#define PSPPK_FMAX 9
#define PSPPK_DEFAULT 10
#define PSPPK_AMIN 11
#define PSPPK_AMAX 12
#define PSPPK_SMOOTH 13
#define PSPPK_SMOOTH_PASS 14

#define PLOT_AM 0
#define PLOT_PH 1

static int psppk_doamph = PLOT_AM;

struct arghdr psppkarg[] = {
	{PSPPK_AM, "AMPLITUDE"	, IHDR, 0, 0, NO, "AMPLITUDE ",3},
	{PSPPK_AM, "AM"	, IHDR, 0, 0, NO, "AMPLITUDE ",-1},
	{PSPPK_PH, "PHASE"	, IHDR, 0, 0, NO, "PHASE", 2},
	{PSPPK_OVERLAY, "OVERLAY" , YHDR, 0, 1, NO, "OVERLAY [ON|OFF] ", 2},
	{PSPPK_SMOOTH , "SMOOTH"  , YHDR, 0, 1, NO, "SMOOTH [ON|OFF]"  , 1},
	{PSPPK_PERPLOT, "PERPLOT" , NHDR, 0, 1, NO, "PERLOT [n|OFF]", 2},
	{PSPPK_XLIN, "XLIN", IHDR, 0, 0, NO, "XLIN", 3},
	{PSPPK_XLOG, "XLOG", IHDR, 0, 0, NO, "XLOG", 3},
	{PSPPK_YLIN, "YLIN", IHDR, 0, 0, NO, "YLIN", 3},
	{PSPPK_YLOG, "YLOG", IHDR, 0, 0, NO, "YLOG", 3},
	{PSPPK_FMIN, "FMIN", RHDR, 0, 1, NO, "FMIN", 3},
	{PSPPK_FMAX, "FMAX", RHDR, 0, 1, NO, "FMAX", 3},
	{PSPPK_AMIN, "AMIN", RHDR, 0, 1, NO, "AMIN", 3},
	{PSPPK_AMAX, "AMAX", RHDR, 0, 1, NO, "AMAX", 3},
	{PSPPK_DEFAULT, "DEFAULT", RHDR, 0, 0, NO, "DEFAULT", 1},
	{0,	""		, IHDR, 0, 0, NO, "",-1}
};

static int psppkperplot = -1;
static int psppk_overlay = NO;
static int psppk_smooth = NO;
static int psppk_yn;
static int psppk_num;
static int psppk_ylin = NO;
static int psppk_xlin = NO;
static float psppk_fmin = -1;
static float psppk_fmax = 1.0e+38;
static float psppk_amin = 0.0;
static float psppk_amax = 1.0e+38;

char *yes_no_str [] = {
	"NO ", "YES" 
};


/* commands for interactive picking */
#define NEXT 0
#define QUIT 1
#define CONT 2
#define PREVIOUS 3
extern struct plt_ctl *pmap ;

void gsac_set_param_psppk(int ncmd, char **cmdstr);
void gsac_exec_psppk(void);
int gsac_show_intsppk(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double fwin,int numperframe, float yh, float yl, float depmin, float depmax);
int sppk_parse_key(float x0, float y0, float xlen, float ylen, int ns, int ne, float ts, float te, float dy, int ntrc, double fwin,int numperframe, float yh, float yl, float cx,float cy,int switchchar,int *lpos, float depmin, float depmax);
void showdec(float xl, float yl, float xh, float yh, int inc);
void do_smooth_5(float *x,int np);

extern int *sortptr;

extern void XviG_Flush();



extern void dogrid(void);

static char ostr1[80];

extern void gsac_plot_fileid(float x0,float y0,float xlen,float ylen, int k, int kmax);
extern void gsac_exec_hold(void);

extern void gsac_psppk(float x0,float y0,float xlen,float ylen,float fe,float fs,int ns,int ne,int ntrc, float depmax, float depmin, int numpsppkframe, float dy, int psppk_doamph, int psppk_overlay, int psppk_xlin, int psppk_ylin, float psppk_amin, float psppk_amax, float *uy, float *ly);
extern int inside(float cx, float cy, float xl, float yl, float xh, float yh);
extern void clearregion(float xl, float yl, float xh, float yh);

extern void traceinfo(float cx,float cy,float x0,float y0, float dy, float xlen,float ylen,int ns, int ne, int ntrc, int *k, float *ux, float *uy , int numperframe);

void psppk_getfreq(float uv, int k, int ns, int ne, int numpsppkframe, int *lpos, float x0);
void psppk_regr(float uv0, float uv1, int k, int ns, int ne, int numpsppkframe, int *lpos, float x0, 
	float y0, float xlen, float ylen,float dy, float depmax, float depmin);




/* these are temporary variables only used here */
float psppk_real[10];

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_psppk(int ncmd, char **cmdstr)
{
	/* the only parameter to be set is MORE 
	 *
	 */
	int i;
	int HasMouse;
	float XminDev, YminDev,
	        XmaxDev, YmaxDev, XminClip,
	        YminClip, XmaxClip, YmaxClip;
	int Color;

	float tmpmx, tmpmn;

	/* initialize graphics */
	if(gsac_control.plotinit == NO){
		if(gsac_control.plotdevice==WIN){
			ginitf("INTEM","GSAC");
			printf("Initializing Interactive Graphics\n");
			gmesg("Initializing Interactive Graphics");
			gsac_control.everinteractive = YES;
			gsac_control.plotinit = YES;
			gsac_control.plotchange = NO;

			ginfo(&HasMouse, &XminDev, &YminDev, 
				&XmaxDev, &YmaxDev, &XminClip, 
				&YminClip, &XmaxClip, &YmaxClip,&Color);
			if(Color >= 4)
				gsac_control.black = 0;
			else
				gsac_control.black = 1;
			gsac_control.kolor = Color%4;
		}
		gsac_control.plotinit = YES;
		gsac_control.plotchange = NO;
	} else {
		gframe(2);
		if(gsac_control.plotchange == YES){
			if(gsac_control.plotdevice==WIN){
				ginitf("INTEM","GSAC");
				gmesg("Initializing Interactive Graphics");
				gsac_control.everinteractive = YES;
				ginfo(&HasMouse, &XminDev, &YminDev, 
					&XmaxDev, &YmaxDev, &XminClip, 
					&YminClip, &XmaxClip, &YmaxClip,&Color);
				gsac_control.XmaxDev = XmaxDev;
				gsac_control.YmaxDev = YmaxDev;
				if(Color >= 4)
					gsac_control.black = 0;
				else
					gsac_control.black = 1;
				gsac_control.kolor = Color%4;
			}


		}
	}
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, psppkarg, YES, YES))
	       	return	;
	/* now determine which was set */
	for(i=0 ; psppkarg[i].key[0] != '\0' ; i++){
		if(psppkarg[i].used > 0){	
			if(psppkarg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, psppkarg[i].key, 
					psppkarg[i].mfit,psppkarg[i].narg, &psppk_yn );
			} else if(psppkarg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, psppkarg[i].key, 
					psppkarg[i].mfit,psppkarg[i].narg, &psppk_num );
			} else if(psppkarg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, psppkarg[i].key, 
					psppkarg[i].mfit,psppkarg[i].narg, psppk_real );
			}
			psppkarg[i].used = 0;
			switch(psppkarg[i].id){
				case PSPPK_PERPLOT:
					psppkperplot = psppk_num;
					break;
				case PSPPK_AM:
					psppk_doamph = PLOT_AM;
					break;
				case PSPPK_PH:
					psppk_doamph = PLOT_PH;
					break;
				case PSPPK_SMOOTH:
					if(psppk_yn == NO){
						psppk_smooth = NO;
					} else  if(psppk_yn == YES) {
						psppk_smooth = YES;
					}
					break;
				case PSPPK_OVERLAY:
					if(psppk_yn == NO){
						psppk_overlay = NO;
					} else if(psppk_yn == YES) {
						psppk_overlay = YES;
					}
					break;
				case PSPPK_XLOG:
					psppk_xlin = NO;
					break;
				case PSPPK_XLIN:
					psppk_xlin = YES;
					break;
				case PSPPK_YLOG:
					psppk_ylin = NO;
					break;
				case PSPPK_YLIN:
					psppk_ylin = YES;
					break;
				case PSPPK_FMIN:
					psppk_fmin = psppk_real[0];;
					break;
				case PSPPK_FMAX:
					psppk_fmax = psppk_real[0];;
					break;
				case PSPPK_AMIN:
					psppk_amin = psppk_real[0];;
					break;
				case PSPPK_AMAX:
					psppk_amax = psppk_real[0];;
					break;
				case PSPPK_DEFAULT:
					psppkperplot = -1;
					psppk_overlay = NO;
					psppk_smooth = NO;
					psppk_ylin = NO;
					psppk_xlin = NO;
					psppk_fmin = -1;
					psppk_fmax = 1.0e+38;
					psppk_amin = 0.0;
					psppk_amax = 1.0e+38;
					break;
				default:
					break;
			}
		}
	}
	/* safety check on fmax, fmin, amax and amin */
	tmpmx = MAX(psppk_fmax, psppk_fmin);
	tmpmn = MIN(psppk_fmax, psppk_fmin);
	psppk_fmax =tmpmx;
	psppk_fmin =tmpmn;
	tmpmx = MAX(psppk_amax, psppk_amin);
	tmpmn = MIN(psppk_amax, psppk_amin);
	psppk_amax =tmpmx;
	psppk_amin =tmpmn;

}

void gsac_exec_psppk(void)
{
int i,j,k,kkk;
int ntrc;
float x0, y0, xlen, ylen, dy;
float depmax, depmin;
int numpsppkframe;
int n2;
float temp;
float fs, fe;
double tr, ti;
int iret;

        double fwin;

        float uy, ly;


	/* if there are no traces return */
	ntrc = gsac_control.number_otraces;
	gsac_exec_hold();
	if(ntrc < 1)
		return;

	if(gsac_control.fft == NO){
		printf("Execute FFT first before plot\n");
		return;
	}
	/* initialize */
	/* must have interactive graphics */
	if(gsac_control.plotdevice!=WIN){
		printf("PLOTSPPK requires interactive graphics\n");
		return;
	}
	if(gsac_control.plotdevice==WIN){
		if(gsac_control.hold == NO){
			gframe(2);
		} else if(gsac_control.hold == 1){
			gframe(2);
			gsac_control.hold++ ;
		}
	}
	if(pmap == (struct plt_ctl *)NULL)
		pmap = (struct plt_ctl *)calloc(ntrc, sizeof(struct plt_ctl));
	else
		pmap = (struct plt_ctl *) realloc(pmap,ntrc*sizeof(struct plt_ctl));
	gsac_control.uxmul = 1.0;
	gsac_control.uxcen = 0.5;

	if(psppkperplot > 0)
		if(psppkperplot > ntrc)
			numpsppkframe = ntrc;
		else
			numpsppkframe = psppkperplot;
	else
		numpsppkframe = ntrc;

	xlen = gsac_control.xlen ;
	ylen = gsac_control.ylen ;
	x0   = gsac_control.x0 ;
	y0   = gsac_control.y0 ;

	if(psppk_overlay == YES){
		dy = ylen;
	} else {
		dy = ylen / numpsppkframe;
	}

	/* temporary */
	/* since all are plotted to the same frequency scale we must
	 * define the lower and upper limits. As a kludge, the minimum
	 * frequency plotted with log x-axis is the DF of first trace 
	 * in memory */
	/* get global fmax */
	fe = MIN(gsac_control.fmax, psppk_fmax);
	if(psppk_xlin == YES){
		fs = MAX (psppk_fmin, 0.0);
	} else {
		fs = MAX (psppk_fmin, sacdata[0].df);
	}

		save_depmin = (float *)realloc(save_depmin,(ntrc)*sizeof(float));
		save_depmax = (float *)realloc(save_depmax,(ntrc)*sizeof(float));

	/* get extreme amplitudes for the case of an overlay */
	if(psppk_overlay == YES ){
		/* never overlay phase spectra */
		if(psppk_doamph == PLOT_AM){
			depmax = 0.0;
			for ( kkk=0 ; kkk < ntrc ; kkk++){
				k = sortptr[kkk];
				n2 = sacdata[k].npow2/2.0  ;
				for(i=0 , j= 0; i <= n2 ; i++){
					tr = sacdata[k].sac_spectra[j++];
					ti = sacdata[k].sac_spectra[j++];
					temp = (float)sqrt(tr*tr + ti*ti);
					if(temp > depmax)
						depmax = temp;
				}
			}
			/* adjust amplitude for AMIN AMAX */
			if(psppk_ylin == YES){
				depmax = MIN(depmax * 1.2, psppk_amax);
				depmin = psppk_amin ;
			} else {
				if(psppk_amax < 1.0e+37)
					depmax = psppk_amax;
				else
					depmax *= 1.5;
				if(psppk_amin > 0.0)
					depmin = psppk_amin;
				else
					depmin = depmax/10000.0 ;
			}
		}
	}

        iret = NEXT;
	for ( kkk=0 ; kkk < ntrc && iret == NEXT ; kkk+=numpsppkframe){
		/* put up the traces */
		gsac_control.uxmul = 1.0;
		gsac_control.uxcen = 0.5;
		gframe(2);

		gsac_psppk(x0,y0,xlen,ylen,fe,fs,kkk,MIN(kkk+numpsppkframe,ntrc),
				ntrc,depmax,depmin,numpsppkframe, dy, 
				psppk_doamph,  psppk_overlay,  psppk_xlin,  psppk_ylin,  psppk_amin,  psppk_amax, &uy, &ly);
		/* do the interactive picking */
                fwin = fe - fs;
		iret = gsac_show_intsppk(x0, y0, xlen, ylen, kkk, 
				MIN(kkk+numpsppkframe,ntrc), fs, fe, dy, ntrc, 
				fwin, numpsppkframe,uy,ly,depmin,depmax);
		/* special case to handle previous */
		if(iret == PREVIOUS){
			kkk = MAX(-numpsppkframe, kkk-numpsppkframe-numpsppkframe);
			iret = NEXT ;
		}
	}	
	/* clean up */
	gclip("off", x0, y0, x0+xlen, y0+dy);
	gcursor("XORArrow");
}

int gsac_show_intsppk(float x0, float y0, float xlen, float ylen, int ns, int ne, float fs, float fe, float dy, int ntrc, double fwin,int numperframe, float yh, float yl, float depmin, float depmax)
{
float cx, cy;
char ch[2];
int switchchar;
int doloop;
int iret;
int lpos;

	/* work interactively with the current traces */
	/* set clip region to the actual plot region not the entire box */
	/* cross hair cursors within the trace plot region , default
	 * 	arrow outside this region */
	gcursor("Cross");
	gclip("on", x0, yl , x0+xlen, yh);
	doloop = YES;
	/* infinite loop */
	ch[1] = '\0';
	lpos = 0;
	while(doloop){
		gsac_control.ylim_rsc = YES;
		curaxy(&cx,&cy,ch);
	/* we must decide if the cursor was in the button region
	 * or if the cursor was in the trace region */
		switchchar = ch[0];
		switchchar = toupper(switchchar);
		if(inside(cx,cy,x0,y0,x0+xlen,y0+ylen)){
			iret = sppk_parse_key(x0, y0, xlen, ylen, ns, ne, 
				fs, fe, dy, ntrc, fwin,numperframe, yh, yl,
				cx,cy,switchchar,&lpos, depmin, depmax);
			if(iret == NEXT)
				return (NEXT);
			else if(iret == PREVIOUS)
				return (PREVIOUS);
			else if(iret == QUIT)
				return (QUIT);
		}
	}
	gsac_control.ylim_rsc = NO;
	return(CONT);
}

int sppk_parse_key(float x0, float y0, float xlen, float ylen, int ns, int ne, float fs, float fe, float dy, int ntrc, double fwin,int numpsppkframe, float yh, float yl, float cx,float cy,int switchchar,int *lpos, float depmin, float depmax) 
{
/* remember that switchchar is upper case from the line before the call */
int k;
char ch[2];
float ux, uy;
float ux0, ux1;
float uv;
float cx_X, twin_X;	/* for use with the X X to define window */
	traceinfo(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,ntrc,&k,&ux,&uy ,numpsppkframe);
	cx_X = -1;
	if(k >= 0)
	switch(switchchar){
		case '-':
		case '_':
			/* compress frequency scale */
			gsac_control.uxcen += (ux - 0.5)/gsac_control.uxmul;
			gsac_control.uxmul /= 2.0 ;
			gsac_control.uxmul = MAX(1.0, gsac_control.uxmul/2.0);
			if(gsac_control.uxmul == 1.0)
				gsac_control.uxcen = 0.5;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			gframe(2);
			clearregion(0.0,1.0,10.0,8.0);
		gsac_psppk(x0,y0,xlen,ylen,fe,fs,ns,ne,
				ntrc,depmax,depmin,numpsppkframe, dy, 
				psppk_doamph,  psppk_overlay,  psppk_xlin,  psppk_ylin,  psppk_amin,  psppk_amax, &yh, &yl);
			gclip("on", x0, yl , x0+xlen, yh);
/*
*/
			break;
		case '+':
		case '=':
			/* expand frequency scale */
			gsac_control.uxcen += (ux - 0.5)/gsac_control.uxmul;
			gsac_control.uxmul *= 2.0 ;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			gframe(2);
			clearregion(0.0,1.0,10.0,8.0);
		gsac_psppk(x0,y0,xlen,ylen,fe,fs,ns,ne,
				ntrc,depmax,depmin,numpsppkframe, dy, 
				psppk_doamph,  psppk_overlay,  psppk_xlin,  psppk_ylin,  psppk_amin,  psppk_amax, &yh, &yl);
			gclip("on", x0, yl , x0+xlen, yh);
/*
*/
			break;
		case ' ':
			/* center cursor point */
			gsac_control.uxcen += (ux - 0.5)/gsac_control.uxmul;
			*lpos = 0;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			gframe(2);
			clearregion(0.0,1.0,10.0,8.0);
		gsac_psppk(x0,y0,xlen,ylen,fe,fs,ns,ne,
				ntrc,depmax,depmin,numpsppkframe, dy, 
				psppk_doamph,  psppk_overlay,  psppk_xlin,  psppk_ylin,  psppk_amin,  psppk_amax, &yh, &yl);
			gclip("on", x0, yl , x0+xlen, yh);
/*
*/
			break;
			break;
		case 'B': /* display previous page  */
			  return PREVIOUS;
			  break;
		case 'L':        /* give cursor location and amplitude
		 			note use a small unclip window */
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			/* note:  ux is the position in the window. psppk_getfreq uses this position
				since it also knows the minimum and maximum frequency of the window
				from the pmap[].tfirst and pmap[].tlast which are set by the plotting
				routine 
				On the other hand the uv value computed for 'space' and 'X'-'X' refers to 
				the original window and not the modified window 
			*/
			  psppk_getfreq(ux,k,ns,ne,numpsppkframe,lpos,x0);
			gclip("on", x0, yl , x0+xlen, yh);
			  break;
		case 'N':  /* next plot of PERPLOT is set  */
			  return NEXT;
			  break;
		case 'O': /* previous plot window  */
			  /* this requires a history of windows
			   * and a replot which requires turning of
			   * clipping while a redraw is done */
			*lpos = 0;
			gsac_control.uxmul = 1.0;
			gsac_control.uxcen = 0.5;
			gclip("off", x0, y0, x0+xlen, y0+ylen);
			gframe(2);
			clearregion(0.0,1.0,10.0,8.0);
			gsac_psppk(x0,y0,xlen,ylen,fe,fs,ns,ne,
				ntrc,depmax,depmin,numpsppkframe, dy, 
				psppk_doamph,  psppk_overlay,  psppk_xlin,  psppk_ylin,  
				psppk_amin,  psppk_amax, &yh, &yl);
			gclip("on", x0, yl , x0+xlen, yh);
			break;
		case 'Q': /* end interactive  */

			return QUIT;
		case 'X': 
			cx_X = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
			/* now invoke the cursor again and only use
			      * of we see a second X */
			curaxy(&cx,&cy,ch);
			switchchar = ch[0];
			switchchar = toupper(switchchar);
			if(inside(cx,cy,x0,y0,x0+xlen,y0+ylen)){
				traceinfo(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					ntrc,&k,&ux,&uy ,numpsppkframe);
			if(switchchar == 'X'){
				uv = gsac_control.uxcen + (ux - 0.5)/gsac_control.uxmul;
				/* re map */
				gsac_control.uxcen = (uv + cx_X)/2.0;
				/* redisplay */
				twin_X=ABS(uv - cx_X)*fwin;
				gsac_control.uxmul = ABS(fwin/twin_X);
				*lpos = 0;
				gclip("off", x0, y0, x0+xlen, y0+ylen);
				clearregion(0.0,0.0,10.0,8.0);
				gsac_psppk(x0,y0,xlen,ylen,fe,fs,ns,ne,
						ntrc,depmax,depmin,numpsppkframe, dy, 
						psppk_doamph,  psppk_overlay,  psppk_xlin,  
						psppk_ylin,  psppk_amin,  psppk_amax, &yh, &yl);
				gclip("on", x0, yl , x0+xlen, yh);
				}
			}
			break; /* invoke twice to define new window */
		case 'Z':        /* give cursor location and amplitude
		 			note use a small unclip window */
			ux0 = ux;
			/* now get the second frequency */
			curaxy(&cx,&cy,ch);
			switchchar = ch[0];
			switchchar = toupper(switchchar);
			if(inside(cx,cy,x0,y0,x0+xlen,y0+ylen)){
				traceinfo(cx,cy,x0,y0,dy,xlen,ylen,ns,ne,
					ntrc,&k,&ux1,&uy ,numpsppkframe);
				if(switchchar == 'Z'){
					/* now we have the two points - do the regression */
					/* permit writing at the top ofthw window */
					gclip("off", x0, y0, x0+xlen, y0+ylen);
			  		psppk_regr(ux0,ux1,k,ns,ne,numpsppkframe,lpos,x0,y0,xlen,ylen,dy,depmax, depmin);
					/* turn off writing at the top of the window */
					gclip("on", x0, y0 , x0+xlen, y0+ylen);
				}
			}
			break;
	}
	return (-1);
}

void gsac_psppk(float x0,float y0,float xlen,float ylen,float fe,float fs,int ns,int ne,int ntrc, float depmax, float depmin, int numpspframe, float dy, int psppk_doamph, int psppk_overlay, int psppk_xlin, int psppk_ylin, float psppk_amin, float psppk_amax, float *uy, float *ly)
{
int kk, kkk, n2;
int i, j, k, is;
double tr, ti;
float depmen;
	int indmax, indmin;
float df;
float yy0, yyy0;
float xx, yy;
float ht;
float freq;
float uxcen, uxmul, fwin;
float fwinl, fwinh;
	/* do the bottom axis */

	ht = 0.10;
	yy0 = 1.0;
	gclip("off", x0, yy0, x0+xlen, yy0+dy);
	uxcen = gsac_control.uxcen;
	uxmul = gsac_control.uxmul;
	*ly =  1.0e+38;
	*uy = -1.0e+38;
	/* ssafety */
	/* linear mapping */
	if(psppk_xlin == YES){ 
		fwin = fe -fs; 
		fwinl = fs + fwin*(uxcen - 0.5/uxmul);
		fwinh = fs + fwin*(uxcen + 0.5/uxmul);
	} else {
		fwinl = fs * pow(fe/fs,(uxcen - 0.5/uxmul) );
		fwinh = fs * pow(fe/fs,(uxcen + 0.5/uxmul) );
	}
	fs = fwinl;
	fe = fwinh;

	/* set the background */
	if(gsac_control.background == YES &&
		gsac_control.background_color >= 0){
		newpen(gsac_control.background_color);
		if(psppk_overlay == YES){
			shader(x0,y0,x0+xlen,y0+ylen,0,0,0.01,0.01);
		} else {
			shader(x0,y0+ylen-(ne-ns)*dy,x0+xlen,y0+ylen,0,0,0.01,0.01);
		}
		newpen(1);
	} else {
		newpen(1);
	}
	for(kkk=ns;kkk < ne;kkk++){
		k = sortptr[kkk];
		kk = kkk%numpspframe;

		n2 = sacdata[k].npow2/2  ;
		y = (float *)realloc(y,(n2+1)*sizeof(float));
		/* now fill array with amplitude or phase spectra */
		for(i=0 , j= 0; i <= n2 ; i++){
			tr = sacdata[k].sac_spectra[j++];
			ti = sacdata[k].sac_spectra[j++];
			if(psppk_doamph == PLOT_AM){
				y[i] = sqrt(tr*tr + ti*ti);
			} else if (psppk_doamph == PLOT_PH){
				y[i] = atan2(ti,tr);
			}
		}
		if(psppk_doamph == PLOT_AM && psppk_overlay == NO){
			getmxmn(y, n2+1,&depmax, &depmin, &depmen,&indmax,&indmin);
			/* adjust amplitude for AMIN AMAX */
			if(psppk_ylin == YES){
				depmax = MIN(depmax * 1.1, psppk_amax);
				depmin = psppk_amin ;
			} else {
				if(psppk_amax < 1.0e+37)
					depmax = psppk_amax;
				else
					depmax *= 1.5;
				if(psppk_amin > 0.0)
					depmin = psppk_amin;
				else
					depmin = depmax/10000.0 ;
			}
		} else if (psppk_doamph == PLOT_PH){
			depmax =  3.1415927;
			depmin = -3.1415927;
			depmen = 0.0;
		}

		if(psppk_smooth == YES)
			do_smooth_5(y,n2);

		df  = sacdata[k].df;

		/* use depmax depmin */
		if(psppk_overlay == YES){
			yy0 = y0;
			yyy0 = yy0 -kk*1.5*ht ;
		} else {
			yy0 = y0 + (numpspframe -1  - kk )*dy;
			yyy0 = yy0;
		}
		if(psppk_xlin == YES){
			if(gsac_control.xgrid == YES)
				dolnxgrid(x0,yy0,yy0+dy,xlen,fe,fs,0.10, YES, 
					gsac_control.xgrid_color, 
					gsac_control.xgrid_type,
					gsac_control.xgrid_minor);
			if(gsac_control.background == YES &&
				gsac_control.background_color >= 0){
				newpen(0);
				dolinx(x0,yy0+dy,xlen,fe,fs,0.10, NO, NO, NO, 0, " ");
				newpen(1);
			} else {
				dolinx(x0,yy0+dy,xlen,fe,fs,0.10, NO, NO, NO, 0, " ");
			}
		} else {
			if(gsac_control.xgrid == YES){
				dologxgrid(x0,yy0,yy0+dy,xlen,fe,fs,0.10,
					gsac_control.xgrid_color,
					gsac_control.xgrid_type,
					gsac_control.xgrid_minor);
			}
			if(gsac_control.background == YES &&
				gsac_control.background_color >= 0){
				newpen(0);
				dologx(x0,yy0+dy,xlen,fe,fs,0.10, NO, NO, NO, 0, " ");
				newpen(1);
			} else {
				dologx(x0,yy0+dy,xlen,fe,fs,0.10, NO, NO, NO, 0, " ");
			}
		}
		if(kkk == ne -1 ){
			if(psppk_xlin == YES) {
				dolinx(x0,yy0,xlen,fe,fs,0.10, YES, NO, YES, 14, "Frequency (Hz)");
				if(gsac_control.background == YES &&
					gsac_control.background_color >= 0){
					newpen(0);
					dolinx(x0,yy0,xlen,fe,fs,0.10, YES, NO, NO, 0, "");
					newpen(1);
				}
			} else {
				dologx(x0,yy0,xlen,fe,fs,0.10,YES,NO,YES,14,"Frequency (Hz)");
				if(gsac_control.background == YES &&
					gsac_control.background_color >= 0){
					newpen(0);
					dologx(x0,yy0,xlen,fe,fs,0.10,YES,NO,NO,0,"");
					newpen(1);
				}
			}
		}
		if(psppk_doamph == PLOT_AM){
			if(psppk_ylin == YES){
				if(gsac_control.ygrid == YES)
					dolnygrid(x0,x0+xlen,yy0,dy,depmax,
						depmin, 0.10, YES, 
						gsac_control.ygrid_color, 
						gsac_control.ygrid_type,
						gsac_control.ygrid_minor);
				doliny(x0,yy0,dy,depmax,depmin,0.10, NO, YES, YES, 1, " ");
				if(gsac_control.background == YES &&
					gsac_control.background_color >= 0){
					newpen(0);
					doliny(x0     ,yy0,dy,depmax,depmin,0.10, NO, YES, NO, 1, " ");
					doliny(x0+xlen,yy0,dy,depmax,depmin,0.10, YES, NO, NO, 0, " ");
					newpen(1);
				}
			} else {
				if(gsac_control.ygrid == YES){
					dologygrid(x0,x0+xlen,yy0,dy,depmax,
						depmin,0.10,
						gsac_control.ygrid_color,
						gsac_control.ygrid_type,
						gsac_control.ygrid_minor);
				}
				dology(x0,yy0,dy,depmax,depmin,0.10,NO,YES,YES,1," ");
				if(gsac_control.background == YES &&
					gsac_control.background_color >= 0){
					newpen(0);
				        dology(x0     ,yy0,dy,depmax,depmin,0.10,NO,YES,NO,1," ");
					dology(x0+xlen,yy0,dy,depmax,depmin,0.10,YES,NO,NO,0," ");
					newpen(1);
				}
			}
		} else {
			if(gsac_control.ygrid == YES){
				dolnygrid(x0,x0+xlen,yy0,dy,depmax,depmin,0.10, 
					YES, gsac_control.ygrid_color, 
					gsac_control.ygrid_type,
					gsac_control.ygrid_minor);
			}
			doliny(x0,yy0,dy,depmax,depmin,0.10, NO, YES, YES, 1, " ");
			doliny(x0+xlen,yy0,dy,depmax,depmin,0.10, YES, NO, NO, 0, " ");
		}
		gbox(x0, yy0, x0+xlen, yy0+dy);
		gclip("on", x0, yy0, x0+xlen, yy0+dy);
		gsac_setcolor(YES, kkk, ntrc);
		gsac_plot_fileid(x0,yyy0,xlen,dy, k, ne - ns + 1);
		/* plot the spectra */
		pmap[kk].ifirst = -1;
		pmap[kk].ilast  = -1;
		for(i=0, is=0 ; i <= n2 ; i++){
			freq = i * df;
			if(freq >= fs && freq <= fe){
				if(pmap[kk].ifirst < 0 )
					pmap[kk].ifirst = i;
				pmap[kk].ilast  = i;
				if(psppk_doamph == PLOT_AM){
					if(psppk_ylin == YES){
				yy = yy0 +dy*(y[i] - depmin)/(depmax - depmin);
					} else {
				yy = yy0 + dy*log10(y[i] / depmin)/log10(depmax / depmin);
					}
				} else {
				yy = yy0 +dy*(y[i] - depmin)/(depmax - depmin);
				}
				if(psppk_xlin == YES){
/*
				u = (freq-fs)/(fe-fs);
				xx = x0 + xlen*(uxmul*(u-uxcen) +0.5 );
*/
				xx = x0 + xlen*(freq-fs)/(fe-fs);
				} else {
/*
				u = log10(freq/fs)/log10(fe/fs);
				xx = x0 + xlen*(uxmul*(u-uxcen) +0.5 );
*/
				xx = x0 + xlen*log10(freq/fs)/log10(fe/fs);
				}
	
		
				if(is == 0 ){
					plot(xx,yy,3);
					is = 1;
				} else {
					plot(xx,yy,2);
				}
			}
		}
save_depmin[k] = depmin;
save_depmax[k] = depmax;
		gsac_setcolor(NO, kkk, ntrc);
		gclip("off", x0, yy0, x0+xlen, yy0+dy);
		/* set up the plot bounds to define subwindows */
		pmap[kk].abstime = YES;
		pmap[kk].xl = x0;
		pmap[kk].xh = x0 +xlen;;
		pmap[kk].yl = yy0;
		pmap[kk].yh = yy0 +dy;
		pmap[kk].k  = k;
		pmap[kk].n  = kk;
		pmap[kk].tfirst = fs;
		pmap[kk].tlast  = fe;
		pmap[kk].npts  = n2-1;
		pmap[kk].xlen = xlen;
		pmap[kk].delta = df;
		/* not used but also float = double */
		pmap[kk].ymult = 1.0;
                if(*uy < pmap[kk].yh)
                        *uy = pmap[kk].yh ;
                if(*ly > pmap[kk].yl)
                        *ly = pmap[kk].yl ;
		XviG_Flush();
	}
		/* annotate with the plot titles
			added 29 MAY 2009 */
}

void psppk_getfreq(float uv, int kk, int ns, int ne, int numpsppkframe, int *lpos, float x0){
/* uv - fraction of window where frequency is picked varies [0,1]
   kk - reference to the frame on the page kk = [0,numpsppkframe-1] except for
	possibly the last page
   ns - begin  of traces on this page
   ne - end index
   numpsppkframe = number of spectra per page
   lpos - controls output of amplitude information at top of page
   x0  - position of lwoer left corner
*/
	/* k is the fraction of the window displayed [0,1] 
	   kk is the reference to the particular frame shown
	   k is the original trace read, before sorted
	*/
	   
	float dfw, fs, fe, df, amp;
	float freq ;
	double tr, ti;
	int k, npow2,i,j,kkks,kkke,kkk;

	/* first set the limits for output */
	/* given the frequency, make into array index */
	if(psppk_overlay == YES ){
		kkks = ns ;
		kkke = ne ;
	} else {
		kkks = kk;
		kkke = kk+1;
	}
	/* originally I just computed the index and then accessed the
		raw spectra. However, since smoothing was implemented, it is
	  	appropriate to compute the entire array. This takes time but, this
		also serves as a prototype for the t* estimation code
	*/
	for(kkk=kkks; kkk < kkke ; kkk++){
		k = pmap[kkk].k ;
		npow2 = sacdata[k].npow2 ;
		y = (float *)realloc(y,(npow2/2+1)*sizeof(float));
		/* now fill array with amplitude or phase spectra */
		for(i=0 , j= 0; i <= npow2/2 ; i++){
			tr = sacdata[k].sac_spectra[j++];
			ti = sacdata[k].sac_spectra[j++];
			y[i] = (float)sqrt(tr*tr + ti*ti);
		}
		if(psppk_smooth == YES)
			do_smooth_5(y,npow2/2);
		df = sacdata[k].df ;
		fs = pmap[kk].tfirst ;
		fe = pmap[kk].tlast  ;
        	dfw = fe - fs;
		if( psppk_xlin == YES )
			freq = fs + uv*dfw ;
		else
			freq = fs * pow( (fe/fs), uv);
		/* freq = j *df where j = [0,npow2] */
		j = (int)(0.5 + freq / df );	
		if( j < 0 )
			j = 0;
		if(j <= npow2/2) {
			amp = y[j];
			/* beware ostr is of length 80 or 79 + 1 there is no limit on the file name */
			sprintf(ostr1,"%f %11.4e %s SMOOTH %s",freq,amp,
				basename(sacdata[k].sac_ifile_name),yes_no_str[psppk_smooth]) ;
			printf("%s\n",ostr1);
	  		gleft(x0,7.8-(*lpos)*0.15,0.10,ostr1,0.0);
			*lpos = *lpos + 1;
		}
	}
}


void psppk_regr(float uv0, float uv1, int kk, int ns, int ne, int numpsppkframe, int *lpos, float x0, 
	float y0, float xlen, float ylen,float dy,float depmax, float depmin)
{
	/* this is essentially psppk_getfreq code except that
		we actually do a regression */
/* do a linear regression of the amplitude spectrum of the form
	ln y = a + b f
*/
/* uv0 - fraction of window where first  frequency is picked varies [0,1]
   uv1 - fraction of window where second frequency is picked varies [0,1]
   kk - reference to the frame on the page kk = [0,numpsppkframe-1] except for
	possibly the last page
   ns - begin  of traces on this page
   ne - end index
   numpsppkframe = number of spectra per page
   lpos - controls output of amplitude information at top of page
   x0  - position of lower left corner
*/
/* k is the fraction of the window displayed [0,1] 
   kk is the reference to the particular frame shown
   k is the original trace read, before sorted
*/
	   
	float dfw, fs, fe, df, tmp;
	float freq, freq0, freq1 ;
	float amp ;
	float yy0 ;
	float vx0, vy0;
	double tr, ti;
	float pval;
	int k,npow2,i,j,j0, j1,kkks,kkke,kkk;
	double sumn, sumx, sumxx, sumy, sumyy, sumxy;
	double aval, bval, det;
	double tstar;
	int ndots;


	/* insure that uv0 < uv1 */
	if(uv1 < uv0){
		tmp = uv0 ; uv0 = uv1 ; uv1 = tmp ;
	}

	/* first set the limits for output */
	/* given the frequency, make into array index */
	if(psppk_overlay == YES ){
		yy0 = y0;
		kkks = ns ;
		kkke = ne ;
	} else {
		yy0 = y0 + (numpsppkframe -1  - kk )*dy;
		kkks = kk;
		kkke = kk+1;
	}
	/* originally I just computed the index and then accessed the
		raw spectra. However, since smoothing was implemented, it is
	  	appropriate to compute the entire array. This takes time but, this
		also serves as a prototype for the t* estimation code
	*/
	for(kkk=kkks; kkk < kkke ; kkk++){
		gsac_setcolor(YES, kkk,  gsac_control.number_otraces );
		k = pmap[kkk].k ;
		npow2 = sacdata[k].npow2 ;
		y = (float *)realloc(y,(npow2/2+1)*sizeof(float));
		/* now fill array with amplitude or phase spectra */
		for(i=0 , j= 0; i <= npow2/2 ; i++){
			tr = sacdata[k].sac_spectra[j++];
			ti = sacdata[k].sac_spectra[j++];
			y[i] = (float)sqrt(tr*tr + ti*ti);
		}
		if(psppk_smooth == YES)
			do_smooth_5(y,npow2/2);
		df = sacdata[k].df ;
		fs = pmap[kk].tfirst ;
		fe = pmap[kk].tlast  ;
        	dfw = fe - fs;
		/* get indices into the array */
		if( psppk_xlin == YES )
			freq0 = fs + uv0*dfw ;
		else
			freq0 = fs * pow( (fe/fs), uv0);
		/* freq = j *df where j = [0,npow2] */
		j0 = (int)(0.5 + freq0 / df );	
		if( psppk_xlin == YES )
			freq1 = fs + uv1*dfw ;
		else
			freq1 = fs * pow( (fe/fs), uv1);
		/* freq = j *df where j = [0,npow2] */
		j1 = (int)(0.5 + freq1 / df );	
		/* initialize */
		sumn = 0 ;
		sumx = 0.0 ;
		sumxx = 0.0 ;
		sumy = 0.0 ;
		sumyy = 0.0 ;
		sumxy = 0.0 ;
		/* safety */

		if( j0 < 0 ) j0 = 0;
		if( j1 < 0 ) j1 = 0;
		if( j0 <= npow2/2 && j1 <= npow2/2 ){
			/* get the partial sums for the regression
			   also replot the selected line segment with a wider line
			   to indicate the line used 
				There is a log of redundant mapping code here
			gwidth(0.05);
			*/
			for(i=j0 ; i <= j1 ; i++){
				if(ABS(y[i] > 0.0) ){
					freq = i*df;
					sumn++ ;
					sumx  += freq;
					sumxx += freq*freq;
					sumxy += freq*log(y[i]);
					sumyy += log(y[i])*log(y[i]);
					sumy  += log(y[i]);
				}
			/*
				if(psppk_xlin == YES){
					vx0 = x0 + xlen*(freq-fs)/(fe-fs);
				} else {
					vx0 = x0 + xlen*log10(freq/fs)/log10(fe/fs);
				}
				if(psppk_ylin == YES){
					vy0 = yy0 +dy*(y[i] - depmin)/(depmax - depmin);
				} else {
					vy0 = yy0 + dy*log10(y[i] / depmin)/log10(depmax / depmin);
				}
				if(i == j0) {
					plot(vx0,vy0,3);
				} else {
					plot(vx0,vy0,2);
				}
			*/
			}
			if(sumn > 0.0){
				det = sumn*sumxx - sumx * sumx ;
				aval = ( sumxx * sumy - sumx * sumxy)/det;
				bval = (-sumx  * sumy + sumn * sumxy)/det;
				tstar = - bval/3.1415927;
				sprintf(ostr1,"lnA = a + bf : a=%f b=%f t*=%f  %s SMOOTH %s",aval,bval,tstar,
				basename(sacdata[k].sac_ifile_name),yes_no_str[psppk_smooth]) ;
				printf("%s\n",ostr1);
				sprintf(ostr1,"t*=%f  %s SMOOTH %s",tstar,
				basename(sacdata[k].sac_ifile_name),yes_no_str[psppk_smooth]) ;
				/* permit writing at top of frame */
				gclip("off", x0, y0, x0+xlen, y0+ylen);
	  			gleft(x0,7.8-(*lpos)*0.15,0.10,ostr1,0.0);
				*lpos = *lpos + 1;
				/* now make the predictions */
				/* now plot in the current frame */
				/* do not permit writing at top of frame */
				gclip("on",x0,pmap[kk].yl, x0+xlen,pmap[kk].yh)  ;
				depmax = save_depmax[k] ;
				depmin = save_depmin[k] ;
				/* plot the regression line. Even though it takes two points to
				draw the line with a linear x axis, for a log x-axis we
				need more points */
				ndots = MAX(2, ABS( xlen*(uv1-uv0)/(2 * 0.10)) );
fprintf(stderr, "xlen %f uv1-uv0 %f ndots %d %f\n",xlen,uv1-uv0,ndots,xlen*(uv1-uv0));
				
				for ( i = 0 ; i < ndots ; i++){
					/* define the frequency */
						pval = (float)(i)/(float)(ndots -1);
					if(psppk_xlin == YES){
						freq = (1.-pval)*freq0 + pval*freq1 ;
						vx0 = x0 + xlen*(freq-fs)/(fe-fs);
					} else {
						freq = freq0*pow(10.0,pval*log10(freq1/freq0)) ;
						vx0 = x0 + xlen*log10(freq/fs)/log10(fe/fs);
					}
						amp = exp( aval + bval * freq ) ;
					if(psppk_ylin == YES){
						vy0 = yy0 +dy*(amp - depmin)/(depmax - depmin);
					} else {
						vy0 = yy0 + dy*log10(amp / depmin)/log10(depmax / depmin);
					}
					fillit("CI",0.05,vx0,vy0);
					newpen(gsac_control.background_color);
/* think carefully about color of outline
   if background = white (0) use 1
   if background = black (1) use 0
   else use 0
*/
if(gsac_control.background_color == 0)
	newpen(0);
else if(gsac_control.background_color == 1)
	newpen(0);
else
	newpen(1);
					curvit("CI",0.06,vx0,vy0);
					gsac_setcolor(YES, kkk,  gsac_control.number_otraces );
				}
			}
			gsac_setcolor(NO, kkk,  gsac_control.number_otraces );
		}
		newpen(1);
	}
}

void do_smooth_5(float *y,  int np){
/* apply a five point smoother */
int m;
	x = (float *)realloc(y,(np+1)*sizeof(float));
	for(m=0;m < np ; m++)
		x[m] = y[m];
        /* special case first 5 points */
        for(m=0;m<2;m++){
                y[m] = (x[m] + x[m+1] + x[m+2] )/3.0;
        }
        /* smoothing interior points */
        for(m=2; m<np -2 ;m++){
                y[m] = (x[m-2] + x[m-1] + x[m] + x[m+1] + x[m+2])/5.0;
        }
        /* special case last 5 points */
        for(m=np-2; m< np;m++){
                y[m] = (x[m-2] +x[m-1] + x[m])/3.0;
        }
}
