/* to do
 * 		}
	for MATCH
	Add, Delete
		Add is a fictitious that must be removed on return
		Delete can delete a previous pick, but also a fictittious
	Mode should highlight the current selection - e.g., XOR
	The zoom box should always plot the proper period at the bottom,
		or should just zoom on the period only not on velocity
	If zooom stat set when tomath reset

	for zoom and interp
		use zoom to define the period limit
	go through the list and select, sort observations then interpolate
	use pointer? since many more -- now for match only have
	one mode to work with so can easily define simple output array



  03 16 00 - mode select implemented - now for do pick we should
	use the mode value and not the original amplitude value
	also select on from already picked instead of from 
        nearest neighbor - also note in the matchprep we use
        instantaneous period not actual period

  10 10 00 - new output STAcomp.ods

  12 27 00 - minor changes to colors
  
  01 30 02 - add AZ to the *.dsp output

  08 23 03 - changed name to mft96.DSP from mft96.dsp for easier copy later
  09 25 03 - fixed improper placement of free(period) and free(mftdsp)
  01 05 04 - error in line 1465 forintf pathname before defined delete 1465
  06 30 05 - changed name to mft96.disp from mft96.DSP for easier copy later 
		under caseless file systems (CYGWIN under windows)
  06 10 06 - permit edited dispersion to be called fname.dsp to facilitate
		interstation green's functions: do_mft -G
  01 14 12 - changed position of the Yes/No Save as dialog
  10 13 13 - added phase velocity from noise cross-correlation
  11 05 13 - put many declarations into do_mft.h and common routines into
             mft_subs.c
*/

#include	<stdio.h>
#include	<ctype.h>
#include	<stdlib.h>
#include	<string.h>
#include	"calplot.h"
#include	<math.h>
#include	"grphsubc.h"
#include	"nmenu.h"
#include	"do_mft.h"

#include	<libgen.h>



FILE *cerr;



struct pos mftpos[3];
struct pos phvpos[3];


#define	ZOOM	1
#define	UNZOOM	2
#define	AUTO	3
#define	MODE	4
#define	PICK	5
#define	RESTART	6
#define	MATCH	7
#define	EXIT	8
#define INTERP	9
#define ADD	10
#define RETURN	11
#define DELETE	12
#define OVERLAYTOMOGPV 13
#define AUTOFEED 14
#define PHASEVELOCITY 15
#define OVERLAYTOMOPHV 16

static struct menu m[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Zoom\0"     , ZOOM   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "UnZoom\0"   , UNZOOM ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Mode\0"     , MODE   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Auto\0"     , AUTO   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Pick\0"     , PICK   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Restart\0"  , RESTART,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Match\0"    , MATCH  ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Exit\0"     , EXIT   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Autofeed\0" , AUTOFEED   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "PhVel\0"    , PHASEVELOCITY   ,-1, 0, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Tomo\0"     , OVERLAYTOMOGPV   ,-1, 0, 1},
};

static struct menu mm[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Zoom\0"   , ZOOM   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "UnZoom\0" , UNZOOM ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Mode\0"   , MODE   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Interp\0" , INTERP ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Add\0"    , ADD    ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Delete\0" , DELETE ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Return\0" , RETURN, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Exit\0"   , EXIT   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Match\0"  , MATCH  ,-1, 1, 1},
};

static struct menu mp[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Zoom\0"     , ZOOM   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "UnZoom\0"   , UNZOOM ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Mode\0"     , MODE   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Auto\0"     , AUTO   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Pick\0"     , PICK   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Restart\0"  , RESTART,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Return\0" , RETURN, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Save\0"     , EXIT   ,-1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "Tomo\0"     , OVERLAYTOMOPHV   ,-1, 0, 1},
};


static struct menu my[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Yes\0" , 1, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "No \0" , 2, -1, 1, 1},
};

#define PHVRESTART 100
/* GLOBAL VARIABLE DECLARATIONS */
	/* display information */
extern int HasMouse; 
extern float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
extern int Color;
extern int doshade;

extern int Units;
extern float border;
extern float title;
extern int black, kolor;
extern int Button_Color;
extern int Button_Color_Light;
extern int Button_Color_Dark;
extern int Button_Color_Fore;
extern int Button_Color_Back;

extern int interstationgreen ;
extern int Cursor;
extern int overlaytomography ;
extern int do_auto_feed;
extern int do_phvel;
extern int AutoMode;
extern int Mode  ;
extern int iMode ;	/* pointer into mode array */

int InPick = OFF;
int InLine = OFF;
int XAxisPeriod = NO;
int ModeSelect  = -1;	/* for use in interpolation selection */
char ostr[280];
char ostr1[280];
struct disp *mftdsp;
struct disp *phvdsp;
int num_mft_disp;
int num_phv_disp;
/* arrays for amplitude and dispersion plots */
/* we will place 3 graphs in the page :
	[0]  spectra, [1] group velocity, [2] trace
	The positions of the original plot file will be mapped
	on to the following boxes:
*/
/* screenplacement */
/*
float  	pxl[3]= {0.7, 5.3, 9.4},
	pyl[3]= {3.0, 3.0, 3.0},  
	pxh[3]= {4.7, 9.3, 9.9},  
	pyh[3]= {7.0, 7.0, 7.0};	
float  	pxl[3]= { 0.7, 4.3, 9.4},
	pxh[3]= { 3.7, 9.3, 9.9},
	pyl[3]= { 2.0, 2.0, 2.0},
	pyh[3]= { 7.0, 7.0, 7.0};
*/
/* placement for spectral amplitude, group velocity and trace
	for MFT nd PMF  */
float  	upxl[3]= { 0.5, 4.1, 9.4},
	upxh[3]= { 3.5, 9.1, 9.9},
	upyl[3]= { 2.0, 2.0, 2.0},
	upyh[3]= { 7.0, 7.0, 7.0};
float   uscl[3], usclx[3], uscly[3];
/* placement for group velocity and phase velocity 
	for PHV
	only the first two entries are used */
float  	ppxl[3]= { 0.6, 5.5, 9.4},
	ppxh[3]= { 4.9, 9.8, 9.9},
	ppyl[3]= { 2.0, 2.0, 2.0},
	ppyh[3]= { 7.0, 7.0, 7.0};
float   pscl[3], psclx[3], pscly[3];

float pxl[3], pxh[3], pyl[3], pyh[3];
float apxl[3], apyl[3], apxh[3], apyh[3];	/* user values of corners */

float opxl[3], opyl[3], opxh[3], opyh[3];	/* original cut corners */
float sclx[3], scly[3], scl[3]		;	/* scale factors */

int xlnlg[3], ylnlg[3];
int uxlnlg[3], uylnlg[3];
int pxlnlg[3], pylnlg[3];

float *period;
int num_period;
int Pick_Color;




/* FUNCTION DECLARATIONS */
extern void gsolid(float x, float y, float h, int symb);
extern char *my_pathfind( char *path, const char *name, const char *mode);
void do_gpv_pick(void);
void get_mft_output(void );
void get_phv_output(void );
void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
float *xup, float *yup, int butrev, int lstrmax) ;
void do_scale(float per,float vel,float *xval,float *yval,int i);
void do_invscale(float *per,float *vel,float xval,float yval,int i);
int get_action_mft(int nm, char *fname);
int get_action_pmf(int nm, char *fname);
int get_action_phv(int nm, char *fname);
void do_pick( float xv, float yv, int lev, int num_mft_disp, struct disp *mftdsp, int porg);
void do_line( float xv, float yv, int num_mft_disp, struct disp *mftdsp, int porg);
void do_zoom(int num_disp, struct disp *dsp, struct pos *pos, int porg);
void do_zoomm( void);
void do_clear( int fig);
void do_unzoom(int num_disp, struct disp *dsp,  struct pos *pos, int porg);
void do_unzoomm( void);
void do_resetzoom(void);
int do_mode( void);
void do_showauto();
void do_showpick1(int num_mft_disp, struct disp *mftdsp, int porg);
void do_showmode1(int mode);
void do_showpk();
void do_showmd(int mode);
void do_showperscl(int i);
void do_restart_mft(int nm);
void do_box(float xl,float yl,float xh,float yh);
void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
void redisplay_pick(int num_disp, struct disp *dsp, int porg);
void do_start(int *nm, int menuon);
	float pinterp(float x,float x1,float x2);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
int do_matchprep(char *fname);
int do_phvelprep(char *fname);
void doout(int *ret, char *fname);
void dooutphv(int *ret, char *fname);
void doout1(int mode);
void do_reset_mft(char *fname, int *nm);

void do_axes(int i, int j);
void do_time(void);
void do_clearscale(int i);
void mgwrtxt(float x0, float y0, char *str, int cmd, int color);
void cleararea(float xl, float yl, float xh, float yh, int kolor);

#define		RED	2
#define		WHITE	0
#define		BLACK	1


char *wavetype;


int do_page4(char *fname, char *type)
{
	int nm, ret;
	gmesg(fname);
	if(doshade)
		Pick_Color = WHITE;
	else
		Pick_Color = RED;
	do_start(&nm, ON);
	wavetype = type;
	ret = get_action_mft(nm,fname);
	newpen(1);
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	gframe(1);
	return (ret);
}




void do_unzoom(int num_disp, struct disp *dsp, struct pos *pos, int porg)
{
	float xfrac, yfrac;
	/* turn off XOR */
	newpen(3000);
	do_clear(1);
	/* reset the pxl and pyl */
	/*
		pxl[1] = mftpos[1].xl ;
		pxh[1] = mftpos[1].xh ;
		pyl[1] = mftpos[1].yl ;
		pyh[1] = mftpos[1].yh ;
		*/
	apxl[1] = pos[1].axl ;
	apxh[1] = pos[1].axh ;
	apyl[1] = pos[1].ayl ;
	apyh[1] = pos[1].ayh ;
	opxl[1] = pos[1].xl ;
	opxh[1] = pos[1].xh ;
	opyl[1] = pos[1].yl ;
	opyh[1] = pos[1].yh ;
	xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
	yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
	sclx[1] = 1.0/xfrac;
	scly[1] = 1.0/yfrac;
	scl[1] = MIN(sclx[1],scly[1]);
	do_clearscale(1 );
	gread(pos[1].fname,pxl[1],pyl[1],
	    opxl[1],opyl[1],opxh[1],opyh[1],
	    1, sclx[1],scly[1]);
	do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
	/* redisplay picks */
	redisplay_pick(num_disp, dsp, porg );
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
}

void get_mft_output(void )
{
	/* get the output of mft96	*/
	FILE *mft96ctl;
	FILE *mft96dsp;
	float per;
	float fval;
	int i, j;
	char instr[180];


	/* GET PLOT CONTROL INFORMATION */
	if(( mft96ctl = fopen("mft96.ctl","r")) == NULL) return;
	fscanf(mft96ctl,"%d",&num_mft_disp);
	fscanf(mft96ctl,"%s",instr);
	if(strncmp("XAXIS-FREQUENCY",instr,15) == 0)
		XAxisPeriod = NO ;
	else
		XAxisPeriod = YES;

	for (i=0;i<3;i++){
		fscanf(mft96ctl,"%f %f %f %f %e %e %e %e %s %s %s",
		    &mftpos[i].xl, &mftpos[i].yl,
		    &mftpos[i].xh, &mftpos[i].yh,
		    &mftpos[i].axl, &mftpos[i].ayl,
		    &mftpos[i].axh, &mftpos[i].ayh,
		    mftpos[i].xlnlg, mftpos[i].ylnlg,
		    mftpos[i].fname);
	}
	/* get the number of filter periods */
	fscanf(mft96ctl,"%d",&num_period);
	/* get the filter periods */
#ifdef DEBUG
	fprintf(stderr,"NUM_PERIOD %d\n",num_period);
	fprintf(stderr,"NUM_MFT_DISP   %d\n",num_mft_disp  );
#endif
	period = (float *)calloc(num_period,sizeof(float));
	if(period == (float *)NULL)
		fprintf(stderr,"PERIOD ALLOC FAIL\n");
	for(i=0; i < num_period; i++){
		fscanf(mft96ctl,"%f",&fval);
		period[i] = fval;
	}
	fclose(mft96ctl);
	/* GET DISPERSION INFORMATION */
	if(( mft96dsp = fopen("mft96.disp","r")) == NULL) return;
	/* allocate the array to hold the dispersion information */
	mftdsp = (struct disp *)calloc(num_mft_disp+20,sizeof(struct disp));
	if(mftdsp == (struct disp *)NULL)
		fprintf(stderr,"MFTDSP ALLOC FAIL\n");
	for(j=0;j<num_mft_disp;j++){
		fscanf(mft96dsp,"%s %s %s %d %f %f %f %f %f %e %f %f %f %f %d %d %f %f %s %s %s %d %d %d %d",
		    mftdsp[j].type,
		    mftdsp[j].lvry, mftdsp[j].cug,
		    &mftdsp[j].mode,
		             &per, &mftdsp[j].vel, 
		    &mftdsp[j].evel,
		    &mftdsp[j].dist, &mftdsp[j].az, &mftdsp[j].amp, 
		    &mftdsp[j].lat1, &mftdsp[j].lon1, 
		    &mftdsp[j].lat2, &mftdsp[j].lon2,
		    &mftdsp[j].ictl, &mftdsp[j].symb,
		    &mftdsp[j].instper,
		    &mftdsp[j].alpha,
		    mftdsp[j].comment,
                    mftdsp[j].kstnm,
                    mftdsp[j].kcmpnm,
                    &mftdsp[j].nzyear,
                    &mftdsp[j].nzjday,
                    &mftdsp[j].nzhour,
                    &mftdsp[j].nzmin);
	if(!XAxisPeriod)
		mftdsp[j].per = 1.0/per;
	else
		mftdsp[j].per = per;

	}
	fclose(mft96dsp);
}


void redisplay_pick(int num_disp, struct disp *dsp, int porg)
{
	float per, vel;
	float xval, yval;
	int i;
	/* redisplay picks */
	/* first show the velocities */
	/* clear the scale */
	if(porg == MFT) {
		do_axes(1,1);
		do_time();
	} else if(porg == PHV) {
		do_axes(1,2);
	}
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	/* turn on XOR */
	newpen(3001);
	newpen(Pick_Color);
	for(i=0; i< num_disp ; i++){
		per = dsp[i].per;
		vel = dsp[i].vel;
		if(dsp[i].mode >= 0){
			do_scale(per,vel,&xval,&yval,1);
			if(inside(xval,yval,pxl[1],pyl[1],pxh[1],pyh[1])){
				/* put up the velocity  point */
				gsolid(xval,yval,0.02*scl[1],dsp[i].symb);
			}
		}
	}
	/* turn off XOR */
	newpen(3000);

}

void do_start(int *nm, int menuon)
{
	int i;
	float xfrac, yfrac;
	for (i=0;i<3;i++){
		pxl[i] = upxl[i];
		pxh[i] = upxh[i];
		pyl[i] = upyl[i];
		pyh[i] = upyh[i];
	}
	get_mft_output();
	for(i=0;i<3;i++){
		if(strncmp(mftpos[i].xlnlg,"lin",3) == 0)
			xlnlg[i] = LIN;
		else
			xlnlg[i] = LOG;
		if(strncmp(mftpos[i].ylnlg,"lin",3) == 0)
			ylnlg[i] = LIN;
		else
			ylnlg[i] = LOG;
		uxlnlg[i] = xlnlg[i];
		uylnlg[i] = ylnlg[i];
		apxl[i] = mftpos[i].axl;
		apyl[i] = mftpos[i].ayl;
		apxh[i] = mftpos[i].axh;
		apyh[i] = mftpos[i].ayh;
		opxl[i] = mftpos[i].xl ;
		opyl[i] = mftpos[i].yl ;
		opxh[i] = mftpos[i].xh ;
		opyh[i] = mftpos[i].yh ;
		/* use this perception to get scale factor */
		xfrac = (opxh[i] - opxl[i])/(pxh[i] - pxl[i]);
		yfrac = (opyh[i] - opyl[i])/(pyh[i] - pyl[i]);
		/* now read in the figure */
		sclx[i] = 1.0/xfrac;
		scly[i] = 1.0/yfrac;
		scl[i] = MIN(sclx[i],scly[i]);
		usclx[i] = sclx[i];
		uscly[i] = scly[i];
	}
	/* initialize plot scale parameters */
	/*
		printf("%d\n",num_mft_disp);
		for(i=0;i<3;i++)
		printf("%f %f %f %f %e %e %e %e %s %s %s\n",
			mftpos[i].xl, mftpos[i].yl,
			mftpos[i].xh, mftpos[i].yh,
			mftpos[i].axl, mftpos[i].ayl,
			mftpos[i].axh, mftpos[i].ayh,
			mftpos[i].xlnlg, mftpos[i].ylnlg,
			mftpos[i].fname);
	
		for(i=0;i<num_mft_disp;i++)
		printf("%s %s %s %d %f %f %f %f %f %e %f %f %f %f %d %d %f %s %s %s %d %d %d %d %f\n",
			mftdsp[i].type, 
			mftdsp[i].lvry, mftdsp[i].cug,
			mftdsp[i].mode,
			mftdsp[i].per, mftdsp[i].vel, 
			mftdsp[i].evel, 
			mftdsp[i].dist, mftdsp[j].az, mftdsp[i].amp, 
			mftdsp[i].lat1, mftdsp[i].lon1, 
			mftdsp[i].lat2, mftdsp[i].lon2,
			mftdsp[i].ictl, mftdsp[i].symb,
			mftdsp[i].instper,
			mftdsp[j].comment,
			mftdsp[j].kstnm,
			mftdsp[j].kcmpnm,
			mftdsp[j].nzyear,
			mftdsp[j].nzjday,
			mftdsp[j].nzhour,
			mftdsp[j].nzmin,
			mftdsp[j].phase);
	*/

	/* put up the plots */
	for(i=0;i<3;i++){
		gread(mftpos[i].fname,pxl[i],pyl[i],
		    mftpos[i].xl,mftpos[i].yl,
		    mftpos[i].xh,mftpos[i].yh,
		    1, sclx[i],scly[i]);
		do_box( pxl[i], pyl[i], pxh[i], pyh[i]);
	}
	/* now put up the plot scales */
	do_axes(0,0);
	do_axes(1,1);
	do_time();
	if(overlaytomography){
		for ( i = 0 ; i < sizeof(m)/sizeof(struct menu) ; i++ ){
			if(m[i].action == OVERLAYTOMOGPV )
				m[i].visible = 1;
		}
	}
	if(do_phvel){
		for ( i = 0 ; i < sizeof(m)/sizeof(struct menu) ; i++ ){
			if(m[i].action == PHASEVELOCITY )
				m[i].visible = 1;
		}
	}
	if(menuon)show_menu(0.5, 1.1, m, sizeof(m), nm);
	if(do_auto_feed == ON){
		Mode = AutoMode;
		do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"o_auto_feed == ON Mode = %d do_Start tdo_mft4 546\n",Mode);
#endif
	} else {
		Mode = -1;
	}
}

void do_restart_mft(int nm)
{
	/* clear the graphics area and do not redisplay the menu */
	gclip("off", XL, YL, XH, YH);
	/* turn off XOR and erase entire display area */
	newpen(3000);
	newpen(0);
	shader(XL, YL, XH, YH, 0, 0, 0.01, 0.01);
	if(do_auto_feed == ON){
		Mode = AutoMode;
	} else {
		Mode = -1;
	}
	iMode = -1;
	InPick = OFF;
	InLine = OFF;
	do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 569 iMode %d\n",Mode,iMode);
fprintf(stderr,"o_auto_feed == ON Mode = %d do_restart\n",Mode);
#endif
	do_showpk();
	do_showauto();
	mgwrtxt(0.3,7.7,wavetype,1,2);
	newpen(1);
	do_start(&nm, OFF);
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	/* should mean reset all picking */
}


void doout(int *ret, char *fname)
{
	int j, myret;
	int i, ii;
	float xv, yv;
	char c[2];
	int cmd, nmy;
	float per;
	FILE *mft96odsp;
	char ofname[100];
	if(interstationgreen != 1){
		sprintf(ofname,"%s%s.dsp",mftdsp[0].kstnm,mftdsp[0].kcmpnm);
	} else {
		/* strip off everything before the file name */
		sprintf(ofname,"%s.dsp",basename(fname));
	}
	
	fprintf(stderr,"%s\n",ofname);
	/* ask whether to also save the file as STACOMP.dsp */
	newpen(3000);
	mgwrtxt(0.1,0.50,"Save as:",1,1);
	mgwrtxt(1.1,0.50,ofname,1,1);
	show_menu(7.1, 0.70, my, sizeof(my), &nmy);
	/* place current mode at top of page */
	myret = -1 ;
	cmd = -1;
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmy ; i++){
			if(inside(xv,yv,my[i].xl,my[i].yl,my[i].xh,my[i].yh))
			{
				cmd = my[i].action;
				gmesg(my[i].str);
				break;
			}
		}
		if(cmd > 0){
			if(cmd == 1)
				myret = 1;
			else
				myret = -1;
		}
	}
	/* if myret == 1 create the new file else use the mft96.ods */
	if(myret > 0){
		if(( mft96odsp = fopen(ofname     ,"w+")) == NULL) {fprintf(stderr,"Cannot open %s \n",ofname); *ret =  -10 ;return;}
		if(*ret != 2)
			*ret = myret;
	} else {
		if(( mft96odsp = fopen("mft96.ods","w+")) == NULL) {fprintf(stderr,"Cannot open mft96.ods \n"); *ret = -10 ;return;}
	}
	for(j=0;j<num_mft_disp;j++){
		if(mftdsp[j].mode >= 0){
			if(!XAxisPeriod) {
				per = 1.0/mftdsp[j].per;
			} else {
				per = mftdsp[j].per;
			}
		fprintf(mft96odsp,"%s %1s %1s %2d %11.4g %11.5f %11.5f %10.4f %6.1f %11.4e %f %f %f %f %d %d %f %s %s %s %d %d %d %d\n",
		    mftdsp[j].type,
		    mftdsp[j].lvry, mftdsp[j].cug,
		    mftdsp[j].mode,
		              per, mftdsp[j].vel, 
		    mftdsp[j].evel,
		    mftdsp[j].dist, mftdsp[j].az, mftdsp[j].amp, 
		    mftdsp[j].lat1, mftdsp[j].lon1, 
		    mftdsp[j].lat2, mftdsp[j].lon2,
		    mftdsp[j].ictl, mftdsp[j].symb,
		    mftdsp[j].instper,
		    mftdsp[j].comment,
                    mftdsp[j].kstnm,
                    mftdsp[j].kcmpnm,
                    mftdsp[j].nzyear,
                    mftdsp[j].nzjday,
                    mftdsp[j].nzhour,
                    mftdsp[j].nzmin);

		}
	}
	fclose(mft96odsp);
}

int do_matchprep(char *fname)
{
	/* interpolate the instantaneous period and the chosen spectra, according to
		the period array. Also select the mode for output
	*/
	int nm;
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	ModeSelect = -1;
	newpen(1);
	do_axes(0,0);
	do_axes(1,1);
	/* now put up a menu and then an event loop */
/*
	for(i=0;i<2;i++){
		do_box( pxl[i], pyl[i], pxh[i], pyh[i]);
	}
*/
	newpen(1);
	do_showpick1(num_mft_disp, mftdsp, MFT);
	/* now show the desired periods */
	do_showperscl(0);
	do_showperscl(1);
	show_menu(0.5, 1.1, mm, sizeof(mm), &nm);
	return(get_action_pmf(nm,fname));
}

int do_phvelprep(char *fname)
{
	/* interpolate the instantaneous period and the chosen spectra, according to
		the period array. Also select the mode for output
	*/
	int i;
	int nm;
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	ModeSelect = -1;
	newpen(1);
	do_axes(0,1);
	do_axes(1,2);
	/* now put up a menu and then an event loop */
	for(i=0;i<2;i++){
		gread(phvpos[i].fname,pxl[i],pyl[i],
		    phvpos[i].xl,phvpos[i].yl,
		    phvpos[i].xh,phvpos[i].yh,
		    1, sclx[i],scly[i]);
		do_box( pxl[i], pyl[i], pxh[i], pyh[i]);
	}
	if(overlaytomography){
		for ( i = 0 ; i < sizeof(mp)/sizeof(struct menu) ; i++ ){
			if(mp[i].action == OVERLAYTOMOPHV ){
				mp[i].visible = 1;
			}
		}
	}
	newpen(1);
	do_showpick1(num_phv_disp, phvdsp, PHV);
/*
*/
	/* now show the desired periods */
	show_menu(0.5, 1.1, mp, sizeof(mp), &nm);
	return(get_action_phv(nm,fname));
}

static char *pathname;
int get_action_pmf(int nm, char *fname)
{
	int cmd, i, ii;
	int ret;
	float xv, yv;
	float xl, yl, xh, yh;
	char c[2];
	int InDel = OFF;
	int InAdd = OFF;
	int oModeSelect ;
	/* annotate top */
	oModeSelect = 1;
	mgwrtxt(2.0,7.4,"None\0",2, 1);
	Mode = -1;
	do_showmd(ModeSelect);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 745\n",ModeSelect);
fprintf(stderr,"o_auto_feed == ON ModeSelect = %d get_action_pmf\n",ModeSelect);
#endif
	mgwrtxt(0.5,7.7,wavetype,1,2);
	/* prohibit writing at the bottom */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	if(Cursor == 0)
		gcursor("XORArrow");
	else 
		gcursor("Cross");
	for(; ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
		/* loop on commands */
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,mm[i].xl,mm[i].yl,mm[i].xh,mm[i].yh))
			{
				cmd = mm[i].action;
				gmesg(mm[i].str);
				ii = i;
				break;
			}
		}
		if(cmd > 0){
			switch(cmd){
			case ZOOM:	/* Zoom */
				draw_button(0.5+ii*1.125,1.1 ,
				    mm[ii].str,&xl,&yl,&xh,&yh,ON ,mm[ii].lstrmx);
				if(Cursor == 0)
					gcursor("XORArrow");
				else 
					gcursor("Cross");
				do_zoomm();
				if(Cursor == 0)
					gcursor("XORArrow");
				else 
					gcursor("Cross");
				draw_button(0.5+ii*1.125,1.1 ,
				    mm[ii].str,&xl,&yl,&xh,&yh,OFF,mm[ii].lstrmx);
				do_clear(1);
				do_clearscale(1);
				do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
				do_axes(1,1);
				do_showpick1(num_mft_disp, mftdsp, MFT);
				do_showmode1(ModeSelect);
				do_showperscl(1);
				break;
			case UNZOOM:	/* UnZoom */
				InPick = OFF;
				InLine = OFF;
				do_unzoomm();
				do_clear(1);
				do_clearscale(1);
				do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
				do_axes(1,1);
				do_showpick1(num_mft_disp, mftdsp, MFT);
				do_showmode1(ModeSelect);
				do_showperscl(1);
				break;
			case MODE:	/* Mode */
				oModeSelect = ModeSelect;
				ModeSelect = do_mode();
                                AutoMode = ModeSelect;
				do_showmode1(oModeSelect);
				do_showmode1(ModeSelect);
				do_showmd(ModeSelect);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 811\n",ModeSelect);
fprintf(stderr,"o_auto_feed == ON ModeSelect = %d get_action_pmf MODE\n",ModeSelect);
#endif
				break;
			case INTERP:
				printf("Interpolate Periods within current box\n");
				break;
			case ADD:
				if(InAdd == ON) {
					InAdd = OFF ;
					InDel = OFF;
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,OFF,mm[ii].lstrmx);
				} else {
					InAdd = ON;
					InDel = OFF;
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,ON ,mm[ii].lstrmx);
				}
				printf("Add a fictitious point for Match \n");
				break;
			case DELETE:
				if(InDel == ON) {
					InDel = OFF ;
					InAdd = OFF ;
					InPick = OFF;
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,OFF,mm[ii].lstrmx);
				} else {
					draw_button(0.5+ii*1.125,1.1 ,
				    		mm[ii].str,&xl,&yl,&xh,&yh,ON ,mm[ii].lstrmx);
					InDel = ON;
					InAdd = OFF ;
					InPick = ON;
				}
				printf("Delete a real/fictitious dispersion for Match \n");
				break;
			case RETURN:	/* Return */
				gcursor("XORArrow");
				gframe(1);
				return(0);
			case MATCH:	/* Match */
				while(ModeSelect < 0 ){
					gmesg("No mode chosen");
					oModeSelect = ModeSelect;
					ModeSelect = do_mode();
					do_showmode1(oModeSelect);
					do_showmode1(ModeSelect);
					do_showmd(ModeSelect);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 859\n",ModeSelect);
fprintf(stderr,"o_auto_feed == ON ModeSelect = %d get_action_pmf MATCH\n",ModeSelect);
#endif
				}
				doout1(ModeSelect);

#ifdef MSDOS
sprintf(ostr,"sacmat96.exe -F %s -D disp.d -AUTO", fname);
#else
        pathname = (char *)my_pathfind(getenv("PATH"), "sacmat96", "rx");
#endif

fprintf(stderr,"ModeSelect: %d pathname: %s\n",ModeSelect, pathname);
sprintf(ostr,"%s -F %s -D disp.d -AUTO",pathname, fname);

fprintf(stderr,"ostr: %s\n",ostr);
					if(strcmp(pathname, " ") != 0){
						printf("run sacmat96 Mode: %d\n",ModeSelect);
							ret = system(ostr);
							if(ret >= 0 ){
								/* success */
								return(2);
							}
					}
				break;
			case EXIT:
				printf("Exit\n");
				gcursor("XORArrow");
				gframe(1);
				return(-1);
			default:
				break;

			}
		} else {
			if(InPick == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_pick(xv, yv, 0, num_mft_disp, mftdsp, MFT);
				}
			}
		}
	}
} 

void do_zoomm(void)
{
	/* only used by PMF - Match */
	char c[2];
	float xv1, xv2, yv1, yv2;
	float xfrac, yfrac;
	float x1,x2,y1,y2;
	float tx1,tx2,ty1,ty2;
	float per1, vel1;
	float per2, vel2;
	float px1,px2,py1,py2;
	/* turn off XOR */
	newpen(3000);
	curaxy(&xv1, &yv1, c);
	if( !inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]))return;
	gcursor("Box");
	curaxy(&xv2, &yv2, c);
	if(Cursor == 0)
		gcursor("XORArrow");
	else 
		gcursor("Cross");
	if(xv1==xv2 && yv1 == yv2) return;  /* need two points for a line */
	/* 
	We work with three coordinate systems, each consisting of two
	opposite corners of a view rectangle

	(pxl,pyl) -> (pxh,pyl)   is display screen which does not change
	(opxl,opyl) -> (opxh,opyh) is the mapped portion of the original plot
	(apxl,apyl) -> (apxh,apyl) are user coordinates of the original space

	To do the mapping, use a simple interpolation

	V(p) = (1-p)V1 + pV2   where  0 <= p <= 1
	*/
	if( inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1])
	    && inside(xv2,yv2, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]) ) {
		/* display all current picks to turn them off */
		/* since the point may extend across clip region ? */
/*
		redisplay_pick();
*/
		do_clear(1);
		/* these are screen coordinates */
		x1 = MIN(xv1,xv2) ;
		x2 = MAX(xv1,xv2) ;
		y1 = MIN(yv1,yv2) ;
		y2 = MAX(yv1,yv2) ;
		px1 = pinterp(x1,pxl[1],pxh[1]);
		px2 = pinterp(x2,pxl[1],pxh[1]);
		py1 = pinterp(y1,pyl[1],pyh[1]);
		py2 = pinterp(y2,pyl[1],pyh[1]);
		/* now estimate the coordinates in the
					original plot space */
		tx1 = (1.0 - px1) * opxl[1] + px1 * opxh[1];
		tx2 = (1.0 - px2) * opxl[1] + px2 * opxh[1];
		ty1 = (1.0 - py1) * opyl[1] + py1 * opyh[1];
		ty2 = (1.0 - py2) * opyl[1] + py2 * opyh[1];
		/* update our perception of the view of original space */
		opxl[1] = tx1;
		opxh[1] = tx2;
		opyl[1] = ty1;
		opyh[1] = ty2;
		/* use this perception to get scale factor */
		xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
		yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
		sclx[1] = 1.0/xfrac;
		scly[1] = 1.0/yfrac;
		scl[1] = MIN(sclx[1],scly[1]);
		/* reset the plot values at the new corners */
		do_invscale(&per1,&vel1,x1,y1,1);
		do_invscale(&per2,&vel2,x2,y2,1);
		apxl[1] = per1;
		apyl[1] = vel1;
		apxh[1] = per2;
		apyh[1] = vel2;
	}
}

void do_unzoomm(void)
{
	/* only used by PMF - Match */
	float xfrac, yfrac;
	/* turn off XOR */
	newpen(3000);
	do_clear(1);
	/* reset the pxl and pyl */
	/*
		pxl[1] = mftpos[1].xl ;
		pxh[1] = mftpos[1].xh ;
		pyl[1] = mftpos[1].yl ;
		pyh[1] = mftpos[1].yh ;
		*/
	apxl[1] = mftpos[1].axl ;
	apxh[1] = mftpos[1].axh ;
	apyl[1] = mftpos[1].ayl ;
	apyh[1] = mftpos[1].ayh ;
	opxl[1] = mftpos[1].xl ;
	opxh[1] = mftpos[1].xh ;
	opyl[1] = mftpos[1].yl ;
	opyh[1] = mftpos[1].yh ;
	xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
	yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
	sclx[1] = 1.0/xfrac;
	scly[1] = 1.0/yfrac;
	scl[1] = MIN(sclx[1],scly[1]);
	/* redisplay picks */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
}

/* display a particular mode - this will permit XOR of
   switching the entire mode on/off and also the interactive picking
*/
void do_showmode1(int mode){
	int j;
	float xval, yval, vel, amp, per;
	if(mode < 0)return;
	for(j=0;j<num_mft_disp;j++){
		if(mftdsp[j].mode == mode  ){
			if(!XAxisPeriod)
				per = 1.0/mftdsp[j].instper;
			else
				per = mftdsp[j].instper;

			vel = mftdsp[j].vel;
			amp = mftdsp[j].amp;
			/* display velocity value */
			/* turn on XOR */
			newpen(3001);
			newpen(1);
			do_scale(per,vel,&xval,&yval,1);
			if(inside(xval,yval,pxl[1],pyl[1],pxh[1],pyh[1])){
				gsolid(xval,yval,0.035*scl[1],mftdsp[j].mode);
			}
			/* display spectral amplitude value */
			do_scale(per,amp,&xval,&yval,0);
			if(inside(xval,yval,pxl[0],pyl[0],pxh[0],pyh[0])){
				gsolid(xval,yval,0.035*scl[0],mftdsp[j].mode);
			}
			/* turn off XOR */
			newpen(3000);
		}
	}
}


void doout1(int mode)
{
	int j;
	FILE *mft96mat;
	FILE *mft96matp;
	if(( mft96mat = fopen("disp.d","w")) == NULL) return;
	if(( mft96matp = fopen("disp.dp","w")) == NULL) return;
	for(j=0;j<num_mft_disp;j++){
		if(mftdsp[j].mode == mode){
		fprintf(mft96mat,"SURF96 %1s %1s T %2d %11.4g %11.4g %11.4g \n",
			mftdsp[j].lvry, 
			mftdsp[j].cug,
			mftdsp[j].mode,
			mftdsp[j].per, 
			mftdsp[j].vel,
			mftdsp[j].evel); 
		fprintf(mft96matp,"SURF96 %1s %1s T %2d %11.4g %11.4g %11.4g \n",
			mftdsp[j].lvry, 
			mftdsp[j].cug,
			mftdsp[j].mode,
			mftdsp[j].instper, 
			mftdsp[j].vel,
			mftdsp[j].evel); 

		}
	}
	fclose(mft96mat);
	fclose(mft96matp);
}
int get_action_mft(int nm ,char *fname)
{
	int cmd, i, ii;
	int ret;
	float xv, yv;
	float xv1, yv1;
	float xl, yl, xh, yh;
	char c[2];
	/* annotate top */
/*
	mgwrtxt(2.0,7.4,"None\0",2, 1);
*/
	do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 1091\n",Mode);
fprintf(stderr,"o_auto_feed == %d Mode = %d get_Action_mft\n",do_auto_feed,Mode);
#endif
	do_showpk();
	do_showauto();
	mgwrtxt(0.5,7.7,wavetype,1,2);
	
	/* prohibit writing at the bottom */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	if(Cursor == 0)
		gcursor("XORArrow");
	else 
		gcursor("Cross");
	for(; ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
		/* loop on commands */
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,m[i].xl,m[i].yl,m[i].xh,m[i].yh))
			{
				cmd = m[i].action;
				gmesg(m[i].str);
				ii = i;
				break;
			}
		}
		if(cmd > 0){
			switch(cmd){
			case ZOOM:	/* Zoom */
				draw_button(0.5+ii*1.125,1.1 ,
				    m[ii].str,&xl,&yl,&xh,&yh,ON ,m[ii].lstrmx);
				if(Cursor == 0)
					gcursor("XORArrow");
				else 	
					gcursor("Cross");
				do_zoom(num_mft_disp, mftdsp, mftpos, MFT);
				if(Cursor == 0)	
					gcursor("XORArrow");
				else 	
					gcursor("Cross");
				draw_button(0.5+ii*1.125,1.1 ,
				    m[ii].str,&xl,&yl,&xh,&yh,OFF,m[ii].lstrmx);
				break;
			case UNZOOM:	/* UnZoom */
				InPick = OFF;
				InLine = OFF;
				do_showpk();
				do_unzoom(num_mft_disp, mftdsp,mftpos, MFT);
				break;
			case AUTO:	/* Automatic Pick from nearest line */
				if(Mode <0){
					gmesg("Choose a Mode First");
					Mode = do_mode();
					AutoMode = Mode;
					do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 1146\n",Mode);
fprintf(stderr,"o_auto_feed == %d Mode = %d get_Action_mft AUTO\n",do_auto_feed,Mode);
#endif
				}
				InPick = OFF;
				InLine = ON;
				do_showpk();
				curaxy(&xv1, &yv1, c);
				do_line(xv1, yv1, num_mft_disp, mftdsp, MFT);
				break;
			case MODE:	/* Mode */
				Mode = do_mode();
				AutoMode = Mode;
				do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 1158\n",Mode);
fprintf(stderr,"o_auto_feed == %d Mode = %d get_Action_mft MODE \n",do_auto_feed,Mode);
#endif
				break;
			case PICK:	/* Pick */
				if(Mode <0){
					gmesg("Choose a Mode First");
					Mode = do_mode();
					do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 1166\n",Mode);
fprintf(stderr,"o_auto_feed == %d Mode = %d get_Action_mft PICK \n",do_auto_feed,Mode);
#endif
					AutoMode = Mode;
				}
				if(InPick == ON) {
					InPick = OFF ;
				} else {
					InPick = ON;
				}
				InLine = OFF;
				do_showpk();
				break;
			case RESTART:	/* Restart */
				do_restart_mft(nm);
				InPick = OFF;
				InLine = OFF;
				do_showpk();
				break;
			case AUTOFEED:
					/* toggle the auto feed */
				if(do_auto_feed == ON)
					do_auto_feed = OFF;
				else
					do_auto_feed = ON;
				do_showauto() ;
				break;
			case OVERLAYTOMOGPV:
				if(overlaytomography){
				/* turn off the overlay button for this page 
 					the MFT96.PLT has been changed */
					for ( i = 0 ; i < sizeof(m)/sizeof(struct menu) ; i++ ){
						if(m[i].action == OVERLAYTOMOGPV ){
							m[i].visible = 0;
							/* clear the Tomo button */
							gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
							cleararea(m[i].xl,m[i].yl,m[i].xh,m[i].yh,0);
							gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
						}
					}
        				pathname = (char *)my_pathfind(getenv("PATH"), "MFTDOOVERLAY", "rx");
					if(strcmp(pathname, " ") != 0){
						sprintf(ostr,"%s \n" ,pathname); 
						/* put up message */

						ret = system(ostr);
						InPick = OFF;
						InLine = OFF;
						do_showpk();
						do_unzoom(num_mft_disp, mftdsp, mftpos, MFT);
					}
				}
				break;
			case MATCH:	/* Match */
				gframe(1);
				gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
				/* turn off XOR */
				newpen(3000);
				/* reset the zoom parameters */
				do_resetzoom();
				ret = do_matchprep(fname);
				gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
				if(ret > 0){
					goto end;
				} else {
                                        do_reset_mft(fname,&nm);
				}
				break;
			case PHASEVELOCITY:	/* inter-station phase velocity */
				ret = PHVRESTART ;
				while(ret == PHVRESTART){
					gframe(1);
					gcursor("XORArrow");
					/* turn off XOR */
					newpen(3000);
					gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
					/* reset the zoom parameters */
					do_resetzoom();
					get_phv_output();
					ret = do_phvelprep(fname);
				}
				gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
                                do_reset_mft(fname,&nm);
				break;
			case EXIT:
				gcursor("XORArrow");
				ret = 1;
				goto end;
			default:
				break;

			}
		} else {
			if(InPick == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_pick(xv, yv, 0, num_mft_disp, mftdsp, MFT);
				}
			}
			if(InLine == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_line(xv, yv, num_mft_disp, mftdsp, MFT);
				}
			}
		}
	}
end:
	/* turn off XOR */
	newpen(3000);
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	gmesg("Save U/Amp File?");
	/* interactively ask if the selected dispersion should be saved */
	doout(&ret,fname);
	/* free period array */
	free(period);
	free(mftdsp);
	return(ret);
} 

void do_zoom(int num_disp, struct disp *dsp,struct pos *pos, int porg)
{
	char c[2];
	float xv1, xv2, yv1, yv2;
	float xfrac, yfrac;
	float x1,x2,y1,y2;
	float tx1,tx2,ty1,ty2;
	float per1, vel1;
	float per2, vel2;
	float px1,px2,py1,py2;
	/* turn off XOR */
	newpen(3000);
	curaxy(&xv1, &yv1, c);
	if( !inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]))return;
	gcursor("Box");
	curaxy(&xv2, &yv2, c);
	if(Cursor == 0)
		gcursor("XORArrow");
	else 
		gcursor("Cross");
	if(xv1==xv2 && yv1 == yv2) return;  /* need two points for a line */
	/* 
	We work with three coordinate systems, each consisting of two
	opposite corners of a view rectangle

	(pxl,pyl) -> (pxh,pyl)   is display screen which does not change
	(opxl,opyl) -> (opxh,opyh) is the mapped portion of the original plot
	(apxl,apyl) -> (apxh,apyl) are user coordinates of the original space

	To do the mapping, use a simple interpolation

	V(p) = (1-p)V1 + pV2   where  0 <= p <= 1
	*/
	if( inside(xv1,yv1, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1])
	    && inside(xv2,yv2, pxl[1]-0.1, pyl[1],
	    pxh[1]+0.1, pyh[1]) ) {
		/* display all current picks to turn them off */
		/* since the point may extend across clip region ? */
		redisplay_pick(num_disp, dsp, porg);
		do_clear(1);
		/* these are screen coordinates */
		x1 = MIN(xv1,xv2) ;
		x2 = MAX(xv1,xv2) ;
		y1 = MIN(yv1,yv2) ;
		y2 = MAX(yv1,yv2) ;
		px1 = pinterp(x1,pxl[1],pxh[1]);
		px2 = pinterp(x2,pxl[1],pxh[1]);
		py1 = pinterp(y1,pyl[1],pyh[1]);
		py2 = pinterp(y2,pyl[1],pyh[1]);
		/* now estimate the coordinates in the
					original plot space */
		tx1 = (1.0 - px1) * opxl[1] + px1 * opxh[1];
		tx2 = (1.0 - px2) * opxl[1] + px2 * opxh[1];
		ty1 = (1.0 - py1) * opyl[1] + py1 * opyh[1];
		ty2 = (1.0 - py2) * opyl[1] + py2 * opyh[1];
		/* update our perception of the view of original space */
		opxl[1] = tx1;
		opxh[1] = tx2;
		opyl[1] = ty1;
		opyh[1] = ty2;
		/* use this perception to get scale factor */
		xfrac = (opxh[1] - opxl[1])/(pxh[1] - pxl[1]);
		yfrac = (opyh[1] - opyl[1])/(pyh[1] - pyl[1]);
		/* now read in the figure */
		sclx[1] = 1.0/xfrac;
		scly[1] = 1.0/yfrac;
		scl[1] = MIN(sclx[1],scly[1]);
		do_clearscale(1 );
		do_clearscale(2 );
		gread(pos[1].fname,pxl[1],pyl[1],
		    opxl[1],opyl[1],opxh[1],opyh[1],
		    1, sclx[1],scly[1]);
		do_box( pxl[1], pyl[1], pxh[1], pyh[1]);
		/* reset the plot values at the new corners */
		do_invscale(&per1,&vel1,x1,y1,1);
		do_invscale(&per2,&vel2,x2,y2,1);
		apxl[1] = per1;
		apyl[1] = vel1;
		apxh[1] = per2;
		apyh[1] = vel2;
		redisplay_pick(num_disp, dsp, porg);
	}
}

void do_reset_mft(char *fname, int *pnm)
{
        int i,j;
	float per, vel, amp, xval, yval;
	for (i=0;i<3;i++){
		pxl[i] = upxl[i];
		pxh[i] = upxh[i];
		pyl[i] = upyl[i];
		pyh[i] = upyh[i];
		opxl[i] = mftpos[i].xl ;
		opyl[i] = mftpos[i].yl ;
		opxh[i] = mftpos[i].xh ;
		opyh[i] = mftpos[i].yh ;
		apxl[i] = mftpos[i].axl;
		apyl[i] = mftpos[i].ayl;
		apxh[i] = mftpos[i].axh;
		apyh[i] = mftpos[i].ayh;
		xlnlg[i] = uxlnlg[i];
		ylnlg[i] = uylnlg[i];
		sclx[i] = usclx[i] ;
		scly[i] = uscly[i] ;
	}
	gmesg(fname);
	show_menu(0.5, 1.1, m, sizeof(m), pnm);
	InPick = OFF;
	InLine = ON;
	if(do_auto_feed == ON){
		Mode = AutoMode;
	} else {
		if(iMode >= 0)
			Mode = iMode ;
		else
			Mode = -1;
	}
	/* redisplay the plots */
	do_resetzoom();
	for(j=0;j<3;j++){
		gread(mftpos[j].fname,pxl[j],pyl[j],
			opxl[j],opyl[j],opxh[j],opyh[j],
			1, sclx[j],scly[j]);
	}
	do_axes(0,0);
	do_axes(1,1);
	do_time();
	/* redisplay the current picked 
	values of both 
	amplitude and velocity in  XOR */
	do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft4 1414 iMode %d \n",Mode,iMode);
fprintf(stderr,"o_auto_feed = %d Mode = %d doresetmft\n",do_auto_feed,Mode);
#endif
	do_showpk();
	mgwrtxt(0.5,7.7,wavetype,1,2);
	/* turn on XOR */
	newpen(3001);
	for(j=0;j<num_mft_disp;j++){
		if(mftdsp[j].mode >= 0){
			per = mftdsp[j].per;
			vel = mftdsp[j].vel;
			amp = mftdsp[j].amp;
			/* display velocity value */
			do_scale(per,vel,&xval,&yval,1);
			newpen(Pick_Color);
			gsolid(xval,yval,0.02*scl[1],mftdsp[j].symb);
			/* display spectral amplitude value */
			do_scale(per,amp,&xval,&yval,0);
			newpen(1);
			gsolid(xval,yval,0.02*scl[0],mftdsp[j].symb);
		}
	}
	/* turn off XOR */
	newpen(3000);
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
}

void get_phv_output(void )
{
	/* get the output of mft96	*/
	FILE *phv96ctl;
	FILE *phv96dsp;
	float per;
	float fval;
	int i, j;
	char instr[180];
	float xfrac, yfrac;


	/* define placement of graphics */
	/* GET PLOT CONTROL INFORMATION */
	if(( phv96ctl = fopen("phv96.ctl","r")) == NULL) return;
	fscanf(phv96ctl,"%d",&num_phv_disp);
	fscanf(phv96ctl,"%s",instr);
	if(strncmp("XAXIS-FREQUENCY",instr,15) == 0)
		XAxisPeriod = NO ;
	else
		XAxisPeriod = YES;

	/* there is only one line of entry */
	for (i=0;i<2;i++){
		fscanf(phv96ctl,"%f %f %f %f %e %e %e %e %s %s %s",
		    &phvpos[i].xl, &phvpos[i].yl,
		    &phvpos[i].xh, &phvpos[i].yh,
		    &phvpos[i].axl, &phvpos[i].ayl,
		    &phvpos[i].axh, &phvpos[i].ayh,
		    phvpos[i].xlnlg, phvpos[i].ylnlg,
		    phvpos[i].fname);
		if(strncmp(phvpos[i].xlnlg,"lin",3) == 0)
			xlnlg[i] = LIN;
		else
			xlnlg[i] = LOG;
		if(strncmp(phvpos[i].ylnlg,"lin",3) == 0)
			ylnlg[i] = LIN;
		else
			ylnlg[i] = LOG;
	  	pxlnlg[i] = xlnlg[i] ;
		pylnlg[i] = ylnlg[i] ;
		pxl[i] = ppxl[i];
		pxh[i] = ppxh[i];
		pyl[i] = ppyl[i];
		pyh[i] = ppyh[i];
		apxl[i] = phvpos[i].axl;
		apyl[i] = phvpos[i].ayl;
		apxh[i] = phvpos[i].axh;
		apyh[i] = phvpos[i].ayh;
		opxl[i] = phvpos[i].xl ;
		opxh[i] = phvpos[i].xh ;
		opyl[i] = phvpos[i].yl ;
		opyh[i] = phvpos[i].yh ;
		xfrac = (opxh[i] - opxl[i])/(pxh[i] - pxl[i]);
		yfrac = (opyh[i] - opyl[i])/(pyh[i] - pyl[i]);
		sclx[i] = 1.0/xfrac;
		scly[i] = 1.0/yfrac;
	}
	/* get the number of filter periods */
	fscanf(phv96ctl,"%d",&num_period);
	/* get the filter periods */
/*
	fprintf(stderr,"NUM_PERIOD %d\n",num_period);
	fprintf(stderr,"NUM_PHV_DISP   %d\n",num_phv_disp  );
	period = (float *)calloc(num_period,sizeof(float));
NEVER ALLOCATE PERIOD HERE SINCE THIS WAS ALREADY ALLOCATED BY
READING mft96.ctl - sacmft96 has same periods in each
	if(period == (float *)NULL)
		fprintf(stderr,"PERIOD ALLOC FAIL\n");
*/
	for(i=0; i < num_period; i++){
		fscanf(phv96ctl,"%f",&fval);
		period[i] = fval;
	}
	fclose(phv96ctl);
	/* GET DISPERSION INFORMATION */
	if(( phv96dsp = fopen("phv96.disp","r")) == NULL) return;
	/* allocate the array to hold the dispersion information */
	phvdsp = (struct disp *)calloc(num_phv_disp+20,sizeof(struct disp));
	if(phvdsp == (struct disp *)NULL)
		fprintf(stderr,"PHVDSP ALLOC FAIL\n");
	for(j=0;j<num_phv_disp;j++){
		fscanf(phv96dsp,"%s %s %s %d %f %f %f %f %f %e %f %f %f %f %d %d %f %f %s %s %s %d %d %d %d %f %f %d",
		    phvdsp[j].type,
		    phvdsp[j].lvry, phvdsp[j].cug,
		    &phvdsp[j].mode,
		             &per, &phvdsp[j].vel, 
		    &phvdsp[j].evel,
		    &phvdsp[j].dist, &phvdsp[j].az, &phvdsp[j].amp, 
		    &phvdsp[j].lat1, &phvdsp[j].lon1, 
		    &phvdsp[j].lat2, &phvdsp[j].lon2,
		    &phvdsp[j].ictl, &phvdsp[j].symb,
		    &phvdsp[j].instper,
		    &phvdsp[j].alpha,
		    phvdsp[j].comment,
                    phvdsp[j].kstnm,
                    phvdsp[j].kcmpnm,
                    &phvdsp[j].nzyear,
                    &phvdsp[j].nzjday,
                    &phvdsp[j].nzhour,
                    &phvdsp[j].nzmin,
		    &phvdsp[j].phase,
		    &phvdsp[j].uvel,
		    &phvdsp[j].nn);
	if(!XAxisPeriod)
		phvdsp[j].per = 1.0/per;
	else
		phvdsp[j].per = per;

	}
	fclose(phv96dsp);
}

int get_action_phv(int nm, char *fname)
{
	int cmd, i, ii;
	int ret;
	float xv, yv;
	float xv1, yv1;
	float xl, yl, xh, yh;
	char c[2];
	/* annotate top */
/*
	mgwrtxt(2.0,7.4,"None\0",2, 1);
*/
	do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft 1565\n",Mode);
fprintf(stderr,"o_auto_feed == ON Mode = %d\n",Mode);
#endif
	do_showpk();
	do_showauto();
	mgwrtxt(0.5,7.7,wavetype,1,2);
	
	/* prohibit writing at the bottom */
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	if(Cursor == 0)
		gcursor("XORArrow");
	else 
		gcursor("Cross");
	for(; ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
		/* loop on commands */
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,mp[i].xl,mp[i].yl,mp[i].xh,mp[i].yh))
			{
				cmd = mp[i].action;
				gmesg(mp[i].str);
				ii = i;
				break;
			}
		}
		if(cmd > 0){
			switch(cmd){
			case ZOOM:	/* Zoom */
				draw_button(0.5+ii*1.125,1.1 ,
				    mp[ii].str,&xl,&yl,&xh,&yh,ON ,mp[ii].lstrmx);
				if(Cursor == 0)
					gcursor("XORArrow");
				else 
					gcursor("Cross");
				do_zoom(num_phv_disp, phvdsp, phvpos, PHV);
				do_showpick1(num_phv_disp, phvdsp, PHV);
				if(Cursor == 0)
					gcursor("XORArrow");
				else 
					gcursor("Cross");
				draw_button(0.5+ii*1.125,1.1 ,
				    mp[ii].str,&xl,&yl,&xh,&yh,OFF,mp[ii].lstrmx);
				break;
			case UNZOOM:	/* UnZoom */
				InPick = OFF;
				InLine = OFF;
				do_showpk();
				do_showpick1(num_phv_disp, phvdsp, PHV);
				do_unzoom(num_phv_disp, phvdsp,phvpos, PHV);
				break;
			case AUTO:	/* Automatic Pick from nearest line */
				if(Mode <0){
					gmesg("Choose a Mode First");
					Mode = do_mode();
					AutoMode = Mode;
					do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft 1623\n",Mode);
fprintf(stderr,"o_auto_feed == ON Mode = %d getACtionphv AUTO\n",Mode);
#endif
				}
				InPick = OFF;
				InLine = ON;
				do_showpk();
				curaxy(&xv1, &yv1, c);
				do_line(xv1, yv1, num_phv_disp, phvdsp, PHV);
				break;
			case MODE:	/* Mode */
				Mode = do_mode();
				AutoMode = Mode;
				do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft 1636\n",Mode);
fprintf(stderr,"o_auto_feed == ON Mode = %d getActionphv MODE\n",Mode);
#endif
				break;
			case PICK:	/* Pick */
				if(Mode <0){
					gmesg("Choose a Mode First");
					Mode = do_mode();
					AutoMode = Mode;
					do_showmd(Mode);
#ifdef DEBUG
fprintf(stderr,"Mode = %d  tdo_mft 1645\n",Mode);
fprintf(stderr,"o_auto_feed == ON Mode = %d getActionphv PICK\n",Mode);
#endif
				}
				if(InPick == ON) {
					InPick = OFF ;
				} else {
					InPick = ON;
				}
				InLine = OFF;
				do_showpk();
				break;
			case RESTART:	/* Restart */
				return(PHVRESTART);
				break;
			case OVERLAYTOMOPHV:
				if(overlaytomography){
				/* turn off the overlay button for this page 
 					the PHV96.PLT has been changed */
					for ( i = 0 ; i < sizeof(mp)/sizeof(struct menu) ; i++ ){
						if(mp[i].action == OVERLAYTOMOPHV ){
							/* clear the Tomo button */
							gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
							cleararea(mp[i].xl,mp[i].yl,mp[i].xh,mp[i].yh,0);
							gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
						}
					}
        				pathname = (char *)my_pathfind(getenv("PATH"), "PHVDOOVERLAY", "rx");
					if(strcmp(pathname, " ") != 0){
						sprintf(ostr,"%s \n" ,pathname); 
						/* put up message */

						ret = system(ostr);
						InPick = OFF;
						InLine = OFF;
						do_showpk();
						do_unzoom(num_phv_disp, phvdsp, phvpos, PHV);
					}
				}
				break;
			case RETURN:	/* Return */
				gcursor("XORArrow");
				gframe(1);
				return(0);
			case EXIT:
				gcursor("XORArrow");
				ret = 1;
				goto end;
			default:
				break;

			}
		} else {
			if(InPick == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_pick(xv, yv, 0, num_phv_disp, phvdsp, PHV);
				}
			}
			if(InLine == ON){
				if(inside(xv,yv,pxl[1]-0.1,pyl[1],pxh[1]+0.1,pyh[1]))
				{
					do_line(xv, yv, num_phv_disp, phvdsp, PHV);
				}
			}
		}
	}
end:
	/* turn off XOR */
	newpen(3000);
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	gmesg("Save C File?");
	/* interactively ask if the selected dispersion should be saved */
	dooutphv(&ret,fname);
	/* clean the screen  */
	gframe(1);
	/* free phvdsp array */
	free(phvdsp);
	return(ret);
} 

void dooutphv(int *ret, char *fname)
{
	int j, myret;
	int i, ii;
	float xv, yv;
	char c[2];
	int cmd, nmy;
	float per;
	FILE *phv96odsp;
	char ofname[100];
	if(interstationgreen != 1){
		sprintf(ofname,"%s%s.phv",phvdsp[0].kstnm,phvdsp[0].kcmpnm);
	} else {
		/* strip off everything before the file name */
		sprintf(ofname,"%s.phv",basename(fname));
	}
	
	fprintf(stderr,"%s\n",ofname);
	/* ask whether to also save the file as STACOMP.dsp */
	newpen(3000);


	mgwrtxt(0.1,0.50,"Save as:",1,1);
	mgwrtxt(1.1,0.50,ofname,1,1);
	show_menu(7.1, 0.70, my, sizeof(my), &nmy);
	/* place current mode at top of page */
	myret = -1 ;
	cmd = -1;
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmy ; i++){
			if(inside(xv,yv,my[i].xl,my[i].yl,my[i].xh,my[i].yh))
			{
				cmd = my[i].action;
				gmesg(my[i].str);
				break;
			}
		}
		if(cmd > 0){
			if(cmd == 1)
				myret = 1;
			else
				myret = -1;
		}
	}
	/* if myret == 1 create the new file else use the phv96.ods */
	if(myret > 0){
		if(( phv96odsp = fopen(ofname     ,"w+")) == NULL) {fprintf(stderr,"Cannot open %s \n",ofname); *ret =  -10 ;return;}
		if(*ret != 2)
			*ret = myret;
	} else {
		if(( phv96odsp = fopen("phv96.ods","w+")) == NULL) {fprintf(stderr,"Cannot open phv96.ods \n"); *ret = -10 ;return;}
	}
	for(j=0;j<num_phv_disp;j++){
		if(phvdsp[j].mode >= 0){
			if(!XAxisPeriod) {
				per = 1.0/phvdsp[j].per;
			} else {
				per = phvdsp[j].per;
			}
		fprintf(phv96odsp,"%s %1s %1s %2d %11.4g %11.5f %11.5f %10.4f %6.1f %11.4e %f %f %f %f %d %d %f %s %s %s %d %d %d %d %f %f %d\n",
		    phvdsp[j].type,
		    phvdsp[j].lvry, phvdsp[j].cug,
		    phvdsp[j].mode,
		              per, phvdsp[j].vel, 
		    phvdsp[j].evel,
		    phvdsp[j].dist, phvdsp[j].az, phvdsp[j].amp, 
		    phvdsp[j].lat1, phvdsp[j].lon1, 
		    phvdsp[j].lat2, phvdsp[j].lon2,
		    phvdsp[j].ictl, phvdsp[j].symb,
		    phvdsp[j].instper,
		    phvdsp[j].comment,
                    phvdsp[j].kstnm,
                    phvdsp[j].kcmpnm,
                    phvdsp[j].nzyear,
                    phvdsp[j].nzjday,
                    phvdsp[j].nzhour,
                    phvdsp[j].nzmin,
                    phvdsp[j].phase,
                    phvdsp[j].uvel,
                    phvdsp[j].nn);

		}
	}
	fclose(phv96odsp);
}
