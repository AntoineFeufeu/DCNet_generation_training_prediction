#include <stdio.h>
#include <math.h>
#include "do_mft.h"
#include "nmenu.h"
#include "calplot.h"
#include "grphsubc.h"

static struct menu md[] = {
		{  -1.0, -1.0, -1.0, -1.0, "Fund\0" , 0, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "1 st\0" , 1, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "2 nd\0" , 2, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "3 rd\0" , 3, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "4 th\0" , 4, -1, 1, 1},
		{  -1.0, -1.0, -1.0, -1.0, "5 th\0" , 5, -1, 1, 1},
};

float pinterp(float x,float x1,float x2);
void do_box(float xl,float yl,float xh,float yh);
void cleararea(float xl, float yl, float xh, float yh, int kolor);
void do_showperscl(int i);
void do_showauto();
void do_showpk();
void do_showmd(int mode);
void do_time();
void do_clearscale(int i );
void do_axes(int i, int j );
int do_mode(void);
void do_invscale(float *per,float *vel,float xval,float yval,int i);
void do_scale(float per,float vel,float *xval,float *yval,int i);
void do_clear(int i );
void do_pick( float xv1, float yv1, int lev, int num_disp, struct disp *dsp, int porg);
void do_line(float xv1, float yv1, int num_disp, struct disp *dsp, int porg);
void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv,
        float xlb, float ylb, float xhb, float yhb);
void mgwrtxt(float x0, float y0, char *str, int cmd, int color);
void gsolid(float x, float y, float h, int symb);
void do_showpick1(int num_disp, struct disp *dsp, int porg);
void do_resetzoom(void);

/* external definitions that are required */
extern float   pxl[3], pyl[3], pxh[3], pyh[3]; 
extern float scl[3], sclx[3], scly[3];
extern float   apxl[3], apyl[3], apxh[3], apyh[3];
extern float opxl[3], opyl[3], opxh[3], opyh[3];
extern int xlnlg[3], ylnlg[3];
extern struct disp *mftdsp;
extern struct disp *phvdsp;

extern struct pos mftpos[3];
extern struct pos phvpos[3];
extern float *period;
extern int num_period;

extern int XAxisPeriod ;
extern int Mode, iMode, ModeSelect, Units, Pick_Color, InLine, InPick;
extern int Cursor;
extern int Units;

extern int do_auto_feed ;
extern int do_auto_curfile;
extern int do_auto_nfiles;
extern int do_auto_curpage;
extern int do_auto_npages;
extern int do_auto_i_next ;
extern int do_auto_i_max ;
extern int do_auto_ihdr11 ;
extern int AutoMode;

void cleararea(float xl, float yl, float xh, float yh, int kolor)
{
	/* clear an aread of the screen. */
	float x[4], y[4];
	x[0] = xl ;
	y[0] = yl;
	x[1] = xh ;
	y[1] = yl;
	x[2] = xh ;
	y[2] = yh;
	x[3] = xl ;
	y[3] = yh;
	newpen(kolor);
	shadep(4,x,y);
	newpen(1);
}

void do_showperscl(int i)
{
	/* This puts in the individual periods - this is useful for the
	   log scale plots so that individual periods used are seen. However
	   for the period axis of the amplitude scale, this may be a little
	   too compressed since that scale often has a wider range of periods
	   than for velocity since sacmft96 worries about 
	   the USER1 and USER2 fields
	*/
	int j;
	float per, vel, xval,yval;
	/*
	float amp;
	*/
	newpen(4);
	for(j=0 ; j < num_period; j++){
		if(!XAxisPeriod)
			per = 1.0/period[j];
		else
			per = period[j];
		vel = mftdsp[0].vel;
		/*
		amp = mftdsp[0].amp;
		*/
		do_scale(per,vel,&xval,&yval,i);
		if(xval >= pxl[i] && xval <= pxh[i]){
			plot(xval,pyl[i]    ,3);
			plot(xval,pyl[i]+0.2,2);
			plot(xval,pyh[i]    ,3);
			plot(xval,pyh[i]-0.2,2);
		}
	}
	newpen(1);
}

void do_showauto()
{
	/* show the status of the autofeed at the top of the figure */
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	mgwrtxt(4.2,7.7,"  Autofeed  ",1,2);
	if(do_auto_feed == ON)
		mgwrtxt(4.2,7.4,"Autofeed On ",1,1);
	else
		mgwrtxt(4.2,7.4,"Autofeed Off",1,1);
}

void do_showpk()
{
	char ostr[20];
	/* display the picking method at the top of the figure */
	gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
	mgwrtxt(2.8,7.7,"   Action   ",1,2);
	if(InPick == ON)
		mgwrtxt(2.8,7.4," Picking On ",1,1);
	else if(InLine == ON)
		mgwrtxt(2.8,7.4,"Auto Picking",1,1);
	else
		mgwrtxt(2.8,7.4," Picking Off",1,1);
	if(do_auto_feed == ON ){
		do_auto_curfile = (do_auto_curpage-1)*10 + (do_auto_i_next-3);
		sprintf(ostr,"Trace %d/%d",do_auto_curfile,do_auto_nfiles);
		mgwrtxt(6.5,7.7,ostr,1,1);
		sprintf(ostr,"ihdr11=%d",do_auto_ihdr11);
		mgwrtxt(8.0,7.7,ostr,1,1);

		sprintf(ostr,"Page %d/%d",do_auto_curpage,do_auto_npages);
		mgwrtxt(6.5,7.4,ostr,1,1);
		sprintf(ostr,"File %d/%d",do_auto_i_next-3,do_auto_i_max-2);
		mgwrtxt(8.0,7.4,ostr,1,1);
	}
}

void do_showmd(int mode)
{
	/* display the current mode value at the top of the figure */
	mgwrtxt(1.8,7.7,"Mode\0",1,2);
	if(mode < 0)
		mgwrtxt(1.8,7.4,"None\0",1,1);
	else
		mgwrtxt(1.8,7.4,md[mode].str,1,1);
}

/* display picked values */
void do_showpick1(int num_disp, struct disp *dsp, int porg)
{
	int j;
	float xval, yval, vel, amp, per;
	for(j=0;j<num_disp;j++){
		if(dsp[j].mode >= 0){
			if(!XAxisPeriod)
				per = 1.0/dsp[j].instper;
			else
				per = dsp[j].instper;

			vel = dsp[j].vel;
			if(porg == MFT) {
				amp = dsp[j].amp;
			} else if(porg == PHV) {
				amp = dsp[j].uvel;
 /*newpen(2);*/
			} else {
				amp = dsp[j].amp;
			}
			/* display velocity value */
			do_scale(per,vel,&xval,&yval,1);
			if(inside(xval,yval,pxl[1],pyl[1],pxh[1],pyh[1])){
				gsolid(xval,yval,0.02*scl[1],dsp[j].mode);
			}
			/* display spectral amplitude value */
			do_scale(per,amp,&xval,&yval,0);
			if(inside(xval,yval,pxl[0],pyl[0],pxh[0],pyh[0])){
				gsolid(xval,yval,0.02*scl[0],dsp[j].mode);
			}
		}
	}
}

void do_resetzoom(void)
{
float xfrac, yfrac;
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
}

#define	NUMVEL	4
void do_time()
{
	/* draw time indicators to the trace 
	   These red lines indicate the relation between
	   group velocity and the time on the trace
	*/	
	float t;
	float y1, y2, dy;
	float py1, py2;
	float y;
	float per, vel;
	int i;
	y2 = pyh[1];	/*     top of figure */
	y1 = pyl[1];	/*  bottom of figure */
	dy = (y2 - y1)/(NUMVEL );
	do_clearscale(2);
	do_box( pxl[2], pyl[2], pxh[2], pyh[2]);
	newpen(2);
	for(i=0 ; i <= NUMVEL ; i++){
		y = y1 + i*dy;
		do_invscale(&per, &vel, pxh[1], y, 1);
		t = mftdsp[0].dist/vel;
		py1 = pinterp(t,apyl[2],apyh[2]);
		py2 = (1.0 - py1)*pyl[2] + py1*pyh[2];
		if(py2 >= pyl[2] && py2 <= pyh[2]){
			plot(pxh[1]+0.05,y,3);
			plot(pxl[2]-0.15,py2,2);
			plot(pxl[2]-0.05,py2,2);

		}
		
	}
	newpen(1);
}

void do_clearscale(int i )
{
	/* i selects the image - 
		for mft96.ctl  0  - amplitude spectra
		               1  - group velocity
		               2  - time series
		for phv96.ctl  1  - phase velocity
	*/
	float xl, yl, xh, yh;
	/* clear y -axis on left */
	if(i == 0)
		xl = 0.0;
	else
		xl = pxh[i-1];
	xl = xl + 0.05;
	yl = pyl[i] - 0.10;
	xh = pxl[i];
	yh = pyh[i] + 0.10;
	newpen(0);
	shader(xl, yl, xh, yh, 0, 0, 0.01, 0.01);
	/* clear y-axis on right */
	xl = pxh[i] +0.0 ;
	yl = pyl[i] - 0.10;
	xh = pxh[i] + 0.15;
	yh = pyh[i] + 0.10;
	newpen(0);
	shader(xl, yl, xh, yh, 0, 0, 0.01, 0.01);
	/* clear x-axis */
	xl = pxl[i];
	xh = pxh[i] + 0.05;
	yl = pyl[i] - 0.40;
	yh = pyl[i];
	shader(xl, yl, xh, yh, 0, 0, 0.01, 0.01);
	newpen(1);
}

void do_axes(int i, int j )
{
	/* i selects the image in array - j selects the caption
		since image and caption differ
                               i
		for mft96.ctl  0  - amplitude spectra
		               1  - group velocity
		               2  - time series
		for phv96.ctl  0  - group velocity
			       1  - phase velocity
             j=0 label with A, j=1 label with U j=2 label with C
`	*/
	float xx, yy, xlen, ylen, xlow, xhgh, ylow, yhgh;
	xx = pxl[i];
	yy = pyl[i];
	xlen = pxh[i] - pxl[i];
	ylen = pyh[i] - pyl[i];
	xlow = apxl[i];
	ylow = apyl[i];
	xhgh = apxh[i];
	yhgh = apyh[i];
	newpen(1);
	if(xlnlg[i] == LIN){
		if(XAxisPeriod){
		dolinx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,10,"Period (s)");
		dolinx(xx,yy+ylen,xlen,xhgh,xlow,0.10,OFF,ON,OFF,10,"Period (s)");
		} else {
		dolinx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,14,"Frequency (Hz)");
		dolinx(xx,yy+ylen,xlen,xhgh,xlow,0.10,OFF,ON,OFF,14,"Frequency (Hz)");
		}
	} else {
		if(XAxisPeriod){
		dologx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,10,"Period (s)");
		dologx(xx,yy+ylen,xlen,xhgh,xlow,0.10,OFF,ON,OFF,10,"Period (s)");
		} else {
		dologx(xx,yy,xlen,xhgh,xlow,0.10,OFF,OFF,ON,14,"Frequency (Hz)");
		dologx(xx,yy+ylen,xlen,xhgh,xlow,0.10,OFF,ON,OFF,14,"Frequency (Hz)");
		}
	}
	if(ylnlg[i] == LIN){
		doliny(xx , yy ,ylen,yhgh,ylow,0.07,ON,ON,ON,0," ");
		doliny(xx+xlen,yy,ylen,yhgh,ylow,0.07,OFF,ON,OFF,0," ");
	} else {
		dology(xx , yy ,ylen,yhgh,ylow,0.07,ON,ON,ON,0," ");
		dology(xx+xlen,yy,ylen,yhgh,ylow,0.07,OFF,OFF,OFF,0," ");
	}
	if(j == 0){
		if(Units == 0)
			gcent(pxl[i]+0.5*xlen,pyh[i]+0.1,0.15,"A(count-sec)",0.0);
		else
			gcent(pxl[i]+0.5*xlen,pyh[i]+0.1,0.15,"A(cm-sec)",0.0);
	} else if ( j == 1){
		gcent(pxl[i]+0.5*xlen,pyh[i]+0.1,0.15,"U(km/sec)",0.0);
	} else if ( j == 2){
		gcent(pxl[i]+0.5*xlen,pyh[i]+0.1,0.15,"C(km/sec)",0.0);
	}
}

void do_box(float xl,float yl,float xh,float yh)
{
	/* draw a box */
	plot(xl,yl,3);
	plot(xh,yl,2);
	plot(xh,yh,2);
	plot(xl,yh,2);
	plot(xl,yl,2);
}

float pinterp(float x,float x1,float x2)
{
	/* linear interpolation parameter */
	/* get the value of p in   x = (1-p)*x1 + p x2 */
	if(x1 == x2){
		return (0.0);
	} else {
		return (  (x - x1) / ( x2 - x1) );
	}
}

int do_mode(void)
{
	int mode;
	int i, cmd, nmd, ii ;
	char c[2];
	float xv, yv;
	/* output a mode menu for mode selection */
	show_menu(1.0, 0.25, md, sizeof(md), &nmd);
	/* place current mode at top of page */
	mode = -1 ;
	cmd = -1;
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	for(; cmd < 0 ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nmd ; i++){
			if(inside(xv,yv,md[i].xl,md[i].yl,md[i].xh,md[i].yh))
			{
				cmd = md[i].action;
				gmesg(md[i].str);
				ii = i;
				iMode = ii;
				break;
			}
		}
		if(cmd >= 0){
                        mode = iMode ;
			AutoMode = iMode;
		}
	}
	/* clear the window */
	gclip("off", XL, YL, XH, YH);
	newpen(0);
	shader(md[0].xl,md[0].yl,
	    md[nmd-1].xh,md[nmd-1].yh, 0, 0, 0.01, 0.01);
	/* put the mode title at the top */
	do_showmd(mode);

#ifdef DEBUG
fprintf(stderr,"Mode = %d  tmft_subs 409 iMode %d\n",Mode,iMode);
#endif
	gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
	return (mode);
}

void do_invscale(float *per,float *vel,float xval,float yval,int i)
{
	/* do the mapping from screen coordinates to user coordinates */
	float xfac, yfac;
	if(xlnlg[i] == LIN){
		xfac = (pxh[i] - pxl[i])
		    / (apxh[i] - apxl[i]);
		*per = apxl[i] + (xval - pxl[i])/xfac ;
	} else{
		xfac = (pxh[i] - pxl[i])
		    / log(apxh[i] / apxl[i]);
		*per = apxl[i] * 
		    exp((xval - pxl[i])/xfac);
	}
	if(ylnlg[i] == LIN){
		yfac = (pyh[i] - pyl[i])
		    / (apyh[i] - apyl[i]);
		*vel = apyl[i] + (yval - pyl[i])/yfac ;
	} else{
		yfac = (pyh[i] - pyl[i])
		    / log(apyh[i] / apyl[i]);
		*vel = apyl[i] * 
		    exp((yval - pyl[i])/yfac);
	}
}

void do_scale(float per,float vel,float *xval,float *yval,int i)
{
/*
	per	FLOAT	x-axis value (either period or frequency)
	vel	FLOAT	y-axis velocity value or spectral amplitude value
	*xval	FLOAT	returned screen x-coordinate
	*yval	FLOAT	returned screen y-coordinate
	i	INT	figure index, e.g., 0 for amp vs per, 1 for vel vs per
*/
	/* do the mapping between user coordinates and screen coordinates */
	float xfac, yfac;
	if(xlnlg[i] == LIN){
		xfac = (pxh[i] - pxl[i])
		    / (apxh[i] - apxl[i]);
		*xval = pxl[i] + (per - apxl[i])
		    * xfac;
	} else{
		xfac = (pxh[i] - pxl[i])
		    / log(apxh[i] / apxl[i]);
		*xval = pxl[i] + log(per / apxl[i])
		    * xfac;
	}
	if(ylnlg[i] == LIN){
		yfac = (pyh[i] - pyl[i])
		    / (apyh[i] - apyl[i]);
		*yval = pyl[i] + (vel - apyl[i])
		    * yfac ;
	} else{
		yfac = (pyh[i] - pyl[i])
		    / log(apyh[i] / apyl[i]);
		*yval = pyl[i] + log(vel / apyl[i])
		    * yfac;
	}
}

void do_clear(int i )
{
	/* clear the graphic for the i'th figure */
	gclip("off", pxl[i],pyl[i],pxh[i],pyh[i]);
	newpen(0);
	shader(pxl[i], pyl[i], pxh[i], pyh[i], 0, 0, 0.01, 0.01);
	newpen(1);
}

void do_line(float xv1, float yv1, int num_disp, struct disp *dsp, int porg)
{
	char c[2];
	float xv2, yv2;
	float x1,x2,y1,y2;
	float ypred;
	float per, vel, amp;
	float dataper;
	float xval, yval;
	float xxval;
	int i,j;
	int ii;
	float sum_min_x, sum_min_y;
	float sum_x, sum_y;
	float xs, ys, as, ps; 
	/*
	float vs;
	*/

	gcursor("Rubber");
	curaxy(&xv2, &yv2, c);
	if(Cursor == 0)gcursor("XORArrow");else gcursor("Cross");
	x1 = MIN(xv1,xv2) ;
	x2 = MAX(xv1,xv2) ;
	y1 = MIN(yv1,yv2) ;
	y2 = MAX(yv1,yv2) ;
	if(x1==x2 && y1 == y2) return;  /* need two points for a line */
	if( inside(x1,y1,pxl[1]-0.1,pyl[1],
	    pxh[1]+0.1, pyh[1])
	    && inside(x2,y2,pxl[1]+0.1,pyl[1],
	    pxh[1]+0.1, pyh[1]) ) {
		newpen(3001);
		/* scan through dispersion file and line nearest to line */
		/* since we are plotting according to the filter periods */
		for(j = 0 ; j < num_period ; j++){
			/* rest goodness of fit */
			sum_min_x = 1.0e+38; 
			sum_min_y = 1.0e+38; 
			if(!XAxisPeriod)
				per = 1.0/period[j];
			else
				per = period[j];
			vel = 1.0;
			do_scale(per,vel,&xxval,&yval,1);
			if(xxval >= x1 && xxval <= x2){
				ii = -1;
				/* this period is within the cursor range */
				/* find the data period closest to this value */
				for(i=0; i< num_disp ; i++){
					sum_x = ABS(per - dsp[i].per);
					if(sum_x < sum_min_x){
						sum_min_x = sum_x;
						dataper = dsp[i].per;
					}
				}

				/* now make a line prediction */
				/* beware of divide by zero */
				ypred = yv1 + (xxval -xv1)*(yv2-yv1)/(xv2-xv1); 
				/*we now find the point closest to this period*/
				for(i=0; i< num_disp ; i++){
					per = dsp[i].per;
					vel = dsp[i].vel;
					if(porg == MFT) {
						amp = dsp[i].amp;
					} else if(porg == PHV) {
						amp = dsp[i].uvel;
					} else {
						amp = dsp[i].amp;
					}
					if(ABS(per - dataper) < 0.001*dataper){
					do_scale(per,vel,&xval,&yval,1);
					if(xval >= x1 && xval <= x2){
						sum_y = ABS(yval-ypred); 
						if(sum_y < sum_min_y){
							sum_min_y = sum_y ;
							ii = i;
							xs = xval; 
							ys = yval;
							ps = per; 
							/*
							vs = vel; 
							*/
							as = amp;
						}
					}
					}
/*
*/
				}
				/* now plot only if we have a hit */
				if(ii >= 0 && ii < num_disp && sum_min_y < 0.10){
					if(dsp[ii].mode < 0){
					/* put up the velocity point on the right figure */
					gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
					if(porg == PHV) {
						newpen(Pick_Color);
 /*newpen(2);*/
					} else {
						newpen(Pick_Color);
					}
					gsolid(xs,ys,0.02*scl[1],dsp[ii].symb);
					gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
					gclip("on", pxl[0],pyl[0],pxh[0],pyh[0]);
					/* put up the amplitude point on left figure for MFT PMF and 
						group velocity for PHV */
					if(porg == PHV) {
						newpen(0);
					} else {
						newpen(2);
					}
					do_scale(ps,as,&xval,&yval,0);
					gsolid(xval,yval,0.02*scl[0],dsp[ii].symb);
					gclip("off", pxl[0],pyl[0],pxh[0],pyh[0]);
					gclip("on", pxl[1],pyl[1],pxh[1],pyh[1]);
					}
					dsp[ii].mode = Mode;
				}
			}
		}
		/* turn off XOR */
		newpen(3000);
	}
}

void do_pick( float xv1, float yv1, int lev, int num_disp, struct disp *dsp, int porg)
{
	/* from the screen coordinates, (xv1, yv1) 
		obtain the nearest point in the dispersion table 
	   lev	= 0 do with original list
		= 1 work with the already picked list and inst period
	*/
	float xval, yval;
	int i;
	float sum_min;
	float sum;
	int ii = -1;
	float per, vel, amp;
	float xs = 0.0, ys = 0.0;
	float ps = 0.0; 
	/* float vs = 0.0; */
	float as = 1.0;
	/* convert screen coordinates to the (period,velocity) pair */
	sum_min = 1.0e+38; 

	/* scan through dispersion file and line nearest to line */
	/* find the closet x-value first, then the closest y-value */
	for(i=0; i< num_disp ; i++){

#ifdef DEBUG
fprintf(stderr,"%s %s %s %f%f %f %d \n",dsp[i].type,dsp[i].lvry,dsp[i].cug,dsp[i].vel,dsp[i].phase,dsp[i].uvel,dsp[i].nn);
#endif
		if(lev == 0) {
			per = dsp[i].per;
			vel = dsp[i].vel;
			if(porg == MFT) {
				amp = dsp[i].amp;
			} else if(porg == PHV) {
				amp = dsp[i].uvel;
			} else {
				amp = dsp[i].amp;
			}
		} else {
                        if(!XAxisPeriod)
                                per = 1.0/dsp[i].instper;
                        else
                                per = dsp[i].instper;
			vel = dsp[i].vel;
			if(porg == MFT) {
				amp = dsp[i].amp;
			} else if(porg == PHV) {
				amp = dsp[i].uvel;
			} else {
				amp = dsp[i].amp;
			}
			if( dsp[i].mode != ModeSelect)
				per = -1.0;
		}
		if(per >= 0.0){
			do_scale(per,vel,&xval,&yval,1);
			sum = (xval-xv1)*(xval -xv1)+(yval-yv1)*(yval-yv1); 
			if(sum < sum_min){
				sum_min = sum ;
				ii = i;
				xs = xval; 
				ys = yval;
				ps = per; 
				/* vs = vel; */
				as = amp;
			}
		}

	}
#ifdef DEBUG
fprintf(stderr,"[%f,%f] [%f,%f]\n",pxl[1],pyl[1],pxh[1],pyh[1]);
fprintf(stderr,"[%f,%f] [%f,%f]\n",apxl[1],apyl[1],apxh[1],apyh[1]);
#endif
	if(ii >= 0 && ii < num_disp){
		/* plot a new XOR symbol there */
		newpen(3001);
		/* put up the velocity  point on right figure */
		if(porg == PHV){ 
			newpen(Pick_Color);
 /*newpen(2);*/
			gsolid(xs,ys,0.02*scl[1],dsp[ii].symb);
		} else {
			newpen(Pick_Color);
			gsolid(xs,ys,0.02*scl[1],dsp[ii].symb);
		}
		/* put up the amplitude point */
		gclip("off", pxl[1],pyl[1],pxh[1],pyh[1]);
		gclip("on" , pxl[0],pyl[0],pxh[0],pyh[0]);
		/* put up the amplitude point on left figure for MFT PMF and 
			group velociuty for PHV */
		if(porg == PHV)
			newpen(0);
		else
			newpen(2);
		do_scale(ps,as,&xval,&yval,0);
		gsolid(xval,yval,0.02*scl[0],dsp[ii].symb);
		gclip("off", pxl[0],pyl[0],pxh[0],pyh[0]);
		newpen(2);
		/* turn off XOR */
		newpen(3000);
		/* toggle the assignment */
		if(dsp[ii].mode > -1)
			dsp[ii].mode = -1;
		else
			dsp[ii].mode = Mode;
		gclip("on" , pxl[1],pyl[1],pxh[1],pyh[1]);
	}

}
