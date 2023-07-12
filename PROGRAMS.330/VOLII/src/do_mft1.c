
#include	"nfmenu.h"
#include	"nmenu.h"
#include	"do_mft.h"
#include	"calplot.h"
#include	<stdio.h>
#include	<string.h>

void show_menu (float x0, float y0, struct menu *men, int size, int *nm);
void noshow_menu (float x0, float y0, struct menu *men, int size, int *nm);
int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);
void clearregion(float xl, float yl, float xh, float yh);
int do_page2(char *fname);

#define		MENU_FILE_QUIT	-1
#define		MENU_FILE_NEXT	-2
#define		MENU_FILE_PREV	-3
#define		MENU_FILE_NUMB	10

static char strout[80];
extern fmenu **file_menu;
fmenu *q;
extern int ndfiles;
char ostrs[80]; 
char ostrr[80]; 

extern int HasMouse; 
extern float XminDev, YminDev, 
	XmaxDev, YmaxDev, XminClip, 
	YminClip, XmaxClip, YmaxClip;
extern int Color;
extern int numperpage;

struct menu menu_p1[] = {
	{  -1.0, -1.0, -1.0, -1.0, "Quit\0" , MENU_FILE_QUIT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " <  \0" , MENU_FILE_PREV, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, " >  \0" , MENU_FILE_NEXT, -1, 1, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1},
	{  -1.0, -1.0, -1.0, -1.0, "    \0" , 10, -1, 0, -1}
};

extern int do_auto_feed ;
extern int do_auto_curfile;
extern int do_auto_nfiles;
extern int do_auto_curpage;
extern int do_auto_npages;
extern int do_auto_i_next ;
extern int do_auto_i_max ;
extern int do_auto_ihdr11 ;

int do_page1(int npage, int *curpage)
{
	int i, cmd, nmd ;
	char c[2];
	float xv, yv;
	float xl, yl, xh, yh;
	int cur_page;
	int ret;
	/* output a mode menu for mode selection */
	cur_page = *curpage;  
/* of auto need no plot but still need coordinate array */
	if(do_auto_feed == YES){
		noshow_menu(1.0, 0.5, menu_p1,sizeof(menu_p1),&nmd);
	} else {
		show_menu(1.0, 0.5, menu_p1,sizeof(menu_p1),&nmd);
	}
	/* place current mode at top of page */
	cmd = -1;
	sprintf(strout,"Page %d of %d",cur_page+1, npage+1);
	gmesg(strout);
	for(; cmd < 0 ;){
                if (  do_auto_feed == OFF ) {
			curaxy(&xv, &yv, c);
		} else {
			do_auto_curpage = cur_page + 1 ;
			do_auto_npages = npage + 1 ;
			if(do_auto_curpage < do_auto_npages)
				do_auto_i_max = nmd -3 + 1 ;
			else
				if ( do_auto_nfiles%10 == 0)
					do_auto_i_max = nmd -3 + 1 ;
				else
					do_auto_i_max = do_auto_nfiles%10  + 2 ;
			/* attempt to implement autofeed by emulating curaxy() */
			if(do_auto_i_next < 0)
				do_auto_i_next = 3;
			i = do_auto_i_next ;
			xl = menu_p1[i].xl;
			yl = menu_p1[i].yl;
			xh = menu_p1[i].xh;
			yh = menu_p1[i].yh;
			xv = (xl + xh ) /2;
			yv = (yl + yh ) /2;
/* 30 MAR 2017 - when we get to the end of a menu e.g., and last page does not fill
       all 10 items, then the next coordinate will not have a valid box and
       the program would stop. We check this here so that we wrap around to the 
       first page. This is why we had wrap around for 10 items on last page but
       program termination for less than 10  */
			if(xl < 0.0){
				do_auto_i_next = -1;
				cur_page++;
				if(cur_page >npage )cur_page= 0;
				*curpage = cur_page;
				return(1);
			}
		}

		cmd = -1;
		/* determine mouse position */
		for(i=0 ; i < nmd ; i++){
			xl = menu_p1[i].xl;
			yl = menu_p1[i].yl;
			xh = menu_p1[i].xh;
			yh = menu_p1[i].yh;
			if(inside(xv,yv,xl,yl,xh,yh)) {
				cmd = menu_p1[i].action;
				if(cmd == MENU_FILE_QUIT){
					gframe(1);
					return(-1);
				} else if(cmd == MENU_FILE_PREV){
					cur_page--;
					if(cur_page < 0)cur_page=npage ;
					*curpage = cur_page;
					gframe(2);
					ginfo(&HasMouse, &XminDev, &YminDev, 
						&XmaxDev, &YmaxDev, &XminClip, 
						&YminClip, &XmaxClip,&YmaxClip,&Color);
					return (1);
				} else if(cmd == MENU_FILE_NEXT){
					cur_page++;
					if(cur_page >npage )cur_page= 0;
					*curpage = cur_page;
					gframe(2);
					ginfo(&HasMouse, &XminDev, &YminDev, 
						&XmaxDev, &YmaxDev, &XminClip, 
						&YminClip, &XmaxClip,&YmaxClip,&Color);
					return (1);
				} else if(cmd >= 1 && cmd <= 10){
					*curpage = cur_page;
					gframe(2);
					ginfo(&HasMouse, &XminDev, &YminDev, 
						&XmaxDev, &YmaxDev, &XminClip, 
						&YminClip, &XmaxClip,&YmaxClip,&Color);
					q = file_menu[i];
					do_auto_i_next = i + 1 ;

					if (do_auto_i_next >= nmd){
						/* reset */
						do_auto_i_next = -1;
						cur_page++;
						if(cur_page >npage )cur_page= 0;
						*curpage = cur_page;
						gframe(2);
						ginfo(&HasMouse, &XminDev, &YminDev, 
							&XmaxDev, &YmaxDev, &XminClip, 
							&YminClip, &XmaxClip,&YmaxClip,&Color);
						return (1);
					}

					do_auto_ihdr11 = q->ihdr11 ;
					ret=do_page2(q->str); 
					q->used = ret;
					if(ret > 1){
					/* insert nodes for the matched and residual files */
						strcpy(ostrs,file_menu[i]->str);
						strcpy(ostrr,file_menu[i]->str);
						strcat(ostrs,"s");   
						strcat(ostrr,"r");   
						insertnode(q, xl, yl, xh, yh,
        						ostrs,  q->action,  q->lstrmx,
							q->type,  q->line,  q->fsize,
        						q->nsamp, q->kstnm, q->kcmpnm, 
							q->datetime, q->page, q->used, q->dist, q->az, q->baz, q->ihdr11) ;
						insertnode(q, xl, yl, xh, yh,
        						ostrr,  q->action,  q->lstrmx,
							q->type,  q->line,  q->fsize,
        						q->nsamp, q->kstnm, q->kcmpnm, 
							q->datetime, q->page, q->used, q->dist, q->az, q->baz, q->ihdr11) ;
					}
					gframe(1);
					return(1);
				}
				break;
			}
		}
	}
	/* never get here but we must have a return */
	return (-1);
}
