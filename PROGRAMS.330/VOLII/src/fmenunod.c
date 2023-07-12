#include "nfmenu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
fmenu *Start, *End;

/* linked list code
	ammeraal, l,
	programs and data structures in c, second edition,
		based on ansi c and c++,
	john wiley & sons, chichester
	1992.
	272 pp

	pages 111-112

Changes:
	10 JUN 2014 - I previously defined
	char kstnm[9], kcmonm[9], datetime[24]
        and char *str in nfmenu.h
        These are all char*
        I use strdup to assign the sizes of these in
        appendnode and insertnode.  gcc 4.8 and 4.9 did not
        like the old code of   strcpy(p->kstnm, kstnm)
        of the memory allocation. The idea of strdup is
	from Kernighan and Ritchie, The C Programming Language, 2nd edition,
	1988, Prentice Hall, pages 139-143
      
*/

void *getmemory(int n)
{
	void *p=malloc(n);
	if(p == NULL){
		fprintf(stderr,"not enough memory\n");
		exit(1);
	}
	return p;
}

fmenu *getnode(void)
{
	return (fmenu *)getmemory(sizeof (fmenu));
}

void deletenode(fmenu *p)
{
	fmenu *q;
	free(p->str);
	q = p->next;
	if(q == End)
		End = p;
	else
		*p = *q;
	free(q);
}

void insertnode(fmenu *p, float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz, int ihdr11)
{
        kstnm[8] = '\0' ;
        kcmpnm[8] = '\0' ;
        datetime[23] = '\0' ;
	fmenu *q;
	q = getnode();
	if(p == End)
		End = q;
	else
		*q = *p;
	p->next = q;
	p->xl = xl;
	p->yl = yl;
	p->xh = xh;
	p->yh = yh;
	p->str = strdup(str);
	p->kstnm = strdup(kstnm);
	p->kcmpnm = strdup(kcmpnm);
	p->datetime = strdup(datetime);
	p->action = action;
	p->lstrmx = lstrmx;
	p->type   = type  ;
	p->line   = line  ;
	p->fsize  = fsize ;
	p->nsamp  = nsamp ;
	p->page   = page  ;
	p->used   = used  ;
	p->dist = dist;
	p->az = az;
	p->baz = baz;
	p->ihdr11 = ihdr11;
}


void appendnode(float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz, int ihdr11)
{
        kstnm[8] = '\0' ;
        kcmpnm[8] = '\0' ;
        datetime[23] = '\0' ;
	fmenu *p=End;
	End = getnode();
	p->next = End;
	p->xl = xl;
	p->yl = yl;
	p->xh = xh;
	p->yh = yh;
	p->str = strdup(str);
	p->kstnm = strdup(kstnm);
	p->kcmpnm = strdup(kcmpnm);
	p->datetime = strdup(datetime);
	p->action = action;
	p->lstrmx = lstrmx;
	p->type   = type  ;
	p->line   = line  ;
	p->fsize  = fsize ;
	p->nsamp  = nsamp ;
	p->page  = page ;
	p->used   = used  ;
	p->dist = dist;
	p->az = az;
	p->baz = baz;
	p->ihdr11 = ihdr11;
}



