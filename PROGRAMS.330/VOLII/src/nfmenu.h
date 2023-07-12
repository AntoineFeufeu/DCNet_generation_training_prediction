/* Changes
	10 JUN 2014 - char kstnm[9], kcmpnm[9], datetime[24]
		replaced by char*
*/
#ifndef _NFMENU_H
#define _NFMENU_H
typedef struct Fmenu { 
	float xl;
	float yl;
	float xh ;
	float yh ;
	char *str ;
	int action;
	int lstrmx;
	int type;
	int line;
	int fsize;
	int nsamp;
	char *kstnm;
	char *kcmpnm;
	char *datetime;
	int page;
	int used;
	float dist;
	float az;
	float baz;
	int ihdr11;
	struct Fmenu *next;
} fmenu;

void *getmemory(int h);
fmenu *getnode(void);
void deletenode(fmenu *p);
void insertnode(fmenu *p, float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz, int ihdr11);
void appendnode(float xl, float yl, float xh, float yh,
	char *str, int action, int lstrmx, int type, int line, int fsize,
	int nsamp, char *kstnm, char *kcmpnm, char *datetime,
	int page, int used, float dist, float az, float baz, int ihdr11);
#endif
