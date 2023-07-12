/* STRUCTURE DEFINITIONS */

struct disp {
	char 	type[7]	;	/* identification of type 	*/
	char	lvry[2]	;	/* A undef L Love R Rayleigh	*/
	char	cug[2]	;	/* U = group, C= phase, G= gamma */
	int	mode	;	/* -1 not defined, 0 = Fund, 1 = first */
	float	per	;	/* period seconds */
	float	vel	;	/* velocity value group for mft96.disp
					phase velocity for phv96.disp	*/
	float	evel	;	/* error velocity value	*/
	float	amp	;	/* spectral amplitude */
	float	dist	;	/* distance */
	float	az	;	/* source->receiver azimuth */
	float	lat1	;	/* event latitude */
	float	lon1	;	/* event longitude */
	float	lat2	;	/* station latitude */
	float	lon2	;	/* station longitude */
	int	ictl	;	/* control flag */
	int	symb	;	/* symbol used in plot */
	float	instper	;	/* instantaneous period */
	float	alpha	;	/* Gaussian filter parameter */
	char	comment[10];	/* comment field */
	char	kstnm[9];	/* station name */
	char	kcmpnm[9];	/* componet name */
	int	nzyear	;	/* year */
	int	nzjday	;	/* julian day */
	int	nzhour	;	/* hour */
	int	nzmin	;	/* minute */
	float	phase	;	/* phv96.disp phase at peak */
	float   uvel    ;	/* phv96.disp group velocity */
	int	nn	;       /* phv96.disp multiple of 2 PI */
};

struct pos {
	float xl;
	float xh;
	float yl;
	float yh;
	float axl;
	float axh;
	float ayl;
	float ayh;
	char xlnlg[4];
	char ylnlg[4];
	char fname[20];
} ;

/* DEFINES */
#define ON      1
#define OFF     0
#define YES     1
#define NO      0
/* plot window for MFT96 graphs
           The plot will be placed between (XL,YL) and (XH,YH)
           The MFT96.PLT will be shifted by adding (XOFF,YOFF) */
#define XL      0.0
#define XH      10.00
#define YL      1.5
#define YH      8.0
#define XOFF    0.0
#define YOFF    -0.5;
#define MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (a) > (b) ? (a):(b) )
#define ABS(a)   ( (a) > (0) ? (a): -(a) )
#define LIN 0
#define LOG 1
#define MFT 0
#define PHV 1
#define PMF 2

