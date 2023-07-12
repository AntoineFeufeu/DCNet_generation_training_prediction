/* Changes
	01 AUG 2007 - corrected error in setting map coordinates 
              when the Station of Epicenter is not plotted, e.g.,
		st off or ep off
	30 OCT 2022 - treat coordinates as doubles because of NVHDR = 7
*/
#include	<stdio.h>
#include	"gsac_docommand.h"
#include        "gsac.h"
#include        "gsac_plot.h"
#include        "gsac_sac.h"
#include        "gsac_arg.h"
#include        "gsac_sachdr.h"

extern struct sacfile_ *sacdata;
extern int *sortptr;


#define	MAP5_DFLT	0
#define	MAP5_NORTH	1
#define	MAP5_EAST	2
#define	MAP5_SOUTH	3
#define	MAP5_WEST	4
#define	MAP5_GLOBAL	5
#define	MAP5_TOPO	6
#define	MAP5_STA		7
#define	MAP5_EPI		8
#define MAP5_RAY		9
#define MAP5_KSTNM	10

static int map5_userlatlon = NO ;
static float map5_minlat = -12345. ;
static float map5_maxlat = -12345. ;
static float map5_minlon = -12345. ;
static float map5_maxlon = -12345. ;
static int map5_dotopo = NO ;
static int map5_dosta = YES ;
static int map5_doepi = YES ;
static int map5_doglobal = NO ;
static int map5_doray = NO ;
static int map5_dokstnm = NO ;

static FILE *gmtscript;

struct arghdr map5arg[] = {
	{MAP5_DFLT, "DEFAULT" , IHDR, NO, 0, NO, "", 1},
	{MAP5_NORTH , "NORTH" , RHDR, NO, 1, NO, "North maxlat", 1},
	{MAP5_EAST  , "EAST"  , RHDR, NO, 1, NO, "East  maxlon", 2},
	{MAP5_EAST  , "E"     , RHDR, NO, 1, NO, "East  maxlon", -1},
	{MAP5_WEST  , "WEST"  , RHDR, NO, 1, NO, "West  minlon", 1},
	{MAP5_SOUTH , "SOUTH" , RHDR, NO, 1, NO, "South minlat", 2},
	{MAP5_SOUTH , "S"     , RHDR, NO, 1, NO, "South minlat", -1},
	{MAP5_TOPO  , "TOPO"  , YHDR, NO, 1, NO, "Topo [ON|OFF] ", 1},
	{MAP5_STA   , "STA"   , YHDR, NO, 1, NO, "STa [ON|OFF] ", 2},
	{MAP5_EPI   , "EPI"   , YHDR, NO, 1, NO, "EPi [ON|OFF] ", 2},
	{MAP5_GLOBAL, "GLO"   , YHDR, NO, 1, NO, "Global [ON|OFF] ", 1},
	{MAP5_RAY   , "RAY"   , YHDR, NO, 1, NO, "Raypath [ON|OFF] ", 1},
	{MAP5_KSTNM , "KSTNM" , YHDR, NO, 1, NO, "KStname [ON|OFF] ", 1},
	{0,	""	     , IHDR, NO, 0, NO, "", -1}
};

static void map5_bound(float *minlat, float *maxlat, float *minlon, float *maxlon);
static void map5_tics(float minlat, float maxlat, float minlon, float maxlon, float *lattic, float *lontic, char *res);
static void map5_header(void);
static void map5_preamble(float map5_minlon,float map5_maxlon,
	float map5_minlat,float map5_maxlat, char res,
	float lontic,float lattic);
static void map5_topo(void);
static void map5_coast(int do_overlay, int do_continue);
static void map5_kstnm(void);
static void map5_epicenter(int do_continue);
static void map5_raypath(void);
static void map5_station(int do_continue);
static void map5_end(void);

/* these are temporary variables only used here */
static float map5_real[10];
static int   map5_int [10];
static int   map5_yn;
static int   map5_num;

/* these are prototypes for global variables to be used by the routine */

void gsac_set_param_map5(int ncmd, char **cmdstr)
{
	float t_n, t_s, t_e, t_w;
	float u_n, u_s, u_e, u_w;
	int i;
	int usercoord;
	map5_userlatlon = NO ;
	usercoord = NO;
	/* initial debug */
	for(i=1; i < ncmd; i++)
		printf("%s ",cmdstr[i]);
	printf("\n");
	/* parsing code here */
	if(ncmd == 1)
		return;
	if(testarg(ncmd, cmdstr, map5arg, NO, YES))
		return;
	/* parse commands */
	for(i=0 ; map5arg[i].key[0] != '\0' ; i++){
		if(map5arg[i].used > 0){
			if(map5arg[i].ricell == RHDR){
				getargr(ncmd, cmdstr, map5arg[i].key, 
					map5arg[i].mfit,map5arg[i].narg, map5_real);
			} else if(map5arg[i].ricell == IHDR){
				getargi(ncmd, cmdstr, map5arg[i].key, 
					map5arg[i].mfit,map5arg[i].narg, map5_int );
			} else if(map5arg[i].ricell == YHDR){
				getargyn(ncmd, cmdstr, map5arg[i].key, 
					map5arg[i].mfit,map5arg[i].narg, &map5_yn );
			} else if(map5arg[i].ricell == NHDR){
				getargn(ncmd, cmdstr, map5arg[i].key, 
					map5arg[i].mfit,map5arg[i].narg, &map5_num );
			}
			switch(map5arg[i].id){
				case MAP5_DFLT:
					map5_userlatlon = NO ;
					map5_dotopo = NO ;
					map5_dosta = YES ;
					map5_doepi = YES ;
					map5_dokstnm = NO ;
					map5_doglobal = NO ;
					break;
				case MAP5_TOPO:
					map5_dotopo   = map5_yn ;
					break;
				case MAP5_STA:
					map5_dosta    = map5_yn ;
					break;
				case MAP5_EPI:
					map5_doepi    = map5_yn ;
					break;
				case MAP5_RAY:
					map5_doray    = map5_yn ;
					break;
				case MAP5_KSTNM:
					map5_dokstnm    = map5_yn ;
					break;
				case MAP5_GLOBAL:
					map5_doglobal = map5_yn ;
					break;
				case MAP5_NORTH:
					t_n = map5_real[0];
					usercoord = YES;
					break;
				case MAP5_SOUTH:
					t_s = map5_real[0];
					usercoord = YES;
					break;
				case MAP5_EAST:
					t_e = map5_real[0];
					usercoord = YES;
					break;
				case MAP5_WEST:
					t_w = map5_real[0];
					usercoord = YES;
					break;

			}
		}
	}
	/* safety */
	if(usercoord == YES){
		u_n = MAX(t_n,t_s);
		u_s = MIN(t_n,t_s);
		u_e = MAX(t_e,t_w);
		u_w = MIN(t_e,t_w);
		if (u_n >= -90.0 && u_n <= 90.0
			&& u_s >= -90.0 && u_s <= 90.0
			&& u_e >= -180.0 && u_e <= 180.0
			&& u_w >= -180.0 && u_w <= 180.0 ){
			map5_minlat = u_s ;
			map5_maxlat = u_n ;
			map5_minlon = u_w ;
			map5_maxlon = u_e ;
			map5_userlatlon = YES ;
		} else {
			map5_userlatlon = NO ;
		}
	} else {
		 map5_minlat = -12345. ;
		 map5_maxlat = -12345. ;
		 map5_minlon = -12345. ;
		 map5_maxlon = -12345. ;
	}
			
		
}


void gsac_exec_map5(void)
{
	/* we always create a map5.sh file 

	Because of the way that GMT works we must be careful to
	indicate when overlay occurs, implying both a -O flag and the
		redirection >> map5.eps
 	and whether there is more to be done which then requires a -K flag
	in GMT
	*/
	int ntrchdr;
	float lattic, lontic;
	char res;
	int do_continue, do_overlay;
	/* if there are no traces return */
	ntrchdr = gsac_control.number_iheaders ;
	if(ntrchdr < 1)
		return;

	if( (gmtscript = fopen("map5.sh", "w+")) == NULL){
		fprintf(stderr,"Error in map5: cannot open the file map5.sh\n");
		fprintf(stderr,"Check file/directory permissions\n");
		return;
	}
	/* POSIX CYGWIN UNIX LINUX */
/*
	chmod("map5.sh", S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | 
		S_IXGRP | S_IROTH | S_IXOTH );
*/

	/* pass 1 - bounds */
	/* get the bounds */
	map5_bound(&map5_minlat, &map5_maxlat, &map5_minlon, &map5_maxlon);

	/* define coastline resolution and tic intervals */
	map5_tics(map5_minlat, map5_maxlat, map5_minlon, map5_maxlon, 
		&lattic, &lontic, &res);

	/* output the shell script header */
	map5_header();

	/* output the preamble */

	map5_preamble(map5_minlon,map5_maxlon,map5_minlat,map5_maxlat,
		res,lontic,lattic);
	
	/* do the topography */
	if(map5_dotopo == YES){
		map5_topo();
		do_overlay = YES;
	} else {
		do_overlay = NO;
	}

	/* do the coastlines */
	if(map5_dosta == YES || map5_doepi == YES)
		do_continue = YES;
	else
		do_continue = NO;
	map5_coast(do_overlay, do_continue);

	/* do the raypaths if epicenter and stations are plotted */
	if(map5_doepi == YES && map5_dosta == YES && map5_doray == YES)
		map5_raypath();

	/* do the  epicenters */
	if(map5_doepi == YES){
		if(map5_dosta == YES)
			do_continue = YES;
		else
			do_continue = NO;
		map5_epicenter(do_continue);
	}

	/* do the  stations */
	if(map5_dosta == YES){
		if(map5_dokstnm == YES)
			do_continue = YES;
		else
			do_continue = NO;
		map5_station(do_continue);
	}
	/* do the  station names */
	if(map5_dokstnm == YES)
		map5_kstnm();
	/* final output */
	map5_end();
	/* close the file */
	fclose(gmtscript);
	printf("Execute using the command: sh map5.sh\n");
}

/* determine the lat, lon bounds */
static void map5_bound(float *minlat, float *maxlat, 
	float *minlon, float *maxlon)
{
	double stla, stlo, evla, evlo ;
	float map5_latmax, map5_latmin, map5_lonmax, map5_lonmin ;
	float lattic, lontic;
	char res;
	int k, ntrchdr;
	ntrchdr = gsac_control.number_iheaders ;
	if(map5_userlatlon == NO){
		map5_latmax =  -90.0 ;
		map5_latmin =   90.0 ;
		map5_lonmax = -180.0 ;
		map5_lonmin =  180.0 ;
		for ( k=0 ; k < ntrchdr ; k ++){
			stla = sacdata[k].sachdr.rhdr[H_STLA];
			stlo = sacdata[k].sachdr.rhdr[H_STLO];
			evla = sacdata[k].sachdr.rhdr[H_EVLA];
			evlo = sacdata[k].sachdr.rhdr[H_EVLO];
			if(map5_dosta == YES){
				if(stlo != -12345.0){
					map5_lonmax = MAX(map5_lonmax,stlo);
					map5_lonmin = MIN(map5_lonmin,stlo);
				}
				if(stla != -12345.0){
					map5_latmax = MAX(map5_latmax,stla);
					map5_latmin = MIN(map5_latmin,stla);
				}
			}
			if(map5_doepi == YES){
				if(evla != -12345.0){
					map5_latmax = MAX(map5_latmax,evla);
					map5_latmin = MIN(map5_latmin,evla);
				}
				if(evlo != -12345.0){
					map5_lonmax = MAX(map5_lonmax,evlo);
					map5_lonmin = MIN(map5_lonmin,evlo);
				}
			}
		}
		/* make a minor adjust to extend slightly beyond the
			precise limits - we will use the tic routine 
			to assist */
			map5_tics(map5_latmin, map5_latmax, map5_lonmin, 
				map5_lonmax, &lattic, &lontic, &res);
		*maxlat = map5_latmax + 0.5*lattic ;
		*minlat = map5_latmin - 0.5*lattic ;
		*maxlon = map5_lonmax + 0.5*lontic ;
		*minlon = map5_lonmin - 0.5*lontic ;
		/* safety */
		if( *maxlat >  90.0) *maxlat =  90.0 ;
		if( *minlat < -90.0) *minlat = -90.0 ;
		if( *maxlon >  180.0) *maxlon =  180.0 ;
		if( *minlon < -180.0) *minlon = -180.0 ;
	}
}

static void map5_tics(float minlat, float maxlat, float minlon, float maxlon, float *lattic, float *lontic, char *res)
{
	/* this is hard to do for aethetics */
	float dflat, dflon;
	dflat = maxlat - minlat ;
	dflon = maxlon - minlon ;
	/* define resolution for coast line 
		Using example in GMT_Tutorial.pdf Section K.4 
			The Five Resolutions */
	if( ABS(dflat)+ABS(dflon) > 90.0)
		*res = 'c' ;
	else if( ABS(dflat)+ABS(dflon) > 30.0)
		*res = 'l' ;
	else if( ABS(dflat)+ABS(dflon) > 5.0)
		*res = 'i' ;
	else if( ABS(dflat)+ABS(dflon) > 1.0)
		*res = 'h' ;
	else 
		*res = 'f' ;

	/* now get the ticmarks - make integer multiples of possible */
	if( dflat > 90.0)
		*lattic = 15.0 ;
	else if( dflat > 45.0)
		*lattic = 10.0 ;
	else if( dflat > 10.0)
		*lattic = 5.0 ;
	else if( dflat > 5.0)
		*lattic = 2.0 ;
	else if( dflat > 1.0)
		*lattic = 1.0 ;
	else 
		*lattic = 0.1 ;

	if( dflon > 90.0)
		*lontic = 30.0 ;
	else if( dflon > 45.0)
		*lontic = 15.0 ;
	else if( dflon > 10.0)
		*lontic = 5.0 ;
	else if( dflon > 5.0)
		*lontic = 2.0 ;
	else if( dflon > 1.0)
		*lontic = 1.0 ;
	else 
		*lontic = 0.1 ;
	
}

static char *map5_proto[] = {
"#!/bin/sh\n",
"\n",
"\n",
"#####\n",
"#       define the color palete for elevation\n",
"#       MIN_ELV R G B MAX_ELEV R G B\n",
"#       (from Chuck Ammon)\n",
"#####\n",
"cat > map.cpt << EOF\n",
"-1000   100     200     255     -500    100     200     255\n",
"-500    150     225     255     0       150     225     255\n",
"0       100     150     100     30      100     150     100\n",
"30      125     175     125     60      125     175     125\n",
"60      150     200     150     122     150     200     150\n",
"122     175     225     175     183     175     225     175\n",
"183     200     255     200     244     200     255     200\n",
"244     212     255     212     305     212     255     212\n",
"305     255     255     225     457     255     255     225\n",
"457     255     225     175     610     255     225     175\n",
"610     255     225     125     702     255     225     125\n",
"702     255     175     75      914     255     175     75\n",
"914     200     150     50      1219    200     150     50\n",
"1219    175     125     50      1450    175     125     50\n",
"1450    150     100     50      1700    150     100     50\n",
"1700    150     125     100     1981    150     125     100\n",
"1981    125     125     125     2134    125     125     125\n",
"2134    150     150     150     2438    150     150     150\n",
"2438    175     175     175     2743    175     175     175\n",
"2743    200     200     200     3048    200     200     200\n",
"3048    233     233     233     9000    233     233     233\n",
"B       100     200     255\n",
"F       100     200     255\n",
"EOF\n",
"\n",
"#####\n",
"#    set GMT defaults\n",
"#####\n",
"# Output format in longitude [-180, 180] \n",
"# with two decimal points. For more precision used ddd.xxx,\n",
"# for less ddd\n",
"# If you do not get a good conversion to PNG, the EPS BoundingBox is wrong.\n",
"# Is is set for a4 paper. change the PS_MEDIA or PS_SCALE_X PS_SCALE_Y\n",
"gmt gmtset MAP_FRAME_TYPE FANCY FORMAT_GEO_MAP ddd.xx  PROJ_LENGTH_UNIT c PS_MEDIA a4 PS_SCALE_X 0.8 PS_SCALE_Y 0.8 \n",
"\n",
"#####\n",
"#    Define the default name of the PostScript output\n",
"#####\n",
"FNAME=\"map5.eps\"\n"
"\n",
""
};

static void map5_header(void){
	/* write the initial map header values */
	int i;
	for(i = 0 ; strlen(map5_proto[i]) > 0 ; i++ )
		fprintf(gmtscript,"%s",map5_proto[i]);
}

static void map5_preamble(float map5_minlon,float map5_maxlon,
	float map5_minlat,float map5_maxlat, char res,
	float lontic,float lattic)
{
	/* output the Shell variables that control the plot 
		Note this is a little complicated since 
		effort is made to assure that the default plot is
		nice */
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define map bounds: MINLON/MAXLON/MINLAT/MAXLAT \n");
fprintf(gmtscript,"#####\n");
if(map5_doglobal == YES){
	fprintf(gmtscript,"LATLON=\"-157/203/-80/80\"\n");
	lattic = 30.0 ; lontic = 30.0 ;
	res = 'c';
} else {
	fprintf(gmtscript,"LATLON=\"%f/%f/%f/%f\"\n",
		map5_minlon,map5_maxlon,MAX(map5_minlat,-80),MIN(map5_maxlat,80));
}
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define raster for topography:  \n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"GRDRAS=\"1\"\n");
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define Mercalli projection: Center_lon/Center_Lat/Plot_Width   \n");
fprintf(gmtscript,"#####\n");
if(map5_doglobal == YES){
	fprintf(gmtscript,"PROJ=\"x0.040\"\n");
} else {
	fprintf(gmtscript,"PROJ=\"M%f/%f/15c\"\n",0.5*(map5_minlon+map5_maxlon),
		0.5*(map5_minlat+map5_maxlat));
}
fprintf(gmtscript,"rm -f ${FNAME}\n");
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define Coastline resolution: one of fhilc \n");
fprintf(gmtscript,"#    (f)ull, (h)igh, (i)ntermediate, (l)ow, and (c)rude)\n");
fprintf(gmtscript,"#  . The resolution drops off by 80%% between data sets.\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"RESCOAST=\"%c\"\n",res);
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define Ticmark interval \n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"TICS=\"a%fg0/a%fg0WSne\"\n",lontic,lattic);
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define epicenter symbol and size  \n");
fprintf(gmtscript,"#####\n");
	switch(res){
		case 'c':
			fprintf(gmtscript,"EPISYM=\"A0.20\"\n");
			break;
		case 'l':
			fprintf(gmtscript,"EPISYM=\"A0.20\"\n");
			break;
		case 'i':
			fprintf(gmtscript,"EPISYM=\"A0.50c\"\n");
			break;
		case 'h':
			fprintf(gmtscript,"EPISYM=\"A0.50c\"\n");
			break;
		case 'f':
			fprintf(gmtscript,"EPISYM=\"A0.50c\"\n");
			break;
		default:
			fprintf(gmtscript,"EPISYM=\"A0.50c\"\n");
	}
fprintf(gmtscript,"EPICOLOR=\"255/255/0\"\n");

fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define station symbol and size  \n");
fprintf(gmtscript,"#####\n");
	switch(res){
		case 'c':
			fprintf(gmtscript,"STASYM=\"c0.15\"\n");
			break;
		case 'l':
			fprintf(gmtscript,"STASYM=\"c0.15\"\n");
			break;
		case 'i':
			fprintf(gmtscript,"STASYM=\"c0.37c\"\n");
			break;
		case 'h':
			fprintf(gmtscript,"STASYM=\"c0.50c\"\n");
			break;
		case 'f':
			fprintf(gmtscript,"STASYM=\"c0.50c\"\n");
			break;
		default:
			fprintf(gmtscript,"STASYM=\"c0.50c\"\n");
	}
fprintf(gmtscript,"STACOLOR=\"255/0/0\"\n");
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define boundaries for pscoast \n");
fprintf(gmtscript,"#    1 = National boundaries \n");
fprintf(gmtscript,"#    2 = State boundaries within the Americas \n");
fprintf(gmtscript,"#    3 = Marine boundaries \n");
fprintf(gmtscript,"#    a = All boundaries (1-3) pscoast \n");
fprintf(gmtscript,"#####\n");
if(map5_doglobal == YES){
	fprintf(gmtscript,"BDRYS=\"1\"\n");
} else {
	fprintf(gmtscript,"BDRYS=\"a\"\n");
}
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    Define resolution for output of grdsample \n");
fprintf(gmtscript,"#####\n");
/* do a quick here using the resolutions inherent in the RES character, e.g.,
	the logic here is that we want something to look good on a screen
	which is nominally 1000 pixels, so */
	switch(res){
		case 'c':
			fprintf(gmtscript,"RESOUT=\"20m\"\n");
			break;
		case 'l':
			fprintf(gmtscript,"RESOUT=\"10m\"\n");
			break;
		case 'i':
			fprintf(gmtscript,"RESOUT=\"2m\"\n");
			break;
		case 'h':
			fprintf(gmtscript,"RESOUT=\"0.5m\"\n");
			break;
		case 'f':
			fprintf(gmtscript,"RESOUT=\"0.25m\"\n");
			break;
		default:
			fprintf(gmtscript,"RESOUT=\"5m\"\n");
	}
fprintf(gmtscript,"\n");
}

static void map5_epicenter(int do_continue)
{
	int k, ntrchdr;
	float evla, evlo;
	ntrchdr = gsac_control.number_iheaders ;
	fprintf(gmtscript,"\n");
	fprintf(gmtscript,"#####\n");
	fprintf(gmtscript,"#    PLOT EPICENTER LOCATIONS\n");
	fprintf(gmtscript,"#####\n");
	if(do_continue){
		fprintf(gmtscript,"gmt psxy -P -J${PROJ} -R${LATLON} -O -: -S${EPISYM} -W0.8 -G${EPICOLOR} -V -K   >> ${FNAME} << EOF\n");
	} else {
		fprintf(gmtscript,"gmt psxy -P -J${PROJ} -R${LATLON} -O -: -S${EPISYM} -W0.8 -G${EPICOLOR} -V      >> ${FNAME} << EOF\n");
	}

	for ( k=0 ; k < ntrchdr ; k ++){
		evla = sacdata[k].sachdr.rhdr[H_EVLA];
		evlo = sacdata[k].sachdr.rhdr[H_EVLO];
		if(evla >= map5_minlat && evla <= map5_maxlat
			&& evlo >= map5_minlon && evlo <= map5_maxlon){
			fprintf(gmtscript,"%f %f\n",evla,evlo);
		}
	}
	fprintf(gmtscript,"EOF\n");
	fprintf(gmtscript,"\n");
}

static void map5_station(int do_continue)
{
	int k, ntrchdr;
	float stla, stlo;
	ntrchdr = gsac_control.number_iheaders ;
	fprintf(gmtscript,"\n");
	fprintf(gmtscript,"#####\n");
	fprintf(gmtscript,"#    PLOT STATION LOCATIONS\n");
	fprintf(gmtscript,"#####\n");
	if(do_continue){
		fprintf(gmtscript,"gmt psxy -P -J${PROJ} -R${LATLON} -O -: -S${STASYM} -W0.8 -G${STACOLOR} -V -K >> ${FNAME} << EOF\n");
	} else { 
		fprintf(gmtscript,"gmt psxy -P -J${PROJ} -R${LATLON} -O -: -S${STASYM} -W0.8 -G${STACOLOR} -V   >> ${FNAME} << EOF\n"); 
	}
	for ( k=0 ; k < ntrchdr ; k ++){
		stla = sacdata[k].sachdr.rhdr[H_STLA];
		stlo = sacdata[k].sachdr.rhdr[H_STLO];
		if(stla >= map5_minlat && stla <= map5_maxlat
			&& stlo >= map5_minlon && stlo <= map5_maxlon){
			fprintf(gmtscript,"%f %f\n",stla,stlo);
		}
	}
	fprintf(gmtscript,"EOF\n");
	fprintf(gmtscript,"\n");
}

static void map5_coast(int do_overlay, int do_continue)
{
	/* note: do_overlay = YES means that the map5.eps has already
		been created with a raster */
	fprintf(gmtscript,"gmt pscoast -P -W  -J${PROJ} -R${LATLON} -B${TICS}");
	if(do_overlay == YES)
		fprintf(gmtscript," -O ");
	else
		fprintf(gmtscript," -G240 ");
	if(do_continue == YES)
		fprintf(gmtscript," -K ");
	fprintf(gmtscript," -N${BDRYS} -D${RESCOAST} -A2500  -V ");
	if(do_overlay == YES)
		fprintf(gmtscript," >> ${FNAME}\n");
	else
		fprintf(gmtscript,"  > ${FNAME}\n");
}

static void map5_end(void)
{
	fprintf(gmtscript,"\n");
	fprintf(gmtscript,"######\n");
	fprintf(gmtscript,"#     Cleanup\n");
	fprintf(gmtscript,"######\n");
	fprintf(gmtscript,"rm -f map.grd cmap.grd MA.grd map.cpt\n");
}

static void map5_topo(void)
{
fprintf(gmtscript,"gmt grdraster ${GRDRAS}  -R${LATLON} -Gcmap.grd -V\n");
fprintf(gmtscript,"gmt grdsample cmap.grd -Gmap.grd -I${RESOUT} -R${LATLON} -V\n");
fprintf(gmtscript,"gmt grdgradient map.grd -A135 -GMA.grd -Nt -V \n");
fprintf(gmtscript,"gmt grdimage -P map.grd -X2.5i -Y1.5i -J${PROJ} -R${LATLON} -Cmap.cpt -K -IMA.grd  -V > ${FNAME} \n");
}

static void map5_raypath(void)
{
	/* trace the raypath between the epicenter and station
		note this only works for sac files having
		stla, stlo, evla, evlo set */
	int k, ntrchdr;
	float stla, stlo, evla, evlo ;
	ntrchdr = gsac_control.number_iheaders ;
fprintf(gmtscript,"\n");
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#    MAP GREAT CIRCLE RAY PATHS USING PROJECT\n");
fprintf(gmtscript,"#####\n");
		for ( k=0 ; k < ntrchdr ; k ++){
			stla = sacdata[k].sachdr.rhdr[H_STLA];
			stlo = sacdata[k].sachdr.rhdr[H_STLO];
			evla = sacdata[k].sachdr.rhdr[H_EVLA];
			evlo = sacdata[k].sachdr.rhdr[H_EVLO];
			if(stla != -12345. && stlo  != -12345.
				&& evla != -12345. && evlo != -12345. ){
	
fprintf(gmtscript,"gmt psxy -P -J${PROJ} -R${LATLON} -W2.0 -O -: -V -K >> ${FNAME} << EOF\n");
fprintf(gmtscript,"%f %f \n", evla,evlo);
fprintf(gmtscript,"%f %f \n", stla,stlo);
fprintf(gmtscript,"EOF\n");
		}
	}
}

static void map5_kstnm(void)
{
	int k, ntrchdr;
	float stla, stlo;
	ntrchdr = gsac_control.number_iheaders ;
fprintf(gmtscript,"#####\n");
fprintf(gmtscript,"#  This will draw a box in black -Wblack\n");
fprintf(gmtscript,"#  The box will be filled with white -G255/255/25\n");
fprintf(gmtscript,"#  the text is 10 point, Helvetica, black, aligned Left Bottom and angle 0\n");
fprintf(gmtscript,"#  To just draw the text remove -G255/255/255 -W0/0/0 from the command line\n");
fprintf(gmtscript,"#####\n");
	fprintf(gmtscript,"gmt pstext -P -J${PROJ} -R${LATLON} -O -:  -V -Dj0.1i/0.1i -G255/255/255 -F+f10p,Helvetica,black+jLB+a0 -Wblack  >> ${FNAME} << EOF\n");

	for ( k=0 ; k < ntrchdr ; k ++){
		stla = sacdata[k].sachdr.rhdr[H_STLA];
		stlo = sacdata[k].sachdr.rhdr[H_STLO];
		if(stla >= map5_minlat && stla <= map5_maxlat
			&& stlo >= map5_minlon && stlo <= map5_maxlon){
			fprintf(gmtscript,"%f %f %s\n",stla,stlo,sacdata[k].schdr[H_KSTNM]);
		}
	}
	fprintf(gmtscript,"EOF\n");
	fprintf(gmtscript,"\n");
}
