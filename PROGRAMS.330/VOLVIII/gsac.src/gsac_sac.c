/* This requires more work so that
 * 1. We test to see if it is a valid SAC file
 * 2. We test to see if it is a BYTE swapped valid SAC file
 */
/* CHANGES
 09 JAN 2005 - added  nerr = -1; in bwsach
 03 MAR 2009 - added logic in brsac and brsach to avoid reading non-sac files
		note this logic is safely duplicated since I usually do a
		gsac_valid_sacfile which calls brsach before one is 
		actually does a brsac 
		Eventually get the file size from the system and check
		file size versus what is estimated from header
 10 SEP 2009 - corrected valid sac file test to include a spectra file, e.g., IXY
		Change
		if(has12345 == 0 || sachdr.ihdr[H_NVHDR] != 6 || sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME){
to
		if(has12345 == 0 || sachdr.ihdr[H_NVHDR] != 6 || (sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME && sachdr.ihdr[H_IFTYPE] !=ENUM_IXY)){
 05 AUG 2022    Implement NVHDR=7 for extended precision for some variables
		NHVDR=7 adds 22 doubble precision floats after the data section. The variables are
		DELTA, B, E, O, A, T0 ... T9, F, EVLO, EVLA, STLO, STLA, SB, and SDELTA
		0      1  2  3  4  5      14 15  16    17    18    19    20       21
		
*/

#include	<stdio.h>
#include	 "gsac.h"
#include	 "gsac_sac.h"
#include	"csstim.h"

static float *x = (float *)NULL;

/* these are absolutely necessary because the addition to cut in
	absoltue time */
#define CUT_GMT		-20
#define CUT_CAL		-21
#define CUT_VEL		-22
#define CUT_P		-23

/* get extended header values for NVHDR = 7 note since the struct sachdr_rhdr is 
	double these values will jsut replace the float values on input
*/

struct dsachdr_ {
	double dhdr[22];
};

/* this is the original sac header For internal computations the float rhdr is now 
   double rhdr */


int brsac(char *fname,struct  sachdr_ *sachdrret, float **data)
{
	/* update 20 OCT 2022  - to implement NVHDR=7 correctly
		read in everything, and then implement the cut
	*/
	int maxpts, nread;
	int newnpts;		/* for cut new number of points */
	double newb, newe;	/* new values for b and e if cut */
	int newbeg, newend;		/* index of first point in new series */
	double oldb, olde;
	double delta;
	double t0, t1, tv;
	double tzref;
	int i;
	int i1, i2;
	int j1, j2, j3;
	int nerr;
	int has12345;
	int rhas12345, ihas12345, chas12345;
        int nvhdr;
	FILE *fptr;
	struct orig_sachdr_ sachdr;
	struct sachdr_ newsachdr;
	struct dsachdr_ dsachdr;
	/*
		determine if the file exists for reading
	*/
	if((fptr=fopen(fname,"rb")) == NULL){
		sachdr.ihdr[H_NPTS] = 0;
		nerr = -1;
		perror("fopen error in brsac:");
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	/*
		file exists get header information only
	*/
	
	nerr = 0;
	if(fread(&sachdr,sizeof( struct orig_sachdr_ ),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	maxpts = sachdr.ihdr[H_NPTS] ;
			nerr = 0 ;

	/* there must be at least one -12345 for this to be a sac file!! */
	has12345 = 0;
	rhas12345 = 0;
	ihas12345 = 0;
	chas12345 = 0;
	for (i=0;i < 70 ; i++){
		if(sachdr.rhdr[i] == -12345.0) {
			has12345++ ; 
			rhas12345++;
		}
	}
	for (i=0;i < 40 ; i++){
		if(sachdr.ihdr[i] == -12345){
			has12345++; 
			ihas12345++;
		}
	}
	for (i=0;i < 24 ; i++){
		if(strncmp(sachdr.chdr[i],  "-12345  ",8)==0){
			has12345++;
			chas12345++;
		}
	}
	if(has12345 == 0 || (sachdr.ihdr[H_NVHDR] != 6 && sachdr.ihdr[H_NVHDR] != 7) || 
		(sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME && sachdr.ihdr[H_IFTYPE] !=ENUM_IXY)){
		if(has12345 == 0){
			fprintf(stderr,"Problem: %s is not a Sac file \n",fname);
		} else {
			fprintf(stderr,"Problem: %s may be a byte-swapped Sac file. Separately run 'saccvt -I'\n",fname);
			nerr = -3;
			fclose(fptr);
			return(nerr);
		}
	}
	/* since the rhdr returned is double and the sachdr_rhdr here is
		float, we must carefully copy the values over */
	for(i=0;i< 70 ; i++)
		newsachdr.rhdr[i] = (double)sachdr.rhdr[i];
	for(i=0 ; i < 40 ; i++)
		newsachdr.ihdr[i] = sachdr.ihdr[i];
	for(i=0; i < 24 ; i++)
		strncpy(newsachdr.chdr[i],sachdr.chdr[i],8);

	/* Now get the wave form and allocate the temporary array to be read in */
	if(x == (float *)NULL)
		x = (float *)calloc(maxpts,sizeof(float));
	else
		x = (float *)realloc(x,maxpts*sizeof(float));
	nread = fread(x,sizeof(float),maxpts,fptr);


	/* at this point if NVHDR = 7 read in extended header values */
	if(newsachdr.ihdr[H_NVHDR] == 7 ){
		if(fread(&dsachdr,sizeof(struct dsachdr_),1,fptr)!=1){
			nerr = -3 ;
		} else {
			newsachdr.rhdr[H_DELTA]  = dsachdr.dhdr[0];
			newsachdr.rhdr[H_B]      = dsachdr.dhdr[1];
			newsachdr.rhdr[H_E]      = dsachdr.dhdr[2];
			newsachdr.rhdr[H_O]      = dsachdr.dhdr[3];
			newsachdr.rhdr[H_A]      = dsachdr.dhdr[4];
			newsachdr.rhdr[H_T0]     = dsachdr.dhdr[5];
			newsachdr.rhdr[H_T1]     = dsachdr.dhdr[6];
			newsachdr.rhdr[H_T2]     = dsachdr.dhdr[7];
			newsachdr.rhdr[H_T3]     = dsachdr.dhdr[8];
			newsachdr.rhdr[H_T4]     = dsachdr.dhdr[9];
			newsachdr.rhdr[H_T5]     = dsachdr.dhdr[10];
			newsachdr.rhdr[H_T6]     = dsachdr.dhdr[11];
			newsachdr.rhdr[H_T7]     = dsachdr.dhdr[12];
			newsachdr.rhdr[H_T8]     = dsachdr.dhdr[13];
			newsachdr.rhdr[H_T9]     = dsachdr.dhdr[14];
			newsachdr.rhdr[H_F]      = dsachdr.dhdr[15];
			newsachdr.rhdr[H_EVLO]   = dsachdr.dhdr[16];
			newsachdr.rhdr[H_EVLA]   = dsachdr.dhdr[17];
			newsachdr.rhdr[H_STLO]   = dsachdr.dhdr[18];
			newsachdr.rhdr[H_STLA]   = dsachdr.dhdr[19];
			newsachdr.rhdr[H_SB]     = dsachdr.dhdr[20];
			newsachdr.rhdr[H_SDELTA] = dsachdr.dhdr[21];
		} 
	}
	/* Implement the cut on read if required */
	/* get tzref - needed to cut CAL or cut GMT */
	
	if(newsachdr.ihdr[H_NZYEAR] == -12345 ||
		newsachdr.ihdr[H_NZJDAY] == -12345 ||
		newsachdr.ihdr[H_NZHOUR] == -12345 ||
		newsachdr.ihdr[H_NZMIN] == -12345 ||
		newsachdr.ihdr[H_NZSEC] == -12345 ||
		newsachdr.ihdr[H_NZMSEC] == -12345){
		/*
		printf("Using 1970 000 00 00 00 000\n");
		*/
		tzref = 0.0;
	} else {
		htoe1(newsachdr.ihdr[H_NZYEAR], 
			newsachdr.ihdr[H_NZJDAY], 
			newsachdr.ihdr[H_NZHOUR],
			newsachdr.ihdr[H_NZMIN],
			newsachdr.ihdr[H_NZSEC],
			newsachdr.ihdr[H_NZMSEC],
			&tzref);
	}

	if(gsac_control.docut && newsachdr.rhdr[gsac_control.cutint[0]]!= -12345.
			&&newsachdr.rhdr[gsac_control.cutint[1]]!= -12345.){

		oldb = newsachdr.rhdr[H_B];
		olde = newsachdr.rhdr[H_E];
		for(i=0;i<2;i++){
			if(gsac_control.cutint[i] == CUT_GMT){
				tv = (float)(gsac_control.cutepoch[i]-tzref);
			} else if(gsac_control.cutint[i] == CUT_CAL){
				tv = (float)(gsac_control.cutepoch[i]-tzref);
			} else if(gsac_control.cutint[i] == CUT_VEL || gsac_control.cutint[i] == CUT_P){
				if ( newsachdr.rhdr[H_O] != -12345. && newsachdr.rhdr[H_DIST] != -12345. ){
					tv = newsachdr.rhdr[H_O] + newsachdr.rhdr[H_DIST]*gsac_control.cutslow[i] 
						+ gsac_control.cutoff[i] - oldb;
				} else {
					tv = 0 ;
				}
			} else {
				tv = newsachdr.rhdr[gsac_control.cutint[i]]
					+gsac_control.cutoff[i] - oldb;
			}
			if(i==0)
				t0 = tv;
			else
				t1 = tv;
		}
		delta = newsachdr.rhdr[H_DELTA];

		newb = MIN(t0,t1);
		newe = MAX(t0,t1);
		newnpts = 1 + (int)((newe - newb)/delta + 0.49);
		if(newb > 0)
			newbeg =  (int)((newb)/delta +0.49);
		else
			newbeg =  (int)((newb)/delta -0.49);
		newb = oldb + newbeg*delta;
		newe = newb + (newnpts -1)*delta;
		newend = newbeg + newnpts -1 ;
	} else {
		oldb = newsachdr.rhdr[H_B];
		olde = newsachdr.rhdr[H_E];
		newb = oldb;
		newe = olde;
		newbeg = 0;
		newnpts = newsachdr.ihdr[H_NPTS] ;
		newend = newbeg + newnpts  ;
	}	
	/* now update the header values */
	newsachdr.rhdr[H_B] = newb;
	newsachdr.rhdr[H_E] = newe;
	newsachdr.ihdr[H_NPTS] = newnpts;
	/* now read into the data array with card, but we need an
	 * offset in one case to skip a few */
	

	/* there are several possibilities here 
	 * if the data are within the original window we cant to go from
	 *   OFF1 to OFF2 where 0 <=OFF1 and OFF2 <= maxpts
	 * if we read off the end then we find that we can only read to 
	 * 	maxpts
	 * if OFF1 < 0 then we fill the data array starting at maxbeg, or
	 * like
	 * 	fread(&(*data)[maxbeg-,sizeof(float),MIN(maxpts,newnpts)
	 * 	et cetera et cetera et cetera */



	i1 = newbeg;
	i2 = newend;
	j1 = MAX(0,-i1);
	j2 = MAX(0, i1);
	j3 = MIN(i2-j2+1,maxpts);
	/* allocate the data values setting them to zero */
	if(*data == NULL){
		if((*data = (float *)calloc(newnpts,sizeof(float))) == NULL){
			nerr = -2;
			fclose(fptr);
			return(nerr);
		}
	} else {
		/* set to zero */
		if((*data = (float *)realloc(*data,newnpts*sizeof(float))) == NULL){
			nerr = -2;
			fclose(fptr);
			return(nerr);
		}
	}
	/* initialize to zero since realloc does not do this */
	for(i=0;i<newnpts;i++)
		(*data)[i] = 0.0;
	/* at this point we have read the file
		next modify the time series and header values is a cut
		command was previously invoked */
	/* ignore error return for fread */
	i = 0;
	while(i < j3){
		if((j2+i) >= 0 && (j2+i)< maxpts){
			/* KLUDGE TO MAKE THIS WORK FOR e -1 e 0.05 */
		(*data)[j1 + i] = x[j2 + i ];

		}
		i++;
	}
	*sachdrret = newsachdr;
	fclose (fptr) ;
	return(nerr);
}

int brsach(char *fname,struct  sachdr_ *sachdrret)
{
	int i, maxpts;
	int nerr;
	int has12345;
	int rhas12345, ihas12345, chas12345;
	FILE *fptr;
	struct orig_sachdr_ sachdr;
	struct sachdr_ newsachdr;
	struct dsachdr_ dsachdr;
	/*
		determine if the file exists for reading
	*/
	if((fptr=fopen(fname,"rb")) == NULL){
		sachdr.ihdr[H_NPTS] = 0;
		nerr = -1;
		perror("fopen error in brsac:");
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	/*
		file exists get header information only
	*/
	
	nerr = 0;
	if(fread(&sachdr,sizeof( struct orig_sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
			nerr = 0 ;

	/* there must be at least one -12345 for this to be a sac file!! */
	has12345 = 0;
	rhas12345 = 0;
	ihas12345 = 0;
	chas12345 = 0;
	for (i=0;i < 70 ; i++){
		if(sachdr.rhdr[i] == -12345.0) {
			has12345++ ; 
			rhas12345++;
		}
	}
	for (i=0;i < 40 ; i++){
		if(sachdr.ihdr[i] == -12345){
			has12345++; 
			ihas12345++;
		}
	}
	for (i=0;i < 24 ; i++){
		if(strncmp(sachdr.chdr[i],  "-12345  ",8)==0){
			has12345++;
			chas12345++;
		}
	}
	/* File 
		must have at least one -12345 sequence, 
		must have the correct * header version, and 
		must be equally spaced time series 
	*/
	if(has12345 == 0 || (sachdr.ihdr[H_NVHDR] == 6 || sachdr.ihdr[H_NVHDR] == 7) == 0 ||
		(sachdr.ihdr[H_IFTYPE] !=ENUM_ITIME && sachdr.ihdr[H_IFTYPE] !=ENUM_IXY)){
		if(has12345 == 0){
			fprintf(stderr,"Problem: %s is not a Sac file \n",fname);
		} else {
			fprintf(stderr,"Problem: %s may be a byte-swapped Sac file. Separately run 'saccvt -I'\n",fname);
		}
		nerr = -3;
		fclose(fptr);
		return(nerr);
	}
	/* since the rhdr returned is double and the sachdr_rhdr here is
		float, we must carefully copy the values over */
	for(i=0;i< 70 ; i++)
		newsachdr.rhdr[i] = (double)sachdr.rhdr[i];
	for(i=0 ; i < 40 ; i++)
		newsachdr.ihdr[i] = sachdr.ihdr[i];
	for(i=0; i < 24 ; i++)
		strncpy(newsachdr.chdr[i],sachdr.chdr[i],8);
	/* at this point if NVHDR = 7 read in extended header values */
	if(sachdr.ihdr[H_NVHDR] == 7 ){
		maxpts = sachdr.ihdr[H_NPTS];
		fseek(fptr, (long)maxpts*4, SEEK_CUR);
		if(fread(&dsachdr,sizeof(struct dsachdr_),1,fptr)!=1){
			nerr = -3 ;
		} else {
			newsachdr.rhdr[H_DELTA]  = dsachdr.dhdr[0];
			newsachdr.rhdr[H_B]      = dsachdr.dhdr[1];
			newsachdr.rhdr[H_E]      = dsachdr.dhdr[2];
			newsachdr.rhdr[H_O]      = dsachdr.dhdr[3];
			newsachdr.rhdr[H_A]      = dsachdr.dhdr[4];
			newsachdr.rhdr[H_T0]     = dsachdr.dhdr[5];
			newsachdr.rhdr[H_T1]     = dsachdr.dhdr[6];
			newsachdr.rhdr[H_T2]     = dsachdr.dhdr[7];
			newsachdr.rhdr[H_T3]     = dsachdr.dhdr[8];
			newsachdr.rhdr[H_T4]     = dsachdr.dhdr[9];
			newsachdr.rhdr[H_T5]     = dsachdr.dhdr[10];
			newsachdr.rhdr[H_T6]     = dsachdr.dhdr[11];
			newsachdr.rhdr[H_T7]     = dsachdr.dhdr[12];
			newsachdr.rhdr[H_T8]     = dsachdr.dhdr[13];
			newsachdr.rhdr[H_T9]     = dsachdr.dhdr[14];
			newsachdr.rhdr[H_F]      = dsachdr.dhdr[15];
			newsachdr.rhdr[H_EVLO]   = dsachdr.dhdr[16];
			newsachdr.rhdr[H_EVLA]   = dsachdr.dhdr[17];
			newsachdr.rhdr[H_STLO]   = dsachdr.dhdr[18];
			newsachdr.rhdr[H_STLA]   = dsachdr.dhdr[19];
			newsachdr.rhdr[H_SB]     = dsachdr.dhdr[20];
			newsachdr.rhdr[H_SDELTA] = dsachdr.dhdr[21];
		}
	}
	*sachdrret = newsachdr;
	
	fclose (fptr) ;
	*sachdrret = newsachdr;
	return(nerr);
}

int bwsac(char *fname,struct  sachdr_ sachdr, float *data){
	FILE *fptr;
	int i, maxpts, nwrite, nerr;
	struct orig_sachdr_ outsachdr;
	struct dsachdr_ dsachdr;

	if((fptr=fopen(fname,"wb")) == NULL){
		nerr = -1;
		perror("fopen error in bwsac:");
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	fseek(fptr, 0L, SEEK_SET);
	for(i=0;i< 70 ; i++)
		outsachdr.rhdr[i] = (float)sachdr.rhdr[i];
	for(i=0 ; i < 40 ; i++)
		outsachdr.ihdr[i] = sachdr.ihdr[i];
	for(i=0; i < 24 ; i++)
		strncpy(outsachdr.chdr[i],sachdr.chdr[i],8);
	if(fwrite(&outsachdr,sizeof( struct orig_sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	maxpts = sachdr.ihdr[H_NPTS] ;
	nwrite = fwrite(data,sizeof(float),maxpts,fptr);
	/* ignore error return for fwrite */
	if(sachdr.ihdr[H_NVHDR] == 7 ){
		/* copy internal variables to the double header */
		dsachdr.dhdr[0] = sachdr.rhdr[H_DELTA]   ;
		dsachdr.dhdr[1] = sachdr.rhdr[H_B]       ;
		dsachdr.dhdr[2] = sachdr.rhdr[H_E]       ;
		dsachdr.dhdr[3] = sachdr.rhdr[H_O]       ;
		dsachdr.dhdr[4] = sachdr.rhdr[H_A]       ;
		dsachdr.dhdr[5] = sachdr.rhdr[H_T0]      ;
		dsachdr.dhdr[6] = sachdr.rhdr[H_T1]      ;
		dsachdr.dhdr[7] = sachdr.rhdr[H_T2]      ;
		dsachdr.dhdr[8] = sachdr.rhdr[H_T3]      ;
		dsachdr.dhdr[9] = sachdr.rhdr[H_T4]      ;
		dsachdr.dhdr[10] = sachdr.rhdr[H_T5]      ;
		dsachdr.dhdr[11] = sachdr.rhdr[H_T6]      ;
		dsachdr.dhdr[12] = sachdr.rhdr[H_T7]      ;
		dsachdr.dhdr[13] = sachdr.rhdr[H_T8]      ;
		dsachdr.dhdr[14] = sachdr.rhdr[H_T9]      ;
		dsachdr.dhdr[15] = sachdr.rhdr[H_F]       ;
		dsachdr.dhdr[16] = sachdr.rhdr[H_EVLO]    ;
		dsachdr.dhdr[17] = sachdr.rhdr[H_EVLA]    ;
		dsachdr.dhdr[18] = sachdr.rhdr[H_STLO]    ;
		dsachdr.dhdr[19] = sachdr.rhdr[H_STLA]    ;
		dsachdr.dhdr[20] = sachdr.rhdr[H_SDELTA]  ;
		dsachdr.dhdr[21] = sachdr.rhdr[H_SB]      ;
		nwrite = fwrite(&dsachdr,sizeof(struct dsachdr_),1,fptr);
	}
	fclose(fptr);
	return(0);
}

int bwsach(char *fname,struct  sachdr_ sachdr){
	int nerr;
	FILE *fptr;
	struct orig_sachdr_ outsachdr;
	struct dsachdr_ dsachdr;
        int i, maxpts, nwrite;
	if((fptr=fopen(fname,"r+b")) == NULL){
		perror("fopen error in bwsach:");
		nerr = -1;
		return(nerr);
	}
#ifdef WIN32
	setmode(fileno(fptr), O_BINARY);
#endif
	fseek(fptr, 0L, SEEK_SET);
	for(i=0;i< 70 ; i++)
		outsachdr.rhdr[i] = (float)sachdr.rhdr[i];
	for(i=0 ; i < 40 ; i++)
		outsachdr.ihdr[i] = sachdr.ihdr[i];
	for(i=0; i < 24 ; i++)
		strncpy(outsachdr.chdr[i],sachdr.chdr[i],8);
	if(fwrite(&outsachdr,sizeof( struct orig_sachdr_),1,fptr)!=1){
		nerr=-3;
		fclose(fptr);
		return(nerr);
	}
	maxpts = sachdr.ihdr[H_NPTS] ;
	/* ignore error return for fwrite */
	if(sachdr.ihdr[H_NVHDR] == 7 ){
		fseek(fptr, (long)maxpts*4, SEEK_CUR);
		/* copy internal variables to the double header */
		dsachdr.dhdr[0] = sachdr.rhdr[H_DELTA]   ;
		dsachdr.dhdr[1] = sachdr.rhdr[H_B]       ;
		dsachdr.dhdr[2] = sachdr.rhdr[H_E]       ;
		dsachdr.dhdr[3] = sachdr.rhdr[H_O]       ;
		dsachdr.dhdr[4] = sachdr.rhdr[H_A]       ;
		dsachdr.dhdr[5] = sachdr.rhdr[H_T0]      ;
		dsachdr.dhdr[6] = sachdr.rhdr[H_T1]      ;
		dsachdr.dhdr[7] = sachdr.rhdr[H_T2]      ;
		dsachdr.dhdr[8] = sachdr.rhdr[H_T3]      ;
		dsachdr.dhdr[9] = sachdr.rhdr[H_T4]      ;
		dsachdr.dhdr[10] = sachdr.rhdr[H_T5]      ;
		dsachdr.dhdr[11] = sachdr.rhdr[H_T6]      ;
		dsachdr.dhdr[12] = sachdr.rhdr[H_T7]      ;
		dsachdr.dhdr[13] = sachdr.rhdr[H_T8]      ;
		dsachdr.dhdr[14] = sachdr.rhdr[H_T9]      ;
		dsachdr.dhdr[15] = sachdr.rhdr[H_F]       ;
		dsachdr.dhdr[16] = sachdr.rhdr[H_EVLO]    ;
		dsachdr.dhdr[17] = sachdr.rhdr[H_EVLA]    ;
		dsachdr.dhdr[18] = sachdr.rhdr[H_STLO]    ;
		dsachdr.dhdr[19] = sachdr.rhdr[H_STLA]    ;
		dsachdr.dhdr[20] = sachdr.rhdr[H_SDELTA]  ;
		dsachdr.dhdr[21] = sachdr.rhdr[H_SB]      ;
		nwrite = fwrite(&dsachdr,sizeof(struct dsachdr_),1,fptr);
	}
	fclose(fptr);
	return(0);
}


int gsac_valid_sacfile(char *name)
{
        struct sacfile_ tsacdata ;
        int iret;
        /* examine the trace header to determine if the
         * file is a valid SAC file
         * return YES or NO
         * The check is performed by opening the file, reading the
         * 3 headers, and determining the following:
         * Is the Version correct, is there one -12345. or -12345
         * in the real/integer header. Also determine if the
         * byte order is reversed
         * */
        strcpy(tsacdata.sac_ifile_name , name);
        strcpy(tsacdata.sac_ofile_name , name);
        iret = brsach(tsacdata.sac_ifile_name,&tsacdata.sachdr);
        if(iret < 0)
                return NO;
        else
                return YES;
}

struct sacfile_ *sacdata ;
int *sortptr;
float *sortfloat;


void gsac_alloc_trace(int oldmax ){
/* allocate a data structure */
/* we need counter here for the number of current successes */
int k;
	if(oldmax == 0 && gsac_control.max_number_traces == 0 ){
		sacdata = (struct sacfile_ *)calloc(1, sizeof(struct sacfile_));
		sortptr = (int *)calloc(1, sizeof(int));
		sortfloat = (float *)calloc(1, sizeof(float));
		sortptr[0] = 0;
		sacdata[0].sac_data = (float *)NULL;
		sacdata[0].sac_spectra = (float *)NULL;
		gsac_control.max_number_traces++;
		sortptr[0] = 0;
	} else {
		if(gsac_control.max_number_traces >(oldmax-1))
		sacdata = realloc(sacdata, sizeof(struct sacfile_)*(gsac_control.max_number_traces +1));
		sortptr = realloc(sortptr, sizeof(int)*(gsac_control.max_number_traces +1));
		sortfloat = realloc(sortfloat, sizeof(float)*(gsac_control.max_number_traces +1));
		sacdata[gsac_control.max_number_traces].sac_data = (float *)NULL;
		sacdata[gsac_control.max_number_traces].sac_spectra = (float *)NULL;
		gsac_control.max_number_traces++;
		for(k=0;k<gsac_control.max_number_traces;k++){
			sortptr[k] = k;
		}
	}
}

