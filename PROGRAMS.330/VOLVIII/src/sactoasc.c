/* converted from asctosac.f */
/*
	Changes
	30 OCT 2022 - changed to support extended sac files with NVHDR =7 in addition ot NVHDR = 6
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <libgen.h>
#include "sacsubc.h"
#include "calplot.h"

#define LN 1000


/* prototypes */
void usage(void);
void gcmdln(int argc, char **argv, char **fin, char **fout);

extern struct sachdr_ sachdr;
extern struct dsachdr_ dsachdr;

void main(int argc, char **argv)
{
	char *fin, *fout;
	float *x;
	int nerr;
	int npts;
	gcmdln(argc, argv,  &fin,  &fout);
	brsac(LN,fin,&x,&nerr);
	npts = sachdr.ihdr[9];
	if(nerr != -1)
		awsac (npts,fout,x);
	free(x);

}

void gcmdln(int argc, char **argv, char **fin, char **fout)
{
	if(argc !=3)
		usage();
	*fin = argv[1];
	*fout = argv[2];
}

void usage(void)
{
        fprintf(stderr,"Usage: sactoasc  SAC_BINARY_FILE SAC_ASCII_FILE \n");
        fprintf(stderr," Convert SAC BINARY FILE TO SAC ASCII\n");
	exit(0);
}
