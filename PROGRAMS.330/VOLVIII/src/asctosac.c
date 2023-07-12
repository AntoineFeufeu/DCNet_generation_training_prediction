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

extern struct sachdr_ sachdr;
extern struct dsachdr_ dsachdr;

/* prototypes */
void usage(void);
void gcmdln(int argc, char **argv, char **fin, char **fout);

void main(int argc, char **argv)
{
	char *fin, *fout;
	float *x;
	int nerr;
	int npts;
	gcmdln(argc, argv,  &fin,  &fout);
	npts = sachdr.ihdr[9];
	arsac(LN,fin,&x,&nerr);
	if(nerr != -1)
		bwsac (npts,fout,x);
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
        fprintf(stderr,"Usage: asctosac  SAC_ASCII_FILE SAC_BINARY_FILE\n");
        fprintf(stderr," Convert SAC ASCII FILE TO SAC BINARY\n");
	exit(0);
}
