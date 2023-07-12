/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: TRTTEK                                                c
c                                                                     c
c      COPYRIGHT (C)  1996 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
	TeraTerm Return Tek
	06 MAY 2014 - return for main and incldue stdio.h (Larry Baker, USGS Menlo Park)
*/

#include <stdio.h>

int main()
{
	putchar('\033');
	putchar('\014');
	putchar('\035');
	putchar('\061');
	putchar('\015');
	putchar('\033');
	putchar('\037');
	putchar('\030');
	return 0;
}
