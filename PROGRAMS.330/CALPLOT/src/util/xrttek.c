/*
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME I                                                       c
c                                                                     c
c      PROGRAM: XRTTEK                                                c
c                                                                     c
c      COPYRIGHT (C)        1991 R. B. Herrmann                       c
c                                                                     c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63103                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
CHANGES:
	06 MAY 2014 main is int with return (Larry Baker, USGS Menlo Park)
*/

#include <stdio.h>

int main()
{
	putchar('\033');
	putchar('\003');
	return 0;
}
