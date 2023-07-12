#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ctype.h>


char strlst(string)
char *string;
{
	while (*string != '\0') { if (*(string+1) == '\0') return(*string); string++;}
	return('\0');
}

/* Trim leading and trailing whitespace from string s (in place) */

void trim( char *s ) {

  char *c;


  if ( s != NULL )
  {

    /* Right trim */
    for ( c = s; *c != '\0'; ++c ) {}		/* c->EOS                    */
    for ( --c; c >= s; --c )			/* c->last non-whitespace    */
      if ( !isspace( *c ) )			/*   character in s, else    */
        break;					/*   s-1 if s is "" or all   */
        					/*   whitespace              */
    *(++c) = '\0';				/* Set EOS at *(c+1)         */

    /* Left trim */
    for ( c = s; isspace( *c ); ++c ) {}	/* c->first non-whitespace   */
    						/*   character in s          */
    if ( c != s )				/* If there were any...      */
      while ( ( *s++ = *c++ ) != '\0' ) {}	/*   ...slide s to the left  */

  }
}

/* ---------------------------------------------------------------------- */
int str_isnum(s)
char *s;

{
	while(*s)
		if (!isdigit(*s++))
			return 0;

	return 1;

}

/* ---------------------------------------------------------------------- */
void pad_it(char *s, int l)
{

        int sl;

        int ii;

	sl = strlen(s);

        if (sl < l)
                for (ii = 0; ii < l - sl; ii++)
                        strcat(s, "0");
}


/* ---------------------------------------------------------------------- */


