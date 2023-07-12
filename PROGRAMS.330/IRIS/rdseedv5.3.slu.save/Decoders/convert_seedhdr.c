/*===========================================================================*/
/* SEED reader     |             convert_seedhdr           |    process data */
/*===========================================================================*/
/*
	Name:		convert_seedhdr
	Purpose:	build a minimum data header from the SEED fixed section data hdr
				in the same form as tables for easy comparison
	Usage:		void convert_seedhdr ();
				char *input_data_ptr;
				struct input_data_hdr *input_data_hdr;
				input_data_ptr = input.data;
				input_data_hdr = (struct input_data_hdr *) input_data_ptr;
				convert_seedhdr (input_data_hdr);
	Input:		input_data_hdr = a structure containing the SEED fixed section
					of data header information
	Output:		data_hdr = a structure containing minimal header information,
					as follows:
					char *station                    = station name
					char *location                   = station location
					char *channel                    = channel name
					struct time time                 =  time
					int nsamples                     = number of samples
					double sample_rate               = sample rate per second
					char activity_flags              = activity flags
					char io_flags                    = i/o flags
					char data_quality_flags          = data quality flags
					long int number_time_corrections = # .0001s time corrections
	Externals:	input.data - pointer to beginning of input data block
				*data_hdr - defined in structures.h, allocated in globals.h
	Warnings:	none
	Errors:		none
	Called by:	process_data
	Calls to:	none
	Algorithm:	mostly direct cast-and-copy of values from the SEED structure
				"input_data_hdr" to a more easily-used form.  Sample rate is
				calculated according to the rules in the SEED description of the
				fixed section of data header.
	Notes:		none
	Problems:	none known
	References:	Halbert et al, 1988; see main routine
	Language:	Sun OS 3.2 C
	Revisions:	07/15/88  Dennis O'Neill  Preliminary initial release 0.9
				11/13/88  Dennis O'Neill  Added smart byte-swapping
				11/21/88  Dennis O'Neill  Production release 1.0
				02/10/89  Dennis O'Neill  fixed the statement
					data_hdr->number_time_corrections += (long int) swap_4byte;
					data_hdr->number_time_corrections += (long int) temp_4byte;
				05/01/92  Allen Nnace     fixed the statement
					data_hdr->number_time_corrections += (long int) temp_4byte; to
					data_hdr->number_time_corrections = (long int) temp_4byte;
*/

#include "rdseed.h"

void convert_seedhdr (input_data_hdr)
struct input_data_hdr *input_data_hdr;
{
	int i;	/* counter */
	double sample_rate;	/* sample rate calculation */
	int sample_rate_multiplier;	/* sample rate calculation */
	struct input_time temptime;	/* temp for byte-swapping */

/*                 +=======================================+                 */
/*=================|   create minimum useful data header   |=================*/
/*                 +=======================================+                 */

	/* copy, strip trailing blanks from text strings, put in structure */
	/* all of these fields are fixed-length, left justified, alpha-numeric characters only, blank padded */

	/* copy station name */
	for ( i = 0; ( i < sizeof( input_data_hdr->station ) ) && ( input_data_hdr->station[i] != ' ' ); ++i )
		data_hdr->station[i] = input_data_hdr->station[i];
	data_hdr->station[i] = '\0';

	/* copy location code */
	for ( i = 0; ( i < sizeof( input_data_hdr->location ) ) && ( input_data_hdr->location[i] != ' ' ); ++i )
		data_hdr->location[i] = input_data_hdr->location[i];
	data_hdr->location[i] = '\0';

	/* copy channel name */
	for ( i = 0; ( i < sizeof( input_data_hdr->channel ) ) && ( input_data_hdr->channel[i] != ' ' ); ++i )
		data_hdr->channel[i] = input_data_hdr->channel[i];
	data_hdr->channel[i] = '\0';

	/* copy network name */
	for ( i = 0; ( i < sizeof( input_data_hdr->network ) ) && ( input_data_hdr->network[i] != ' ' ); ++i )
		data_hdr->network[i] = input_data_hdr->network[i];
	data_hdr->network[i] = '\0';

	/* copy time information */
	data_hdr->time.year 	= input_data_hdr->time.year;
	data_hdr->time.day     	= input_data_hdr->time.day;
	data_hdr->time.hour    	= input_data_hdr->time.hour;
	data_hdr->time.minute  	= input_data_hdr->time.minute;
	data_hdr->time.second  	= input_data_hdr->time.second;
	data_hdr->time.fracsec 	= input_data_hdr->time.fracsec;

	/* copy numerical information */
	data_hdr->nsamples = input_data_hdr->nsamples;

	data_hdr->number_time_corrections = 
			input_data_hdr->number_time_corrections;

	/* recover sample rate information */
	sample_rate            = (double) input_data_hdr->sample_rate;
	sample_rate_multiplier = (int) input_data_hdr->sample_rate_multiplier;

	sample_rate = input_data_hdr->sample_rate;
 
	sample_rate_multiplier = input_data_hdr->sample_rate_multiplier;

	/* adjust sample rate according to SEED definitions */
	if (sample_rate < 0) sample_rate = 1 / (-sample_rate);
	if (sample_rate_multiplier > 0) 
		sample_rate = sample_rate * sample_rate_multiplier;
	else if (sample_rate_multiplier < 0) 
		sample_rate = sample_rate / (-sample_rate_multiplier);

	/* put sample rate into the structure */
	data_hdr->sample_rate = sample_rate;

	/* copy activity, io, and data quality flags */
	data_hdr->activity_flags     = input_data_hdr->activity_flags;
	data_hdr->io_flags           = input_data_hdr->io_flags;
	data_hdr->data_quality_flags = input_data_hdr->data_quality_flags;

	data_hdr->number_blockettes = (int) input_data_hdr->number_blockettes;

	/* starting byte number of data */
	data_hdr->bod = input_data_hdr->bod;

	/* starting byte number of blockettes */
	data_hdr->bofb = input_data_hdr->bofb;

}
