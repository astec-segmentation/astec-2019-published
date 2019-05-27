/*************************************************************************
 * iobmp.c - I procedures for BMP raw images
 *
 * $Id: iobmp.c,v 1.1 1999/10/06 15:55:45 greg Exp $
 *
 * AUTHOR:
 * 
 * CREATION DATE: 
 *
 * ADDITIONS, CHANGES
 * adapted from test.c
 *         from bmp.zip, see the url http://www.ddj.com/ftp/1995/1995.03/
 *         author Dr. Dobb's
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include "bmptypes.h"
#include "bmpendian.h"
#include "readbmp.h"

#include <iobmp.h>

static int _VERBOSE_ = 1;



void *_readBmpImage( char *name, 
		     int *dimx, int *dimy, int *dimz )
{
  char *proc="_readBmpImage";
  void *buf = (void*)NULL;
  unsigned char *myBuf;

  FILE *fp;
  RGB **argbs;
  char **xorMasks, **andMasks;
  UINT32 *heights, *widths, row, col;
  UINT16 fileType;
  long filePos;
  int numImages, i;
  int rc;
    
    fp = fopen(name, "rb");
    if (fp == NULL) {
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: error in opening %s\n", proc, name );
      return( (void*)NULL );
    }


    /*
     * Read the first two bytes as little-endian to determine the file type.
     * Preserve the file position.
     */
    filePos = ftell(fp);
    rc = readUINT16little(fp, &fileType);
    if (rc != 0) {
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: error in getting file type %s\n", proc, name );
      return( (void*)NULL );
    }
    fseek(fp, filePos, SEEK_SET);

    /*
     * Read the images.
     */
    switch (fileType) {
    case TYPE_ARRAY:
	/*
	 * If this is an array of images, read them.  All the arrays we need
	 * will be allocated by the reader function.
	 */
	rc = readMultipleImage(fp, &argbs, &xorMasks, &andMasks, &heights,
			       &widths, &numImages);
	break;
    case TYPE_BMP:
    case TYPE_ICO:
    case TYPE_ICO_COLOR:
    case TYPE_PTR:
    case TYPE_PTR_COLOR:
	/*
	 * If this is a single-image file, we've a little more work.  In order
	 * to make the output part of this test program easy to write, we're
	 * going to allocate dummy arrays that represent what
	 * readMultipleImage would have allocated.  We'll read the data into
	 * those arrays.
	 */
	argbs = (RGB **)calloc(1, sizeof(RGB *));
	if (argbs == NULL)
	{
	    rc = 1005;
	    break;
	}
	xorMasks = (char **)calloc(1, sizeof(char *));
	if (xorMasks == NULL)
	{
	    vtfree(argbs);
	    rc = 1005;
	    break;
	}
	andMasks = (char **)calloc(1, sizeof(char *));
	if (andMasks == NULL)
	{
	    vtfree(argbs);
	    vtfree(xorMasks);
	    rc = 1005;
	    break;
	}
	heights = (UINT32 *)calloc(1, sizeof(UINT32));
	if (heights == NULL)
	{
	    vtfree(argbs);
	    vtfree(xorMasks);
	    vtfree(andMasks);
	    rc = 1005;
	    break;
	}
	widths = (UINT32 *)calloc(1, sizeof(UINT32));
	if (widths == NULL)
	{
	    vtfree(argbs);
	    vtfree(xorMasks);
	    vtfree(andMasks);
	    vtfree(heights);
	    rc = 1005;
	    break;
	}
	numImages = 1;

	/*
	 * Now that we have our arrays allocted, read the image into them.
	 */
	switch (fileType) {
	case TYPE_BMP:
	    rc = readSingleImageBMP(fp, argbs, widths, heights);
	    break;
	case TYPE_ICO:
	case TYPE_PTR:
	    rc = readSingleImageICOPTR(fp, xorMasks, andMasks, widths,
				       heights);
	    break;
	case TYPE_ICO_COLOR:
	case TYPE_PTR_COLOR:
	    rc = readSingleImageColorICOPTR(fp, argbs, xorMasks, andMasks,
					    widths, heights);
	    break;
	}
	break;
    default:
	rc = 1000;
    }
    fclose(fp);


  
    /*
     * At this point, everything's been read.  Display status messages based
     * on the return values.
     */
    switch (rc) {
    case 1000:
    case 1006:
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: File is not a valid bitmap file\n", proc );
      break;
    case 1001:
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: Illegal information in an image\n", proc );
      break;
    case 1002:
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: Legal information that I can't handle yet in an image\n", proc );
      break;
    case 1003:
    case 1004:
    case 1005:
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: Ran out of memory\n", proc );
      break;
    case 0:
      if ( _VERBOSE_ > 1 ) 
	fprintf( stderr, "%s: Got good data from file, writing results\n", proc );
      break;
    default:
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: Error reading file rc=%d\n", proc,rc );
      break;
    }

    /*
     * If the return value wasn't 0, something went wrong.
     */
    if (rc != 0)
    {
	if (rc != 1000 && rc != 1005)
	{
	    for (i=0; i<numImages; i++)
	    {
		if (argbs[i] != NULL)
		    vtfree(argbs[i]);
		if (andMasks[i] != NULL)
		    vtfree(andMasks[i]);
		if (xorMasks[i] != NULL)
		    vtfree(xorMasks[i]);
	    }
	    vtfree(argbs);
	    vtfree(andMasks);
	    vtfree(xorMasks);
	    vtfree(widths);
	    vtfree(heights);
	}
	return( (void*)NULL );
    }
    
    
    /*
     * Dump the images.
     */
    if ( _VERBOSE_ > 1 ) 
      fprintf (stderr, "%s: There are %d images in the file\n", proc, numImages);



    /*
     * my stuff:
     * just reading one bmp image
     */
    if ( (numImages > 0) &&
	 (argbs[0] != NULL) ) {

      buf = (void*)vtmalloc( widths[0]*heights[0]*3 * sizeof( unsigned char ) );
      if ( buf == (void*)NULL ) {
	if ( _VERBOSE_ ) 
	  fprintf( stderr, "%s: error in allocating data buffer for %s\n", proc, name );

	for (i=0; i<numImages; i++) {
	  if (argbs[i] != NULL)
	    vtfree(argbs[i]);
	  if (andMasks[i] != NULL)
	    vtfree(andMasks[i]);
	  if (xorMasks[i] != NULL)
	    vtfree(xorMasks[i]);
	}
	free(argbs);
	free(andMasks);
	free(xorMasks);
	free(widths);
	free(heights);

	return( (void*)NULL );
      }

      myBuf = (unsigned char*)buf;
      i = 0;
      for (row = 0; row < heights[0]; row++)
      for (col = 0; col < widths[0]; col++,i++) {
	myBuf[i]                          = argbs[0][row * widths[0] + col].red;
	myBuf[heights[0]*widths[0] + i]   = argbs[0][row * widths[0] + col].green;
	myBuf[2*heights[0]*widths[0] + i] = argbs[0][row * widths[0] + col].blue;
      }
      
      *dimx = widths[0];
      *dimy = heights[0];
      *dimz = 3;

    } else {
      if ( _VERBOSE_ ) 
	fprintf( stderr, "%s: no image or null image\n", proc );
      
      for (i=0; i<numImages; i++) {
	if (argbs[i] != NULL)
	  vtfree(argbs[i]);
	if (andMasks[i] != NULL)
	  vtfree(andMasks[i]);
	if (xorMasks[i] != NULL)
	  vtfree(xorMasks[i]);
      }
      vtfree(argbs);
      vtfree(andMasks);
      vtfree(xorMasks);
      vtfree(widths);
      vtfree(heights);
      
      return( (void*)NULL );
    }



    for (i=0; i<numImages; i++) {
      if (argbs[i] != NULL)
	free(argbs[i]);
      if (andMasks[i] != NULL)
	free(andMasks[i]);
      if (xorMasks[i] != NULL)
	free(xorMasks[i]);
    }
    vtfree(argbs);
    vtfree(andMasks);
    vtfree(xorMasks);
    vtfree(widths);
    vtfree(heights);
    
    return( buf );





    /*
     * old stuff from test.c
     */

    for (i=0; i<numImages; i++)
    {
	/*
	 * Loop through all the images that were returned.
	 */
      if ( _VERBOSE_ ) {
	fprintf (stderr, "%s: Doing image number %d\n\n", proc, i+1);
	fprintf (stderr, "%s: Image dimensions: (%ld,%ld)\n", proc, widths[i], heights[i]);
      }

      if (argbs[i] != NULL) {
	/*
	 * If the image has colors, dump them (BMP, color ICO and color
	 * PTR files
	 */
	fprintf(stderr, "Colors");
	for (row = 0; row < heights[i]; row++)
	  {
	    fprintf (stderr, "\n\nRow %ld pixels (R,G,B), hex values:\n",
		     row);
	    for (col = 0; col < widths[i]; col++)
	      {
		fprintf (stderr, "(%2.2x,%2.2x,%2.2x)",
			 argbs[i][row * widths[i] + col].red,
			 argbs[i][row * widths[i] + col].green,
			 argbs[i][row * widths[i] + col].blue);
	      }
	  }
      } else {
	/*
	 * If there image has no colors, say so.  (monochrome ICO and PTR
	 * files) 
	 */
	fprintf (stderr, "No color image\n");
      }


      
      if (xorMasks[i] != NULL) {
	/*
	 * If the image has an xor mask, dump it.  (ICO and PTR files)
	 */
	fprintf (stderr, "\nXOR mask\n");
	for (row = 0; row < heights[i]; row++)
	  {
	    for (col = 0; col < widths[i]; col++)
	      {
		fprintf (stderr, "%c",
			 xorMasks[i][row * widths[i] + col] ? '@' : '.');
	      }
	    fprintf (stderr, "\n");
	  }
      } else {
	/*
	 * If the image has no xor mask, say so.  (BMP files).
	 */
	fprintf (stderr, "No xor mask\n");
      }
      


      if (andMasks[i] != NULL) {
	/*
	 * If the image has an and mask, dump it.  (ICO and PTR files)
	 */
	fprintf (stderr, "\nAND mask\n");
	for (row = 0; row < heights[i]; row++)
	  {
	    for (col = 0; col < widths[i]; col++)
	      {
		fprintf (stderr, "%c",
			 andMasks[i][row * widths[i] + col] ? '@' : '.');
	      }
	    fprintf (stderr, "\n");
	  }
      } else {
	/*
	 * If the image has noand mask, say so.  (BMP files)
	 */
	fprintf (stderr, "No and mask\n");
      }
      

      if (i != numImages-1)
	fprintf (stderr, "\n------------------------------------------\n\n");
      
    }
    
    /*
     * Dumping is complete.  Free all the arrays and quit
     */
    for (i=0; i<numImages; i++)
    {
	if (argbs[i] != NULL)
	    vtfree(argbs[i]);
	if (andMasks[i] != NULL)
	    vtfree(andMasks[i]);
	if (xorMasks[i] != NULL)
	    vtfree(xorMasks[i]);
    }
    vtfree(argbs);
    vtfree(andMasks);
    vtfree(xorMasks);
    vtfree(widths);
    vtfree(heights);
    
    return( buf );
}








void IoBmp_verbose ( )
{
  if ( _VERBOSE_ <= 0 )
    _VERBOSE_ = 1;
  else 
    _VERBOSE_ += 1;
}

void IoBmp_noverbose ( )
{
  _VERBOSE_ = 0;
}



