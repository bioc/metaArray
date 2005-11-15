#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <assert.h>
/* #include <Rmath.h> */
#include "meta.h"
/* #include "meta2.h" */

/************************************************

This function implements Durstenfeld's random 
permutation algorithm, for a vector of integers,
as described in The Art of Computer Programming: 
Volume 2 (Seminumerical Algorithms), by D. E. 
Knuth, pg. 145.

This implementation is Copyright (C) 2003 by 
Andrew Eckford, and is released into the public
domain as long as this entire comment is 
preserved in the source.

************************************************/

void perm(int *in, int *out, int t)  {

/* "in" is a pointer to the vector to be permuted, and the
permutation is written to the vector "out".  Both vectors
are of length "t". */

	int c;
	int j,k,temp;
	double U;

	/* Map "in" into "out"; "out" now becomes our scratch pad */

	for (c = 0; c < t; c++) {
		*(out+c) = *(in+c);
	}


	/* Initialize ... (note: j begins at t-1, not t as in Knuth because
	vectors in C go from 0 .. t-1 */

	/* stop at j=1 since no exchange is possible at j=0 */

	for (j = t-1; j > 0; --j)  {

		/* Generate U, a random variable uniformly
		distributed between zero and one */
	
		U = ((double)rand()) / ((double)RAND_MAX);
	
	
		/* Exchange: First set k, which is now a random variable
		between 0 and j-1 (note: drop the +1 from Knuth because
		vectors in C begin at index 0 */
	
		k = (int)(j*U);
	
		/* Now exchange out[k] and out[j] */
	
		temp = *(out+j);
		*(out+j) = *(out+k);
		*(out+k) = temp;

	}

}
