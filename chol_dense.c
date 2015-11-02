#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* for Fortran - j*n + i */
//#define IDX(i, j, n)	(((i) * (n)) + (j))
#define IDX(i, j, n) (((j)+ (i)*(n)))

int
chol(double *A, unsigned int n)
{
	unsigned int i;
	unsigned int j;
	unsigned int k;

	for (j = 0; j < n; j++) {	
		for (i = j; i < n; i++) {
			for (k = 0; k < j; ++k) {
				A[IDX(i, j, n)] -= A[IDX(i, k, n)] *
				    A[IDX(j, k, n)];
			}
		}

		if (A[IDX(j, j, n)] < 0.0) {
			return (1);
		}

		A[IDX(j, j, n)] = sqrt(A[IDX(j, j, n)]);
		for (i = j + 1; i < n; i++) 
			A[IDX(i, j, n)] /= A[IDX(j, j, n)];
	}

	return (0);
}

int
main()
{
	double *A;
	int i, j, n, ret;

	n = 3;
	A = calloc(n*n, sizeof(double));
	assert(A != NULL);

	A[IDX(0, 0, n)] = 4.0;   A[IDX(0, 1, n)] = 12.0;  A[IDX(0, 2, n)] = -16.0;
	A[IDX(1, 0, n)] = 12.0;  A[IDX(1, 1, n)] = 37.0;  A[IDX(1, 2, n)] = -43.0;
	A[IDX(2, 0, n)] = -16.0; A[IDX(2, 1, n)] = -43.0; A[IDX(2, 2, n)] = 98.0;

	if (chol(A, 3)) {
		fprintf(stderr, "Error: matrix is either not symmetric or not positive definite.\n");
	} else {
		fprintf(stdout, "Tri(L) = \n");
		for (i = 0; i < 3; i++) {
			for (j = 0; j <= i; j++)
				printf("%2.8lf\t", A[IDX(i, j, n)]);
			printf("\n");
		}
	}

	free(A);
	return 0;
}

