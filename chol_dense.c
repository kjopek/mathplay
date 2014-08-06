#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

void chol_complex(complex double **A, complex double **L, unsigned int n)
{
// TODO
}

int chol(double **A, double **L, unsigned int n)
{
    unsigned int i=0;
    unsigned int j=0;
    unsigned int k=0;

    for (j=0; j<n; ++j) {
        L[j][j] = A[j][j];

        for (k=0; k<j; ++k) {
            L[j][j] -= L[j][k]*L[j][k];
        }

        if (L[j][j] < 0) {
            return 1;
        }

        L[j][j] = sqrt(L[j][j]);

        for (i=j+1; i<n; ++i) {
            L[i][j] = A[i][j];
            for (k=0; k<j; ++k) {
                L[i][j] -= L[i][k] * L[j][k];
            }
            L[i][j] /= L[j][j];
        }
    }
    return 0;
}

int main()
{
    double ** A = (double**) malloc(3*sizeof(double*));
    A[0] = (double*) malloc(3*3*sizeof(double));
    A[1] = A[0]+3;
    A[2] = A[0]+6;
    A[0][0] = 4;   A[0][1] = 12;  A[0][2] = -16;
    A[1][0] = 12;  A[1][1] = 37;  A[1][2] = -43;
    A[2][0] = -16; A[2][1] = -43; A[2][2] = 98;

    double ** L = (double**) malloc(3*sizeof(double*));
    L[0] = (double*) malloc(3*3*sizeof(double));
    L[1] = L[0]+3;
    L[2] = L[0]+6;
    chol(A, L, 3);

    if ((unsigned long) (A) % sizeof(double*) != 0 ||
        (unsigned long) (A[0]) % sizeof(double) != 0) {
        printf("A is not aligned\n");
    }

    if ((unsigned long) (L) % sizeof(double*) != 0 || 
        (unsigned long) (L[0]) % sizeof(double) != 0) {
        printf("L is not aligned\n");
    }

    int i, j;

    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            printf("%lf ", L[i][j]);
        }
        printf("\n");
    }

    free(A[0]); free(A);
    free(L[0]); free(L);

    return 0;
}

