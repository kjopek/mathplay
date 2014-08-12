#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <sys/time.h>

#define CHUNKSIZE 5

typedef enum {
    COL_IDX,
    ROW_IDX
} idx;

typedef struct
{
  // vector of values
  double *x;
  // not zero elements
  unsigned int nz;
  // PTR and IDX values
  unsigned int *ptr;
  unsigned int *ind;
  // size of matrix
  unsigned int n;
  // Compressed Row Storage or Compressed Column Storage?
  idx indexing;
} sm;

typedef struct
{
    unsigned int row;
    unsigned int col;
    double val;
} in_data;

sm * alloc_matrix(int nz, int n, idx indexing)
{
    sm *mat = (sm*) calloc(1, sizeof(sm));
    mat->x = (double*) calloc(nz, sizeof(double));
    mat->ptr = (unsigned int *) calloc(n+1, sizeof(mat->ptr));
    mat->ind = (unsigned int *) calloc(nz, sizeof(mat->ind));
    mat->n = n;
    mat->indexing = indexing;
    mat->nz = nz;
    
    if (!mat || !mat->x || !mat->ptr || !mat->ind) {
        free(mat);
        free(mat->x);
        free(mat->ptr);
        free(mat->ind);
        return NULL;
    }

    return mat;
}

void free_matrix(sm *mat)
{
    free(mat->x);
    free(mat->ptr);
    free(mat->ind);
    free(mat);
}

void fill_matrix(sm *matrix, in_data ** input_data, int nz)
{
    int n = 0;
    matrix->nz = nz;
    
    // we assume lower-triangular matrix on input!
    
    if (matrix->indexing == COL_IDX) {
        for (n=0; n<nz; ++n) {
            matrix->x[n] = input_data[n]->val;
            matrix->ind[n] = input_data[n]->col;
            if (n>0) {
                if (input_data[n]->row != input_data[n-1]->row) {
                    matrix->ptr[input_data[n]->row] = n;
                }
            } else {
                matrix->ptr[input_data[n]->row] = 0;
            }
        }
        matrix->ptr[input_data[nz-1]->row+1] = nz;
    } else {
        for (n=0; n<nz; ++n) {
            matrix->x[n] = input_data[n]->val;
            matrix->ind[n] = input_data[n]->row;
            if (n>0) {
                if (input_data[n]->col != input_data[n-1]->col) {
                    matrix->ptr[input_data[n]->col] = n;
                }
            } else {
                matrix->ptr[input_data[n]->col] = 0;
            }
        }
        matrix->ptr[input_data[nz-1]->col+1] = nz;
    }
}

int *elim_tree (const sm *A)
{
    int i, k, p, inext, *w, *parent, *ancestor, *prev ;
    unsigned n, *Ap, *Ai;
    n = A->n ; Ap = A->ptr ; Ai = A->ind ;
    parent = (int *) calloc (n, sizeof (int)) ;
    /* allocate result */
    w = (int *) calloc(n, sizeof (int)) ;
    ancestor = w ; prev = w + n ;

    for (k = 0 ; k < n ; k++) {
        parent [k] = -1 ;
        /* node k has no parent yet */
        ancestor [k] = -1 ;
        /* nor does k have an ancestor */
        for (p = Ap [k] ; p < Ap [k+1] ; p++) {
            i = Ai [p] ;
            for ( ; i != -1 && i < k ; i = inext) {
                /* traverse from i to k */
                inext = ancestor [i] ;
                /* inext = ancestor of i */
                ancestor [i] = k ;
                /* path compression */
                if (inext == -1) parent [i] = k ;
                /* no anc., parent is k */
            }
        }
    }
    free(w);
    return parent;
}

int bin_search(sm *mat, unsigned int r, unsigned int c)
{
    int imin, imax;
    int ikey = mat->indexing == COL_IDX ? c : r;

    if (mat->indexing == COL_IDX) {
        imin = mat->ptr[r];
        imax = mat->ptr[r+1];
    } else {
        imin = mat->ptr[c];
        imax = mat->ptr[c+1];
    }

    while (imax >= imin) {
        int mid = (imax+imin)/2;
        if (mat->ind[mid] == ikey) {
            return mid;
        }
        else if (mat->ind[mid] < ikey) {
            imin = mid + 1;
        }
        else {
            imax = mid-1;
        }
    }
    return -1;
}

void usage(char *progname)
{
	printf("Usage: %s <data file>\n", progname);
}

double get_time(struct timeval *t1, struct timeval *t2)
{
    return (t2->tv_sec - t1->tv_sec) + (t2->tv_usec - t1->tv_usec)/1000000.0;
}

void print_time(struct timeval *t1, struct timeval *t2)
{
    printf("%lf [s]\n", get_time(t1, t2));
}

int indata_cmp_row(const void *p1, const void *p2)
{
    in_data *d1 = (in_data*) p1;
    in_data *d2 = (in_data*) p2;
    if (d1->row < d2->row) return -1;
    if (d1->row > d2->row) return 1;
    return d1->col - d2->col;
}

int indata_cmp_col(const void *p1, const void *p2)
{
    in_data *d1 = (in_data*) p1;
    in_data *d2 = (in_data*) p2;
    if (d1->col < d2->col) return -1;
    if (d1->col > d2->col) return 1;
    return d1->row - d2->row;
}

void debug_matrix(sm *matrix)
{
    int i;
    printf("x = [");
    for (i=0;i<matrix->nz; ++i) {
        printf("%f ", matrix->x[i]);
    }
    printf("]\n");
    printf("idx = [");
    for (i=0;i<matrix->nz; ++i) {
        printf("%d ", matrix->ind[i]);
    }
    printf("]\n");
    printf("ptr = [");
    for (i=0;i<matrix->n+1; ++i) {
        printf("%d ", matrix->ptr[i]);
    }
    printf("]\n");

}

void print_data(in_data ** data, int nz)
{
    int i;
    for (i=0; i<nz; ++i) {
        printf("%d %d %e\n", data[i]->row, data[i]->col, data[i]->val);
    }
}

void chol(sm *mat, sm *mat_col)
{
    int i, j, k, ki;

    for (i=0; i<mat->n; ++i) {

        double x;
        double y;
        // assuming we have lower-triangular matrix L we can do this:
        int p = mat->ptr[i+1]-1;
        int q;
        
        x = mat->x[p];
        for (j=mat->ptr[i]; j<mat->ptr[i+1]-1; ++j) {
            x -= mat->x[j]*mat->x[j];
        }
        mat->x[p] = sqrt(x);

        #pragma omp parallel for private(q,y,k,ki) schedule(dynamic, CHUNKSIZE)
        for (j=mat_col->ptr[i]+1; j<mat_col->ptr[i+1]; ++j) {
            // current hotspot - Tim Davies has answer :)
            q = bin_search(mat, mat_col->ind[j], i);

            y = mat->x[q];
            for (k=mat->ptr[mat_col->ind[j]], ki=mat->ptr[i]; 
                 k<mat->ptr[mat_col->ind[j]+1] && mat->ind[k] < i; 
                 ++k) {
                // below we perform simple incrementation, maybe we should perform binary search?
                while (mat->ind[ki] < mat->ind[k] && ki<mat->ptr[i+1]) ++ki;
                if (mat->ind[ki] == mat->ind[k])
                    y -= mat->x[k] * mat->x[ki] * (mat->ind[ki] == mat->ind[k]);
            }
            mat->x[q] = y/mat->x[p];
        }
    }
}

int main(int argc, char ** argv)
{
    unsigned int nz = 0;
    unsigned int n = 0;
    unsigned int i, j;
    double val;

    struct timeval t1, t2;

    in_data ** input_data;
    sm *matrix;
    sm *matrix_col;

    if (argc != 2) {
        usage(argv[0]);
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL) {
        usage(argv[0]);
  	    return 2;
    }

    while (fscanf(fp, "%u %u  %lg", &i, &j, &val) != EOF) {
        ++nz;
        n = j > n ? j : n;
        n = i > n ? i : n;
    }
    fseek(fp, 0, SEEK_SET);
    ++n;;
    printf("Thr/proc #: %d / %d \n", omp_get_max_threads(), omp_get_num_procs());
    printf("Mat stat: nz=%u n=%u\n", nz, n);
    matrix = alloc_matrix(nz, n, COL_IDX);
    matrix_col = alloc_matrix(nz, n, ROW_IDX);

    input_data = (in_data**) calloc(nz, sizeof(in_data*));
    input_data[0] = (in_data*) calloc(nz, sizeof(in_data));
    n = 0;
    while (fscanf(fp, "%u %u %lg", &(input_data[n]->row), &(input_data[n]->col), &(input_data[n]->val)) != EOF) {
        ++n;
        input_data[n] = input_data[0]+n;
    }

    qsort((void*) input_data[0], nz, sizeof(in_data), indata_cmp_col);

    fill_matrix(matrix_col, input_data, nz);

    if (matrix->indexing == COL_IDX)
        qsort((void*) input_data[0], nz, sizeof(in_data), indata_cmp_row);
    else
        qsort((void*) input_data[0], nz, sizeof(in_data), indata_cmp_col);
    fill_matrix(matrix, input_data, nz);

    free(input_data[0]);
    free(input_data);
    fclose(fp);

    gettimeofday(&t1, NULL);
    chol(matrix, matrix_col);
    gettimeofday(&t2, NULL);
    printf("SOLVE time: ");
    print_time(&t1, &t2);

    //for (i=0; i<matrix->n; ++i) {
    //    for (j=matrix->ptr[i]; j<matrix->ptr[i+1]; ++j) {
    //        printf("%d %d %E\n", i, matrix->ind[j], matrix->x[j]);
    //    }
    //}

    free_matrix(matrix);
    free_matrix(matrix_col);
    return 0;

}
