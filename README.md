Mathematical games
==================

Cholesky factorization
----------------------

As far, there are two algorithms for dense and sparse Cholesky 
factorization. Dense version is rather easy to understand, while
sparse factorization code requires some explanation. Firstly, it
requires symbolic factorization that calculates fill-in required for 
this algorithm. Ordering of symbolic factorization should follow
ordering of matrix.

A total complexity of the algorithm is O(n) = n^3 and memory requirement
is O(n) = 4*nz + 2*n + 2 where nz is the number of non-zero elements
and n is the size of matrix. Memory consumption can be reduced to
O(n) = 3*nz + 2*n + 2.

