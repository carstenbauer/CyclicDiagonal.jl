# CyclicDiagonal is a Julia package to find approximate factorizations of square matrices M as products of diagonal and circulant matrices in the form M ≈ D1 * C * D2, where D1 and D2 are diagonal matrices and C is a circulant matrix

## Use
Factorization approximation is implemented in the function iterDCD. You can call it as follows

  D1, C, D2 = CyclicDiagonal.iterDCD( M )

The function also takes to additional parameters: the maximum number of iterations and the tolerance used to determine when to stop the iteration.

## Details
Factorization uses an algorithm developed by the author alternating the optimization from left and right sides, i.e.

1. Initialize D1, C, and D2
2. Find optimal D1 & C keeping D2 fixed
3. Find optimal C & D2 keeping D1 fixed
4. Return to step 2 if solution not within tolerance
