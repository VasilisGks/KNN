# KNN
Serial and Parallel(using MPI communication protocol) implementation of KNN, distributed algorithm. Running on multiple cores (up to 64).

We will assume that we have two set of points Q and C  came from uniform distribution (inside the three dimensional unit cube)  containing Nq and Nc points respectively.

Searching for nearest neighbour, for a given point q (set Q) and a distance measure D: Find a point c (set C) such that D(q,c) is minimum.

<br />Run program with the following command:  ./KNNmpi n m k  Nc Nq
<br />
<br />n , m , k - grid dimensions: 2<sup>[12:16]</sup> 
<br />Nc, Nq - Size for each point set: 2<sup>[17:22]</sup>

