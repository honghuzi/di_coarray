# di_coarray
QTMC implemented with coarray Fortran

With coarrays, we can use very large arrays (5,000,000 * 10 * 16 cores * 10 nodes) in CNGrid, which extraordinarily increase the ensembles in the CTMC or QTMC model.

Coarrays can be declared as shared arrays locally or as distributed arrays with MPI.
