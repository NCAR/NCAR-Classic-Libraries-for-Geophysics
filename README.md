# NCAR-Classic-Libraries-for-Geophysics

FFTPACK 5.1 - A Fortran77 library of fast Fourier transforms 

Abstract
Library FFTPACK 5.1 contains 1D, 2D, and multiple fast Fourier subroutines, 
written in Fortran 77, for transforming real and complex data, real even and 
odd wave data, and real even and odd quarter-wave data. 

The FFTPACK 5.1 routines are grouped in triplets e.g. {CFFT1I, CFFT1F, CFFT1B} 
where suffix I denotes initialize, F forward (as in forward transform) and B backward. 
In an application program, before calling B or F routines for the first time, or before 
calling them with a different length, users must initialize an array by calling the 
I routine of the appropriate pair or triplet. 

Initialization routines need not be called each time before a B or F routine is called.

All of the transform routines in FFTPACK 5.1 are normalized.

Error messages are written to unit 6 by routine XERFFT. The standard version of XERFFT 
issues an error message and halts execution, so that no FFTPACK 5.1 routine will return 
to the calling program with error return IER different than zero. Users may consider 
modifying the STOP statement in order to call system-specific 
exception-handling facilities.

Caveat
The internal implementation of the FFTPACK 5.1 routines does not always conform to the Fortran standard.

References
(1) Vectorizing the Fast Fourier Transforms, by Paul Swarztrauber, Parallel Computations, G. Rodrigue, 
     ed., Academic Press, New York 1982.
(2) Fast Fourier Transforms Algorithms for Vector Computers, 
    by Paul Swarztrauber, Parallel Computing, (1984) pp.45-63.
