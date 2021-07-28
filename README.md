# NCAR-Classic-Libraries-for-Geophysics

Several mathematical libraries developed in the years 1970-1990 remain popular in the geophysics community.

These libraries are:

**FFTPACK**: A library of fast Fourier transforms

**FISHPACK**: Fortran subprograms for solving separable elliptic partial differential equations (PDEs)

**FISHPACK 90**: FISHPACK subprograms with a Fortran 90 interface

**MUDPACK**: Multigrid Fortran subprograms for solving separable and non-separable elliptic PDEs

**SPHEREPACK**: A Fortran library for modeling geophysical processes


All of these library routines are written primarily in Fortran 77. Their internal implementation does not 
always conform to the Fortran Standard. FISHPACK90 provides a Fortran 90 interface to the FISHPACK routines. 
Only MUDPACK is written with parallelism in mind; it uses OpenMP directives for shared-memory parallelism. 
The other libraries were designed to run on a single processor.

These libraries represent many person-years of development, and though they are no longer under development, 
NCAR continues to make them available to the public at no cost under a software licensing agreement. The 
libraries are best suited to Linux and UNIX environments and require a directory structure, tar, and gmake 
commands.

Paul N. Swarztrauber (deceased August 2011), John C. Adams, and Roland A. Sweet (both retired) spearheaded 
the development of these libraries for the NCAR Scientific Computing Division (SCD), which is now known as 
the Computational & Information Systems Laboratory (CISL).

-------------------------------------------------------------------------------------

### FFTPACK 5.1 - A Fortran77 library of fast Fourier transforms**

#### Abstract:
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
