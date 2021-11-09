# NCAR-Classic-Libraries-for-Geophysics

Several mathematical libraries developed at the National Center for Atmospheric Research (NCAR) in the years 1970-1990 remain popular in the geophysics community.

Those libraries are:

[**FFTPACK**](https://github.com/NCAR/NCAR-Classic-Libraries-for-Geophysics/tree/main/FFTPack): A library of fast Fourier transforms

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

###### FFTPACK5, Copyright (C) 2004-2011, Computational Information Systems Laboratory, University Corporation for Atmospheric Research
-------------------------------------------------------------------------------------

### FFTPACK 5.1 - A Fortran77 library of fast Fourier transforms

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


### FISHPACK - Efficient FORTRAN Subprograms for the Solution of Separable Elliptic Partial Differential Equations

#### Abstract:
FISHPACK contains a collection of Fortran77 subroutines that solve second- and fourth-order finite difference approximations to separable elliptic Partial Differential Equations (PDEs). These include Helmholtz equations in cartesian, polar, cylindrical, and spherical coordinates, as well as more general separable elliptic equations. The solvers use the cyclic reduction algorithm. When the problem is singular, a least-squares solution is computed. Singularities induced by the coordinate system are handled, including at the origin r=0 in cylindrical coordinates, and at the poles in spherical coordinates.

Test programs are provided for the 19 solvers. Each serves two purposes: as a template to guide you in writing your own codes utilizing the FISHPACK solvers, and as a demonstration on your computer that you can correctly produce FISHPACK executables.

The FISHPACK library and programs are intended to be installed on your computer using the Makefile provided when you download the files in this distribution. The Makefile builds the library and driver executables under the compiler you specify when you run "make".

If your application requires solution of nonseparable elliptic PDEs, or a mix of separable and nonseparable ones, consider using the MUDPACK library instead of FISHPACK. MUDPACK uses multigrid iteration to approximate separable and nonseparable elliptic PDEs. The software is available on NCAR's web pages. If you are solving separable elliptic PDEs only, and prefer Fortran90 syntax, then you may want to use FISHPACK90, also available on NCAR's web pages. Both FISHPACK and FISHPACK90 have the same functionality.

CAVEAT: FISHPACK source code is known to break the Fortran77 Standard in various ways. In particular, some of the routines pass arguments of one type and use them as another. We have not compiled a comprehensive list of FISHPACK infractions. Prospective users who require complete adherence to the standard for their applications are advised that this package is not compliant.

###### FISHPACK, Copyright (C) 2004-2011, Computational Information Systems Laboratory, University Corporation for Atmospheric Research
-------------------------------------------------------------------------------------


### FISHPACK90 - Efficient FORTRAN Subprograms for the Solution of Separable Elliptic Partial Differential Equations

#### Abstract:
FISHPACK90 is a modernization of the original FISHPACK, employing Fortran90 to slightly simplify and standardize the interface to some of the routines. It is a collection of Fortran programs and subroutines that solve second- and fourth-order finite difference approximations to separable elliptic Partial Differential Equations (PDEs). These include Helmholtz equations in cartesian, polar, cylindrical, and spherical coordinates, as well as more general separable elliptic equations. The solvers use the cyclic reduction algorithm. When the problem is singular, a least-squares solution is computed. Singularities induced by the coordinate system are handled, including at the origin r=0 in cylindrical coordinates, and at the poles in spherical coordinates.

Test programs are provided for the 19 solvers. Each serves two purposes: as a template to guide you in writing your own codes utilizing the FISHPACK90 solvers, and as a demonstration on your computer that you can correctly produce FISHPACK90 executables.

The FISHPACK90 library and programs are intended to be installed on your computer using the Makefile provided when you download the files in this distribution. The Makefile builds the library and driver executables under the compiler you specify when you run "make".

If your application requires solution of nonseparable elliptic PDEs, or a mix of separable and nonseparable ones, consider using the MUDPACK library instead of FISHPACK90. MUDPACK uses multigrid iteration to approximate separable and nonseparable elliptic PDEs. The software is available on NCAR's web pages. If you are solving separable elliptic PDEs only, and prefer Fortran77 syntax, then you may want to use FISHPACK, also available on NCAR's web pages. Both FISHPACK and FISHPACK90 have the same functionality, though their calling sequences are slightly different and the packages must not be used interchangably without making the appropriate syntax changes.

CAVEAT:

FISHPACK90 source code is known to break the Fortran Standard in various ways. In particular, the subsidiary FFT routines sometimes pass arguments of one type and use them as another. We have not compiled a comprehensive list of FISHPACK90 infractions. Prospective users who require complete adherence to the standard for their applications are advised that this package is not compliant.

###### FISHPACK90, Copyright (C) 2004-2011, Computational Information Systems Laboratory, University Corporation for Atmospheric Research
-------------------------------------------------------------------------------------


### MUDPACK: Multigrid Software for Elliptic Partial Differential Equations

#### Abstract:
MUDPACK was first released in 1990, and has remained at version 5.0.1 the past five years. It is a collection of portable, mostly Fortran 77 subprograms (the code also employs a few Fortran 90 extensions). Its purpose is the efficient solving of linear elliptic Partial Differential Equations (PDEs) -- both separable and nonseparable -- using multigrid iteration. MUDPACK solvers can achieve parallel speedup via OpenMP directives in the 5.0.1 subroutines, but the user will need to activate the directives by providing appropriate OpenMP compiler options when the library and application are built. Speedup is dependent on problem size and characteristics of the shared multi-processor platform where the code is compiled and run.

###### MUDPACK, Copyright (C) 2004-2011, Computational Information Systems Laboratory, University Corporation for Atmospheric Research
-------------------------------------------------------------------------------------


### SPHEREPACK 3.2 - A Package for Modeling Geophysical Processes

#### Abstract:
SPHEREPACK 3.2 is a collection of FORTRAN77 programs and subroutines facilitating computer modeling of geophysical processes. The package contains subroutines for computing common differential operators including divergence, vorticity, latitudinal derivatives, gradients, the Laplacian of both scalar and vector functions, and the inverses of these operators. For example, given divergence and vorticity, the package can be used to compute velocity components, then the Laplacian inverse can be used to solve the scalar and vector Poisson equations. The package also contains routines for computing the associated Legendre functions, Gauss points and weights, multiple fast Fourier transforms, and for converting scalar and vector fields between geophysical and mathematical spherical coordinates.

Example programs are provided for solving these equations on the full sphere:

advection
Helmholz
shallow-water
Each program serves two purposes: as a template to guide you in writing your own codes utilizing the SPHEREPACK routines, and as a demonstration on your computer that you can correctly produce SPHEREPACK executables.

The SPHEREPACK library and programs are intended to be installed on your computer using the Makefile provided when you download the files in this distribution. The Makefile builds the library and driver executables under the compiler you specify when you run "make".

This work was partially supported by the Computer Hardware, Advanced Mathematics, and Model Physics (CHAMMP) Program which is administered by the Office of Energy Research under the Office of Health and Environmental Research in the U.S. Department of Energy, Environmental Sciences Division.

###### SPHEREPACK, Copyright (C) 2004-2011, Computational Information Systems Laboratory, University Corporation for Atmospheric Research
