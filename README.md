# momprojectfive
The FORTRAN f90 code to run my method of moments (MoM) wire scatterer of an electromagnetic plane wave simulation. The code contains two main files:
MoMWireScattering.f90 - the main file containting references to the FORTRAN std_lib, a dependency, as well as calling all subroutines needed to perform the method of moments and solve for the surface current sources, and compute the radar cross section (RCS) of a cylindrical wire. 

subroutinesMoM.f90 - the module contianing all necessary subroutines to run MoM algorithm

The code was written in FortMata (thanks to Christopher Erickson!) which includes a plotting module based on matplotlib, however this is not a native FORTRAN feature, and so the code includes commented out code that is only compatible with fortmata capabilities. I have added a pyplot-fortran external library to allow plotting in fortran.

## Dependencies:

std_lib - the standard fortran library which despite the name is not a native library adn should be built with cmake from the source code git repository. This includes the linear algebra solvers used to compute theh currents via matrix inversion, whihc is not a native FORTRAN feature. It also includes a quadrature feature which I use for the numerical integrations.

pyplot-fortran : an analogue to pyplot from matplotlib as seen in python. This allows simple plotting capability without main core matplotlib features but enough to get the job done without having to use python externally, which was good enough for me. 

## Building the Solver

The solver integreates the external libaries via the fortran package manager which should be installed to run this project. To build the project, because FORTRAN line length rules were not strictly followed, you must pass the line (fpm build --flag "-ffree-line-length-none") which allows the compiler to ignore the line length limit that F90 imposes automatically at 130 some odd lines. 

To run the solver you must then again write (fpm run --flag "-ffree-line-length-none"), which should work just fine if the code compiled properly.