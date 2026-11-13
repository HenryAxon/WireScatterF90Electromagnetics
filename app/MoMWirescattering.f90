program MoMWirescattering
! An algorithm for solving the wire scattering problem via the method of moments (MoM). Written for use as part of 
! ECE 541 at Colorado State university, and as such, is not especially user friendly. The ability to sweep lengths for 
! computing monostatic RCS across different electrical lengths is included via commented out code. Easy to change 

! Plotting options currently are to plot the current distribution along the wire induced by the incident electromagnetic
! wave. Commented out is a routine to plot the monostatic RCS as a function of length via a python matplotlib wrapper.
! good luck! gitignore should work now

! written in F90 (Fortran 1990) for gfortran compiler
    use stdlib_linalg, only: solve, det, eye
    use stdlib_quadrature, only: gauss_legendre
    use pyplot_module
    use subroutines 
    use pyplot_module, only: pyplot, wp => pyplot_wp

    implicit none
    type(pyplot) :: plt
    
    integer, parameter :: Nsubs = 100
    real(8), parameter :: L = 2
    real(8) :: a 
    !real(8), allocatable :: mesh(,:)
    real(8) :: d 
    real(8), parameter :: f = 300D6
    real(8), parameter :: c = 299792458
    real(8), parameter :: pi = DACOS(-1.D0)
    real(8), parameter :: mu = (4 * pi) * 1.0D-7
    real(8), parameter :: eps = 8.8541878188D-12
    real(8), parameter:: e0 = 1
    real(8), parameter :: beta = 2* pi * f/c
    real(8), parameter :: incidentAngle = 5*pi/12
    real(8), parameter :: omega = 2 * pi * f
    real(8) :: z(Nsubs)
    complex(8) :: In(Nsubs)
    complex(8) :: excitationIn(Nsubs)
    complex(8) :: zMatIn(Nsubs, Nsubs)
    complex(8) :: coeffresult(Nsubs)
    integer, parameter :: Nz = Nsubs +2
    integer :: i
    real(8) :: zplot(Nz)
    complex(8) :: Iplot(Nz)
    
    !real(8) ::  L(100) 
   
    !complex(8) :: RCS(100)
    complex(8) :: RCS
    complex(8) :: current(Nsubs)
    complex(8) :: Q
    !call linspace(0.25D0, 2.0D0, L)
    a = 12.5D-3
    
    !do i = 1, 100
        !a = L(i) / 160
        call linspace(0.0D0, L, zplot)
        

        !call meshing(Nsubs, L, a, mesh, d)
        d = L/ (Nsubs + 1)
        !d = L(i) / (Nsubs + 1)
        
        !call excitationVec(e0, incidentAngle, beta, Nsubs, L(i),excitationIn)
        call excitationVec(e0, incidentAngle, beta, Nsubs, L,excitationIn)
        
        !call fillZ(L(i), Nsubs, beta,a,f, zMatIn)
        call fillZ(L, Nsubs, beta,a,f, zMatIn)
        
        call solveSystem(excitationIn, zMatIn, coeffresult)
        
        call currents(Nsubs, d, coeffresult,zplot, Iplot)

        call QCompute(Nsubs,d,Iplot,beta,incidentAngle, Q)

        !call radarCrossSection(Nsubs,d,Q,E0,incidentAngle, mu, omega,pi, RCS(i))
        call radarCrossSection(Nsubs,d,Q,E0,incidentAngle, mu, omega,pi, RCS)
            
        call plt%initialize(grid=.true.,xlabel='Z (m)',ylabel='I(z) Amps', title='Plot of Surface Currents on Wire Scatterer with Theta = 30',legend=.true.)
        call plt%add_plot(zplot,real(Iplot),label='Currents',linestyle='b-o',markersize=5,linewidth=2)
        call plt%savefig('currents75_fixed_zmat.png', pyfile='currents75_fixed_zmat.py')

    !end do

    !call plt%initialize(grid=.true.,xlabel='Length (m)',ylabel='RCS', title='Plot of Radar Cross Section with Theta = 90',legend=.true.)
    !call plt%add_plot(L, abs(RCS),label='Radar Cross Section',linestyle='b-o',markersize=5,linewidth=2)
    !call plt%savefig('RCSsweepnodb_correctedzmat.png', pyfile='RCSsweepnodb_correctedzmat.py')
    
    
    
end program MoMWirescattering