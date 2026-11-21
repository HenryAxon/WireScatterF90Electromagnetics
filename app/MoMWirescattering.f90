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
    !use fplot

    implicit none
    type(pyplot) :: plt
    
    integer, parameter :: Nsubs = 300
    !real(8), parameter :: L = 2
    real(8) :: a 
    !real(8), allocatable :: mesh(,:)
    real(8) :: d 
    ! select what metrics to compute: 0- current distribution at a given incidence and other
    !inputs, 1 - the RCS swept across a range of lengths.
    


    real(8), parameter :: f = 300D6
    real(8), parameter :: c = 299792458
    real(8), parameter :: pi = DACOS(-1.D0)
    real(8), parameter :: mu = (4 * pi) * 1.0D-7
    real(8), parameter :: eps = 8.8541878188D-12
    real(8), parameter:: e0 = 1
    real(8), parameter :: beta = 2* pi * f/c
    real(8), parameter :: incidentAngle = pi/2
    real(8), parameter :: omega = 2 * pi * f
    real(8) :: z(Nsubs)
    real(8) :: phase(Nsubs+2)
    !complex(8) :: In(Nsubs)
    complex(8),allocatable :: excitationIn(:)
    complex(8), allocatable :: zMatIn(:,: ) ! main difficulty with Fortran - allocatables and recognizing that I cannot define things later 
    complex(8), allocatable :: coeffresult(:)
    integer, parameter :: Nz = Nsubs +2
    integer :: i
    real(8),allocatable :: zplot(:)
    complex(8), allocatable :: Iplot(:)
    complex(8), allocatable :: current(:)
    integer :: x
    complex(8) :: Q
    real(8) ::  Lsingle
    real(8) :: asingle
    ! for computing a hihg density reference solution to compare convergence:
    integer, parameter :: numNs = 10
    integer, parameter :: denseSolnN = 1001
    integer :: ListOfN(numNs)
    real(8) :: currents_midpoint_residuals(numNs)
    complex(8) :: currents_midpoint
    complex(8) :: currents_midpoint_densesolve
    real(8) ::  L(100) 
    complex(8) :: RCS(100)
    integer :: Nzres
    complex(8) :: exciteDense(denseSolnN)
    complex(8) :: zmatdense(denseSolnN, denseSolnN)
    complex(8) :: currentcoeffdense(denseSolnN)
    complex(8) :: currentsdense(denseSolnN+2)
    real(8) :: dDense(denseSolnN)
    integer :: simFlag
    real(8) :: zplotdense(denseSolnN+2)

    simFlag = 1

    if (simFlag == 1) then ! branch of script for sweeping L and computing RCS
        allocate(excitationIn(Nsubs))
        allocate(zMatIn(Nsubs, Nsubs))
        allocate(coeffresult(Nsubs))
        allocate(zplot(Nz))
        allocate(Iplot(Nz))
        call linspace(0.25D0, 2.0D0, L)
        
    !complex(8) :: RCS
        !call linspace(0.25D0, 2.0D0, L)
    !a = 12.5D-3
    
        do i = 1, 100
            a = L(i) / 160
            call linspace(0.0D0, L(i), zplot)
        

            !call meshing(Nsubs, L, a, mesh, d)
            !d = L/ (Nsubs + 1)
            d = L(i) / (Nsubs + 1)
        
            call excitationVec(e0, incidentAngle, beta, Nsubs, L(i),excitationIn)
        !call excitationVec(e0, incidentAngle, beta, Nsubs, L,excitationIn)
        
            call fillZ(L(i), Nsubs, beta,a,f, zMatIn)
        !call fillZ(L, Nsubs, beta,a,f, zMatIn)
        
            call solveSystem(excitationIn, zMatIn, coeffresult)
        
            call currents(Nsubs, d, coeffresult,zplot, Iplot)

            call QCompute(Nsubs,d,Iplot,beta,incidentAngle, Q)

            call radarCrossSection(Nsubs,d,Q,E0,incidentAngle, mu, omega,pi, RCS(i))
        !call radarCrossSection(Nsubs,d,Q,E0,incidentAngle, mu, omega,pi, RCS)
            
        ! call plt%initialize(grid=.true.,xlabel='Z (m)',ylabel='I(z) Amps', title='Plot of Surface Currents on Wire Scatterer with Theta = 30',legend=.true.)
        ! call plt%add_plot(zplot,abs(Iplot),label='ABS(Currents)',linestyle='b--o',markersize=5,linewidth=2)
        ! call plt%add_plot(zplot,real(Iplot),label='RE(Currents)',linestyle='r-*',markersize=5,linewidth=2)
        ! call plt%add_plot(zplot,aimag(Iplot),label='IM(Currents)',linestyle='g-s',markersize=5,linewidth=2)
        ! call plt%add_str('plt.tight_layout()')
        ! call plt%savefig('currents30_mag_imag_real.png', pyfile='currents30_fixed_mag_imag.py')

        end do

        call plt%initialize(grid=.true.,xlabel='Length (m)',ylabel='RCS (m^2)', title='Plot of Radar Cross Section with Theta = 90, f=300 MHz',legend=.true.)
        call plt%add_plot(L, abs(RCS),label='Radar Cross Section',linestyle='b-o',markersize=5,linewidth=2)
        call plt%savefig('RCSsweepnodb_correctedzmattheta90_300Mhz.png', pyfile='RCSsweepnodb_correctedzmat_theta90300MHz.py')
        deallocate(excitationIn,zMatIn,coeffresult,zplot,Iplot)
        return
    else if (simFlag == 2) then
        allocate(excitationin(Nsubs))
        allocate(zMatIn(Nsubs,Nsubs))
        allocate(zplot(Nsubs+2))
        allocate(Iplot(Nsubs+2))
        allocate(coeffresult(Nsubs))
        !complex(8) :: RCS
        Lsingle = 2
        asingle = 12.5D-3
        d = Lsingle/ (Nsubs + 1)
        call linspace(0.0D0, Lsingle, zplot)

        call excitationVec(e0, incidentAngle, beta, Nsubs, Lsingle,excitationIn)


        call fillZ(Lsingle, Nsubs, beta,asingle,f, zMatIn)

        call solveSystem(excitationIn, zMatIn, coeffresult)

        call currents(Nsubs, d, coeffresult,zplot, Iplot)

        !call QCompute(Nsubs,d,Iplot,beta,incidentAngle, Q)

        !call radarCrossSection(Nsubs,d,Q,E0,incidentAngle, mu, omega,pi, RCS)
        ! adding a plotting sequence for phase -    
        phase = atan2(real(Iplot), aimag(Iplot)) * 180.0d0 / pi! I think this one line is all i need.
        
            ! plotting with the adapted matplotlib fortran library - currents
        call plt%initialize(grid=.true.,xlabel='Z (m)',ylabel='I(z) Amps', title='Plot of Surface Currents on Wire Scatterer with Theta = 75',legend=.true.)
        call plt%add_plot(zplot,abs(Iplot),label='ABS(Currents)',linestyle='b--o',markersize=5,linewidth=2)
        call plt%add_plot(zplot,real(Iplot),label='RE(Currents)',linestyle='r-*',markersize=5,linewidth=2)
        call plt%add_plot(zplot,aimag(Iplot),label='IM(Currents)',linestyle='g-s',markersize=5,linewidth=2)
        call plt%add_str('plt.tight_layout()')
        call plt%savefig('currents75_mag_imag_real_new.png', pyfile='currents75_fixed_mag_imag_real_new.py')
        

                    ! plotting with the adapted matplotlib fortran library - phase
        call plt%initialize(grid=.true.,xlabel='Z (m)',ylabel='Phase (Degress)', title='Surface Current Phasing Along Wire, Theta = 30',legend=.true.)
        call plt%add_plot(zplot,phase,label='ABS(Currents)',linestyle='b--o',markersize=5,linewidth=2)

        call plt%add_str('plt.tight_layout()')
        call plt%savefig('Phase30.png', pyfile='Phase30.py')
        deallocate(zmatin, excitationin, Iplot, coeffresult)


        return
    else
            ! must do dense solve OUTSIDE the loop
        Lsingle = 2.0D0
        asingle = 12.5D-3
        d = Lsingle/ (denseSolnN + 1)
        call linspace(0.0D0, Lsingle, zplotdense)
        call excitationVec(e0, incidentAngle, beta, denseSolnN, Lsingle, exciteDense)


        call fillZ(Lsingle, denseSolnN, beta,asingle,f, zmatdense)

        call solveSystem(exciteDense, zmatdense, currentcoeffdense)

        call currents(denseSolnN, d, currentcoeffdense,zplotdense,currentsdense)

        currents_midpoint_densesolve = currentcoeffdense((denseSolnN + 1) / 2)
        ListofN = [5,11,21,31,51,61,101,151,201,301] ! must use odd numbers to get a center properly

        do  x = 1, numNs
            
            !allocate(currents_midpoint_residuals(ListOfN(x)))
            !allocate(In(ListOfN(x)))
            !allocate(z(ListOfN(x)))
            allocate(excitationIn(ListofN(x)))
            allocate(zMatIn(ListofN(x),ListofN(x)))
            allocate(coeffresult(ListofN(x)))
            Nzres = ListofN(x) + 2
            allocate(zplot(Nzres))
            allocate(Iplot(Nzres))
            


            

                
                
                
            d = Lsingle/ (ListOfN(x) + 1)

                

            call excitationVec(e0, incidentAngle, beta, ListofN(x), Lsingle,excitationIn)


            call fillZ(Lsingle, ListofN(x), beta,asingle,f, zMatIn)

            call solveSystem(excitationIn, zMatIn, coeffresult)

            call currents(ListofN(x), d, coeffresult,zplot, Iplot)

            !call QCompute(ListofN(x),d,Iplot,beta,incidentAngle, Q)

            !call radarCrossSection(ListofN(x),d,Q,E0,incidentAngle, mu, omega,pi, RCS)
                ! copute residuals by comparaing excitation returned from one N value to the next
            currents_midpoint = coeffresult((ListOfN(x) + 1) / 2)

            currents_midpoint_residuals(x) = abs(currents_midpoint - currents_midpoint_densesolve) / currents_midpoint_densesolve

            deallocate(coeffresult,zplot,zMatIn,excitationIn,Iplot)
        end do


            

            ! plot the convergence
        call plt%initialize(grid=.true.,xlabel='N',ylabel='Res(I_n)', title='Residual Convergence of Surface Current Coefficients',legend=.true.)
        call plt%add_plot(real(ListofN,wp),currents_midpoint_residuals,label='ABS(Currents)',linestyle='b--o',markersize=5,linewidth=2)! type problem because list of N is an integer array
        ! and x is integer, anyhow plotting doesnt like that. if I turn it to real it owrks just fine!
        call plt%add_str('plt.tight_layout()')
        call plt%savefig('currentconvergence.png', pyfile='currentconvergence.py')
        return
        

    end if


end program MoMWirescattering