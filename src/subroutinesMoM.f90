module subroutines
! use, intrinsic :: iso_fortran_env, only: dp => real64
use stdlib_linalg, only: solve, det, eye !! using the standard fortran library to get a linear algebra solver and 
! gauss quadrature capabilities. build the library with CMAKE, ninja, fypp!!
use stdlib_quadrature, only: gauss_legendre
use pyplot_module, only: pyplot, wp => pyplot_wp
contains 

    subroutine linspace(from, to, array)
        implicit none
        ! we define the inputs, outputs, their types, as well as all other
        ! needed variables
        real(8), intent(in) :: from, to
        real(8), intent(out) :: array(:)
        real(8) :: range
        integer :: n, i

        n = size(array)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from 
            return
        end if

        do i = 1, n
            array(i) = from + range * (i-1) / (n-1)
        end do

    end subroutine linspace  
    
        
    !subroutine meshing(N, L, a, mesh, d)
        ! mesher function for a wire scatterer
        
        ! params:
        ! N - number of nodes
        ! L - length of the wire
        ! a - the radius of the wire
        
        ! Outputs:
        ! mesh - Nx2 array of positions on the wire to compute SIE at
        ! d - the z-direction distances between mesh nodes.
        !implicit none
        ! defines the type of the input to the mesh
        !real(8), intent(in) :: a
        !integer, intent(in) :: N
        !real(8), intent(in) :: L
        ! defines the type of the outputs
        !real(8), intent(out) :: mesh(N,2) ! not sure about this
        !real(8), intent(out) :: d
        !real(8) :: yvals(N)
        !real(8) :: xvals(2)
        !integer :: i
        !integer :: j
        !xvals(1) = 0.0D0
        !xvals(2) = a
        
        !real(dp) :: mesh(N,N)
        
        !d = L / (N + 1)
        !call linspace(0.0D0, L, yvals) 
        
        ! filling the matrix of positions with the x values being the points on the center z-axis
        ! and the y values being the points on wire surface or the z'- axis
        !do i = 1, N
            
        !    mesh(i,1) = yvals(i)
        !    mesh(i,2) = xvals(i)
            
        !end do
           
    !end subroutine meshing
    
    subroutine excitationVec(e0, theta, beta, N,L, excitation)
        ! subroutine to compute the [V_m] excitation vector of the plane wave
        ! incident on the wire. 
        
        ! params:
        ! e0 - initial input E-field vector for E at the bottom of the wire surface (a,0) 
        ! theta - the angle of incidence
        ! beta - the wave number
        ! N - number of subdivisions of the mesh
        
        ! Outputs:
        ! excitation - the vector of excitation values
        
        real(8), intent(in) :: e0
        real(8), intent(in) :: theta
        real(8), intent(in) :: beta
        integer, intent(in) :: N
        real(8), intent(in) :: L
        complex(8), intent(out) :: excitation(N)
        !real(8) :: zprime(N)
        real(8) :: zprime
        integer :: i
        complex(8) :: comp
        real(8) :: d
        d = L/ (N +1)
        comp = (0.0, 1.0)
        !call linspace(0.0D0, L, zprime)
        
        
        
        do i = 1, N
            zprime = i * d
            excitation(i) = - e0 * sin(theta) * exp(comp * beta * cos(theta) * zprime)
        end do
       
        
    end subroutine excitationVec
    
    subroutine sciOne(m, n, d, Nsub, a, beta, integralOne)
        ! sciOne computes the first integral 
        ! params:
        ! (n,m) - the nodes in the mesh on z and z' axes to compute integrals at
        ! d - distance between mesh codes
        ! N - number of subdivisions
        
        ! Outputs:
        ! integralOne - the integral!
  
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(8), intent(in) :: d
        integer, intent(in) :: Nsub
        real(8), intent(in) :: a
        real(8), intent(in) :: beta
        integer :: i
        integer, parameter :: k = 5
        real(8) :: x(k), w(k) ! pre-defining the gauss weigths and nodes
        complex(8), intent(out) :: integralOne
        
        real(8) :: analyticalPart
        complex(8) :: numericalMequalN
        complex(8) :: numericalMdneN
        real(8) :: RMequalN(k)
        real(8) :: RMdneN(k)
        real(8) :: interval(2)
        complex(8) :: comp
        
        comp = (0.0, 1.0)
        
        
        if (m == n) then
            ! perform quadrature integration from -d to d on the z axis
            
            analyticalPart =  2 * cos(beta * a) * LOG((d + SQRT(d**2 + a**2)) / a)
            call gauss_legendre(x, w, [0.0d0, d])
            RMequalN = sqrt(x ** 2 + a**2)
            numericalMequalN = 2*sum(w * (1 - x/d) * ((exp(-comp * beta * RMequalN))/RMequalN - cos(beta*a)/RMequalN))
            integralOne = numericalMequalN+ analyticalPart
            return
        else if (m /= n) then
            ! perform gaussian quadrature, define the  
            call gauss_legendre(x, w, [-d, d])
            RMdneN = sqrt(((n - m) * d - x)**2 + a**2)
            numericalMdneN = sum(w * (1- abs(x)/d) * exp(-comp *beta * RMdneN)/RMdneN)
            
            integralOne = numericalMdneN
            return
        end if
        
       
    end subroutine sciOne 
    
    subroutine sciTwo(m,n, d, a, Nsub,beta, integralTwo)
        ! sciTwo computes the first integral 
        ! params:
        ! (n,m) - the nodes in the mesh on z and z' axes to compute integrals at
        ! d - distance between mesh codes
        ! N - number of subdivisions
        
        ! Outputs:
        ! integralTwo - the integral result!
        integer, intent(in) :: m
        integer, intent(in) :: n
        real(8), intent(in) :: d
        real(8), intent(in) :: a
        integer, intent(in) :: Nsub
        real(8), intent(in) :: beta
        
        complex(8), intent(out) :: integralTwo
        complex(8) :: comp
        
        real(8) :: Rone 
        real(8) :: Rtwo
        real(8) :: R
        
        Rone = sqrt(((n - m - 1) * d) ** 2 + a ** 2)
        Rtwo = sqrt(((n - m + 1) * d) ** 2 + a ** 2)
        R = sqrt(((n - m) * d) ** 2 + a ** 2)
        comp = (0.0, 1.0)
        
        integralTwo = (exp(-comp * beta * Rone) / Rone - 2 * (exp(-comp * beta * R) / R) + exp(-comp * beta * Rtwo) / Rtwo)
        
        
    end subroutine sciTwo
    
    subroutine sciThree(n, d,beta,theta, Nsub, integralthree)
        ! computing the integral needed for RCS computation - another far field radiation integral - analytical!
        ! Params:
        ! n - the cell being examined integer
        ! d - the real valued width of a cell
        ! beta - the wave number which is real valued
        ! theta - the angle of incidence
        ! Nsub - the integer number of subdivisions
        
        ! Outputs:
        ! integralfour - the integral scifour complex number
        integer, intent(in) :: n
        real(8), intent(in) :: d
        real(8), intent(in) :: beta
        real(8), intent(in) :: theta
        integer, intent(in) :: Nsub
        complex(8), intent(out) :: integralthree
        complex(8) :: compUnit

        compUnit = (0.0D0, 1.0D0)
        if (abs(cos(theta)) <= 1.0D-8) then
            integralthree = d
            return
        else if (abs(cos(theta)) > 1.0D-8) then
            integralthree = (exp((n-1)*compUnit * beta * d * cos(theta))) / (compUnit * beta * cos(theta)) * (exp(compUnit*beta*d*cos(theta)) - 1)
            return
        end if

    end subroutine sciThree
    
    subroutine sciFour(n, d,beta,theta, Nsub, integralfour)
        ! computing the integral needed for RCS computation - the radiation integral in the far field!
        ! Params:
        ! n - the cell being examined integer
        ! d - the real valued width of a cell
        ! beta - the wave number which is real valued
        ! theta - the angle of incidence
        ! Nsub - the integer number of subdivisions
        
        ! Outputs:
        ! integralfour - the integral scifour complex number

        integer, intent(in) :: n
        real(8), intent(in) :: d
        real(8), intent(in) :: beta
        real(8), intent(in) :: theta
        integer, intent(in) :: Nsub
        complex(8), intent(out) :: integralfour
        complex(8) :: compUnit
        complex(8) :: intpartone
        complex(8) :: intparttwo
        complex(8) :: intpartthree

        compUnit = (0.0D0, 1.0D0) 
        
        if (abs(cos(theta)) > 1.0D-8) then
            intpartone = (exp(compUnit * beta * d * cos(theta)* (n-1))) 
            intparttwo = compUnit * beta * cos(theta)
            intpartthree = d * (exp(compUnit*beta*d* cos(theta))) - 1.0D0/(intparttwo) * (exp(compUnit*d*beta*cos(theta)) - 1)
            integralfour = intpartone / intparttwo * intpartthree
            return
        else if (abs(cos(theta)) <= 1.0D-8) then
            integralfour = 0.5D0 * d **2
            return
        end if


    end subroutine sciFour
    
    subroutine fillZ(L,Nsub,beta,a,f,zMat)
        ! fill the Z - matrix with the proper values/entries
        ! params:
        ! L - length of the wire
        ! N - the number of subdivisions
        ! beta - the wavenumber of incident wave
        
        ! Ouputs:
        ! zMat - the filled z matrix
        
        real(8), intent(in) :: L
        integer, intent(in) :: Nsub
        real(8), intent(in) :: beta
        real(8), intent(in) :: f
        real(8), intent(in) :: a
        complex(8), intent(out) :: zMat(Nsub,Nsub)
        real(8) :: d 
        real(8), parameter :: pi = DACOS(-1.D0)
        integer :: n
        integer :: m 
        real(8) :: omega
        complex(8) :: comp
        real(8) :: mu
        complex(8) :: intOne
        complex(8) :: intTwo
        
        comp = (0.0, 1.0)
        omega = 2 * pi * f
        mu = 4 * pi * 1.0D-7
        
        d = L / (Nsub + 1)
        
        do n = 1, Nsub
            do m = 1, Nsub
                call sciOne(m,n, d, Nsub, a,beta ,intOne)
                call sciTwo(m,n, d, a, Nsub, beta ,intTwo)
                zMat(n,m) = comp * omega * mu / (4 * pi) * (intOne + intTwo/((beta ** 2) * d))
            end do
        end do
    end subroutine fillZ

    subroutine rooftopBasis(z, zn, d, Tz)
        ! Compute the rooftop basis function values at a particular position
        ! on the z-axis
        
        ! params:
        ! z - the z position where we seek to compute the current
        ! zn - the z position we are offsetting the center to
        ! d - the width of a discretization
        
        ! Outputs:
        ! Tz - the value of the rooftop at a position
        
        real(8), intent(in) :: z
        real(8), intent(in) :: zn
        real(8), intent(in) :: d
        real(8), intent(out) :: Tz
        
        if (abs(z-zn) <= d) then
            Tz = 1 - abs(z-zn) / d    
            return 
        else if (abs(z-zn) > d) then
            Tz = 0 
            return
        end if
    end subroutine rooftopbasis
    
    
    subroutine solveSystem(vVec, zMat, coeff)
        ! Solves the system of equations and spits out the coefficients
        ! to be used to solve for the currents
        
        ! Params:
        ! vVec - the excitation vector
        ! zMat - the 
        
        ! Outputs:
        ! coeff - the coefficients resulting from a solve on the system
        
        complex(8), intent(in) :: vVec(:)
        complex(8), intent(in) :: zMat(:,:)
        complex(8), intent(out) :: coeff(:)
        
        coeff = solve(zMat, vVec)
    end subroutine solveSystem
    
    
    subroutine currents(Nsub, d,coefficients,zeval, In)
        ! compute the currents from the resulting coefficients        
        ! Params: 
        ! N - the number of subdivisions
        ! coefficients
        
        ! Outputs: 
        ! In - the currents
        
        integer, intent(in) :: Nsub
        real(8), intent(in) :: d
        complex(8), intent(in) :: coefficients(:)
        real(8), intent(in) :: zeval(:)
        complex(8), intent(out) :: In(size(zeval))
        
        integer :: i
        integer :: n
        real(8) :: zn
        real(8) :: Tz
        In = (0.0D0, 0.0D0)

        do i = 1, size(zeval)
            do n = 1, Nsub
                zn = n * d
                call rooftopBasis(zeval(i), zn, d, Tz)
                In(i) = In(i) + coefficients(n) * Tz
            end do
        end do
        
    end subroutine currents

    subroutine QCompute(Nsubs,d,current,beta,theta, Q)
    ! A subroutine written with the goal of computing the Q summation which is the primary quantity
    ! used for the computation of RCS!
    ! Params:
    ! Nsubs - integer number of nodes
    ! d - real size of a mesh subdivision
    ! current - the complex valued current array 
    ! beta - the wavenumber real values
    ! theta - the real valued angle of incidence 

    ! Outputs:
    ! Q - the q summation of current contributions used to compute RCS
        integer, intent(in) :: Nsubs
        real(8), intent(in) :: d
        complex(8), intent(in) :: current(:)
        real(8), intent(in) :: beta
        real(8), intent(in) :: theta
        complex(8), intent(out) :: Q
        complex(8) :: intthree
        complex(8) :: intfour

        integer :: n
        !integer :: m 
        complex(8) :: unitComp
        unitComp = (0.0, 1.0)
        Q = 0.0
        do n = 2, Nsubs+1
            call sciThree(n, d,beta,theta, Nsubs, intthree)
            call sciFour(n, d,beta,theta, Nsubs, intfour)
            Q = Q + current(n - 1) * intthree + (current(n) - current(n-1))/d * intfour
        end do

    end subroutine QCompute
    
    subroutine radarCrossSection(Nsubs,d,Q,E0,theta, mu, omega,pi, RCS)
    ! A subroutine for computing radar cross section from the hopefully previously computed Q quantity
    ! Params:
    ! Nsubs - the integer number of nodes in mesh
    ! d - real valued width of mesh subdivision
    ! Q - the sum of current contributions computed in the above subroutine
    ! E0 - the incident electric field!
    ! theta - the incident electric fields angle of attack
    ! mu - the permeability of (in this case) free space
    ! omega -  the angular frequency of the incident wave
    ! pi - the irrational real valued number pi

    ! Outputs:
    ! RCS - radar cross section - monostatic.

        integer, intent(in) :: Nsubs
        real(8), intent(in) :: d

        complex(8), intent(in) :: Q
        real(8), intent(in) :: E0
        real(8), intent(in) :: theta
        real(8), intent(in) :: mu
        real(8), intent(in) :: omega
        real(8), intent(in) :: pi
        complex(8), intent(out) ::RCS


        RCS = (omega * mu * sin(theta)* ABS(Q)) ** 2 / (4 * pi * E0 **2)

    end subroutine radarCrossSection
        
end module subroutines