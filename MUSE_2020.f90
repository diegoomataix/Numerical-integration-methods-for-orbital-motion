module MUSE_2020

        use dislin 
        use Linear_systems
        use Non_Linear_Systems
        use Cauchy_Problem
        use Temporal_Schemes
        use Temporal_scheme_interface 
        use Stability_regions
        !use Stability
        use Temporal_error
        use API_Example_Cauchy_Problem
        use plots 
        
        implicit none 

       real, parameter :: PI = 4 * atan(1d0)
       real            :: mu = 0.0121505856
       
    contains  

!**************************************************************************************************
!subroutine for user to choose the Milestone using the console
!**************************************************************************************************
    
subroutine Orbits_and_Numerical_Methods

    integer :: option

     write(*,*) " Welcome to the milestone menu "
     write(*,*) " Please select an option " 
     write(*,*) " 0. Exit menu"
     write(*,*) " 2. Milestone 2: Kepler Orbit Integrated with a Numerical Scheme"
     write(*,*) " 3. Milestone 3: Stability of a Kepler Orbit"
     write(*,*) " 4. Milestone 4: Error of a Kepler Orbit Integrated with a Numerical Scheme" 
     write(*,*) " 5. Milestone 5: Circular Restricted Three Body Problem"
     write(*,*) " 6. Milestone 6: Circular Restricted Three Body Problem"
     
    read(*,*) option
    
    
select case(option)
    
    case(2)
        call Milestone2
    
    case(3)
        call Milestone3
        
    case(4)
        call Milestone4
        
    case(5)
        call Milestone5
        
    case(6)
        call Milestone6
        
    case default
        write(*,*) "Not implemented"
    
    end select
    
end subroutine

!**************************************************************************************************
! subroutines for linking all the "sub-subroutines" to each Milestone n.
!**************************************************************************************************

! MILESTONE 2
subroutine Milestone2

    integer :: N

    integer :: option
    
     write(*,*) " Please select an option " 
     write(*,*) " 0. Exit menu"
     write(*,*) " 1. Euler scheme"
     write(*,*) " 2. Inverse Euler scheme"
     write(*,*) " 3. Crank-Nicolson scheme" 
     write(*,*) " 4. Runge-Kutta (4) scheme"
     
    read(*,*) option
    
    select case(option)
    
    case(0)
        call Orbits_and_Numerical_Methods
    
    case(1)
        write(*,*) "Euler scheme"
        write(*,*) " Select number of steps (reference:5000)"
        read(*,*) N  
        call Milestone_2_E(N)
        call Orbits_and_Numerical_Methods
    
    case(2)
        write(*,*) "Inverse Euler scheme"
        write(*,*) " Select number of steps (reference:3000)"
        read(*,*) N  
        call Milestone_2_IE(N)
        call Orbits_and_Numerical_Methods
        
    case(3)
        write(*,*) "Crank Nicolson scheme"
        write(*,*) " Select number of steps (reference:160)"
        read(*,*) N  
        call Milestone_2_CN(N)
        call Orbits_and_Numerical_Methods
        
    case(4)
        write(*,*) "Runge Kutta (4) scheme"
        write(*,*) " Select number of steps (reference:150)"
        read(*,*) N  
        call Milestone_2_RK4(N)
        call Orbits_and_Numerical_Methods
        
    case default
        write(*,*) "Not implemented"
    
    end select
    

    ! Initial time-step
    !call Milestone_2_E(N=5000)
    !call Milestone_2_IE(N=3000)
    !call Milestone_2_CN(N=160)
    !call Milestone_2_RK4(N=150)
    
    !call Kepler_orbit
    
end subroutine

! MILESTONE 3 
subroutine Milestone3

     integer :: option

     write(*,*) " Please select an option " 
     write(*,*) " 0. Exit menu"
     write(*,*) " 1. Euler scheme"
     write(*,*) " 2. Inverse Euler scheme"
     write(*,*) " 3. Crank-Nicolson scheme" 
     write(*,*) " 4. Runge-Kutta (4) scheme"
     write(*,*) " 5. GBS scheme"
     write(*,*) " 6. Convergence rates" 
     write(*,*) " 7. Variable step simulation"
     
    read(*,*) option
    
    
select case(option)
    
    case(0)
        call Orbits_and_Numerical_Methods
    
    case(1)
        write(*,*) "Euler scheme"
        call Kepler_orbit_polar(Scheme = Euler)
        call Error_Kepler_orbit(Scheme = Euler)
        call Phase_diagram_Kepler_orbit_polar(Scheme = Euler)
        call Orbits_and_Numerical_Methods
        
    case(2)
        write(*,*) "Inverse Euler scheme"
        call Kepler_orbit_polar(Scheme = Inverse_Euler)
        call Error_Kepler_orbit(Scheme = Inverse_Euler)
        call Phase_diagram_Kepler_orbit_polar(Scheme = Inverse_Euler)
        call Orbits_and_Numerical_Methods
        
    case(3)
        write(*,*) "Crank Nicolson scheme"
        call Kepler_orbit_polar(Scheme = Crank_Nicolson)
        call Error_Kepler_orbit(Scheme = Crank_Nicolson)
        call Phase_diagram_Kepler_orbit_polar(Scheme = Crank_Nicolson)
        call Orbits_and_Numerical_Methods
        
    case(4)
        write(*,*) "Runge Kutta (4) scheme"
        call Kepler_orbit_polar(Scheme = Runge_Kutta4)
        call Error_Kepler_orbit(Scheme = Runge_Kutta4)
        call Phase_diagram_Kepler_orbit_polar(Scheme = Runge_Kutta4)
        call Orbits_and_Numerical_Methods
        
    case(5)
        write(*,*) "GBS scheme"
        call Kepler_orbit_polar(Scheme = GBS_scheme)
        call Error_Kepler_orbit(Scheme = GBS_scheme)
        call Phase_diagram_Kepler_orbit_polar(Scheme = GBS_scheme)
        call Orbits_and_Numerical_Methods
    
    case(6)
         call Convergence_rate_Kepler
    
    case(7)
        call Variable_step_simulation
        
    case(8)
        call Error_Kepler_orbit(Scheme = Adams_Bashforth2)
        call Error_Kepler_orbit(Scheme = ERK_scheme)
        call Error_Kepler_orbit(Scheme = GBS_scheme)
        call Error_Kepler_orbit(Scheme = PC_ABM)
        
    case default
        write(*,*) "Not implemented"
    
    end select
    
    
end subroutine

! MILESTONE 4
subroutine Milestone4

    integer :: N(3)     = [ 1000, 400, 200 ]
    real    :: U0(4)    = [1., 0., 0., 0.7 ]
    real    :: U1(2)    = [1., 0.]
    
    integer :: option

     write(*,*) " Please select an option " 
     write(*,*) " 0. Exit menu"
     write(*,*) " 2. Eigenvalues_Kepler_polar"
     write(*,*) " 3. Stability_region_examples"
     write(*,*) " 4. Orbits_and_schemes" 
     write(*,*) " 5. Phase_space_Kepler_orbit"
     write(*,*) " 6. Phase_space_Kepler_orbit_polar"
     
    read(*,*) option
    
    
select case(option)
    
    case(2)
        call Eigenvalues_Kepler_polar  
    
    case(3)
        call Stability_region_examples
        
    case(4)
        write (*, *) ' Orbits_and_schemes with Oscillator function using Euler scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( Oscillator, Euler, 2, U1, N )
    
        write (*, *) ' Orbits_and_schemes with Oscillator function using Runge_kutta4 scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( Oscillator, Runge_kutta4, 2, U1, [20, 13, 100] )
    
        write (*, *) ' Orbits_and_schemes with Oscillator function using GBS scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( Oscillator, GBS_scheme, 2, U1, [20, 13, 100] )
    
        write (*, *) ' Orbits_and_schemes with F_Kepler function using Euler scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( F_Kepler, Euler,        4, U0, N )
    
        write (*, *) ' Orbits_and_schemes with F_Kepler function using Runge_kutta4 scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( F_Kepler, Runge_kutta4, 4, U0, N )
    
        write (*, *) ' Orbits_and_schemes with F_Kepler function using GBS scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( F_Kepler, GBS_scheme,        4, U0, N ) 
        
        write (*, *) ' Orbits_and_schemes with F_Kepler function using GBS scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( F_Kepler, ERK_scheme,        4, U0, N ) 
        
        write (*, *) ' Orbits_and_schemes with F_Kepler function using GBS scheme '
        write(*,*)  "press enter "; read(*,*)
        call Orbits_and_schemes_( F_Kepler, Adams_Bashforth2,        4, U0, N ) 
        
    case(5)
        call Phase_space_Kepler_orbit
        
    case(6)
        call Phase_space_Kepler_orbit_polar
        
    case default
        write(*,*) "Not implemented"
    
    end select
    
    
end subroutine

! MILESTONE 5
subroutine Milestone5
    call CR3BP_Lagrange_points_and_stability
     
end subroutine

! MILESTONE 6
subroutine Milestone6 
    
integer :: option

     write(*,*) " Welcome to the milestone menu "
     write(*,*) " Please select an option " 
     write(*,*) " 0. Exit menu"
     write(*,*) " 1. Arenstorf_orbit (Comparing Embedded RK: WDOPRI5 vs DOPRI54)"
     write(*,*) " 2. Arenstorf_orbit (Effect of tolerance: Embedded RK: WDOPRI5, DOPRI54)"
     write(*,*) " 3. Arenstorf_orbit (GBS and Embedded GBS)"
     write(*,*) " 4. Arenstorf_orbit (Effect of tolerance: wGBS, GBS)"
     write(*,*) " 5. Arenstorf_orbit (ABM and Embedded ABM)"
     write(*,*) " 6. Arenstorf_orbit (Effect of tolerance: wABM, ABM)"
     write(*,*) " 7. Shampine_Gordon_orbit_RK4"
     write(*,*) " 8. Shampine_Gordon_orbit_GBS"
     write(*,*) " 9. Shampine_Gordon_orbit_(Effect of tolerance: RK4, GBS)"
     write(*,*) " 10. Lyapunov_orbots_L1_L2_L3" 
     write(*,*) " 11. HenonHeiles_system (GBS)"
     write(*,*) " 12. Computational_effort"
     write(*,*) " 13. Computational_effort_comparison"
     
    read(*,*) option
    
    
select case(option)
    
    case(1)
        call Arenstorf_orbit_eRK
    
    case(2)
        call Arenstorf_orbit__eRK_tolerances
    
    case(3)
        call Arenstorf_orbit_GBS
    
    case(4)
        call Arenstorf_orbit_GBS_tolerances
        
    case(5)   
        call Arenstorf_orbit_ABM
        
    case(6)
        call Arenstorf_orbit_ABM_tolerances
    
    case(7)
        call Shampine_Gordon_orbit_RK4
    
    case(8)
        call Shampine_Gordon_orbit_GBS
        
    case(9)   
        call Shampine_Gordon_orbit_RK4_GBS_tolerance
    
    case(10)
        call Lyapunov_orbots_L1_L2_L3
        
    case(11)
        call HenonHeiles_system
        
    case(12)
        call Computational_effort
        
    case(13)
        call Computational_effort_comparison
        
    case default
        write(*,*) "Not implemented"
    
    end select
    
    
end subroutine
    
!**************************************************************************************************
! ******* USEFUL FUNCTIONS & SUBROUTINES *******
!**************************************************************************************************
! system matrix function

function System_matrix_1( F, U0, t ) result(A)
            procedure (ODES) :: F
            real, intent(in) :: U0(:), t
            real :: A( size(U0), size(U0) )
            
        real :: delta( size(U0) ), eps = 1d-6
        integer :: j, N
        
        N = size(U0)
        
        do j=1, N
            
            delta = 0
            delta(j) = eps
            A(:, j) = ( F( U0 + delta, t ) - F( U0 - delta, t ) ) / (2*eps)
            
        end do
        
end function


! Periodic_U function

function Periodic_U(N, t) result(U)

integer, intent(in) :: N
    real, intent(in)    :: t
    real :: U( N )
    
        real :: r(2), drdt(2)
        
        r       = [  cos(t), sin(t) ]
        drdt    = [ -sin(t), sin(t) ]
        
        U = [ r, drdt ]
        
end function

! Oscillator Function

function Oscillator_(U, t) result(F)

    real    :: U(:),t
    real    :: F(size(U))
    
    F = [ U(2), - U(1) ] 
    
end function


! CR3BP function

function CR3BP(U, t) result(F)
    real    :: U(:), t
    real    :: F ( size(U) )
    
        real :: x, y, z, vx, vy, vz, dvxdt, dvydt, dvzdt
        real :: d, r
        
        x   = U(1);  y = U(2);  z = U(3);
        vx  = U(4); vy = U(5); vz = U(6);
        
        d = sqrt( (x+mu)**2 + y**2 + z**2 )
        r = sqrt( (x-1+mu)**2 + y**2 + z**2 )
        
        dvxdt = x + 2 * vy - (1-mu) * ( x + mu )/d**3 - mu * (x-1+mu)/r**3
        dvydt = y - 2 * vx - (1-mu) * y / d**3 - mu * y/r**3
        dvzdt = - (1-mu)*z/d**3 - mu*z/r**3
        
        F = [ vx, vy, vz, dvxdt, dvydt, dvzdt ]
        
end function

! CR3BP 2D function
function CR3BP_2D(U, t) result(F)
    real    :: U(:), t
    real    :: F ( size(U) )
    
        real :: x, y, vx, vy, dxdt, dydt, dvxdt, dvydt
        real :: D1, D2
        
        x   = U(1);  y = U(2);  vx = U(3); vy = U(4)
        
        D1 = sqrt( (x+mu)**2    + y**2 )**3
        D2 = sqrt( (x-(1-mu))**2  + y**2 )**3
        
        dxdt = vx
        dydt = vy
        dvxdt = x + 2 * vy - (1-mu) * ( x + mu )/D1 - mu * (x-(1-mu))/D2
        dvydt = y - 2 * vx - (1-mu) * y /D1 - mu * y/D2
        
        F = [ dxdt, dydt, dvxdt, dvydt ]
        
end function

!**************************************************************************************************
! ******* SUBROUTINES MILESTONE 2 *******
!**************************************************************************************************

!**************************************************************************************************
! Kepler Orbit          EULER
!**************************************************************************************************
subroutine Milestone_2_E(N)
 
    integer, intent(in) :: N
    real :: Time(0:N), U(0:N,4)
    real :: a = 10., b = 28., c=2.66666666
    real :: t0 = 0, tf=6*PI
    integer :: i 
   
        Time = [ (t0 + (tf-t0)*i/real(N), i=0, N) ]
   
        U(0,:) = [1, 0, 0, 1] 
   
   call Cauchy_ProblemS(Time_Domain=Time, Differential_operator=F_Kepler, &
       Solution = U, Scheme = Euler)
   
   write(*, *)'Kepler orbit (Euler Method)'
   write(*,*) "press enter" ; read(*,*)
   call plot_parametrics(U(:,1),U(:,2:2), ["Kepler orbit"], "x", "y")
   call savedata(Filename = "./Milestone_2_E", x = U(:,1:2), y = U(:,2:2),N=1)
   
    contains
   
    function F_Kepler(U, t) result(F)
        real :: U(:),t
        real :: F(size(U))
    
        real :: r(2), drdt(2)
    
        r = U(1:2); drdt = U(3:4);
    
        F = [ drdt, -r / norm2(r)**3 ]
    
    end function 

end subroutine

!**************************************************************************************************
! Kepler Orbit      INVERSE EULER
!*****************************************************************************************
subroutine Milestone_2_IE(N)
 
    integer, intent(in) :: N
    real :: Time(0:N), U(0:N,4)
    real :: a = 10., b = 28., c=2.66666666
    real :: t0 = 0, tf=6*PI
    integer :: i 
   
        Time = [ (t0 + (tf-t0)*i/real(N), i=0, N) ]
   
        U(0,:) = [1, 0, 0, 1] 
   
   call Cauchy_ProblemS(Time_Domain=Time, Differential_operator=F_Kepler, &
       Solution = U, Scheme = Inverse_Euler)
   
   write(*, *)'Kepler orbit (Inverse Euler Method)'
   write(*,*) "press enter" ; read(*,*)
   call plot_parametrics(U(:,1),U(:,2:2), ["Kepler orbit"], "x", "y")
   call savedata(Filename = "./Milestone_2_IE", x = U(:,1:2), y = U(:,2:2), N = 1)
   
    end subroutine
    

!**************************************************************************************************
! Kepler Orbit      CRANK NICOLSON
!*****************************************************************************************
subroutine Milestone_2_CN(N)
 
    integer, intent(in) :: N
    real :: Time(0:N), U(0:N,4)
    real :: a = 10., b = 28., c=2.66666666
    real :: t0 = 0, tf=6*PI
    integer :: i 
   
        Time = [ (t0 + (tf-t0)*i/real(N), i=0, N) ]
   
        U(0,:) = [1, 0, 0, 1] 
   
   call Cauchy_ProblemS(Time_Domain=Time, Differential_operator=F_Kepler, &
       Solution = U, Scheme = Crank_Nicolson)
   
   write(*, *)'Kepler orbit (Crank-Nicolson Method)'
   write(*,*) "press enter" ; read(*,*)
   call plot_parametrics(U(:,1),U(:,2:2), ["Kepler orbit"], "x", "y")
   call savedata(Filename = "./Milestone_2_CN", x = U(:,1:2), y = U(:,2:2), N = 1)

    end subroutine
    
!**************************************************************************************************
! Kepler Orbit      RK4
!*****************************************************************************************
subroutine Milestone_2_RK4(N)
 
    integer, intent(in) :: N
    real :: Time(0:N), U(0:N,4)
    real :: a = 10., b = 28., c=2.66666666
    real :: t0 = 0, tf=6*PI
    integer :: i 
   
        Time = [ (t0 + (tf-t0)*i/real(N), i=0, N) ]
   
        U(0,:) = [1, 0, 0, 1] 
   
   call Cauchy_ProblemS(Time_Domain=Time, Differential_operator=F_Kepler, &
       Solution = U, Scheme = Runge_Kutta4)
   
   write(*, *)'Kepler orbit (RK4 Method)'
   write(*,*) "press enter" ; read(*,*)
   call plot_parametrics(U(:,1),U(:,2:2), ["Kepler orbit"], "x", "y")
   call savedata(Filename = "./Milestone_2_RK4", x = U(:,1:2), y = U(:,2:2), N = 1)

    end subroutine
        

!**************************************************************************************************
! Kepler_orbit
!**************************************************************************************************

subroutine Kepler_orbit
    integer, parameter :: N = 20000
    real    :: Time(0:N), U(0:N,4)
    real    :: t0 = 0, tf = 48*PI
    integer :: i

    integer, parameter :: Np = 10
    real :: x(0:N,0:Np), y(0:N,0:Np)
    real :: dx(0:N,0:Np), dy(0:N,0:Np)
    real :: r(0:N,0:Np), v(0:N,0:Np)
    real :: v0(0:Np)
    real :: vi = -0.3, vf = 0.3
    character(len=10) :: legends(0:Np)
    
    Time = [ (t0 + (tf - t0 ) * i/real(N), i=0, N ) ]
    v0   = [ (vi + (vf - vi ) * i/real(Np), i=0, Np ) ]
    legends = "v0"
    
    do i=0, Np
        U(0,:) = [1., 0., v0(i), 1.+v0(i)]
        
        call Cauchy_Problems ( Time_Domain=Time, Differential_operator = F_Kepler, &
                                !Solution = U, Scheme = Euler )
                                Solution = U, Scheme = Runge_kutta4 )
        
        x(:, i)  = U(:,1)
        y(:, i)  = U(:,2)
        dx(:, i) = U(:,3)
        dy(:, i) = U(:,4)
        
        r(:,i)   = sqrt( x(:,i)**2 + y(:,1)**2 )
        v(:,i)   = ( x(:,i)*dx(:,i) + y(:,i)*dy(:,i) )/r(:,i)
    end do
    
        write (*, *) 'Kepler orbit  '
        write(*,*)  "press enter "; read(*,*)
        call plot_parametrics(x(:,0:Np), y(:,0:Np), legends, "x", "y")
        
        write (*, *) 'Phase diagram: Kepler orbits  '
        write(*,*)  "press enter "; read(*,*)
        call plot_parametrics(r(:,0:Np), v(:,0:Np), legends, "r", "v") 
        
end subroutine
    
!**************************************************************************************************
! ******* SUBROUTINES MILESTONE 3 *******
!**************************************************************************************************

!**************************************************************************************************
! KEPLER ORBIT POLAR
!**************************************************************************************************

subroutine Kepler_orbit_polar(Scheme)

    procedure (Temporal_Scheme), optional :: Scheme 
    integer, parameter  :: N=1000
    real :: Time(0:N), U(0:N,2)
    real :: t0 = 0, tf = 12 * PI
    
    integer :: i
    
    Time = [ (t0 + (tf-t0 ) * i/real(N), i=0, N ) ]
    
     U(0,:) = [1., 0.01]
     
     
     call Cauchy_Problems( Time_Domain=Time,                        &
                           Differential_operator = F_Kepler_polar,  &
                           Solution = U, Scheme= Scheme)
     
     write(*, *) ' Polar coordinates: Kepler orbit '
     write(*,*) "press enter "; read(*,*)
     
     !call plot_parametrics(U(:,1), U(:,2:2), ["r"], "time", "r")
     call plot_parametrics(U(:,1), U(:,2:2), ["r"], "r", "v")
         
end subroutine

!**************** F_Kepler_polar function ****************
function F_Kepler_polar(U, t) result(F)
    real :: U(:), T
    real :: F(size(U))
    
    real :: r, drdt
    
    r = U(1); drdt = U(2);
    
    F = [ drdt, 1/r**2 * ( 1/r - 1 ) ]
    
end function

!**************************************************************************************************
! Error estimation by means of Richardson extrapolation
!**************************************************************************************************
subroutine Error_Kepler_orbit(Scheme)

    procedure (Temporal_Scheme), optional :: Scheme 
    real :: t0 = 0, tf = 16     ! Time domain for Kepler problem I Time steps of the first grid
    integer, parameter :: N = 1600  ! Time steps of the first grid
    integer, parameter :: Nv = 4    ! Number of variables
    real :: Time (0:N)              ! Time domain
    real :: Solution(0:N, Nv)       ! Solution
    real :: Error(0:N, Nv)          ! Error
    integer :: i 
    real :: x_aux(0:N,1), y_aux(0:N,1)

    x_aux(:,1) = Time(:)
    !y_aux(:,1) = Error(:, 4)
    y_aux(:,1) = sqrt(Error(:,1)**2 + Error(:,2)**2)
     
    write(*,*) "Error determination by Richardson extrapolation" 
    write(*,*) "press enter "; read(*,*)
    read(*,*)
    
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N) ] 
    Solution(0,:) =  [1, 0, 0, 1] !Initial conditions for the Kepler equation
    
    call Error_Cauchy_Problem( Time, F_Kepler,      &
                     Scheme, 1, Solution, Error ) 
    
    !call plot_parametrics( Solution(:,1), Solution(:,2:2), ["Kepler orbit"], "$ t $ ", " $x$ " ) 
    call plot_parametrics( Time, sqrt(Error(:,1:1)**2 + Error(:,2:2)**2), ["Kepler error"], "$ t $ ", "$ E $")
    call savedata(Filename = "./Error_Kepler_orbit(Scheme)", x = x_aux(:,1:1), y = y_aux(:, 1:1), N = 1)

end subroutine 

!**************************************************************************************************
! Convergence rate of numerical methods with the number of time steps 
!**************************************************************************************************
subroutine Convergence_rate_Kepler

    real :: t0 = 0, tf = 30         ! Time domain
    integer, parameter :: N = 500   ! Time steps of the first grid
    integer, parameter :: Nv = 4    ! Number of variables  
    integer, parameter :: M = 10    ! number of time grids
    real :: U(0:N,Nv)               ! Solution
    real :: Time (0:N), U0(Nv)      ! Time domain and initial conditions
    integer, parameter :: Np = 6    ! number of numerical schemes to be evaluated
    real :: log_E(M, Np), log_N(M)  ! Error versus time steps
                                    ! first index: grid, second index: scheme 
    
    integer :: i 
    character (len=20) :: names(Np) = [ "Euler", "Inverse_Euler", "Crank_Nicolson", "RK4", "Adams_Bashforth2", "GBS"]
    
    write(*,*) "Convergence rate: Error versus number of time steps" 
    write(*,*) "Temporal scheme : Euler, Inverse Euler, Crank Nicolson, Rk4, Adams_Bashforth2 y GBS " 
    write(*,*) "Press enter " ; read(*,*)
    
    Time = [ (t0+ (tf - t0 ) * i / real(N), i=0, N ) ] 
    U0 = [1,0,0,1] !Initial conditions 
    U(0,:) = U0
    
    call Cauchy_Problems( Time, F_Kepler, U, Euler ) 
    
    call Temporal_convergence_rate( Time, F_Kepler, U0,         &
                                    Euler, 2, log_E(:,1), log_N ) 
    
    call Temporal_convergence_rate( Time, F_Kepler, U0,         &
                                    Inverse_Euler, 2, log_E(:,2), log_N ) 
    
    call Temporal_convergence_rate( Time, F_Kepler, U0,         &
                                    Crank_Nicolson, 4, log_E(:,3), log_N ) 
    
    call Temporal_convergence_rate( Time, F_Kepler, U0,         &
                                    Runge_Kutta4, 4, log_E(:,4), log_N ) 
    
    call Temporal_convergence_rate( Time, F_Kepler, U0,         &
                                    Adams_Bashforth2, 4, log_E(:,5), log_N )
    
    call Temporal_convergence_rate( Time, F_Kepler, U0,         &
                                    GBS_scheme, 4, log_E(:,6), log_N )      
    
    !write(*, *) "Kepler orbit solution (Euler)"
    !write(*,*) "press enter"; read(*,*)
    !call plot_parametrics( U(:,1), U(:,2:2), [" "], "$x$", "$y$")
    
    !write(*, *) "Error versus time steps"
    !write(*,*) "press enter"; read(*,*)
    call plot_parametrics( log_N, log_E, names, "$\log N $", "$\log E$") 

end subroutine 

!**************** F_Kepler function ****************
function F_Kepler(U, t) result(F)
    real :: U(:),t
    real :: F(size(U))
    
    real :: r(2), drdt(2)
    
    r = U(1:2); drdt = U(3:4);
    
    F = [ drdt, -r / norm2(r)**3 ]
    
end function


!**************************************************************************************************
! Kepler_orbit in polar coordinates
!**************************************************************************************************
subroutine Phase_diagram_Kepler_orbit_polar(Scheme)

    procedure (Temporal_Scheme), optional :: Scheme 
    integer, parameter :: N = 1000
    real :: Time(0:N), U(0:N,2)
    real :: t0 = 0, tf = 12* PI
    
    integer, parameter :: Np = 10
    real :: r(0:N,0:Np), v(0:N,0:Np)
    real :: v0(0:Np)
    real :: vi = -0.8, vf = 0.8
    character(len=10) :: legends(0:Np)
    
    integer :: i
    
    Time = [ (t0 + (tf - t0) * i/real(N), i=0, N ) ]
    v0   = [ (vi + (vf - vi) * i/real(Np), i=0, Np ) ]
    legends = "v0"
    write(*,*) "v0 = ", v0
    
    do i=0, Np
        U=0
        U(0,:) = [1., v0(i)]
        
        call Cauchy_Problems( Time_Domain=Time,                         &
                              Differential_operator = F_Kepler_polar,    &
                              Solution = U, Scheme = Scheme)
        
        r(:, i) = U(:,1)
        v(:, i) = U(:,2)
    end do
    
        write(*, *) 'Polar coordinates: Kepler orbit    '
        write(*,*) "press enter "; read(*,*)
        
        call plot_parametrics(r(:,0:Np), v(:,0:Np), legends, "r", "drdt") 

end subroutine


!**************************************************************************************************
! ******* SUBROUTINES MILESTONE 4 *******
!**************************************************************************************************



!**************************************************************************************************
!   Region of absolute stability (examples)
!**************************************************************************************************

subroutine Stability_region_examples
    integer, parameter :: N = 50
    integer :: i, j
    real :: x(0:N), y(0:N), Region(0:N,0:N)
    !real :: x0 = -4d0, xf = 4d0, dx
    !real :: y0 = -4d0, yf = 4d0, dy
    real :: x0 = -10d0, xf = 3d0, dx    ! USE FOR GBS
    real :: y0 = -10d0, yf = 10d0, dy   ! USE FOR GBS

    real :: x_aux(0:N,1), y_aux(0:N,1)
    
    integer, parameter :: N_levels = 9
    real :: levels(0:N_levels)
    
    dx = (xf-x0) / N; dy = (yf-y0) / N
    x = [ ( x0 + dx * i , i=0,N ) ]
    y = [ ( y0 + dy * j , j=0,N ) ]
    
    x_aux(:,1) = x(:)
    y_aux(:,1) = y(:)
    
    write(*, *) "Region of absolute stability: Euler"
    write(*,*) "press enter"; read(*,*)
    
    levels = [ ( j/real(N_levels) , j=0, N_levels )]
    
    call Absolute_Stability_Region(Euler, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_E',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    
    
    write(*, *) "Region of absolute stability: Inverse_Euler"
    write(*,*) "press enter"; read(*,*)
    
    call Absolute_Stability_Region(Inverse_Euler, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_IE',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    
    
    write(*, *) "Region of absolute stability: Runge_Kutta4"
    write(*,*) "press enter"; read(*,*)
    
    call Absolute_Stability_Region(Runge_Kutta4, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_RK4',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    
    
    write(*, *) "Region of absolute stability: Crank_Nicolson"
    write(*,*) "press enter"; read(*,*)
    
    call Absolute_Stability_Region(Crank_Nicolson, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_CN',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    !
    
    write(*, *) "Region of absolute stability: Adams_Bashforth2"
    write(*,*) "press enter"; read(*,*)
    
    call Absolute_Stability_Region(Adams_Bashforth2, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_AB2',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    !    
    write(*, *) "Region of absolute stability: GBS"
    write(*,*) "press enter"; read(*,*)
    
    call Absolute_Stability_Region(GBS_scheme, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_GBS',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    !
     write(*, *) "Region of absolute stability: ERK"
    write(*,*) "press enter"; read(*,*)
    
    call Absolute_Stability_Region(ERK_scheme, x, y, Region)
    
    call plot_contour(x, y, Region, "$\Re(z)$","\Im(z)$",   &
                      levels = levels, graph_type = "isolines")
    call savedata(Filename='./Milestone4_Stability_ERK',x=x_aux(:,1:1),y=y_aux(:,1:1),z=Region(:,:),N=1)
    
end subroutine

!**************************************************************************************************
!  Region of absolute stability (Eigenvalues_Kepler_polar)
!**************************************************************************************************

subroutine Eigenvalues_Kepler_polar

    integer, parameter :: N = 2
    real :: t, A(N, N), U(N)
    integer :: i
    complex :: lambda(N)
    
    U = [ 5, 2]
    t = 0
    
    A = System_matrix_1( F_Kepler_polar, U, t )
    
    do i=1, N
        write(*,'(A, 5f8.2)') "A =", A(i,:)
    end do
    
    call Eigenvalues_QR(A, lambda)
    
    do i=1, N
        write(*,'(A, 2f8.2)') "lambda =", lambda(i)
    end do
    
end subroutine
!**************************************************************************************************
!   Phase_space_Kepler_orbit
!**************************************************************************************************

subroutine Phase_space_Kepler_orbit
    
    integer, parameter  :: N = 10000
    real :: Time(0:N), U(0:N,4)
    real :: t0 = 0, tf = 48*PI
    integer :: i
    
    integer, parameter :: Np = 10
    real    :: x(0:N,0:Np), y(0:N,0:Np)
    real    :: dx(0:N,0:Np), dy(0:N,0:Np)
    real    :: r(0:N,0:Np), v(0:N,0:Np)
    real    :: v0(0:Np)
    real    :: vi = -0.3, vf = 0.3
    
  character(len=10) :: legends(0:Np)
  
    Time = [ (t0 + (tf - t0 ) * i/real(N), i=0, N ) ]
    v0   = [ (vi + (vf - vi ) * i/real(Np), i=0, Np ) ]
    legends = "v0"
    
    do i=0, Np
        U(0,:) = [1., 0., v0(i), 1.+v0(i)]
        
        call Cauchy_Problems( Time_Domain=Time,                         &
                              Differential_operator = F_Kepler,    &
                              Solution = U, Scheme = Runge_kutta4 )
                              !Solution = U, Scheme = Euler )
                              !Solution = U, Scheme = Crank_Nicolson )
        x(:, i) = U(:,1)
        y(:, i) = U(:,2)
        dx(:, i) = U(:,3)
        dy(:, i) = U(:,4)
        
        r(:,i) = sqrt( x(:,i)**2 + y(:,i)**2 )
        v(:,i) = ( x(:,i)*dx(:,i) + y(:,i)*dy(:,i) )/r(:,i)
    end do
    
    write (*, *) 'Kepler orbit  '
    write(*,*) "press enter "; read(*,*)
    call plot_parametrics(x(:,0:Np), y(:,0:Np), legends, "x", "y")
    
    write (*, *) 'Phase diagram: Kepler orbits  '
    write(*,*) "press enter "; read(*,*)
    
  call plot_parametrics(r(:,0:Np), v(:,0:Np), legends, "r", "v")
  
end subroutine
        

!**************************************************************************************************
!   Phase_space_Kepler_orbit_polar
!**************************************************************************************************

subroutine Phase_space_Kepler_orbit_polar

    integer, parameter  :: N = 1000
    real :: Time(0:N), U(0:N,2)
    real :: t0 = 0, tf = 12*PI
    
    integer, parameter  :: Np = 10
    real    :: r(0:N,0:Np), v(0:N,0:Np)
    real    :: v0(0:Np)
    real    :: vi = -0.8, vf = 0.8
    character(len=10)   :: legends(0:Np)
    
    integer :: i
    
    Time = [ (t0 + (tf - t0 ) * i/real(N), i=0, N ) ]
    v0   = [ (vi + (vf - vi ) * i/real(Np), i=0, Np ) ]
    legends = "v0"
    write(*,*) " v0 = ", v0
    
    do i=0, Np
        U = 0
        U(0,:) = [1., v0(i)]
        
        call Cauchy_Problems( Time_Domain=Time,                         &
                              Differential_operator = F_Kepler_polar,   &
                              Solution = U, Scheme = Runge_kutta4 )
                              !Solution = U, Scheme = Euler )
                              !Solution = U, Scheme = Crank_Nicolson )
        
        r(:, i) = U(:,1)
        v(:, i) = U(:,2)
    end do
    
    write(*, *) 'Polar coordinates: Kepler orbit    '
    write(*,*) "press enter " ; read(*,*)
    call plot_parametrics(r(:,0:Np), v(:,0:Np), legends, "r", "drdt")
    
end subroutine

!**************************************************************************************************
!   Orbits and Schemes 
!**************************************************************************************************
subroutine Orbits_and_schemes_(F, Scheme, Nv, U0, N )
    procedure (ODES)    :: F
    procedure (Temporal_Scheme), optional :: Scheme 
    integer, intent(in) :: Nv
    real, intent(in)    :: U0(:)
    integer, intent(in) :: N(:)
    
    real, allocatable   :: Time(:), U(:,:)
    real :: t0 = 0, tf=12*PI
    integer :: i, j
    integer :: Np
    
    Np = size(N)
    call ini_plot(xmax = 2.,  xmin = -2., ymax = 2., ymin = -2. )
    
    do i=1, Np
    
        allocate( Time(0:N(i)), U(0:N(i), Nv) )
        Time = [ (t0 + ( tf - t0 ) * j/real(N(i)), j=0, N(i) ) ]
        U(0,:) = U0(:)
        
        
         call Cauchy_Problems( Time, F, U, Scheme )
         
         call chncrv( 'color')
         call curve( U(:,1), U(:,2), N(i)+1 )
         
         call savedata(Filename = "./orb_n_scheme", x = U(:, 1:1), y = U(:, 1:2), N = 3)
         deallocate( Time, U )
         
    end do
    call disfin
    
end subroutine
    
subroutine ini_plot( xmax, xmin, ymax, ymin )
    real, intent(in)    :: xmax, xmin, ymax, ymin
    
    call metafl("xwin")
    call PAGE (4000, 4000)
    call scrmod("reverse")
    call disini
    call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10)
    
end subroutine


!**************************************************************************************************
! ******* SUBROUTINES MILESTONE 5 *******
!**************************************************************************************************

!**************************************************************************************************
!   CR3BP_Lagrange_points_and_stability
!**************************************************************************************************
subroutine CR3BP_Lagrange_points_and_stability
    
    integer, parameter  :: N = 20000 ! time steps
    integer, parameter  :: M = 6     ! number of variables
    real    :: U(0:N, M)             ! state vector
    real    :: E(0:N, M)             ! error for different time steps
    real    :: Time(0:N)             ! time domain
    integer, parameter  :: NL = 5    ! number of Lagrange points
    real    :: U0(NL, M)             ! Lagrange points
    real    :: eps(M)                ! disturbances around Lagrange points
    real    :: A(M,M)                ! Jacobian
    complex :: lambda(M)             ! Eigenvalues
    
    real    :: t0 = 0
    real    :: tf = 4*PI/0.3         ! final integration time
    integer :: i, j                  ! index
    
    
    Time = [ (t0 + (tf-t0)*i/N, i=0, N ) ]
    
    U0(1,:) = [ 0.8, 0.6, 0., 0., 0., 0. ]
    U0(2,:) = [ 0.8, -0.6, 0., 0., 0., 0. ]
    U0(3,:) = [ -0.1, 0.0, 0., 0., 0., 0. ]
    U0(4,:) = [ 0.1, 0.0, 0., 0., 0., 0. ]
    U0(5,:) = [ 1.1, 0.0, 0., 0., 0., 0. ]

    do i=1, NL
    ! Lagrange points L1, L2, L3, L4, L5
        call Newton( G, U0(i, 1:3) )
        
    ! Jacobian matrix at Lagrange points
        A = System_matrix_1( F = CR3BP, U0 = U0(i,:), t=0.)
        
    ! Stability of the Lagrange points
        call Eigenvalues_QR( A, lambda )
        A = System_matrix_1( F = CR3BP, U0 = U0(i,:), t=0.)
        
        write(*,'(A,6f12.4)') " Lagrange point =", U0(i,:)
        write(*,'(A,6E16.4)') " F in Lagrange point =", CR3BP( U0(i,:), 0.0 )
        
        do j=1, M
            write(*,*) " Eigenvalue =", lambda(j)
        end do
        read(*,*)
        
    ! orbit around the Lagrange point
        call RANDOM_NUMBER(eps)
        eps = 1d-2 * eps
        U(0,:) = U0(i,:) + eps
        
        call Cauchy_Problems( Time, CR3BP, U )
        
        call plot(U, mu, N)
        call qplot( U(:,1), U(:,2), N+1)
        
    ! Error of the orbit
        call Error_Cauchy_Problem( Time, CR3BP, Runge_kutta4, 4, U, E )
        call qplot( Time, E(:,1), N+1)
        
    end do
    
contains
!-------------------------------------------------------------------------
function G(Y)
    real, intent(in)    :: Y(:)
    real    :: G(size(Y))
    
        real    :: X(6), GX(6)
        X(1:3)  = Y
        X(4:6)  = 0.
        
        GX = CR3BP(X, 0.)
        G  = GX(4:6)
        
end function

end subroutine

!-------------------------------------------------------------------------
subroutine plot(U, mu, N)
    real, intent(in) :: U(:,:), mu
    integer, intent(in) :: N
    
    call metafl("xwin")
    CALL PAGE (4000, 4000)
    call scrmod("reverse")
    call disini
    call graf(-1.5, 1.5, -1.5, 0.5, -1.5, 1.5, -1.5, 0.5)
    call curve( U(:,1), U(:,2), N+1)
    
    call color("blue"); call marker(21); call incmrk(-1)
    call curve([-mu, -mu], [0.,0.], 2 )
    call color("red"); call marker(21); call incmrk(-1)
    call curve([1-mu, 1-mu], [0.,0.], 2 )
    
    call disfin
    
end subroutine

!**************************************************************************************************
! ******* SUBROUTINES MILESTONE 6 *******
!**************************************************************************************************

!**************************************************************************************************
!   Lyapunov_orbots_L1_L2_L3
!**************************************************************************************************
subroutine Lyapunov_orbots_L1_L2_L3

	real	:: t0 = 0
	real	:: tf = 10 ! Time domain for Arenstorf orbit
	integer, parameter  :: N = 2000 ! number of time steps
	integer, parameter  :: Nv = 6   ! number of variables
	real	:: Time(0:N), U0(Nv)	! time domain and initial conditions
	integer, parameter  :: Np = 3   ! number of graphs
	real	:: U(0:N, Nv, Np)   	! solutions
                                	! first index: time, second: componenet, third scheme
	real	:: Error(0:N, Np)   	! error associated to different schemes
	integer :: i, j             	! indexes
	character (len=40) :: names(Np) ! name of parametric curves
	character (len=40) :: family(Np)! name of parametric curves
	real	:: tolerance = 1d-8
    
 	!write(*,*) " Lyapunov orbits L1 L2 L3 "
 	!read(*,*)
	 
 	!U(0,:,1) = [0.6089, 0., 0., 0., 0.85246, 0. ]
 	!U(0,:,2) = [1.3220, 0., 0., 0., -0.60859, 0. ]
 	!U(0,:,3) = [-1.9118, 0., 0., 0., 1.70633, 0. ]
    
    U(0,:,1) = [0.5, 0., 0., 0., 0.8, 0. ]
 	U(0,:,2) = [1.2, 0., 0., 0., -0.5, 0. ]
 	U(0,:,3) = [-1.8, 0., 0., 0., 1.6, 0. ]
 	Time = [ (t0 + (tf-t0) * i / real(N), i=0, N ) ]
	 
 	call set_solver( "eRK", "RK87" )
 	call set_tolerance(tolerance)
 	do j=1, Np
     	call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    write(*,*) " Lyapunov orbits L1 L2 L3 [RK4]"
 	read(*,*)
 	call plot_parametrics( U(:, 1, :), U(:, 2, :), ["L1", "L2", "L3"], "$x$", "$y$" )
    call savedata(Filename = "./Lyapunov_orbots_L1_L2_L3_RK4", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
    !call plot_parametrics( U(:, 1, :), U(:, 2, :), ["L1"], "$x$", "$y$" )
    !call savedata(Filename = "./Lyapunov_orbots_L1_L2_L3_RK4", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
    
    call set_solver( "GBS", "GBS" )
 	call set_tolerance(tolerance)
 	do j=1, Np
     	call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    write(*,*) " Lyapunov orbits L1 L2 L3 [GBS]"
 	read(*,*)
 	call plot_parametrics( U(:, 1, :), U(:, 2, :), ["L1", "L2", "L3"], "$x$", "$y$" )
    call savedata(Filename = "./Lyapunov_orbots_L1_L2_L3_GBS", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
	 
end subroutine


!**************************************************************************************************
!  Arenstorf_orbit_eRK
!**************************************************************************************************

subroutine Arenstorf_orbit_eRK

    real :: t0 = 0
    real(kind=8) :: tf = 17.0652165601579625    ! time domain for arenstorf orbit
    integer, parameter  :: N  = 1500    ! number of time teps
    integer, parameter  :: Nv = 6       ! number of variables
    real :: Time(0:N), U0(Nv)           ! time domain and initial conditions
    !integer, parameter  :: Np = 3       ! number of graphs (USAR ESTE SI SE INCLUYE GBS EN SUBROUTINE)
    integer, parameter  :: Np = 2       ! number of graphs
    real :: U(0:N, Nv, Np)              ! solutions
                                        ! first index: time, second: component, third: scheme
    integer :: i, j                     ! indexes
    character (len=40)  :: names(Np)    ! name of parametric curves
    character (len=40)  :: family(Np)   ! name of numerical scheme famuly
    
    real    :: tolerance = 1d-4
    mu = 0.012277471
    
    
    do j=1, Np
        U(0,:, j) = [ 0.994, 0., 0., 0., -2.0015851063790825, 0. ]
    end do
    
    Time = [ (t0 + (tf - t0) * i / real(N), i=0, N) ]
    
    call Cauchy_Problems( Time, CR3BP, U(:, :, 1) )
    
    !write(*, *) "Arenstorf orbit with RK4 scheme"
    !write(*,*) "press enter"; read(*,*)
    !call plot_parametrics( U(:, 1, 1:1), U(:, 2, 1:1), ["RK4"], "$x$", "$y$" )
    
    !family = ["weRK", "eRK", "GBS" ]           ! (USAR ESTE SI SE INCLUYE GBS EN SUBROUTINE)
    !names  = ["WDOP853", "RK87", "GBS" ]       ! (USAR ESTE SI SE INCLUYE GBS EN SUBROUTINE)
    
    family = ["weRK", "eRK"]
    names  = ["WDOP853", "RK87"] 
    
    do j=1, Np
        call set_solver( family(j), names(j) )
        call set_tolerance(tolerance)
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
    !write(*, *) "Arenstorf orbit with WDOP853, RK87 and GBS schemes"    ! (USAR ESTE SI SE INCLUYE GBS EN SUBROUTINE)
    write(*, *) "Arenstorf orbit with WDOP853 and RK87 schemes"
    write(*,*) "press enter"; read(*,*)
    call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
    !call savedata(Filename = "./Arenstorf_orbit_eRK", x =  U(:, 1, :), y=  U(:, 2, :),N=Np)
    call savedata(Filename = "./Arenstorf_orbit_RK4", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine 

!**************************************************************************************************
!  Arenstorf_orbit_ Effect of tolerances    eRK
!**************************************************************************************************

subroutine Arenstorf_orbit__eRK_tolerances

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 1500   ! Time steps of the first grid 
    integer, parameter :: Nv = 4     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 

    character (len=40) :: names(Np) =    & 
                          [ "$\epsilon=10^{-3}$", "$\epsilon=10^{-4}$", "$\epsilon=10^{-5}$"]
    real :: tolerances(Np) = [ 1d-3, 1d-4, 1d-5 ]  
    integer :: No, i, j  
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/ArenstorfWRKa", &
                                        "./doc/chapters/Cauchy_problem/figures/ArenstorfRKb"  ] 
 
     do j=1, Np 
       U(0,:,j) = [0.994, 0., 0., -2.0015851063790825 ]
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ] 
     
     call set_solver(family_name="weRK", scheme_name="WDOPRI5")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     write(*, *) "Arenstorf orbit with WDOP853 scheme, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names,           &
                            "$x$", "$y$", "(a)", path(1) )
     call savedata(Filename = "./Arenstorf_orbit_RK4_tolerances_WDOPRI5", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
     
     call set_solver(family_name="eRK", scheme_name="DOPRI54")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do
     write(*, *) "Arenstorf orbit with DOPRI54 scheme, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(b)", path(2) )
     call savedata(Filename = "./Arenstorf_orbit_RK4_tolerances_DOPRI54", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
end subroutine

!**************************************************************************************************
!  Arenstorf_orbit_GBS
!**************************************************************************************************

subroutine Arenstorf_orbit_GBS

    real :: t0 = 0
    real(kind=8) :: tf = 17.0652165601579625    ! time domain for arenstorf orbit
    integer, parameter  :: N  = 1000    ! number of time teps
    integer, parameter  :: Nv = 6       ! number of variables
    real :: Time(0:N), U0(Nv)           ! time domain and initial conditions
    integer, parameter  :: Np = 2       ! number of graphs
    reaL :: U(0:N, Nv, Np)              ! solutions
                                        ! first index: time, second: component, third: scheme
    integer :: i, j                     ! indexes
    character (len=40)  :: names(Np)    ! name of parametric curves
    character (len=40)  :: family(Np)   ! name of numerical scheme famuly
    
    real    :: tolerance = 1d-4     ! estaría bien poder cambiar esto o como input o poner varias curvas en el mismo gráfico
    
    mu = 0.012277471
    
    
    do j=1, Np
        U(0,:, j) = [ 0.994, 0., 0., 0., -2.0015851063790825, 0. ]
    end do
    
    Time = [ (t0 + (tf - t0) * i / real(N), i=0, N) ]
    
    call Cauchy_Problems( Time, CR3BP, U(:, :, 1) )
    
    !write(*, *) "Arenstorf orbit with RK4 scheme"
    !write(*,*) "press enter"; read(*,*)
    !call plot_parametrics( U(:, 1, 1:1), U(:, 2, 1:1), ["RK4"], "$x$", "$y$" )
    
    family = ["GBS", "wGBS" ]
    names  = ["GBS", "wGBS" ] 
    
    do j=1, Np
        call set_solver( family(j), names(j) )
        call set_tolerance(tolerance)
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
    write(*, *) "Arenstorf orbit with GBS, wGBS"
    write(*,*) "press enter"; read(*,*)
    call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
    call savedata(Filename = "./Arenstorf_orbit_GBS", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine 

!**************************************************************************************************
!  Arenstorf_orbit_ Effect of tolerances GBS
!**************************************************************************************************

subroutine Arenstorf_orbit_GBS_tolerances

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 1500   ! Time steps of the first grid 
    integer, parameter :: Nv = 4     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 

    character (len=40) :: names(Np) =    & 
                          [ "$\epsilon=10^{-3}$", "$\epsilon=10^{-4}$", "$\epsilon=10^{-5}$"]
    real :: tolerances(Np) = [ 1d-3, 1d-4, 1d-5 ]  
    integer :: No, i, j  
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/ArenstorfWRKa", &
                                        "./doc/chapters/Cauchy_problem/figures/ArenstorfRKb"  ] 
 
     do j=1, Np 
       U(0,:,j) = [0.994, 0., 0., -2.0015851063790825 ]
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ] 
     
     call set_solver(family_name="GBS", scheme_name="GBS")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     write(*, *) "Arenstorf orbit with GBS scheme, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names,           &
                            "$x$", "$y$", "(a)", path(1) ) 
    call savedata(Filename = "./Arenstorf_orbit_GBS_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
     
     call set_solver(family_name="wGBS", scheme_name="wGBS")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do
     write(*, *) "Arenstorf orbit with wGBS scheme, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(b)", path(2) )
    call savedata(Filename = "./Arenstorf_orbit_wGBS_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine

!**************************************************************************************************
!  Arenstorf_orbit_ABM
!**************************************************************************************************

subroutine Arenstorf_orbit_ABM
    
    real :: t0 = 0
    real(kind=8) :: tf = 17.065216560157962!5    ! time domain for arenstorf orbit
    integer, parameter  :: N  = 1000    ! number of time teps
    integer, parameter  :: Nv = 6       ! number of variables
    real :: Time(0:N), U0(Nv)           ! time domain and initial conditions
    integer, parameter  :: Np = 1       ! number of graphs
    reaL :: U(0:N, Nv, Np)              ! solutions
                                        ! first index: time, second: component, third: scheme
    integer :: i, j                     ! indexes
    character (len=40)  :: names(Np)    ! name of parametric curves
    character (len=40)  :: family(Np)   ! name of numerical scheme famuly
    
    real    :: tolerance = 1d-4     !
    
    mu = 0.012277471
    
    do j=1, Np
        U(0,:, j) = [ 0.994, 0., 0., 0., -2.0015851063790825, 0. ]
    end do
    
    Time = [ (t0 + (tf - t0) * i / real(N), i=0, N) ]
    
    call Cauchy_Problems( Time, CR3BP, U(:, :, 1) )
    
    !family = ["ABM", "wABM" ]
    !names  = ["ABM", "wABM" ] 
    family = ["wABM" ]
    names  = ["wABM" ] 
    
    do j=1, Np
        call set_solver( family(j), names(j) )
        call set_tolerance(tolerance)
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
    write(*, *) "Arenstorf orbit with ABM, wABM"
    write(*,*) "press enter"; read(*,*)
    call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$")
    call savedata(Filename = "./Arenstorf_orbit_ABM", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine 

!**************************************************************************************************
!  Arenstorf_orbit_ Effect of tolerances ABM
!**************************************************************************************************

subroutine Arenstorf_orbit_ABM_tolerances

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 1500   ! Time steps of the first grid 
    integer, parameter :: Nv = 4     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 

    character (len=40) :: names(Np) =    & 
                          [ "$\epsilon=10^{-3}$", "$\epsilon=10^{-4}$", "$\epsilon=10^{-5}$"]
    real :: tolerances(Np) = [ 1d-3, 1d-4, 1d-5 ]  
    integer :: No, i, j  
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/ArenstorfWRKa", &
                                        "./doc/chapters/Cauchy_problem/figures/ArenstorfRKb"  ] 
 
     do j=1, Np 
       U(0,:,j) = [0.994, 0., 0., -2.0015851063790825 ]
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ] 
     
     call set_solver(family_name="ABM", scheme_name="ABM")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     !write(*, *) "Arenstorf orbit with ABM scheme, effect of chosen tolerance"
     !write(*,*) "press enter"; read(*,*)
     !call plot_parametrics( U(:, 1, :), U(:, 2, :), names,           &
     !                       "$x$", "$y$", "(a)", path(1) ) 
     !call savedata(Filename = "./Arenstorf_orbit_ABM_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
     
    
     call set_solver(family_name="wABM", scheme_name="wABM")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do
     write(*, *) "Arenstorf orbit with wABM scheme, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(b)", path(2) )
    call savedata(Filename = "./Arenstorf_orbit_wABM_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
     
end subroutine

!**************************************************************************************************
!  Shampine_Gordon_orbit_RK4
!**************************************************************************************************
subroutine Shampine_Gordon_orbit_RK4

    real :: t0 = 0
    real(kind=8) :: tf = 17.0652165601579625    ! time domain 
    integer, parameter  :: N  = 2000    ! number of time teps
    integer, parameter  :: Nv = 6       ! number of variables
    real :: Time(0:N), U0(Nv)           ! time domain and initial conditions
    integer, parameter  :: Np = 3       ! number of graphs
    reaL :: U(0:N, Nv, Np)              ! solutions
    real :: Error(0:N, Np)              ! error associated to different schemes 
                                        ! first index: time, second: component, third: scheme
    integer :: i, j                     ! indexes
    character (len=40)  :: names(Np)    ! name of parametric curves
    character (len=40)  :: family(Np)   ! name of numerical scheme famuly
    real    :: tolerance = 1d-5

    !write(*,*) "Shampine Gordon orbit "
    !read(*,*)
    
    do j=1, Np
        U(0,:,j) = [1.2, 0. ,0., 0., -1.049357509830, 0. ] !shampine orbit initial conditions
        !U(0,:, j) = [ 0.394, 0., 0., 0., -1.0015851063790825, 0. ]
    end do
    Time = [ (t0 + (tf - t0 ) * i / real(N), i=0, N ) ]
    
    call Cauchy_Problems( Time , CR3BP, U(:, :, 1) )
    call plot_parametrics( U(:, 1, 1:1), U(:, 2, 1:1), ["RK4"], "$x$", "$y$" )
    
    
    names  = ["WDOP853", "DOPRI54", "RK87" ] 
    family = ["weRK", "eRK", "eRK" ]
    
    do j=1, Np
        call set_solver( family(j), names(j) )
        call set_tolerance(tolerance)
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
    call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
    call savedata(Filename = "./Shampine_Gordon_orbit_RK4", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine 

!**************************************************************************************************
!  Shampine_Gordon_orbit_GBS
!**************************************************************************************************
subroutine Shampine_Gordon_orbit_GBS

    real :: t0 = 0
    real(kind=8) :: tf = 17.0652165601579625    ! time domain 
    integer, parameter  :: N  = 2000    ! number of time teps
    integer, parameter  :: Nv = 6       ! number of variables
    real :: Time(0:N), U0(Nv)           ! time domain and initial conditions
    integer, parameter  :: Np = 2       ! number of graphs
    reaL :: U(0:N, Nv, Np)              ! solutions
    real :: Error(0:N, Np)              ! error associated to different schemes 
                                        ! first index: time, second: component, third: scheme
    integer :: i, j                     ! indexes
    character (len=40)  :: names(Np)    ! name of parametric curves
    character (len=40)  :: family(Np)   ! name of numerical scheme famuly
    real    :: tolerance = 1d-5

    !write(*,*) "Shampine Gordon orbit "
    !read(*,*)
    
    do j=1, Np
        U(0,:,j) = [1.2, 0. ,0., 0., -1.049357509830, 0. ]
    end do
    Time = [ (t0 + (tf - t0 ) * i / real(N), i=0, N ) ]
    
    call Cauchy_Problems( Time , CR3BP, U(:, :, 1) )
    call plot_parametrics( U(:, 1, 1:1), U(:, 2, 1:1), ["GBS"], "$x$", "$y$" )
    
    
    names  = ["GBS", "wGBS" ] 
    family = ["GBS", "wGBS" ]
    
    do j=1, Np
        call set_solver( family(j), names(j) )
        call set_tolerance(tolerance)
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
    call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
   call savedata(Filename = "./Shampine_Gordon_orbit_GBS", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine 


!**************************************************************************************************
!  Shampine_Gordon_orbit_RK4(RK87) & GBS tolerance
!**************************************************************************************************
subroutine Shampine_Gordon_orbit_RK4_GBS_tolerance 
     
    real :: t0 = 0
    real(kind=8) :: tf = 17.0652165601579625    ! time domain 
    integer, parameter  :: N  = 2000    ! number of time teps
    integer, parameter  :: Nv = 6       ! number of variables
    real :: Time(0:N), U0(Nv)           ! time domain and initial conditions
    integer, parameter  :: Np = 3       ! number of graphs
    reaL :: U(0:N, Nv, Np)              ! solutions
    real :: Error(0:N, Np)              ! error associated to different schemes 
    
                                        ! first index: time, second: component, third: scheme
    integer :: i, j                     ! indexes
    
    character (len=40) :: names(Np) =    & 
                          [ "$\epsilon=10^{-3}$", "$\epsilon=10^{-4}$", "$\epsilon=10^{-5}$"]
    real :: tolerances(Np) = [ 1d-3, 1d-4, 1d-5 ]  

    !write(*,*) "Shampine Gordon orbit "
    !read(*,*)
    
    do j=1, Np
        U(0,:,j) = [1.2, 0. ,0., 0., -1.049357509830, 0. ]
    end do
    Time = [ (t0 + (tf - t0 ) * i / real(N), i=0, N ) ]
    
    
    call set_solver(family_name="eRK", scheme_name="RK87")
    do j=1, Np
        call set_tolerance(tolerances(j))
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
     write(*, *) "Shampine_Gordon orbit_RK87, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(a)")
     call savedata(Filename = "./Shampine_Gordon_orbit_RK87_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
     
     !
     
    call set_solver(family_name="weRK", scheme_name="WDOP853")
    do j=1, Np
        call set_tolerance(tolerances(j))
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
     write(*, *) "Shampine_Gordon orbit_WDOP853, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(b)")
     call savedata(Filename = "./Shampine_Gordon_orbit_WDOP853_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
  
!
     
    call set_solver(family_name="eRK", scheme_name="DOPRI54")
    do j=1, Np
        call set_tolerance(tolerances(j))
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
     write(*, *) "Shampine_Gordon orbit_DOPRI54, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(c)" )
     call savedata(Filename = "./Shampine_Gordon_orbit_DOPRI54_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
     
     !
     
     call set_solver(family_name="GBS", scheme_name="GBS")
    do j=1, Np
        call set_tolerance(tolerances(j))
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
     write(*, *) "Shampine_Gordon orbit_GBS, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(d)")
     call savedata(Filename = "./Shampine_Gordon_orbit_GBS_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
     
    !
     
     call set_solver(family_name="wGBS", scheme_name="wGBS")
    do j=1, Np
        call set_tolerance(tolerances(j))
        call Cauchy_Problems( Time, CR3BP, U(:, :, j) )
    end do
    
     write(*, *) "Shampine_Gordon orbit_wGBS, effect of chosen tolerance"
     write(*,*) "press enter"; read(*,*)
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(e)")
     call savedata(Filename = "./Shampine_Gordon_orbit_wGBS_tolerances", x =  U(:, 1, :), y= U(:, 2, :), N=Np)
    
end subroutine 

!**************************************************************************************************
! ******* COMPUTATIONAL EFFORT OF RUNGE-KUTTA SCHEME *******
!**************************************************************************************************

subroutine Computational_effort

    real :: t0 = 0, tf = 30                 ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 300           ! Time steps of the first grid 
    integer, parameter :: Nv = 2            ! Number of variables van de pol
    integer, parameter :: M = 9             ! number of different tolerances           
    real :: Time (0:N)                      ! Time domain 
    real :: U0(Nv)                          ! Initial conditions  
    real :: U(0:N,Nv)                       ! Solution 
    integer, parameter :: Np = 10            ! number of numerical schemes to be evaluated 
    real :: log_mu(M), log_effort(M,Np)      ! log time steps versus  log ( 1/ tolerance)
                                            ! first index: grid, second index: scheme  
    
    
    character (len=20) :: names(Np) = [ "HeunEuler21", "RK21", "BogackiShampine", "CashKarp", "RK87", "RK65", "DOPRI54", "Fehlberg54", "Fehlberg87", "Verner65"] 
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/Tstepsa", &
                                        "./doc/chapters/Cauchy_problem/figures/Tstepsb"  ] 
    integer :: i, j
    
    write(*,*) "Convergence rate: Error versus tolerance "   
    write(*,*) "Temporal schemes :  embedded Runge Kutta "  
    write(*,*)  "press enter "; read(*,*) 
    
    U0 = [3, 4]
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    log_mu = [( i, i=1, M ) ]
       
    do j=1, Np 
      call set_solver(family_name = "eRK", scheme_name = names(j) )
      call Temporal_effort_with_tolerance( Time, F_kepler_polar, U0,&        !VanDerPol_equation             
                                           log_mu, log_effort(:,j)  )
     end do 
    
     call plot_parametrics( log_mu, log_effort(:,1:3), names(1:3),      &
                         "$-\log \epsilon$", "$\log M$ ", "(a)", path(1))

     
     call plot_parametrics( log_mu, log_effort(:,4:7), names(4:7),      &
                         "$-\log \epsilon$", "$\log M$ ", "(b)", path(2))

     
     
end subroutine


!**************************************************************************************************
! ******* COMPUTATIONAL EFFORT OF RUNGE-KUTTA SCHEME *******
!**************************************************************************************************

subroutine Computational_effort_comparison

    real :: t0 = 0, tf = 30                 ! Time domain 
    integer, parameter :: N = 300           ! Time steps of the first grid 
    integer, parameter :: Nv = 2            ! Number of variables van de pol
    integer, parameter :: M = 9             ! number of different tolerances           
    real :: Time (0:N)                      ! Time domain 
    real :: U0(Nv)                          ! Initial conditions  
    real :: U(0:N,Nv)                       ! Solution 
    integer, parameter :: Np = 3            ! number of numerical schemes to be evaluated 
    real :: log_mu(M), log_effort(M,Np)      ! log time steps versus  log ( 1/ tolerance)
                                            ! first index: grid, second index: scheme  
    !
    !real :: x_aux(0:N,1), y_aux(:,1:3)
    real :: x_aux(M,Np)
    !
    
    character (len=20) :: family(Np) = [ "weRK", "wGBS", "GBS"]!, "Euler", "Inverse_Euler", "Crank_Nicolson" ]
    character (len=20) :: names(Np) = [ "WDOP853", "wGBS", "GBS"]!, "Euler", "Inverse_Euler", "Crank_Nicolson" ] 

    integer :: i, j
    
    write(*,*) "Convergence rate: Error versus tolerance "   
    write(*,*) "Temporal schemes :  RK87 vs WDOP853 vs GBS vs wGBS "  
    write(*,*)  "press enter "; read(*,*) 
    
    U0 = [3, 4]
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    log_mu = [( i, i=1, M ) ]
       
    do j=1, Np 
      call set_solver(family_name = family(j), scheme_name = names(j) )
      call Temporal_effort_with_tolerance( Time, F_kepler_polar, U0,&                     
                                           log_mu, log_effort(:,j)  )
     end do 
    
     call plot_parametrics( log_mu, log_effort(:,1:3), names(1:3),      &
                         "$-\log \epsilon$", "$\log M$ ", "(a)")
     
      do j=1, Np
         x_aux(:,j) = log_mu(:)
      end do
      
     call savedata(Filename = "./Computational_effort", x = x_aux(:,:), y = log_effort(:,:), N = Np)
     
     
end subroutine

!**************************************************************************************************
! ******* HENON-HEILES SYSTEM *******
!**************************************************************************************************
! Henon equation funtion !!! already in numerical hub

!function Henon( U, t) result(F)
!        real    :: U(:)
!        real    :: F(size(U))
!        
!        real    :: x, y, px, py, lambda = -1
!        
!        x = U(1) ; y = U(2); px = U(3); py = U(4)
!        
!        F = [ px, py, -x-2* lambda*x*y, -y-lambda*(x**2 - y**2) ]        
!end function

!! METERLE OTROS ESQUEMAS PARA COMPARAR??

subroutine HenonHeiles_system

    integer, parameter  :: N = 1000, Nv = 4, M = 1 ! Time steps
    real, parameter     :: dt = 0.1
    real    :: t0 = 0, tf = dt * N
    real    :: Time (0:N), U(0:N, Nv), H(0:N)
    integer :: i
    
    Time = [ ( t0 + (tf-t0) * i / (1d0 * N), i=0, N ) ]
    
    !U(0,:) = [0., 0., 0.6, 0.]
    !U(0,:) = [0.5, 0.5, 0., 0.]     ! Chaotic behaviour 1
    U(0,:) = [0.2, 0.2, 0.4, 0.]     ! Chaotic behaviour 2
    call set_solver("GBS")
    call set_tolerance(1d-2)
    call Cauchy_Problems( Time, Henon_equation, U )
    
    write (*, *) ' Henon-Heiles system, trayectory of the star '
    write(*,*)  "press enter "; read(*,*)
    
    call plot_parametrics(U(:,1),U(:,2:2), ["Henon-Heiles"], "x", "y")
    call savedata(Filename = "./Henon-Heiles_x_y", x = U(:, 1:2), y = U(:, 2:2), N = 1)
    
    write (*, *) ' Henon-Heiles system, projection (x, dx) of the solution in the phase plane '
    write(*,*)  "press enter "; read(*,*)
    
    call plot_parametrics(U(:,1),U(:,3:3), ["Henon-Heiles"], "x", "px")
    call savedata(Filename = "./Henon-Heiles_x_px", x = U(:, 1:2), y = U(:, 3:3), N = 1)
    
    write (*, *) ' Henon-Heiles system, projection (y, dy) of the solution in the phase plane '
    write(*,*)  "press enter "; read(*,*)
    
    call plot_parametrics(U(:,2),U(:,4:4), ["Henon-Heiles"], "y", "py")
    call savedata(Filename = "./Henon-Heiles_y_py", x = U(:, 2:3), y = U(:, 4:4), N = 1)
    
    write (*, *) ' Henon-Heiles system, projection (dx, dy) of the solution in the phase plane '
    write(*,*)  "press enter "; read(*,*)
    
    call plot_parametrics(U(:,3),U(:,4:4), ["Henon-Heiles"], "px", "py")
    call savedata(Filename = "./Henon-Heiles_px_py", x = U(:, 3:4), y = U(:, 4:4), N = 1)
    
end subroutine



!**************************************************************************************************
! ******* SUBROUTINES FOR SAVING FILES *******
!**************************************************************************************************
subroutine savedata(Filename,x,y,z,N)
 !   
	!character (len=*), intent(in) :: Filename
	!real, dimension(:,:),intent(in) :: x, y
	!real, dimension(:,:),intent(in), optional :: z
	!integer, intent(in), optional :: N
 !   
	!integer :: i, j
 !   
	!write(*,*) "Saving ", Filename, " ..."
 !   
	!if (present(z)) then
 !   	open (unit = 1, file = Filename)
 !       	do j = 1,N
 !           	write(1,*) '%'
 !           	do i=1, size(x(:,j))
 !               	write(1,*) x(i,j), y(i,j), z(i,j)
 !           	end do
 !       	end do
 !   	close(1)
 !   else
 !   	open (unit = 1, file = Filename)
 !       	do j = 1,N
 !           	write(1,*) '%'
 !           	do i=1, size(x(:,j))
 !               	write(1,*) x(i,j), y(i,j)
 !           	end do
 !       	end do
 !   	close(1)
	!end if
 !!   
 !
	!character (len=*), intent(in) :: Filename
	!real, dimension(:,:),intent(in) :: x, y
	!real, dimension(:,:),intent(in), optional :: z
	!integer, intent(in), optional :: N
 !   
 !   integer :: i, j
 !       
	!write(*,*) "Saving ", Filename, " ..."
 !   
	!if (present(z)) then
 !   	open (unit = 1, file = Filename)
 !           write(1,*) '%'
 !           do i=1,size(x(:,1))
 !               write(1,*) x(i,1), y(i,1), z(i,1)
 !           end do 
 !           do j = 2, size(z(:,1))
 !               do i=1, size(x(:,j))
 !                   write(1,*) 'Nan     ', 'Nan     ' ,z(i,j)
 !               end do
 !           end do
 !   	close(1)
 !   else
 !   	open (unit = 1, file = Filename)
 !       	do j = 1,N
 !           	write(1,*) '%'
 !           	do i=1, size(x(:,j))
 !               	write(1,*) x(i,j), y(i,j)
 !           	end do
 !       	end do
 !   	close(1)
	!end if
    	character (len=*), intent(in) :: Filename
	real, dimension(:,:),intent(in) :: x, y
	real, dimension(:,:),intent(in), optional :: z
	integer, intent(in), optional :: N
    
    integer :: i, j
        
	write(*,*) "Saving ", Filename, " ..."
    
	if (present(z)) then
        open (unit = 1, file = Filename)
        	do j = 1,N
            	write(1,*) '%'
            	do i=1, size(x(:,j))
                	write(1,*) x(i,j), y(i,j), z(i,j)
            	end do
        	end do
    	close(1)
    else
    	open (unit = 1, file = Filename)
        	do j = 1,N
            	write(1,*) '%'
            	do i=1, size(x(:,j))
                	write(1,*) x(i,j), y(i,j)
            	end do
        	end do
    	close(1)
    end if

end subroutine

    
end module  