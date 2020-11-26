!***************************************************
!* Book:  How to learn Applied maths
!***************************************************    
program main_NumericalHUB

       
       use API_Example_Systems_of_Equations
       use API_Example_Lagrange_Interpolation
       use API_Example_spectral_series
       use API_Example_Cauchy_Problem
       use API_Example_Finite_Differences
       use API_Example_Boundary_Value_Problem
       use API_Example_Initial_Boundary_Value_Problem
       use API_Example_IBVP_and_BVP
       
       use MUSE_2020 
          
       implicit none 
       integer :: option = 1  
     
           
do while (option>0) 
    
     write(*,*) "Welcome to NumericalHUB" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Systems of equations  "
     write(*,*) " 2. Lagrange interpolation  "
     write(*,*) " 3. Chebyshev interpolation  "
     write(*,*) " 4. ODE Cauchy problems   "
     write(*,*) " 5. Finite difference   "
     write(*,*) " 6. Boundary value problems  "
     write(*,*) " 7. Initial-boundary value problems  "
     write(*,*) " 8. Mixed problems: IBVP+BVP  "
     write(*,*) " 9. Advanced methods ODE methods "
     write(*,*) " 10. MUSE 2020 "
     
     read(*,*) option 
     
     select case(option)
     case(1) 
         call Systems_of_Equations_examples
         
     case(2) 
         call Lagrange_Interpolation_examples 
         
     case(3) 
         call Fourier_and_Chebyshev_examples
         
     case(4)
         call Cauchy_problem_examples
         
     case(5) 
         call Finite_difference_examples
      
     case(6) 
         call BVP_examples
      
     case(7) 
         call IBVP_examples 
      
     case(8) 
         call Nonlinear_Plate_Vibrations
       
     case(9) 
         call Advanced_Cauchy_problem_examples
         
     case(10)
         call Orbits_and_Numerical_Methods
         
         case default
              
     end select 
     
end do
  
 
end program  

