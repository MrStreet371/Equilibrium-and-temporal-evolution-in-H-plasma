!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!           Newton-Raphson Saha solver for H plasma                    C
!                                                                      C
!   Author:   Cayetano Hernández Sánchez                               C
!                                                                      C
!   Version:   23.09.26                                                C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



module plasma_physics 
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307) !To ensure double precission
contains
    ! Two steps gaussian elimination subroutine to solve Jac * dx = F 
    subroutine gauss_elimination(Jac, F, dx, n) 
        integer, intent(in) :: n  !Only read n 
        real(dp), intent(inout) :: Jac(n,n), F(n) !Modified inside
        real(dp), intent(out) :: dx(n) !Only write it out 
        integer :: i, j, k
        real(dp) :: factor, temp

        !Step 1: Forward elimination 
        
        do k = 1, n-1 !Loop initialization walking all pivots except last...  
            do i = k+1, n
                factor = Jac(i,k) / Jac(k,k) !x how much each pivot
                do j = k, n !Now we substract factor times row k from i colum by colum from k to n 
                    Jac(i,j) = Jac(i,j) - factor * Jac(k,j)                     
                end do
                F(i) = F(i) - factor * F(k) !same operation to maintain linear equivalence 
            end do
        end do !Now Jac is triangular superior 

        !Step 2: Backwards subsitution
        !Now solve from last to first unknown...
        dx(n) = F(n) / Jac(n,n)
        do i = n-1, 1, -1
            temp = F(i)
            do j = i+1, n 
                temp = temp - Jac(i,j) * dx(j) !Substract contribution of already known unknown...
            end do
            dx(i) = temp / Jac(i,i) !solve by dividing by the pertinent element 
        end do
    end subroutine gauss_elimination
end module plasma_physics


!------------------------------------------------------------

program plasma_saha
    use plasma_physics
    implicit none

    integer, parameter :: n_eq = 4
    integer, parameter :: n_steps = 100
    
    ! Constantes físicas
    real(dp) :: k_B_erg = 1.380649e-16_dp !Boltzmann cte erg/K 
    real(dp) :: k_B_eV  = 8.617333e-5_dp !Boltzmann cte ev/K
    real(dp) :: m_e     = 9.109383e-28_dp !Mass of the electron g
    real(dp) :: h       = 6.626070e-27_dp !Planck cte erg*s
    real(dp) :: T       = 7000.0_dp !T in Kelvin 
    real(dp) :: pi      = 3.14159265358979323846_dp !pi

    real(dp) :: factor_termico, K1, K2
    real(dp) :: N_total(n_steps), N_ini, N_end
    real(dp) :: x(n_eq), F(n_eq), Jac(n_eq, n_eq), dx(n_eq)
    integer  :: i, iter
    real(dp) :: factor_escala

    ! Set up of Saha constants...
    factor_termico = (2.0_dp * pi * m_e * k_B_erg * T / (h**2))**1.5_dp !factor termico por comodidad 
    K1 = factor_termico * exp(-13.60_dp / (k_B_eV * T)) !Cte for ionization 
    K2 = 4.0_dp * factor_termico * exp(-0.754_dp / (k_B_eV * T)) !Cte for electron capture 

    ! Configuration of the range of densities 
    N_ini = 1.0e18_dp
    N_end = 1.0e8_dp
    do i = 1, n_steps ! Generate 100 points separated logarithmically between N_ini y N_end.
        N_total(i) = N_ini * (N_end/N_ini)**(real(i-1, dp)/real(n_steps-1, dp))
    end do

    open(unit=10, file='densidades_equilibrio.dat', status='replace')
    write(10, *) 'N_total            N_H                N_e                N_H+               N_H-'

    ! Initial guess for the high density limit (N = 10^18)
    ! x = [N_H, N_e, N_H+, N_H-]
    x(1) = N_total(1)
    x(2) = sqrt(K1 * N_total(1))
    x(3) = sqrt(K1 * N_total(1))
    x(4) = (N_total(1) * sqrt(K1 * N_total(1))) / K2

    ! Start of the main resolution loop 
    do i = 1, n_steps
        ! Start of the Newton-Raphson
        do iter = 1, 50

            if (x(2) < 1.0e-30_dp) x(2) = 1.0e-30_dp   ! Limit minima to avoid division by 0

            ! Build F(x)
            F(1) = x(3) - K1 * (x(1) / x(2))
            F(2) = x(4) - (x(1) * x(2)) / K2
            F(3) = x(1) + x(3) + x(4) + x(2) - N_total(i)
            F(4) = x(3) - x(2) - x(4)

            ! Check convergence
            if (maxval(abs(F)) < 1.0e-10_dp) exit

            ! Build the jacobian matrix
            Jac(1,1) = -K1 / x(2)
            Jac(1,2) = K1 * x(1) / (x(2)**2)
            Jac(1,3) = 1.0_dp
            Jac(1,4) = 0.0_dp

            Jac(2,1) = -x(2) / K2
            Jac(2,2) = -x(1) / K2
            Jac(2,3) = 0.0_dp
            Jac(2,4) = 1.0_dp

            Jac(3,1) = 1.0_dp
            Jac(3,2) = 1.0_dp
            Jac(3,3) = 1.0_dp
            Jac(3,4) = 1.0_dp

            Jac(4,1) = 0.0_dp
            Jac(4,2) = -1.0_dp
            Jac(4,3) = 1.0_dp
            Jac(4,4) = -1.0_dp

            ! Solve Jac * dx = -F by calling gaussian elimination subroutine 
            F = -F
            call gauss_elimination(Jac, F, dx, n_eq)

            ! Damping steps // This is to avoid divergence and negative densities 
            if (x(1) + dx(1) <= 0.0_dp) dx(1) = -0.9_dp * x(1)
            if (x(2) + dx(2) <= 0.0_dp) dx(2) = -0.9_dp * x(2)
            if (x(3) + dx(3) <= 0.0_dp) dx(3) = -0.9_dp * x(3)
            if (x(4) + dx(4) <= 0.0_dp) dx(4) = -0.9_dp * x(4)

            x = x + dx
        end do

        ! Save data
        write(10, '(5E19.6)') N_total(i), x(1), x(2), x(3), x(4)

        ! Scale and use as starting point for next iter
        if (i < n_steps) then
            factor_escala = N_total(i+1) / N_total(i)
            x = x * factor_escala
        end if
    end do

    close(10)
    print *, "exito!!"

end program plasma_saha

