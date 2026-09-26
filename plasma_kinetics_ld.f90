!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!           Backward euler for H plasma time evolution (N \sim 10^8)   C
!                                                                      C
!   Author:   Cayetano Hernández Sánchez                               C
!                                                                      C
!   Version:   24.09.26                                                C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


module plasma_physics 
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307) !To ensure double precision
contains
    ! Two steps gaussian elimination subroutine to solve Jac * dx = -F 
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
                do j = k, n !Now we subtract factor times row k from i column by column from k to n 
                    Jac(i,j) = Jac(i,j) - factor * Jac(k,j)                     
                end do
                F(i) = F(i) - factor * F(k) !same operation to maintain linear equivalence 
            end do
        end do !Now Jac is triangular superior 

        !Step 2: Backwards substitution
        !Now solve from last to first unknown...
        dx(n) = F(n) / Jac(n,n)
        do i = n-1, 1, -1
            temp = F(i)
            do j = i+1, n 
                temp = temp - Jac(i,j) * dx(j) !Subtract contribution of already known unknown...
            end do
            dx(i) = temp / Jac(i,i) !solve by dividing by the pertinent element 
        end do
    end subroutine gauss_elimination
end module plasma_physics


program plasma_kinetics
    use plasma_physics
    implicit none

    integer, parameter :: n_eq = 4
    integer, parameter :: n_steps = 1000
    
    ! Constantes físicas
    real(dp) :: k_B_erg = 1.380649e-16_dp !Boltzmann cte erg/K 
    real(dp) :: k_B_eV  = 8.617333e-5_dp !Boltzmann cte eV/K
    real(dp) :: m_e     = 9.109383e-28_dp !Mass of the electron g
    real(dp) :: h       = 6.626070e-27_dp !Planck cte erg*s
    real(dp) :: T       = 7000.0_dp !T in Kelvin 
    real(dp) :: pi      = 3.14159265358979323846_dp !pi

    ! Reaction rates variables
    real(dp) :: alpha1, alpha2, alpha3, alpha4
    real(dp) :: factor_termico, K1, K2, zeta1, zeta2
    real(dp) :: rr, pi_r, ra, pd, mn, ip 
    ! rr (recombination), pi_r (photoionization), ra (association), pd (photodetachment), mn (mutual neutralization), ip (ion pair production)

    ! Time mesh and state variables
    real(dp) :: t_array(n_steps), dt
    real(dp) :: x_old(n_eq), x_new(n_eq)
    real(dp) :: F_rates(n_eq), J_rates(n_eq, n_eq)
    real(dp) :: G(n_eq), Jac_G(n_eq, n_eq), dx(n_eq)
    
    integer :: i, step, iter

    ! Direct rates as for table 2
    alpha1 = 3.6e-12_dp * (T / 300.0_dp)**(-0.75_dp)
    alpha2 = 3.0e-16_dp * (T / 300.0_dp)**(0.95_dp) * exp(-T / 9320.0_dp)
    alpha3 = 4.0e-8_dp  * (T / 300.0_dp)**(-0.50_dp)

    ! We define Saha constants and detailed balance for inverse rates
    factor_termico = (2.0_dp * pi * m_e * k_B_erg * T / (h**2))**1.5_dp !factor termico
    K1 = factor_termico * exp(-13.60_dp / (k_B_eV * T)) !Cte for ionization 
    K2 = 4.0_dp * factor_termico * exp(-0.754_dp / (k_B_eV * T)) !Cte for H-

    zeta1  = alpha1 * K1
    zeta2  = alpha2 * K2
    alpha4 = alpha3 * (K1 / K2)

    ! Stablishing the initial guess 
    ! x = [N_H, N_H+, N_H-, N_e]
    x_old(1) = 1.0e8_dp
    x_old(2) = 2.0e8_dp
    x_old(3) = 1.0e8_dp
    x_old(4) = x_old(2) - x_old(3) !N_e by charge conservation

    ! Configuration of the temporal mesh 
    do i = 1, n_steps ! Generate points logarithmically between 10^-8 and 10^10 s
        t_array(i) = 1.0e-8_dp * (1.0e10_dp / 1.0e-8_dp)**(real(i-1, dp)/real(n_steps-1, dp))
    end do

    open(unit=20, file='evolucion_temporal2.dat', status='replace')
    write(20, *) 'Tiempo(s)          N_H                N_H+               N_H-               N_e'
    write(20, '(5E19.6)') 0.0_dp, x_old(1), x_old(2), x_old(3), x_old(4)

    ! Start of the main resolution loop (Backward Euler)
    do step = 2, n_steps
        dt = t_array(step) - t_array(step-1)
        x_new = x_old
        
        ! Start of the Newton-Raphson iteration...
        do iter = 1, 50
            ! Evaluate the physical processes 
            rr = alpha1 * x_new(2) * x_new(4)
            pi_r = zeta1 * x_new(1)
            ra = alpha2 * x_new(1) * x_new(4)
            pd = zeta2 * x_new(3)
            mn = alpha3 * x_new(2) * x_new(3)
            ip = alpha4 * (x_new(1)**2)

            ! Build the F_rates vector
            F_rates(1) = rr - pi_r - ra + pd + 2.0_dp*mn - 2.0_dp*ip
            F_rates(2) = -rr + pi_r - mn + ip
            F_rates(3) = ra - pd - mn + ip
            F_rates(4) = -rr + pi_r - ra + pd

            ! Build the jacobian matrix of the rates...
            J_rates(1,1) = -zeta1 - alpha2*x_new(4) - 4.0_dp*alpha4*x_new(1)
            J_rates(1,2) = alpha1*x_new(4) + 2.0_dp*alpha3*x_new(3)
            J_rates(1,3) = zeta2 + 2.0_dp*alpha3*x_new(2)
            J_rates(1,4) = alpha1*x_new(2) - alpha2*x_new(1)

            J_rates(2,1) = zeta1 + 2.0_dp*alpha4*x_new(1)
            J_rates(2,2) = -alpha1*x_new(4) - alpha3*x_new(3)
            J_rates(2,3) = -alpha3*x_new(2)
            J_rates(2,4) = -alpha1*x_new(2)

            J_rates(3,1) = alpha2*x_new(4) + 2.0_dp*alpha4*x_new(1)
            J_rates(3,2) = -alpha3*x_new(3)
            J_rates(3,3) = -zeta2 - alpha3*x_new(2)
            J_rates(3,4) = alpha2*x_new(1)

            J_rates(4,1) = zeta1 - alpha2*x_new(4)
            J_rates(4,2) = -alpha1*x_new(4)
            J_rates(4,3) = zeta2
            J_rates(4,4) = -alpha1*x_new(2) - alpha2*x_new(1)

            ! Build the objective function G
            G = x_new - x_old - dt * F_rates
            
            ! Build the jacobian matrix for Newton-Raphson: Jac_G = I - dt*J_rates
            Jac_G = -dt * J_rates
            do i = 1, n_eq
                Jac_G(i,i) = Jac_G(i,i) + 1.0_dp
            end do

            ! Solve Jac_G * dx = -G by calling gaussian elimination subroutine 
            G = -G
            call gauss_elimination(Jac_G, G, dx, n_eq)
            
            ! Update densities 
            x_new = x_new + dx
            
            ! Impose floor limit 
            do i = 1, n_eq
                if (x_new(i) < 1.0e-30_dp) x_new(i) = 1.0e-30_dp
            end do

            ! Check convergence
            if (maxval(abs(dx)/x_new) < 1.0e-6_dp) exit
        end do
        
        ! Use as starting point for next time step
        x_old = x_new
        write(20, '(5E19.6)') t_array(step), x_old(1), x_old(2), x_old(3), x_old(4)
    end do

    close(20)
    print *, "exito!!"

end program plasma_kinetics
