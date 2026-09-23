!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!           Backward euler for H plasma time evolution                 C
!                                                                      C
!   Author:   Cayetano Hernández Sánchez                               C
!                                                                      C
!   Version:   23.09.26                                                C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


module plasma_physics
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
contains
    ! Eliminación Gaussiana reciclada de la Parte A
    subroutine gauss_elimination(Jac, F, dx, n)
        integer, intent(in) :: n
        real(dp), intent(inout) :: Jac(n,n), F(n)
        real(dp), intent(out) :: dx(n)
        integer :: i, j, k
        real(dp) :: factor, temp

        do k = 1, n-1
            do i = k+1, n
                factor = Jac(i,k) / Jac(k,k)
                do j = k, n
                    Jac(i,j) = Jac(i,j) - factor * Jac(k,j)
                end do
                F(i) = F(i) - factor * F(k)
            end do
        end do

        dx(n) = F(n) / Jac(n,n)
        do i = n-1, 1, -1
            temp = F(i)
            do j = i+1, n
                temp = temp - Jac(i,j) * dx(j)
            end do
            dx(i) = temp / Jac(i,i)
        end do
    end subroutine gauss_elimination
end module plasma_physics

program plasma_kinetics
    use plasma_physics
    implicit none

    integer, parameter :: n_eq = 4
    integer, parameter :: n_steps = 500
    
    ! Constantes
    real(dp) :: k_B_erg = 1.380649e-16_dp
    real(dp) :: k_B_eV  = 8.617333e-5_dp
    real(dp) :: m_e     = 9.109383e-28_dp
    real(dp) :: h       = 6.626070e-27_dp
    real(dp) :: T       = 7000.0_dp
    real(dp) :: pi      = 3.14159265358979323846_dp

    ! Tasas de reaccion
    real(dp) :: alpha1, alpha2, alpha3, alpha4
    real(dp) :: factor_termico, K1, K2, zeta1, zeta2
    real(dp) :: rr, pi_r, ra, pd, mn, ip

    ! Variables temporales y de estado
    real(dp) :: t_array(n_steps), dt
    real(dp) :: y_old(n_eq), y_new(n_eq)
    real(dp) :: F_rates(n_eq), J_rates(n_eq, n_eq)
    real(dp) :: G(n_eq), Jac_G(n_eq, n_eq), dy(n_eq)
    
    integer :: i, step, iter

    ! 1. CÁLCULO DE TASAS SEGÚN TABLA 2[cite: 2]
    alpha1 = 3.6e-12_dp * (T / 300.0_dp)**(-0.75_dp)
    alpha2 = 3.0e-16_dp * (T / 300.0_dp)**(0.95_dp) * exp(-T / 9320.0_dp)
    alpha3 = 4.0e-8_dp  * (T / 300.0_dp)**(-0.50_dp)

    ! Constantes de Saha y Balance Detallado para inversas
    factor_termico = (2.0_dp * pi * m_e * k_B_erg * T / (h**2))**1.5_dp
    K1 = factor_termico * exp(-13.60_dp / (k_B_eV * T))
    K2 = 4.0_dp * factor_termico * exp(-0.754_dp / (k_B_eV * T))

    zeta1  = alpha1 * K1
    zeta2  = alpha2 * K2
    alpha4 = alpha3 * (K1 / K2)

	! 2. CONDICIONES INICIALES (Escaladas a N ~ 10^8 cm^-3)
    ! y = [N_H, N_H+, N_H-, N_e]
    y_old(1) = 1.0e8_dp
    y_old(2) = 2.0e8_dp
    y_old(3) = 1.0e8_dp
    y_old(4) = y_old(2) - y_old(3) ! N_e por conservación de carga

	! Malla temporal logarítmica (de 10^-3 s a 10^10 s)
    do i = 1, n_steps
        t_array(i) = 1.0e-3_dp * (1.0e10_dp / 1.0e-3_dp)**(real(i-1, dp)/real(n_steps-1, dp))
    end do
    open(unit=20, file='evolucion_temporal.dat', status='replace')
    write(20, *) 'Tiempo(s)          N_H                N_H+               N_H-               N_e'
    write(20, '(5E19.6)') 0.0_dp, y_old(1), y_old(2), y_old(3), y_old(4)

    ! 3. INTEGRACIÓN IMPLÍCITA TEMPORAL (Backward Euler)
    do step = 2, n_steps
        dt = t_array(step) - t_array(step-1)
        y_new = y_old
        
        ! Iteración Newton-Raphson para el paso temporal implícito
        do iter = 1, 50
            ! Evaluación de procesos y sus derivadas
            rr = alpha1 * y_new(2) * y_new(4)
            pi_r = zeta1 * y_new(1)
            ra = alpha2 * y_new(1) * y_new(4)
            pd = zeta2 * y_new(3)
            mn = alpha3 * y_new(2) * y_new(3)
            ip = alpha4 * (y_new(1)**2)

            ! Sistema de ecuaciones diferenciales F(y)[cite: 2]
            F_rates(1) = rr - pi_r - ra + pd + 2.0_dp*mn - 2.0_dp*ip
            F_rates(2) = -rr + pi_r - mn + ip
            F_rates(3) = ra - pd - mn + ip
            F_rates(4) = -rr + pi_r - ra + pd

            ! Jacobiano Analítico de las tasas de reacción
            J_rates(1,1) = -zeta1 - alpha2*y_new(4) - 4.0_dp*alpha4*y_new(1)
            J_rates(1,2) = alpha1*y_new(4) + 2.0_dp*alpha3*y_new(3)
            J_rates(1,3) = zeta2 + 2.0_dp*alpha3*y_new(2)
            J_rates(1,4) = alpha1*y_new(2) - alpha2*y_new(1)

            J_rates(2,1) = zeta1 + 2.0_dp*alpha4*y_new(1)
            J_rates(2,2) = -alpha1*y_new(4) - alpha3*y_new(3)
            J_rates(2,3) = -alpha3*y_new(2)
            J_rates(2,4) = -alpha1*y_new(2)

            J_rates(3,1) = alpha2*y_new(4) + 2.0_dp*alpha4*y_new(1)
            J_rates(3,2) = -alpha3*y_new(3)
            J_rates(3,3) = -zeta2 - alpha3*y_new(2)
            J_rates(3,4) = alpha2*y_new(1)

            J_rates(4,1) = zeta1 - alpha2*y_new(4)
            J_rates(4,2) = -alpha1*y_new(4)
            J_rates(4,3) = zeta2
            J_rates(4,4) = -alpha1*y_new(2) - alpha2*y_new(1)

            ! Construir la función G = y_new - y_old - dt*F_rates = 0
            G = y_new - y_old - dt * F_rates
            
            ! Construir la matriz del Jacobiano del sistema G: Jac_G = I - dt*J_rates
            Jac_G = -dt * J_rates
            do i = 1, n_eq
                Jac_G(i,i) = Jac_G(i,i) + 1.0_dp
            end do

            ! Resolver Jac_G * dy = -G
            G = -G
            call gauss_elimination(Jac_G, G, dy, n_eq)
            
            ! Actualizar variables
            y_new = y_new + dy
            
            ! Control de positividad
            do i = 1, n_eq
                if (y_new(i) < 1.0e-30_dp) y_new(i) = 1.0e-30_dp
            end do

            ! Criterio de convergencia de Newton
            if (maxval(abs(dy)/y_new) < 1.0e-6_dp) exit
        end do
        
        y_old = y_new
        write(20, '(5E19.6)') t_array(step), y_old(1), y_old(2), y_old(3), y_old(4)
    end do

    close(20)
    print *, "Evolución temporal completada. Datos en 'evolucion_temporal2.dat'"

end program plasma_kinetics
