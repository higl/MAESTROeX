
module initdata_module

  use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor
  use network, only: nspec
  use amrex_fort_module, only : amrex_spacedim, amrex_random
  use amrex_constants_module
  use fundamental_constants_module
  use base_state_geometry_module, only: nr_fine, max_radial_level
  use meth_params_module, only: nscal, rho_comp, rhoh_comp, temp_comp, spec_comp, pi_comp, &
       prob_lo, prob_hi, grav_const
  use probin_module, only: velpert_amplitude, velpert_radius, velpert_scale, & 
       velpert_steep, velpert_min, p0_ext, dtheta0, pottemp
  use extern_probin_module, only: eos_gamma
  use eos_module
  use eos_type_module
  
  implicit none

  private

contains

  subroutine initdata(lev, time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, &
       dx) bind(C, name="initdata")

    integer         , intent(in   ) :: lev, lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), 1:nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), 1:nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)

    ! Local variables
    integer :: i, j, n, k, comp
    double precision :: x, y, r
    double precision :: x0, y0, r0
    double precision :: gammainv, entK, intc, intf
    double precision :: dtheta, GPhi
    double precision :: p, rho, rho0_ext
    double precision :: GasConst
    
    type (eos_t) :: eos_state

    ! random numbers between -1 and 1
    double precision :: alpha(3,3), beta(3,3), gamma(3,3)

    ! random numbers between 0 and 2*pi
    double precision :: phix(3,3), phiy(3,3)
    
    ! a random number
    double precision :: rand

    ! L2 norm of k
    double precision :: norm(3,3)

    ! Local variables
    integer :: iloc, jloc

    ! cos and sin of (2*pi*kx/L + phix), etc
    double precision :: cx(3,3), cy(3,3), cz(3,3)
    double precision :: sx(3,3), sy(3,3), sz(3,3)

    ! radius, or distance, to center of star
    double precision :: rloc

    ! the point we're at
    double precision :: xloc(2)
    
    ! the center
    double precision :: xc(2)

    ! perturbational velocity to add
    double precision :: upert(2)


       ! generate random numbers
       ! random numbers are not currently supported
       ! use functions that result in numbers between (-1,1)
       
          do j=1,3
             do i=1,3
                rand = amrex_random()
                rand = 2.0d0*rand - 1.0d0
                ! rand = (.5)**i * (.7)**j * (.3)**k * (-1.)**i
                alpha(i,j) = rand
                rand = amrex_random()
                rand = 2.0d0*rand - 1.0d0
                ! rand = (.5)**i * (.3)**j * (.7)**k * (-1.)**j
                beta(i,j) = rand
                rand = amrex_random()
                rand = 2.0d0*rand - 1.0d0
                ! rand = (.3)**i * (.5)**j * (.7)**k * (-1.)**k
                gamma(i,j) = rand
                rand = amrex_random()
                ! rand = (.3)**i * (.7)**j * (.5)**k
                rand = 2.0d0*M_PI*rand
                phix(i,j) = rand
                rand = amrex_random()
                ! rand = (.7)**i * (.3)**j * (.5)**k
                rand = 2.0d0*M_PI*rand
                phiy(i,j) = rand
                rand = amrex_random()
             enddo
          enddo
      
          do j=1,3
             do i=1,3
                norm(i,j) = sqrt(dble(i)**2+dble(j)**2)
             enddo
          enddo
    ! initialize the velocity to zero everywhere
    vel = ZERO
    
    ! define where center of star is

    xc(1) = 0.5d0*(prob_lo(1)+prob_hi(1))
    xc(2) = 0.5d0*(prob_lo(2)+prob_hi(2))
    
    ! now do the big loop over all points in the domain
    do k = lo(3), hi(3)
        do jloc = lo(2),hi(2)
           do iloc = lo(1),hi(1)

                 ! set perturbational velocity to zero
                 upert = ZERO
                
                 ! compute where we physically are
                 xloc(1) = prob_lo(1) + (dble(iloc)+0.5d0)*dx(1)
                 xloc(2) = prob_lo(2) + (dble(jloc)+0.5d0)*dx(2)

                 ! compute distance to the center of the star
                 rloc = xloc(2)

                 ! loop over the 9 combinations of fourier components
                 do i=1,3
                    do j=1,3
                          ! compute cosines and sines
                          cx(i,j) = cos(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(i,j))
                          cy(i,j) = cos(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(i,j))
                          sx(i,j) = sin(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(i,j))
                          sy(i,j) = sin(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(i,j))
                    enddo
                 enddo

                 ! loop over the 9 combinations of fourier components
                 do i=1,3
                    do j=1,3
                          ! compute contribution from perturbation velocity from each mode
                          upert(1) = upert(1) + &
                               (-gamma(i,j)*dble(j)*cx(i,j)*sy(i,j) &
                                 +beta(i,j)*dble(j)*cy(i,j)*sx(i,j)) &
                                 / norm(i,j)
                                
                          upert(2) = upert(2) + &
                               (gamma(i,j)*dble(i)*cy(i,j)*sx(i,j) &
                               -alpha(i,j)*dble(j)*cx(i,j)*sy(i,j)) &
                                 / norm(i,j)
                                
                    enddo
                 enddo

                 ! apply the cutoff function to the perturbational velocity
                 do i=1,2
                    upert(i) = velpert_amplitude*upert(i) &
                         *(0.5d0+0.5d0*tanh((velpert_radius-rloc)/velpert_steep))
                 enddo

                 ! add perturbational velocity to background velocity
                 do i=1,2
                    vel(iloc,jloc,k,i) = vel(iloc,jloc,k,i) + upert(i)
                 enddo

           enddo
        enddo
    enddo 
    
    !compute the gas constant R
    GasConst = k_B * n_A
    
    gammainv = 1.0/eos_gamma
    rho0_ext = p0_ext/(pottemp* GasConst)
    entK = p0_ext/rho0_ext**eos_gamma
    intc = p0_ext**(1.0 - gammainv)
    intf = (1.0 - gammainv)/entK**gammainv

    ! Parameters defining the bubble.
    x0 = 5.0e5
    y0 = 2.75e5
    r0  = 2.5e5

    ! initial the domain with the base state
    scal = ZERO

    ! initialize the scalars
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              scal(i,j,k,rho_comp)  = s0_init(lev,j,rho_comp)
              scal(i,j,k,rhoh_comp) = s0_init(lev,j,rhoh_comp)
              scal(i,j,k,temp_comp) = s0_init(lev,j,temp_comp)
              scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                   s0_init(lev,j,spec_comp:spec_comp+nspec-1)
           enddo
        enddo
    enddo   
    ! initial the density perturbations
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            y = prob_lo(2) + (dble(j)+HALF) * dx(2)
           
           
            do i = lo(1), hi(1)
                x = prob_lo(1) + (dble(i)+HALF) * dx(1)
                          
                          
                r   = sqrt(((x-x0)/r0)**2+((y-y0)/r0)**2)

                GPhi = grav_const*y
                p = (intc + intf*GPhi)**(1.0/(1.0-gammainv))

                if (r < 1.0) then
                   dtheta = dtheta0*(cos(0.5*M_PI*r)**2)
                else
                   dtheta = 0.0
                end if

                rho = p0_ext**(1.0-gammainv)*p**gammainv/GasConst/(pottemp+dtheta)

              
                !make the state thermodynamically consistent again
                eos_state%T     = 1d-15 !a temperature guess
                eos_state%rho   = rho
                eos_state%p     = p
                eos_state%xn(:) = scal(i,j,k,spec_comp:spec_comp+nspec-1)/scal(i,j,k,rho_comp)

                ! (rho,p) --> T, h
                call eos(eos_input_rp, eos_state)

                scal(i,j,k,rho_comp) = rho
                scal(i,j,k,rhoh_comp) = rho * eos_state%h
                scal(i,j,k,temp_comp) = eos_state%T
                scal(i,j,k,spec_comp:spec_comp+nspec-1) = &
                    eos_state%xn(:)*rho
                
           enddo
        enddo
    enddo    

  end subroutine initdata


  subroutine initdata_sphr(time, lo, hi, &
       scal, scal_lo, scal_hi, nc_s, &
       vel, vel_lo, vel_hi, nc_v, &
       s0_init, p0_init, &
       dx, r_cc_loc, r_edge_loc, &
       cc_to_r, ccr_lo, ccr_hi) bind(C, name="initdata_sphr")

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: scal_lo(3), scal_hi(3), nc_s
    integer         , intent(in   ) :: vel_lo(3), vel_hi(3), nc_v
    double precision, intent(in   ) :: time
    double precision, intent(inout) :: scal(scal_lo(1):scal_hi(1), &
         scal_lo(2):scal_hi(2), &
         scal_lo(3):scal_hi(3), nc_s)
    double precision, intent(inout) :: vel(vel_lo(1):vel_hi(1), &
         vel_lo(2):vel_hi(2), &
         vel_lo(3):vel_hi(3), nc_v)
    double precision, intent(in   ) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(in   ) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: dx(3)
    double precision, intent(in   ) :: r_cc_loc (0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine)
    integer         , intent(in   ) :: ccr_lo(3), ccr_hi(3)
    double precision, intent(in   ) :: cc_to_r(ccr_lo(1):ccr_hi(1), &
         ccr_lo(2):ccr_hi(2),ccr_lo(3):ccr_hi(3))

  end subroutine initdata_sphr

end module initdata_module
