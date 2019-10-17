module base_state_module
  ! init_base_state is used to initialize the base state arrays from the
  ! model file.  The actual reading of the model file is handled by the
  ! model_parser_module in Util/
  !
  ! Note: The initial base state quantities returned from this routine
  ! are only a temporary base state.  These quantities are mapped onto
  ! the full 2- or 3-d state in initscaldata.f90 and a new base state is
  ! created after initialization by averaging the density and calling
  ! enforce_HSE in initialize.f90.

  use model_parser_module
  use eos_type_module
  use eos_module
  use amrex_constants_module
  use simple_log_module
  use inlet_bc_module
  use fundamental_constants_module
  use amrex_fort_module, only: amrex_spacedim
  use network, only: nspec
  use meth_params_module, only: nscal, model_file, spherical, base_cutoff_density, &
                                do_2d_planar_octant, do_planar_invsq_grav, rho_comp, &
                                rhoh_comp, spec_comp, temp_comp, grav_const, &
                                planar_invsq_mass, print_init_hse_diag, prob_lo
                                
  use base_state_geometry_module, only: nr_fine, dr, nr, max_radial_level
  use probin_module, only:  p0_ext, dtheta0, pottemp
  use extern_probin_module, only: eos_gamma
  
  
  implicit none

  private

contains

  subroutine init_base_state(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init) &
       bind(C, name="init_base_state")
    ! Binds to C function ``init_base_state``

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)

   ! local
    integer         :: i,j,r,n,comp
    double precision:: rloc,rmid
    double precision:: d_ambient,t_ambient,p_ambient,xn_ambient(nspec),e_ambient
    double precision:: GasConst
    double precision:: gammainv, rho0_ext
    double precision:: entK, intc, intf
    double precision:: x0,y0,r0
    double precision:: GPhi
    double precision:: t_guess

    double precision, parameter :: TINY = 1.0d-10

    double precision :: mencl, g, r_l, r_r, dpdr, rhog
    double precision :: max_hse_error
    double precision:: mod_dr

    type (eos_t) :: eos_state

    if (spherical .eq. 1)  then
       call amrex_error("ERROR: rt base_state is not valid for spherical or polar")
    endif
    
887 format(78('-'))
888 format(a60,g18.10)
   do n=0,max_radial_level
        
        !compute the gas constant R
        GasConst = k_B * n_A
         
        ! Parameters defining the atmosphere.
        gammainv = 1.0/eos_gamma
        rho0_ext = p0_ext/(pottemp* GasConst)
        entK = p0_ext/rho0_ext**eos_gamma
        intc = p0_ext**(1.0 - gammainv)
        intf = (1.0 - gammainv)/entK**gammainv

        write(*,*) pottemp, rho0_ext, p0_ext, eos_gamma
        ! Parameters defining the bubble.
        x0 = 5.0e5
        y0 = 2.75e5
        r0  = 2.5e5
        !dtheta0 = 6.6e0
        !dtheta0 = 0.0

        ! set a guess for the temperature for the EOS calls
        t_guess = 1.e-8
        
        ! fill the base state arrays
        do r=0,nr(n)-1

            ! height above the bottom of the domain
            rloc = prob_lo(size(dr)) + (dble(r) + HALF)*dr(n)
            
            GPhi = grav_const*rloc
            
            p_ambient = (intc + intf*GPhi)**(1.0/(1.0-gammainv))
            d_ambient = p0_ext**(1.0-gammainv)*p_ambient**gammainv/GasConst/(pottemp)
            xn_ambient(:) = ONE
            t_ambient = pottemp
            
            ! use the EOS to make the state consistent
            eos_state%T     = t_ambient
            eos_state%rho   = d_ambient
            eos_state%p     = p_ambient
            eos_state%xn(:) = xn_ambient(:)

            ! (rho,p) --> e, h
            call eos(eos_input_rp, eos_state)

            s0_init(n,r, rho_comp) = d_ambient
            s0_init(n,r,rhoh_comp) = d_ambient * eos_state%h
            s0_init(n,r,spec_comp:spec_comp+nspec-1) = xn_ambient(1:nspec) * d_ambient
            p0_init(n,r) = p_ambient
            s0_init(n,r,temp_comp) = eos_state%T 

        end do

        ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
        rho0 = s0_init(:,:,rho_comp)
        rhoh0 = s0_init(:,:,rhoh_comp)
        tempbar = s0_init(:,:,temp_comp)
        tempbar_init = s0_init(:,:,temp_comp)
        p0 = p0_init

        max_hse_error = -1.d30

        do r=1,nr(n)-1


            rloc = prob_lo(size(dr)) + (dble(r) + HALF)*dr(n)

            dpdr = (p0_init(n,r) - p0_init(n,r-1))/dr(n)
            rhog = HALF*(s0_init(n,r,rho_comp) + s0_init(n,r-1,rho_comp))*grav_const

            if (print_init_hse_diag) then
             if ( parallel_IOProcessor() ) then
                print *, 'r, dpdr, rhog, err: ', rloc, p0_init(n,r),s0_init(n,r,rho_comp) , dpdr, rhog, &
                     abs(dpdr - rhog)/abs(rhog)
             endif
            endif

            max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))

        enddo

        if ( parallel_IOProcessor() ) then
           write (*,*) ' '
           write (*,*) 'Maximum HSE Error = ', max_hse_error
           write (*,*) '   (after putting initial model into base state arrays, and'
           write (*,*) '    for density < base_cutoff_density)'
           write (*,887)
           write (*,*) ' '
        endif

       ! initialize any inlet BC parameters
       call set_inlet_bcs()
       
    end do ! end loop over levels

  end subroutine init_base_state


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_base_state_irreg(s0_init,p0_init,rho0,rhoh0,p0,tempbar,tempbar_init, &
                                     r_cc_loc, r_edge_loc) &
       bind(C, name="init_base_state_irreg")
    ! Binds to C function ``init_base_state_irreg``

    double precision, intent(inout) :: s0_init(0:max_radial_level,0:nr_fine-1,1:nscal)
    double precision, intent(inout) :: p0_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::    rho0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::   rhoh0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) ::      p0(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar(0:max_radial_level,0:nr_fine-1)
    double precision, intent(inout) :: tempbar_init(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) ::   r_cc_loc(0:max_radial_level,0:nr_fine-1)
    double precision, intent(in   ) :: r_edge_loc(0:max_radial_level,0:nr_fine  )

   ! local
    integer         :: i,j,r,n,comp
    double precision:: rloc,rmid
    double precision:: d_ambient,t_ambient,p_ambient,xn_ambient(nspec),e_ambient
    double precision:: GasConst
    double precision:: gammainv, rho0_ext
    double precision:: entK, intc, intf
    double precision:: x0,y0,r0
    double precision:: GPhi
    double precision:: t_guess
    double precision:: eos_gamma

    double precision, parameter :: TINY = 1.0d-10

    double precision :: mencl, g, r_l, r_r, dpdr, rhog
    double precision :: max_hse_error
    double precision:: mod_dr

    type (eos_t) :: eos_state

    if (spherical .eq. 1)  then
       call amrex_error("ERROR: rt base_state is not valid for spherical or polar")
    endif
    
887 format(78('-'))
888 format(a60,g18.10)

   do n=0,max_radial_level
        
        eos_gamma = eos_state%gam1
        
        !compute the gas constant R
        GasConst = k_B * n_A
        
        ! Parameters defining the atmosphere.
        gammainv = 1.0/eos_gamma
        rho0_ext = p0_ext/(pottemp* GasConst)
        entK = p0_ext/rho0_ext**eos_gamma
        intc = p0_ext**(1.0 - gammainv)
        intf = (1.0 - gammainv)/entK**gammainv

        ! Parameters defining the bubble.
        x0 = 5.0e5
        y0 = 2.75e5
        r0  = 2.5e5
        !dtheta0 = 6.6e0
        !dtheta0 = 0.0

        ! set a guess for the temperature for the EOS calls
        t_guess = 1.e-8
        
        ! fill the base state arrays
        do r=0,nr(n)-1

            ! height above the bottom of the domain
            rloc = prob_lo(size(dr)) + (dble(r) + HALF)*dr(n)
            
            GPhi = grav_const*rloc
            
            p_ambient = (intc + intf*GPhi)**(1.0/(1.0-gammainv))
            d_ambient = p0_ext**(1.0-gammainv)*p_ambient**gammainv/GasConst/(pottemp)
            xn_ambient(:) = ONE
            t_ambient = pottemp
            
            ! use the EOS to make the state consistent
            eos_state%T     = t_ambient
            eos_state%rho   = d_ambient
            eos_state%p     = p_ambient
            eos_state%xn(:) = xn_ambient(:)

            ! (rho,p) --> e, h
            call eos(eos_input_rp, eos_state)

            s0_init(n,r, rho_comp) = d_ambient
            s0_init(n,r,rhoh_comp) = d_ambient * eos_state%h
            s0_init(n,r,spec_comp:spec_comp+nspec-1) = xn_ambient(1:nspec) * d_ambient
            p0_init(n,r) = p_ambient
            s0_init(n,r,temp_comp) = eos_state%T 

        end do

        ! copy s0_init and p0_init into rho0, rhoh0, p0, and tempbar
        rho0 = s0_init(:,:,rho_comp)
        rhoh0 = s0_init(:,:,rhoh_comp)
        tempbar = s0_init(:,:,temp_comp)
        tempbar_init = s0_init(:,:,temp_comp)
        p0 = p0_init

        max_hse_error = -1.d30

        do r=1,nr(n)-1


            rloc = prob_lo(size(dr)) + (dble(r) + HALF)*dr(n)

            dpdr = (p0_init(n,r) - p0_init(n,r-1))/dr(n)
            rhog = HALF*(s0_init(n,r,rho_comp) + s0_init(n,r-1,rho_comp))*grav_const

            if (print_init_hse_diag) then
             if ( parallel_IOProcessor() ) then
                print *, 'r, dpdr, rhog, err: ', rloc, p0_init(n,r),s0_init(n,r,rho_comp) , dpdr, rhog, &
                     abs(dpdr - rhog)/abs(rhog)
             endif
            endif

            max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))

        enddo

        if ( parallel_IOProcessor() ) then
           write (*,*) ' '
           write (*,*) 'Maximum HSE Error = ', max_hse_error
           write (*,*) '   (after putting initial model into base state arrays, and'
           write (*,*) '    for density < base_cutoff_density)'
           write (*,887)
           write (*,*) ' '
        endif

       ! initialize any inlet BC parameters
       call set_inlet_bcs()
       
    end do ! end loop over levels

  end subroutine init_base_state_irreg
end module base_state_module
