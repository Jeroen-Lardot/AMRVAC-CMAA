!> RadiatHydrodynamics physics module
module mod_rhd_phys

  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: rhd_energy = .true.

  !> Whether thermal conduction is added
  logical, public, protected              :: rhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: rhd_radiative_cooling = .false.

  !> Whether dust is added
  logical, public, protected              :: rhd_dust = .false.

  !> Whether viscosity is added
  logical, public, protected              :: rhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: rhd_gravity = .false.

  !> Whether particles module is added
  logical, public, protected              :: rhd_particles = .false.

  !> Number of tracer species
  integer, public, protected              :: rhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> The adiabatic index
  !> Radiative gasses have a adiabatic index of 4/3, and not 5/3
  double precision, public, protected     :: rhd_gamma = 4.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: rhd_adiab = 1.0d0

  !> The smallest allowed energy
  double precision, protected             :: small_e

  !> The smallest allowed radiation energy
  double precision, public, protected             :: small_r_e = 0.d0

  !> Helium abundance over Hydrogen
  double precision, public, protected     :: He_abundance = 0.1d0

  !> Index of the radiation energy
  integer, public, protected              :: r_e

  !> Formalism to treat radiation
  character(len=8), public :: rhd_radiation_formalism = 'fld'

  !> In the case of no rhd_energy, how to compute pressure
  character(len=8), public :: rhd_pressure = 'Trad'

  !> Treat radiation fld_Rad_force
  logical, public, protected :: rhd_radiation_force = .true.

  !> Treat radiation-gas energy interaction
  logical, public, protected :: rhd_energy_interact = .true.

  !> Treat radiation energy diffusion
  logical, public, protected :: rhd_radiation_diffusion = .true.

  !> Treat radiation advection
  logical, public, protected :: rhd_radiation_advection = .true.

  !> Do a running mean over the radiation pressure when determining dt
  logical, protected :: radio_acoustic_filter = .false.
  integer, protected :: size_ra_filter = 1

  !> kb/(m_p mu)* 1/a_rad**4,
  double precision, public :: kbmpmua4

  !> Use the speed of light to calculate the timestep
  logical :: dt_c = .false.


  ! Public methods
  public :: rhd_phys_init
  public :: rhd_kin_en
  public :: rhd_get_pthermal
  public :: rhd_get_pradiation
  public :: rhd_get_ptot
  public :: rhd_to_conserved
  public :: rhd_to_primitive
  public :: rhd_get_tgas
  public :: rhd_get_trad
  public :: rhd_set_mg_bounds

contains

  !> Read this module's parameters from a file
  subroutine rhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /rhd_list/ rhd_energy, rhd_pressure, rhd_n_tracer, rhd_gamma,&
        rhd_adiab, rhd_dust, rhd_thermal_conduction, rhd_radiative_cooling,&
        rhd_viscosity, rhd_gravity, He_abundance, SI_unit, rhd_particles,&
        rhd_radiation_formalism,rhd_radiation_force, rhd_energy_interact,&
        rhd_radiation_diffusion, rhd_radiation_advection,&
        radio_acoustic_filter, size_ra_filter, dt_c

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, rhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine rhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine rhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = rhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine rhd_write_info

  !> Add fluxes in an angular momentum conserving way
  subroutine rhd_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim)
    use mod_global_parameters
    use mod_dust, only: dust_n_species, dust_mom
    use mod_geometry
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim),  wnew(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    integer, intent(in)                :: idim
    integer                            :: hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
        kxCmin1,kxCmin2,kxCmax1,kxCmax2, iw
    double precision                   :: inv_volume(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    logical isangmom

    ! shifted indexes
    hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
    hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
    ! all the indexes
    kxCmin1=hxOmin1;kxCmin2=hxOmin2;
    kxCmax1=ixOmax1;kxCmax2=ixOmax2;

    inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
       1.0d0/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    select case(coordinate)
    case (cylindrical)
       do iw=1,nwflux
        isangmom = (iw==iw_mom(phi_))
        if (rhd_dust) isangmom = (isangmom .or. any(dust_mom(phi_,&
           1:dust_n_species) == iw))
        if (idim==r_ .and. isangmom) then
          fC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw,idim)= fC(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2,iw,idim)*(x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             r_)+half*block%dx(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idim))
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * (inv_volume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim))
        else
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        endif
      enddo
     case (spherical)
      if (rhd_dust) call mpistop("Error: rhd_angmomfix is not implemented &
      
      
        &with dust and coordinate==sperical")
      do iw=1,nwflux
        if     (idim==r_ .and. (iw==iw_mom(2) .or. iw==iw_mom(phi_))) then
          fC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw,idim)= fC(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2,iw,idim)*(x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             idim)+half*block%dx(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idim))
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * (inv_volume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim))
        elseif (idim==2  .and. iw==iw_mom(phi_)) then
          fC(kxCmin1:kxCmax1,kxCmin2:kxCmax2,iw,idim)=fC(kxCmin1:kxCmax1,&
             kxCmin2:kxCmax2,iw,idim)*sin(x(kxCmin1:kxCmax1,kxCmin2:kxCmax2,&
             idim)+half*block%dx(kxCmin1:kxCmax1,kxCmin2:kxCmax2,idim)) !(x(4,3,1)-x(3,3,1)))
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * (inv_volume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)/sin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)))
        else
          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
             idim)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
             idim)) * inv_volume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        endif
      enddo

    end select

  end subroutine rhd_angmomfix

  !> Initialize the module
  subroutine rhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_dust, only: dust_init
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_fld
    use mod_physics

    integer :: itr, idir

    call rhd_read_params(par_files)

    physics_type = "rhd"
    phys_energy  = rhd_energy
    use_particles = rhd_particles

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (rhd_energy) then
       e_ = var_set_energy()
       p_ = e_
    else
       e_ = -1
       p_ = -1
    end if

    !> set radiation energy
    r_e = var_set_radiation_energy()

    phys_get_dt              => rhd_get_dt
    phys_get_cmax            => rhd_get_cmax
    phys_get_cbounds         => rhd_get_cbounds
    phys_get_flux            => rhd_get_flux
    phys_get_v_idim          => rhd_get_v
    phys_add_source_geom     => rhd_add_source_geom
    phys_add_source          => rhd_add_source
    phys_to_conserved        => rhd_to_conserved
    phys_to_primitive        => rhd_to_primitive
    phys_check_params        => rhd_check_params
    phys_check_w             => rhd_check_w
    phys_get_pthermal        => rhd_get_pthermal
    phys_get_tgas            => rhd_get_tgas
    phys_get_trad            => rhd_get_trad
    phys_write_info          => rhd_write_info
    phys_handle_small_values => rhd_handle_small_values
    phys_angmomfix           => rhd_angmomfix
    phys_set_mg_bounds       => rhd_set_mg_bounds

    ! Whether diagonal ghost cells are required for the physics
    phys_req_diagonal = .false.

    ! derive units from basic units
    call rhd_physical_units()

    if (rhd_dust) call dust_init(rho_, mom(:), e_)

    select case (rhd_radiation_formalism)
    case('fld')
      call fld_init(He_abundance, rhd_radiation_diffusion, rhd_gamma)
    case default
      call mpistop('Radiation formalism unknown')
    end select

    allocate(tracer(rhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, rhd_n_tracer
       tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! initialize thermal conduction module
    if (rhd_thermal_conduction) then
      if (.not. rhd_energy) call mpistop(&
         "thermal conduction needs rhd_energy=T")
      call thermal_conduction_init(rhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (rhd_radiative_cooling) then
      if (.not. rhd_energy) call mpistop(&
         "radiative cooling needs rhd_energy=T")
      call radiative_cooling_init(rhd_gamma,He_abundance)
    end if

    if (rhd_energy_interact) then
      if (.not. rhd_energy) call mpistop(&
         "energy interaction needs rhd_energy=T")
    end if

    ! Initialize viscosity module
    if (rhd_viscosity) call viscosity_init(phys_wider_stencil,&
       phys_req_diagonal)

    ! Initialize gravity module
    if (rhd_gravity) call gravity_init()

    ! Initialize particles module
    if (rhd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    nvector      = 1 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1

    kbmpmua4 = unit_pressure**(-3./4)*unit_density*const_kB/(&
       const_mp*fld_mu)*const_rad_a**(-1.d0/4)

  end subroutine rhd_phys_init

  subroutine rhd_check_params
    use mod_global_parameters
    use mod_dust, only: dust_check_params

    if (.not. rhd_energy .and. rhd_pressure == 'adiabatic') then
       if (rhd_gamma <= 0.0d0) call mpistop ("Error: rhd_gamma <= 0")
       if (rhd_adiab < 0.0d0) call mpistop  ("Error: rhd_adiab < 0")
       small_pressure= rhd_adiab*small_density**rhd_gamma
    elseif (rhd_pressure == 'Tcond') then
      small_pressure = smalldouble
    else
       if (rhd_gamma <= 0.0d0 .or. rhd_gamma == 1.0d0) call mpistop &
          ("Error: rhd_gamma <= 0 or rhd_gamma == 1.0")
       small_e = small_pressure/(rhd_gamma - 1.0d0)
    end if

    small_r_e = small_pressure/(rhd_gamma - 1.0d0)

    if (rhd_dust) call dust_check_params()

    if (rhd_radiation_diffusion .and. .not. use_imex_scheme) call &
       mpistop("Use an IMEX scheme when doing FLD")

    if (use_multigrid) call rhd_set_mg_bounds()

  end subroutine rhd_check_params

  !> Set the boundaries for the diffusion of E
  subroutine rhd_set_mg_bounds
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_usr_methods

    integer :: iB

    ! Set boundary conditions for the multigrid solver
    do iB = 1, 2*ndim
       select case (typeboundary(r_e, iB))
       case ('symm')
          ! d/dx u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! u = 0
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          ! d/dx u = 0
          ! mg%bc(iB, mg_iphi)%bc_type = mg_bc_continuous
          mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(iB, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case ('noinflow')
          call usr_special_mg_bc(iB)
       case ('special')
          call usr_special_mg_bc(iB)
       end select
    end do
  end subroutine rhd_set_mg_bounds

  subroutine rhd_physical_units
    use mod_global_parameters
    double precision :: mp,kB
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
    else
      mp=mp_cgs
      kB=kB_cgs
    end if
    if(unit_velocity==0) then
      !> Set numberdensity, temperature and length
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)&
         *unit_numberdensity*kB*unit_temperature
      unit_velocity=dsqrt(unit_pressure/unit_density)
      unit_time=unit_length/unit_velocity
    else
      !> Set numberdensity, velocity and length
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_time=unit_length/unit_velocity
    end if

      unit_radflux = unit_velocity*unit_pressure
      unit_opacity = one/(unit_density*unit_length)
  end subroutine rhd_physical_units

  !> Returns 0 in argument flag where values are ok
  subroutine rhd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, flag, smallw)
    use mod_global_parameters

    logical, intent(in)          :: primitive
    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    integer, intent(inout)       :: flag(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out) :: smallw(1:nw)
    double precision             :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)

    smallw=1.d0
    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0

    if (rhd_energy) then
       if (primitive) then
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = e_
          if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == e_)) smallw(e_) = &
             minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))
       else
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (rhd_gamma - &
             1.0d0)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) - rhd_kin_en(w,&
              ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
             ixOmax2))
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < small_pressure) &
             flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = e_
          if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == e_)) smallw(e_) = &
             minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))
       endif
    end if

    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == rho_)) smallw(rho_) = &
       minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        r_e) < small_r_e) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = r_e
    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == r_e)) smallw(r_e) = &
       minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e))

  end subroutine rhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine rhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_conserved
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision                :: invgam
    integer                         :: idir, itr

    if (check_small_values .and. small_values_use_primitive) then
      call rhd_handle_small_values(.true., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'rhd_to_conserved')
    end if

    if (rhd_energy) then
       invgam = 1.d0/(rhd_gamma - 1.0d0)
       ! Calculate total energy from pressure and kinetic energy
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, e_) * invgam + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(:))**2, dim=ndim+1) * w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))
    end do

    if (rhd_dust) then
      call dust_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, w, x)
    end if

    if (check_small_values .and. .not. small_values_use_primitive) then
      call rhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'rhd_to_conserved')
    end if
  end subroutine rhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine rhd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    use mod_dust, only: dust_to_primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: itr, idir
    double precision                :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (check_small_values .and. .not. small_values_use_primitive) then
      call rhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'rhd_to_primitive')
    end if

    inv_rho = 1.0d0 / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    if (rhd_energy) then
       ! Compute pressure
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) = (rhd_gamma - 1.0d0) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) - rhd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2, inv_rho))
    end if

    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idir)) * inv_rho
    end do

    ! Convert dust momentum to dust velocity
    if (rhd_dust) then
      call dust_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, w, x)
    end if

    if (check_small_values .and. small_values_use_primitive) then
      call rhd_handle_small_values(.true., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'rhd_to_primitive')
    end if

  end subroutine rhd_to_primitive

  subroutine e_to_rhos(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)

    if (rhd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) = (rhd_gamma - 1.0d0) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)**(1.0d0 - rhd_gamma) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) - rhd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2))
    else
       call mpistop("energy from entropy can not be used with -eos = iso !")
    end if
  end subroutine e_to_rhos

  subroutine rhos_to_e(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision             :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)

    if (rhd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)**(rhd_gamma - 1.0d0) * w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, e_) / (rhd_gamma - 1.0d0) + rhd_kin_en(w, ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    else
       call mpistop("entropy from energy can not be used with -eos = iso !")
    end if
  end subroutine rhos_to_e

  !> Calculate v_i = m_i / rho within ixO^L
  subroutine rhd_get_v(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, idim, v)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2)

    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mom(idim)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
  end subroutine rhd_get_v

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine rhd_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                          :: csound(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision                          :: v(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call rhd_get_v(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, v)
    call rhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,csound)
    csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sqrt(csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))
    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = abs(v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (rhd_dust) then
      call dust_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    end if
  end subroutine rhd_get_cmax

  !> Calculate cmax_idim = csound + abs(v_idim) within ixO^L
  subroutine rhd_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
    use mod_global_parameters
    use mod_dust, only: dust_get_cmax

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative left and right status
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    ! primitive left and right status
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    double precision :: wmean(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: umean,&
        dmean, csoundL, csoundR, tmp1,tmp2,tmp3

    if (typeboundspeed/='cmaxmean') then
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.

      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=1.d0/(sqrt(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_))+sqrt(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

      if(rhd_energy) then
        csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rhd_gamma*wLp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)/wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
        csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rhd_gamma*wRp(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)/wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
      else
        select case (rhd_pressure)
        case ('Trad')
          csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rhd_gamma*kbmpmua4*wLp(&
             ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e)**(1.d0/4)
          csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rhd_gamma*kbmpmua4*wRp(&
             ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e)**(1.d0/4)
        case ('adiabatic')
          csoundL(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rhd_gamma*rhd_adiab*wLp(&
             ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**(rhd_gamma-one)
          csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=rhd_gamma*rhd_adiab*wRp(&
             ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**(rhd_gamma-one)
        end select
      end if

      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundL(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)) * tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) + 0.5d0*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))-wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)))**2

      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(dmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))
      if(present(cmin)) then
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(umean(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if

      if (rhd_dust) then
        wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
        call dust_get_cmax(wmean, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
      end if

    else

      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=wmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)
      call rhd_get_csound2(wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,csoundR)
      csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = sqrt(csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))

      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if

      if (rhd_dust) then
        call dust_get_cmax(wmean, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
      end if
    end if

  end subroutine rhd_get_cbounds

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine rhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,&
       ixImin2:ixImax2)


    if(rhd_energy) then
      call rhd_get_ptot(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(rhd_gamma,&
         4.d0/3.d0)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)
    else
      call rhd_get_ptot(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(rhd_gamma,&
         4.d0/3.d0)*csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)
    end if
  end subroutine rhd_get_csound2

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho) within ixO^L
  subroutine rhd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, pth)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    if (rhd_energy) then
       pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = (rhd_gamma - 1.0d0) * &
          (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) - rhd_kin_en(w, ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2))
    else
      select case (rhd_pressure)
      case ('Trad')
        pth(ixImin1:ixImax1,ixImin2:ixImax2) = (w(ixImin1:ixImax1,&
           ixImin2:ixImax2,r_e)*unit_pressure/const_rad_a)**&
           0.25d0/unit_temperature*w(ixImin1:ixImax1,ixImin2:ixImax2, rho_)
      case ('adiabatic')
        pth(ixImin1:ixImax1,ixImin2:ixImax2) = rhd_adiab * w(ixImin1:ixImax1,&
           ixImin2:ixImax2, rho_)**rhd_gamma
      case ('Tcond')
        pth(ixImin1:ixImax1,ixImin2:ixImax2) = &
           (rhd_gamma-1.d0)*w(ixImin1:ixImax1,ixImin2:ixImax2,r_e)
      case default
        call mpistop('rhd_pressure unknown, use Trad or adiabatic')
      end select
    end if

  end subroutine rhd_get_pthermal

  !> Calculate radiation pressure within ixO^L
  subroutine rhd_get_pradiation(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, prad)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: prad(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        1:ndim, 1:ndim)

    select case (rhd_radiation_formalism)
    case('fld')
      call fld_get_radpress(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, prad)
    case default
      call mpistop('Radiation formalism unknown')
    end select
  end subroutine rhd_get_pradiation

  !> calculates the sum of the gas pressure and max Prad tensor element
  subroutine rhd_get_ptot(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, ptot)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision             :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision             :: prad_tensor(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, 1:ndim, 1:ndim)
    double precision             :: prad_max(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision, intent(out):: ptot(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: ix1,ix2

    call rhd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, pth)
    call rhd_get_pradiation(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, prad_tensor)

    do ix1 = ixOmin1,ixOmax1
    do ix2 = ixOmin2,ixOmax2
      prad_max(ix1,ix2) = maxval(prad_tensor(ix1,ix2,:,:))
    enddo
    enddo

    !> filter cmax
    if (radio_acoustic_filter) then
      call rhd_radio_acoustic_filter(x, ixImin1,ixImin2,ixImax1,ixImax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2, prad_max)
    endif

    ptot(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + prad_max(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

  end subroutine rhd_get_ptot

  !> Filter peaks in cmax due to radiation energy density
  subroutine rhd_radio_acoustic_filter(x, ixImin1,ixImin2,ixImax1,ixImax2,&
      ixOmin1,ixOmin2,ixOmax1,ixOmax2, prad_max)
    use mod_global_parameters

    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)              :: x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout)           :: prad_max(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    double precision :: tmp_prad(ixImin1:ixImax1,ixImin2:ixImax2)
    integer :: ix1,ix2, filter, idim

    if (size_ra_filter .lt. 1) call mpistop(&
       "ra filter of size < 1 makes no sense")
    if (size_ra_filter .gt. nghostcells) call &
       mpistop("ra filter of size < nghostcells makes no sense")

    tmp_prad(ixImin1:ixImax1,ixImin2:ixImax2) = zero
    tmp_prad(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = prad_max(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    do filter = 1,size_ra_filter
      do idim = 1,ndim
        ! {do ix^D = ixOmin^D+filter,ixOmax^D-filter\}
        do ix1 = ixOmin1,ixOmax1
        do ix2 = ixOmin2,ixOmax2
            prad_max(ix1,ix2) = min(tmp_prad(ix1,ix2),&
               tmp_prad(ix1+filter*kr(idim,1),ix2+filter*kr(idim,2)))
            prad_max(ix1,ix2) = min(tmp_prad(ix1,ix2),&
               tmp_prad(ix1-filter*kr(idim,1),ix2-filter*kr(idim,2)))
        enddo
        enddo
      enddo
    enddo
  end subroutine rhd_radio_acoustic_filter

  !> Calculates gas temperature
  subroutine rhd_get_tgas(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, tgas)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision             :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out):: tgas(ixImin1:ixImax1,ixImin2:ixImax2)

    double precision :: mu

    call rhd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, pth)

    tgas(ixImin1:ixImax1,ixImin2:ixImax2) = pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)

  end subroutine rhd_get_tgas

  !> Calculates radiation temperature
  subroutine rhd_get_trad(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, trad)
    use mod_global_parameters
    use mod_constants

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out):: trad(ixImin1:ixImax1,ixImin2:ixImax2)

    trad(ixImin1:ixImax1,ixImin2:ixImax2) = (w(ixImin1:ixImax1,ixImin2:ixImax2,&
       r_e)*unit_pressure/const_rad_a)**(1.d0/4.d0)/unit_temperature

  end subroutine rhd_get_trad

  ! Calculate flux f_idim[iw]
  subroutine rhd_get_flux_cons(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
        v(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir, itr

    call rhd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, pth)
    call rhd_get_v(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idim, v)

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) = v(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir))
    end do

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idim)) = f(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mom(idim)) + pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if(rhd_energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) + pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    end if

    if (rhd_radiation_advection) then
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_e) = v(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e)
    else
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_e) = zero
    endif

    do itr = 1, rhd_n_tracer
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tracer(itr)) = v(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tracer(itr))
    end do

    ! Dust fluxes
    if (rhd_dust) then
      call dust_get_flux(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, f)
    end if

  end subroutine rhd_get_flux_cons

  ! Calculate flux f_idim[iw]
  subroutine rhd_get_flux(wC, w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters
    use mod_dust, only: dust_get_flux_prim
    use mod_viscosity, only: visc_get_flux_prim ! viscInDiv

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative w
    double precision, intent(in)    :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    ! primitive w
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir, itr

    if (rhd_energy) then
       pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,p_)
    else
       call rhd_get_pthermal(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2, pth)
    end if

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_) = w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(idim)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)

    ! Momentum flux is v_i*m_i, +p in direction idim
    do idir = 1, ndir
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim)) * wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mom(idir))
    end do

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idim)) = f(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2, mom(idim)) + pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if(rhd_energy) then
      ! Energy flux is v_i*e + v*p ! Check? m_i/rho*p
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim)) * (wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_))
    end if

    if (rhd_radiation_advection) then
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_e) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim)) * wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e)
    else
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_e) = zero
    endif

    do itr = 1, rhd_n_tracer
       f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, tracer(itr)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(idim)) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           tracer(itr))
    end do

    ! Dust fluxes
    if (rhd_dust) then
      call dust_get_flux_prim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, f)
    end if

    ! Viscosity fluxes - viscInDiv
    if (rhd_viscosity) then
      call visc_get_flux_prim(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2, idim, f, rhd_energy)
    endif

  end subroutine rhd_get_flux

  !> Add geometrical source terms to w
  !>
  !> Notice that the expressions of the geometrical terms depend only on ndir,
  !> not ndim. Eg, they are the same in 2.5D and in 3D, for any geometry.
  !>
  !> Ileyk : to do :
  !>     - address the source term for the dust in case (typeaxial == 'spherical')
  subroutine rhd_add_source_geom(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, wCT, w, x)
    use mod_global_parameters
    use mod_viscosity, only: visc_add_source_geom ! viscInDiv
    use mod_dust, only: dust_n_species, dust_mom, dust_rho, dust_small_to_zero,&
        set_dusttozero, dust_min_rho
    use mod_geometry
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2, 1:nw)
    ! to change and to set as a parameter in the parfile once the possibility to
    ! solve the equations in an angular momentum conserving form has been
    ! implemented (change tvdlf.t eg)
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
        source(ixImin1:ixImax1,ixImin2:ixImax2), minrho
    integer                         :: iw,idir, h1xmin1,h1xmin2,h1xmax1,&
       h1xmax2, h2xmin1,h2xmin2,h2xmax1,h2xmax2
    integer :: mr_,mphi_ ! Polar var. names
    integer :: irho, ifluid, n_fluids

    if (rhd_dust) then
       n_fluids = 1 + dust_n_species
    else
       n_fluids = 1
    end if

    select case (coordinate)
    case (cylindrical)
       do ifluid = 0, n_fluids-1
          ! s[mr]=(pthermal+mphi**2/rho)/radius
          if (ifluid == 0) then
             ! gas
             irho  = rho_
             mr_   = mom(r_)
             if(phi_>0) mphi_ = mom(phi_)
             call rhd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2,&
                 ixOmin1,ixOmin2,ixOmax1,ixOmax2, source)
             minrho = 0.0d0
          else
             ! dust : no pressure
             irho  = dust_rho(ifluid)
             mr_   = dust_mom(r_, ifluid)
             mphi_ = dust_mom(phi_, ifluid)
             source(ixImin1:ixImax1,ixImin2:ixImax2) = zero
             minrho = dust_min_rho
          end if
          if (phi_ > 0) then
             where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho) > minrho)
                source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
                   source(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    mphi_)**2 / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho)
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mr_) = w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2, mr_) + qdt * source(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_)
             end where
             ! s[mphi]=(-mphi*mr/rho)/radius
             if(.not. angmomfix) then
                where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho) > minrho)
                   source(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2) = -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mphi_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mr_) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, irho)
                   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mphi_) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       mphi_) + qdt * source(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                       r_)
                end where
             end if
          else
             ! s[mr]=2pthermal/radius
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mr_) = w(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2, mr_) + qdt * source(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, r_)
          end if
       end do
    case (spherical)
       if (rhd_dust) then
          call mpistop&
("Dust geom source terms not implemented yet with spherical geometries")
       end if
       mr_   = mom(r_)
       if(phi_>0) mphi_ = mom(phi_)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmax1=ixOmax1-kr(1,1);h1xmax2=ixOmax2-kr(1,2)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmax1=ixOmax1-kr(2,1);h2xmax2=ixOmax2-kr(2,2);
       ! s[mr]=((mtheta**2+mphi**2)/rho+2*p)/r
       call rhd_get_pthermal(wCT, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2, pth)
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = pth(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1) - block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,&
           1)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if (ndir > 1) then
         do idir = 2, ndir
           source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) + wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(idir))**2 / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mr_) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mr_) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1)

       
       ! s[mtheta]=-(mr*mtheta/rho)/r+cot(theta)*(mphi**2/rho+p)/r
       source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = pth(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1) * (block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           2) - block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
           2)) / block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if (ndir == 3) then
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) + (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))**2 / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_)) / tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 2))
       end if
       if (.not. angmomfix) then
          source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = source(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) - (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mr_)) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(2)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(2)) + qdt * source(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1)

       if (ndir == 3) then
         ! s[mphi]=-(mphi*mr/rho)/r-cot(theta)*(mtheta*mphi/rho)/r
         if (.not. angmomfix) then
           source(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = -(wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2, mom(3)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mr_)) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_)- (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(2)) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mom(3))) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               rho_) / tan(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 2))
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(3)) = w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2, mom(3)) + qdt * source(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2, 1)
         end if
       end if
      
    end select

    if (rhd_dust .and. dust_small_to_zero) then
       call set_dusttozero(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,  wCT,  w, x)
    end if

    if (rhd_viscosity) call visc_add_source_geom(qdt,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x)

  end subroutine rhd_add_source_geom

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine rhd_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_dust, only: dust_add_source, dust_mom, dust_rho, dust_n_species
    use mod_viscosity, only: viscosity_add_source
    use mod_usr_methods, only: usr_gravity
    use mod_gravity, only: gravity_add_source, grav_split
    use mod_dust, only: dust_small_to_zero, set_dusttozero


    use mod_fld


    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    double precision :: gravity_field(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    integer :: idust, idim
    integer :: i

    if(rhd_dust) then
      call dust_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    end if

    if(rhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    end if

    if(rhd_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,rhd_energy,qsourcesplit,active)
    end if

    if (rhd_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,rhd_energy,qsourcesplit,active)

      if (rhd_dust .and. qsourcesplit .eqv. grav_split) then
         active = .true.

         call usr_gravity(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
            ixOmax1,ixOmax2, wCT, x, gravity_field)
         do idust = 1, dust_n_species
            do idim = 1, ndim
               w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idim,&
                   idust)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, dust_mom(idim,&
                   idust)) + qdt * gravity_field(ixOmin1:ixOmax1,&
                  ixOmin2:ixOmax2, idim) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   dust_rho(idust))
            end do
         end do
         if (dust_small_to_zero) then
            call set_dusttozero(qdt, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
               ixOmin2,ixOmax1,ixOmax2,  wCT,  w, x)
         end if
      end if
    end if
   
    call rhd_add_radiation_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)

  end subroutine rhd_add_source

  subroutine rhd_add_radiation_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    use mod_constants
    use mod_global_parameters
    use mod_usr_methods
    use mod_fld

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active
    double precision :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)

    if (fld_diff_scheme .eq. 'mg') call fld_get_diffcoef_central(w, wCT, x,&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    ! if (fld_diff_scheme .eq. 'mg') call set_mg_bounds(wCT, x, ixI^L, ixO^L)

    select case(rhd_radiation_formalism)
    case('fld')

      !> radiation force
      if (rhd_radiation_force) call get_fld_rad_force(qdt,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,rhd_energy,&
         qsourcesplit,active)

      if (check_small_values .and. small_values_use_primitive) then
        call rhd_handle_small_values(.true., w, x, ixImin1,ixImin2,ixImax1,&
           ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'fld_e_interact')
      end if

      !> photon tiring, heating and cooling
      if (rhd_energy) then
      if (rhd_energy_interact) call get_fld_energy_interact(qdt,ixImin1,&
         ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,&
         rhd_energy,qsourcesplit,active)
      endif

    case default
      call mpistop('Radiation formalism unknown')
    end select

    ! !>  NOT necessary for calculation, just want to know the grid-dependent-timestep
    !call rhd_get_cmax(w, x, ixI^L, ixO^L, 2, cmax)
    !w(ixI^S,i_test) = cmax(ixI^S)

  end subroutine rhd_add_radiation_source

  subroutine rhd_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
    use mod_global_parameters
    use mod_dust, only: dust_get_dt
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt
    use mod_fld, only: fld_radforce_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:2)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(inout) :: dtnew

    dtnew = bigdouble

    if (.not. dt_c) then

      if(rhd_dust) then
        call dust_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
      end if

      if(rhd_radiation_force) then
        call fld_radforce_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
      endif

      if(rhd_radiative_cooling) then
        call cooling_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
      end if

      if(rhd_viscosity) then
        call viscosity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
      end if

      if(rhd_gravity) then
        call gravity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
      end if
    else
      
       dtnew = min(dx1*unit_velocity/const_c,dx2*unit_velocity/const_c)
    endif

  end subroutine rhd_get_dt

  function rhd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)                    :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)           :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw)
    double precision                       :: ke(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
    end if
  end function rhd_kin_en

  function rhd_inv_rho(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2) result(inv_rho)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    ! Can make this more robust
    inv_rho = 1.0d0 / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
  end function rhd_inv_rho

  subroutine rhd_handle_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname)
    use mod_global_parameters
    use mod_small_values
    ! use mod_fld
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: smallw(1:nw)
    integer :: idir, flag(ixImin1:ixImax1,ixImin2:ixImax2)

    if (small_values_method == "ignore") return

    call rhd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, w, flag, smallw)

    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)) then
      select case (small_values_method)
      case ("replace")
        if (small_values_fix_iw(rho_)) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == rho_) &
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_) = small_density
        end if

        ! do idir = 1, ndir
        !   if (small_values_fix_iw(mom(idir))) then
        !     where(flag(ixO^S) == rho_) w(ixO^S, mom(idir)) = 0.0d0
        !   end if
        ! end do

        if (small_values_fix_iw(r_e)) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) == r_e) &
             w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,r_e) = small_r_e
        end if

        if (rhd_energy) then
          if (small_values_fix_iw(e_)) then
            if(primitive) then
              where(flag(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 e_) = small_pressure
            else
              where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)
                ! Add kinetic energy
                w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   e_) = small_e + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2, mom(:))**2,&
                    dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
              end where
            end if
          end if
        end if

      case ("average")
        call small_values_average(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, w, x, flag)
      case default
        call small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2, flag, subname, smallw)
      end select
    end if
  end subroutine rhd_handle_small_values

end module mod_rhd_phys
