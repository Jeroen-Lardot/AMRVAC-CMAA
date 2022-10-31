!> Magneto-hydrodynamics module
module mod_mhd_phys
  use mod_global_parameters, only: std_len
  implicit none
  private

  !> Whether an energy equation is used
  logical, public, protected              :: mhd_energy = .true.

  !> Whether thermal conduction is used
  logical, public, protected              :: mhd_thermal_conduction = .false.

  !> Whether radiative cooling is added
  logical, public, protected              :: mhd_radiative_cooling = .false.

  !> Whether viscosity is added
  logical, public, protected              :: mhd_viscosity = .false.

  !> Whether gravity is added
  logical, public, protected              :: mhd_gravity = .false.

  !> Whether Hall-MHD is used
  logical, public, protected              :: mhd_Hall = .false.

  !> Whether particles module is added
  logical, public, protected              :: mhd_particles = .false.

  !> Whether magnetofriction is added
  logical, public, protected              :: mhd_magnetofriction = .false.

  !> Whether GLM-MHD is used
  logical, public, protected              :: mhd_glm = .false.

  !> Whether auxiliary internal energy is solved
  logical, public, protected              :: mhd_solve_eaux = .false.

  !> Whether internal energy is solved instead of total energy
  logical, public, protected              :: mhd_internal_e = .false.

  !> Whether divB cleaning sources are added splitting from fluid solver
  logical, public, protected              :: source_split_divb = .false.

  !> GLM-MHD parameter: ratio of the diffusive and advective time scales for div b
  !> taking values within [0, 1]
  double precision, public                :: mhd_glm_alpha = 0.5d0

  !> Use Boris approximation
  character(len=20) :: mhd_boris_method = "none"

  integer, parameter :: boris_none           = 0
  integer, parameter :: boris_reduced_force  = 1
  integer, parameter :: boris_simplification = 2
  integer            :: mhd_boris_type       = boris_none

  !> Speed of light for Boris' approximation. If negative, test changes to the
  !> momentum equation with gamma_A = 1
  double precision                        :: mhd_boris_c = 0.0d0

  !> MHD fourth order
  logical, public, protected              :: mhd_4th_order = .false.

  !> Number of tracer species
  integer, public, protected              :: mhd_n_tracer = 0

  !> Index of the density (in the w array)
  integer, public, protected              :: rho_

  !> Indices of the momentum density
  integer, allocatable, public, protected :: mom(:)

  !> Index of the energy density (-1 if not present)
  integer, public, protected              :: e_

  !> Index of the gas pressure (-1 if not present) should equal e_
  integer, public, protected              :: p_

  !> Indices of the magnetic field
  integer, allocatable, public, protected :: mag(:)

  !> Indices of the GLM psi
  integer, public, protected :: psi_

  !> Indices of auxiliary internal energy
  integer, public, protected :: eaux_
  integer, public, protected :: paux_

  !> Indices of the tracers
  integer, allocatable, public, protected :: tracer(:)

  !> The adiabatic index
  double precision, public                :: mhd_gamma = 5.d0/3.0d0

  !> The adiabatic constant
  double precision, public                :: mhd_adiab = 1.0d0

  !> The MHD resistivity
  double precision, public                :: mhd_eta = 0.0d0

  !> The MHD hyper-resistivity
  double precision, public                :: mhd_eta_hyper = 0.0d0

  !> TODO: what is this?
  double precision, public                :: mhd_etah = 0.0d0

  !> The small_est allowed energy
  double precision, protected             :: small_e

  !> The number of waves
  integer :: nwwave=8

  !> Method type to clean divergence of B
  character(len=std_len), public, protected :: typedivbfix  = 'linde'

  !> Method type of constrained transport
  character(len=std_len), public, protected :: type_ct  = 'uct_contact'

  !> Whether divB is computed with a fourth order approximation
  logical, public, protected :: mhd_divb_4thorder = .false.

  !> Method type in a integer for good performance
  integer :: type_divb

  !> Coefficient of diffusive divB cleaning
  double precision :: divbdiff     = 0.8d0

  !> Update all equations due to divB cleaning
  character(len=std_len) ::    typedivbdiff = 'all'

  !> Use a compact way to add resistivity
  logical :: compactres   = .false.

  !> Add divB wave in Roe solver
  logical, public :: divbwave     = .true.

  !> clean initial divB
  logical, public :: clean_initial_divb     = .false.

  !> Helium abundance over Hydrogen
  double precision, public, protected  :: He_abundance=0.1d0

  !> To control divB=0 fix for boundary
  logical, public, protected :: boundary_divbfix(2*2)=.true.

  !> To skip * layer of ghost cells during divB=0 fix for boundary
  integer, public, protected :: boundary_divbfix_skip(2*2)=0

  !> B0 field is force-free
  logical, public, protected :: B0field_forcefree=.true.

  !> Whether an total energy equation is used
  logical :: total_energy = .true.

  !> gamma minus one and its inverse
  double precision :: gamma_1, inv_gamma_1

  ! DivB cleaning methods
  integer, parameter :: divb_none          = 0
  integer, parameter :: divb_multigrid     = -1
  integer, parameter :: divb_glm           = 1
  integer, parameter :: divb_powel         = 2
  integer, parameter :: divb_janhunen      = 3
  integer, parameter :: divb_linde         = 4
  integer, parameter :: divb_lindejanhunen = 5
  integer, parameter :: divb_lindepowel    = 6
  integer, parameter :: divb_lindeglm      = 7
  integer, parameter :: divb_ct            = 8


  ! Public methods
  public :: mhd_phys_init
  public :: mhd_kin_en
  public :: mhd_get_pthermal
  public :: mhd_get_v
  public :: mhd_get_v_idim
  public :: mhd_to_conserved
  public :: mhd_to_primitive
  public :: mhd_get_csound2
  public :: mhd_face_to_center
  public :: get_divb
  public :: get_current
  public :: get_normalized_divb
  public :: b_from_vector_potential
  public :: mhd_mag_en_all
  
  public :: mhd_clean_divb_multigrid
 

contains

  !> Read this module"s parameters from a file
  subroutine mhd_read_params(files)
    use mod_global_parameters
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /mhd_list/ mhd_energy, mhd_n_tracer, mhd_gamma, mhd_adiab,mhd_eta,&
        mhd_eta_hyper, mhd_etah, mhd_glm_alpha, mhd_magnetofriction,&
       mhd_thermal_conduction, mhd_radiative_cooling, mhd_Hall, mhd_gravity,&
       mhd_viscosity, mhd_4th_order, typedivbfix, source_split_divb, divbdiff,&
       typedivbdiff, type_ct, compactres, divbwave, He_abundance, SI_unit,&
        B0field,B0field_forcefree, Bdip, Bquad, Boct, Busr, mhd_particles,&
       boundary_divbfix, boundary_divbfix_skip, mhd_divb_4thorder,&
        mhd_boris_method, mhd_boris_c, clean_initial_divb, mhd_solve_eaux,&
        mhd_internal_e

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, mhd_list, end=111)
111    close(unitpar)
    end do

  end subroutine mhd_read_params

  !> Write this module's parameters to a snapsoht
  subroutine mhd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 1
    double precision                    :: values(n_par)
    character(len=name_len)             :: names(n_par)
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)

    names(1) = "gamma"
    values(1) = mhd_gamma
    call MPI_FILE_WRITE(fh, values, n_par, MPI_DOUBLE_PRECISION, st, er)
    call MPI_FILE_WRITE(fh, names, n_par * name_len, MPI_CHARACTER, st, er)
  end subroutine mhd_write_info

  subroutine mhd_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim)
    use mod_global_parameters
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

    call mpistop("to do")

  end subroutine mhd_angmomfix

  subroutine mhd_phys_init()
    use mod_global_parameters
    use mod_thermal_conduction
    use mod_radiative_cooling
    use mod_viscosity, only: viscosity_init
    use mod_gravity, only: gravity_init
    use mod_particles, only: particles_init
    use mod_magnetofriction, only: magnetofriction_init
    use mod_physics
    
    use mod_multigrid_coupling
   

    integer :: itr, idir

    call mhd_read_params(par_files)

    physics_type = "mhd"
    phys_energy=mhd_energy
    phys_internal_e=mhd_internal_e
    phys_solve_eaux=mhd_solve_eaux

    if(mhd_energy.and..not.mhd_internal_e) then
      total_energy=.true.
    else
      total_energy=.false.
    end if
    phys_total_energy=total_energy

    if(mhd_internal_e.and.mhd_solve_eaux) then
      mhd_solve_eaux=.false.
      if(mype==0) write(*,*)&
          'WARNING: set mhd_solve_eaux=F when mhd_internal_e=T'
    end if

    if(.not. mhd_energy) then
      if(mhd_internal_e) then
        mhd_internal_e=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_internal_e=F when mhd_energy=F'
      end if
      if(mhd_solve_eaux) then
        mhd_solve_eaux=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_solve_eaux=F when mhd_energy=F'
      end if
      if(mhd_thermal_conduction) then
        mhd_thermal_conduction=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_thermal_conduction=F when mhd_energy=F'
      end if
      if(mhd_radiative_cooling) then
        mhd_radiative_cooling=.false.
        if(mype==0) write(*,*)&
            'WARNING: set mhd_radiative_cooling=F when mhd_energy=F'
      end if
    end if

    if(mhd_solve_eaux) prolongprimitive=.true.

    ! set default gamma for polytropic/isothermal process
    if(.not.mhd_energy) mhd_gamma=1.d0
    use_particles=mhd_particles
    if(ndim==1) typedivbfix='none'
    select case (typedivbfix)
    case ('none')
       type_divb = divb_none
    
    case ('multigrid')
       type_divb = divb_multigrid
       use_multigrid = .true.
       mg%operator_type = mg_laplacian
       phys_global_source_after => mhd_clean_divb_multigrid
   
    case ('glm')
      mhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_glm
    case ('powel', 'powell')
      type_divb = divb_powel
    case ('janhunen')
      type_divb = divb_janhunen
    case ('linde')
      type_divb = divb_linde
    case ('lindejanhunen')
      type_divb = divb_lindejanhunen
    case ('lindepowel')
      type_divb = divb_lindepowel
    case ('lindeglm')
      mhd_glm          = .true.
      need_global_cmax = .true.
      type_divb        = divb_lindeglm
    case ('ct')
      type_divb = divb_ct
      stagger_grid = .true.
    case default
      call mpistop('Unknown divB fix')
    end select

    select case (mhd_boris_method)
    case ("none")
      mhd_boris_type = boris_none
    case ("reduced_force")
      mhd_boris_type = boris_reduced_force
    case ("simplification")
      mhd_boris_type = boris_simplification
    case default
      call mpistop&
         ("Unknown mhd_boris_method (none, reduced_force, simplification)")
    end select

    ! Determine flux variables
    rho_ = var_set_rho()

    allocate(mom(ndir))
    mom(:) = var_set_momentum(ndir)

    ! Set index of energy variable
    if (mhd_energy) then
      nwwave = 8
      e_     = var_set_energy() ! energy density
      p_     = e_               ! gas pressure
    else
      nwwave = 7
      e_     = -1
      p_     = -1
    end if

    allocate(mag(ndir))
    mag(:) = var_set_bfield(ndir)

    if (mhd_glm) then
      psi_ = var_set_fluxvar('psi', 'psi', need_bc=.false.)
    else
      psi_ = -1
    end if

    !  set auxiliary internal energy variable
    if(mhd_energy .and. mhd_solve_eaux) then
      eaux_ = var_set_internal_energy()
      paux_ = eaux_
    else
      eaux_ = -1
      paux_ = -1
    end if

    allocate(tracer(mhd_n_tracer))

    ! Set starting index of tracers
    do itr = 1, mhd_n_tracer
      tracer(itr) = var_set_fluxvar("trc", "trp", itr, need_bc=.false.)
    end do

    ! determine number of stagger variables
    if(stagger_grid) nws=ndim

    nvector      = 2 ! No. vector vars
    allocate(iw_vector(nvector))
    iw_vector(1) = mom(1) - 1   ! TODO: why like this?
    iw_vector(2) = mag(1) - 1   ! TODO: why like this?

    ! Check whether custom flux types have been defined
    if (.not. allocated(flux_type)) then
       allocate(flux_type(ndir, nw))
       flux_type = flux_default
    else if (any(shape(flux_type) /= [ndir, nw])) then
       call mpistop("phys_check error: flux_type has wrong shape")
    end if

    if(ndim>1) then
      if(mhd_glm) then
        flux_type(:,psi_)=flux_special
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_special
        end do
      else
        do idir=1,ndir
          flux_type(idir,mag(idir))=flux_tvdlf
        end do
      end if
    end if

    select case (mhd_boris_method)
    case ("none")
      mhd_boris_type = boris_none
    case ("reduced_force")
      mhd_boris_type = boris_reduced_force
    case ("simplification")
      mhd_boris_type = boris_simplification
      do idir = 1, ndir
        phys_iw_methods(mom(idir))%inv_capacity => mhd_gamma2_alfven
      end do
    case default
      call mpistop&
         ("Unknown mhd_boris_method (none, reduced_force, simplification)")
    end select

    phys_get_dt              => mhd_get_dt
    phys_get_cmax            => mhd_get_cmax
    phys_get_a2max           => mhd_get_a2max
    phys_get_tcutoff         => mhd_get_tcutoff
    phys_get_cbounds         => mhd_get_cbounds
    phys_get_flux            => mhd_get_flux
    phys_get_v_idim          => mhd_get_v_idim
    phys_add_source_geom     => mhd_add_source_geom
    phys_add_source          => mhd_add_source
    phys_to_conserved        => mhd_to_conserved
    phys_to_primitive        => mhd_to_primitive
    phys_ei_to_e             => mhd_ei_to_e
    phys_e_to_ei             => mhd_e_to_ei
    phys_check_params        => mhd_check_params
    phys_check_w             => mhd_check_w
    phys_get_pthermal        => mhd_get_pthermal
    phys_write_info          => mhd_write_info
    phys_angmomfix           => mhd_angmomfix
    phys_handle_small_values => mhd_handle_small_values
    phys_energy_synchro      => mhd_energy_synchro

    if(type_divb==divb_glm) then
      phys_modify_wLR => mhd_modify_wLR
    end if

    ! if using ct stagger grid, boundary divb=0 is not done here
    if(stagger_grid) then
      phys_update_faces => mhd_update_faces
      phys_face_to_center => mhd_face_to_center
      phys_modify_wLR => mhd_modify_wLR
    else if(ndim>1) then
      phys_boundary_adjust => mhd_boundary_adjust
    end if

    ! Whether diagonal ghost cells are required for the physics
    if(type_divb < divb_linde) phys_req_diagonal = .false.

    ! derive units from basic units
    call mhd_physical_units()

    if(.not. mhd_energy .and. mhd_thermal_conduction) then
      call mpistop("thermal conduction needs mhd_energy=T")
    end if
    if(.not. mhd_energy .and. mhd_radiative_cooling) then
      call mpistop("radiative cooling needs mhd_energy=T")
    end if

    ! initialize thermal conduction module
    if (mhd_thermal_conduction) then
      phys_req_diagonal = .true.
      call thermal_conduction_init(mhd_gamma)
    end if

    ! Initialize radiative cooling module
    if (mhd_radiative_cooling) then
      call radiative_cooling_init(mhd_gamma,He_abundance)
    end if

    ! Initialize viscosity module
    if (mhd_viscosity) call viscosity_init(phys_wider_stencil,&
       phys_req_diagonal)

    ! Initialize gravity module
    if(mhd_gravity) then
      call gravity_init()
    end if

    ! Initialize particles module
    if(mhd_particles) then
      call particles_init()
      phys_req_diagonal = .true.
    end if

    ! initialize magnetofriction module
    if(mhd_magnetofriction) then
      phys_req_diagonal = .true.
      call magnetofriction_init()
    end if

    ! For Hall, we need one more reconstructed layer since currents are computed
    ! in getflux: assuming one additional ghost layer (two for FOURTHORDER) was
    ! added in nghostcells.
    if (mhd_hall) then
       phys_req_diagonal = .true.
       if (mhd_4th_order) then
          phys_wider_stencil = 2
       else
          phys_wider_stencil = 1
       end if
    end if

  end subroutine mhd_phys_init

  subroutine mhd_check_params
    use mod_global_parameters

    ! after user parameter setting
    gamma_1=mhd_gamma-1.d0
    if (.not. mhd_energy) then
       if (mhd_gamma <= 0.0d0) call mpistop ("Error: mhd_gamma <= 0")
       if (mhd_adiab < 0.0d0) call mpistop ("Error: mhd_adiab < 0")
       small_pressure = mhd_adiab*small_density**mhd_gamma
    else
       if (mhd_gamma <= 0.0d0 .or. mhd_gamma == 1.0d0) call mpistop &
          ("Error: mhd_gamma <= 0 or mhd_gamma == 1")
       inv_gamma_1=1.d0/gamma_1
       small_e = small_pressure * inv_gamma_1
    end if

    if (mhd_boris_type > 0 .and. abs(mhd_boris_c) <= 0.0d0) then
      call mpistop("You have not specified mhd_boris_c")
    end if

  end subroutine mhd_check_params

  subroutine mhd_physical_units()
    use mod_global_parameters
    double precision :: mp,kB,miu0,c_lightspeed
    ! Derive scaling units
    if(SI_unit) then
      mp=mp_SI
      kB=kB_SI
      miu0=miu0_SI
      c_lightspeed=c_SI
    else
      mp=mp_cgs
      kB=kB_cgs
      miu0=4.d0*dpi
      c_lightspeed=const_c
    end if
    if(unit_velocity==0) then
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=(2.d0+3.d0*He_abundance)&
         *unit_numberdensity*kB*unit_temperature
      unit_velocity=sqrt(unit_pressure/unit_density)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    else
      unit_density=(1.d0+4.d0*He_abundance)*mp*unit_numberdensity
      unit_pressure=unit_density*unit_velocity**2
      unit_temperature=unit_pressure/((2.d0+&
         3.d0*He_abundance)*unit_numberdensity*kB)
      unit_magneticfield=sqrt(miu0*unit_pressure)
      unit_time=unit_length/unit_velocity
    end if
    ! Additional units needed for the particles
    c_norm=c_lightspeed/unit_velocity
    unit_charge=unit_magneticfield*unit_length**2/unit_velocity/miu0
    if (.not. SI_unit) unit_charge = unit_charge*const_c
    unit_mass=unit_density*unit_length**3

  end subroutine mhd_physical_units

  subroutine mhd_check_w(primitive,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w,flag,smallw)
    use mod_global_parameters

    logical, intent(in) :: primitive
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    integer, intent(inout) :: flag(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out) :: smallw(1:nw)
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2)

    smallw=1.d0
    flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0
    where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        rho_) < small_density) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = rho_
    if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==rho_)) &
       smallw(rho_)=minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    if(mhd_energy) then
      if(primitive) then
        where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) < small_pressure) flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = e_
        if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==e_)) &
           smallw(e_)=minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))
      else
        if(mhd_internal_e) then
          where(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              e_) < small_pressure*inv_gamma_1) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = e_
          if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==e_)) &
             smallw(e_)=minval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_))
        else
          ! Calculate pressure=(gamma-1)*(e-0.5*(2ek+2eb))
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,&
             ixOmin1,ixOmin2,ixOmax1,ixOmax2)-mhd_mag_en(w,ixImin1,ixImin2,&
             ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2)
          where(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2) < &
             small_pressure*inv_gamma_1) flag(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) = e_
          if(any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)==e_)) &
             smallw(e_)=minval(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
        end if
      end if
    end if

  end subroutine mhd_check_w

  !> Transform primitive variables into conservative ones
  subroutine mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    integer                         :: idir, itr

    if (check_small_values .and. small_values_use_primitive) then
      call mhd_handle_small_values(.true., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'mhd_to_conserved')
    end if

    ! Calculate total energy from pressure, kinetic and magnetic energy
    if(mhd_energy) then
      if(mhd_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)*inv_gamma_1
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,p_)*inv_gamma_1+half*sum(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(:))**2,dim=ndim+1)*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,rho_)+mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2,&
            ixOmin1,ixOmin2,ixOmax1,ixOmax2)
        if(mhd_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,paux_)*inv_gamma_1
      end if
    end if

    ! Convert velocity to momentum
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idir))
    end do

    if (check_small_values .and. .not. small_values_use_primitive) then
      call mhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'mhd_to_conserved')
    end if
  end subroutine mhd_to_conserved

  !> Transform conservative variables into primitive ones
  subroutine mhd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)
    double precision                :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    integer                         :: itr, idir

    if (check_small_values .and. .not. small_values_use_primitive) then
      call mhd_handle_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'mhd_to_primitive')
    end if

    inv_rho=1.0d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    ! Calculate pressure = (gamma-1) * (e-ek-eb)
    if(mhd_energy) then
      if(mhd_internal_e) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)*gamma_1
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)=gamma_1*(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2,inv_rho)-mhd_mag_en(w,ixImin1,&
           ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))
        if(mhd_solve_eaux) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           paux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)*gamma_1
      end if
    end if
    
    ! Convert momentum to velocity
    do idir = 1, ndir
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, mom(idir))*inv_rho
    end do

    if (check_small_values .and. small_values_use_primitive) then
      call mhd_handle_small_values(.true., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'mhd_to_primitive')
    end if
  end subroutine mhd_to_primitive

  !> Transform internal energy to total energy
  subroutine mhd_ei_to_e(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    ! Calculate total energy from internal, kinetic and magnetic energy
    if(total_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2)+mhd_mag_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
    end if

  end subroutine mhd_ei_to_e

  !> Transform total energy to internal energy
  subroutine mhd_e_to_ei(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:ndim)

    ! Calculate ei = e - ek - eb
    if(total_energy) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2)-mhd_mag_en(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
    end if

  end subroutine mhd_e_to_ei

  subroutine mhd_energy_synchro(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: pth1(ixImin1:ixImax1,ixImin2:ixImax2),&
       pth2(ixImin1:ixImax1,ixImin2:ixImax2),alfa(ixImin1:ixImax1,&
       ixImin2:ixImax2),beta(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, parameter :: beta_low=0.005d0,beta_high=0.05d0

    ! add the source of internal energy equation
    call internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,eaux_)
!    double precision :: vtot(ixI^S),cs2(ixI^S),mach(ixI^S)
!    double precision, parameter :: mach_low=20.d0,mach_high=200.d0

    ! get magnetic energy
    alfa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_mag_en(w,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2)
    pth1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gamma_1*(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,e_)-mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2)-alfa(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    pth2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       eaux_)*gamma_1
    ! get plasma beta
    beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(pth1(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2),pth2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))/alfa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    ! whether Mach number should be another criterion ?
!    vtot(ixO^S)=sum(w(ixO^S,mom(:))**2,dim=ndim+1)
!    call mhd_get_csound2(w,x,ixI^L,ixO^L,cs2)
!    mach(ixO^S)=sqrt(vtot(ixO^S)/cs2(ixO^S))/w(ixO^S,rho_)
    where(beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .ge. beta_high)
!    where(beta(ixO^S) .ge. beta_high .and. mach(ixO^S) .le. mach_low)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)=pth1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*inv_gamma_1
    else where(beta(ixOmin1:ixOmax1,ixOmin2:ixOmax2) .le. beta_low)
!    else where(beta(ixO^S) .le. beta_low .or. mach(ixO^S) .ge. mach_high)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)-pth1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*inv_gamma_1+&
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)
    else where
      alfa(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dlog(beta(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/beta_low)/dlog(beta_high/beta_low)
!      alfa(ixO^S)=min(dlog(beta(ixO^S)/beta_low)/dlog(beta_high/beta_low),
!                      dlog(mach_high(ixO^S)/mach(ixO^S))/dlog(mach_high/mach_low))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)=(pth2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*(one-alfa(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))+pth1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*alfa(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*inv_gamma_1
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)-pth1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*inv_gamma_1+&
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)
    end where
  end subroutine mhd_energy_synchro

  subroutine mhd_handle_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname)
    use mod_global_parameters
    use mod_small_values
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

    call mhd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, w, flag, smallw)

    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)) then
      select case (small_values_method)
      case ("replace")
        if (small_values_fix_iw(rho_)) then
          where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,rho_) = small_density
        end if

        do idir = 1, ndir
          if (small_values_fix_iw(mom(idir))) then
            where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0) &
               w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(idir)) = 0.0d0
          end if
        end do

        if (mhd_energy) then
          if (small_values_fix_iw(e_)) then
            if(mhd_solve_eaux) then
              if(primitive) then
                where(flag(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   p_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,paux_)
              else
                where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)
                  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2,eaux_) + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2, mom(:))**2,&
                      dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      rho_) + mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2,&
                      ixOmin1,ixOmin2,ixOmax1,ixOmax2)
                end where
              end if
            else
              if(primitive) then
                where(flag(ixOmin1:ixOmax1,&
                   ixOmin2:ixOmax2) /= 0) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                   p_) = small_pressure
              else
                where(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)
                  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                     e_) = small_e + 0.5d0 * sum(w(ixOmin1:ixOmax1,&
                     ixOmin2:ixOmax2, mom(:))**2,&
                      dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      rho_) + mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2,&
                      ixOmin1,ixOmin2,ixOmax1,ixOmax2)
                end where
              end if
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
  end subroutine mhd_handle_small_values

  !> Convert energy to entropy
  subroutine e_to_rhos(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision,intent(inout)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if(mhd_energy) then
      if(.not.mhd_internal_e) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_)-mhd_kin_en(w, ixImin1,&
         ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
         ixOmax2) -mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_)=gamma_1*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, rho_)**(1.0d0 - mhd_gamma)*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, e_)    
    else
      call mpistop("e_to_rhos can not be used without energy equation!")
    end if
  end subroutine e_to_rhos

  !> Convert entropy to energy
  subroutine rhos_to_e(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,w,x)
    use mod_global_parameters
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    if(mhd_energy) then
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_)=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, rho_)**gamma_1 * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_) * inv_gamma_1
       if(.not.mhd_internal_e) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, e_)+mhd_kin_en(w, ixImin1,&
          ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
          ixOmax2) + mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2)
    else
       call mpistop("rhos_to_e can not be used without energy equation!")
    end if
  end subroutine rhos_to_e

  !> Calculate v vector
  subroutine mhd_get_v(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2,ndir)

    integer :: idir

    do idir=1,ndir
      v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, mom(idir)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    end do

  end subroutine mhd_get_v

  !> Calculate v component
  subroutine mhd_get_v_idim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,v)
    use mod_global_parameters

    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out) :: v(ixImin1:ixImax1,ixImin2:ixImax2)

    v(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mom(idim)) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

  end subroutine mhd_get_v_idim

  !> Calculate cmax_idim=csound+abs(v_idim) within ixO^L
  subroutine mhd_get_cmax(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idim,cmax)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_get_csound(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idim,cmax)

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mom(idim))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_))+cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

  end subroutine mhd_get_cmax

  subroutine mhd_get_a2max(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,a2max)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: a2max(ndim)
    double precision :: a2(ixImin1:ixImax1,ixImin2:ixImax2,ndim,nw)
    integer :: gxOmin1,gxOmin2,gxOmax1,gxOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
       jxOmin1,jxOmin2,jxOmax1,jxOmax2,kxOmin1,kxOmin2,kxOmax1,kxOmax2,i,j

    a2=zero
    do i = 1,ndim
      !> 4th order
      hxOmin1=ixOmin1-kr(i,1);hxOmin2=ixOmin2-kr(i,2);hxOmax1=ixOmax1-kr(i,1)
      hxOmax2=ixOmax2-kr(i,2);
      gxOmin1=hxOmin1-kr(i,1);gxOmin2=hxOmin2-kr(i,2);gxOmax1=hxOmax1-kr(i,1)
      gxOmax2=hxOmax2-kr(i,2);
      jxOmin1=ixOmin1+kr(i,1);jxOmin2=ixOmin2+kr(i,2);jxOmax1=ixOmax1+kr(i,1)
      jxOmax2=ixOmax2+kr(i,2);
      kxOmin1=jxOmin1+kr(i,1);kxOmin2=jxOmin2+kr(i,2);kxOmax1=jxOmax1+kr(i,1)
      kxOmax2=jxOmax2+kr(i,2);
      a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i,1:nw)=abs(-w(kxOmin1:kxOmax1,&
         kxOmin2:kxOmax2,1:nw)+16.d0*w(jxOmin1:jxOmax1,jxOmin2:jxOmax2,&
         1:nw)-30.d0*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nw)+16.d0*w(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nw)-w(gxOmin1:gxOmax1,&
         gxOmin2:gxOmax2,1:nw))
      a2max(i)=maxval(a2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i,&
         1:nw))/12.d0/dxlevel(i)**2
    end do
  end subroutine mhd_get_a2max

  !> get adaptive cutoff temperature for TRAC (Johnston 2019 ApJL, 873, L22)
  subroutine mhd_get_tcutoff(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,tco_local,Tmax_local)
    use mod_global_parameters
    use mod_geometry
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
       w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(out) :: tco_local, Tmax_local

    double precision, parameter :: delta=0.5d0
    double precision :: tmp1(ixImin1:ixImax1,ixImin2:ixImax2),&
       Te(ixImin1:ixImax1,ixImin2:ixImax2),lts(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir) :: bunitvec
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim) :: gradT
    integer :: idims
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2)

    if(mhd_internal_e) then
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)
    else
      tmp1(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
         e_)-0.5d0*(sum(w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mom(:))**2,&
         dim=ndim+1)/w(ixImin1:ixImax1,ixImin2:ixImax2,&
         rho_)+sum(w(ixImin1:ixImax1,ixImin2:ixImax2,iw_mag(:))**2,&
         dim=ndim+1))
    end if
    Te(ixImin1:ixImax1,ixImin2:ixImax2)=tmp1(ixImin1:ixImax1,&
       ixImin2:ixImax2)/w(ixImin1:ixImax1,ixImin2:ixImax2,&
       rho_)*(mhd_gamma-1.d0)

    Tmax_local=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2))

    ! temperature gradient at cell centers
    do idims=1,ndim
      call gradient(Te,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
         ixOmax2,idims,tmp1)
      gradT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)=tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    end do
    ! B vector
    if(B0field) then
      bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw_mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,&
         0)
    else
      bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,iw_mag(:))
    end if
    ! |B|
    tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsqrt(sum(bunitvec(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,:)**2,dim=ndim+1))
    where(tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)/=0.d0)
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=1.d0/tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    elsewhere
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=bigdouble
    end where
    ! b unit vector: magnetic field direction vector
    do idims=1,ndim
      bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idims)=bunitvec(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idims)*tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do
    ! temperature length scale inversed
    lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(sum(gradT(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,1:ndim)*bunitvec(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       1:ndim),dim=ndim+1))/Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    ! fraction of cells size to temperature length scale
    lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=minval(dxlevel)*lts(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    lrlt=.false.
    where(lts(ixOmin1:ixOmax1,ixOmin2:ixOmax2) > delta)
      lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=.true.
    end where
    block%special_values(1)=0.d0
    if(any(lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2))) then
      block%special_values(1)=maxval(Te(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
          mask=lrlt(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
    end if

  end subroutine mhd_get_tcutoff

  !> Estimating bounds for the minimum and maximum signal velocities
  subroutine mhd_get_cbounds(wLC,wRC,wLp,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,cmax,cmin)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2, nw)
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
    integer                            :: idimE,idimN

    if (typeboundspeed=='cmaxmean') then
      wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)=0.5d0*(wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:nwflux)+wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux))
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=wmean(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))/wmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)
      call mhd_get_csound(wmean,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,csoundR)
      if(present(cmin)) then
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=max(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
        cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=min(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
      else
        cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(tmp1(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2))+csoundR(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
    else
      ! This implements formula (10.52) from "Riemann Solvers and Numerical
      ! Methods for Fluid Dynamics" by Toro.
      tmp1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(abs(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)))
      tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(abs(wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)))
      tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=1.d0/(tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))
      umean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idim))*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))*tmp3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      call mhd_get_csound_prim(wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,csoundL)
      call mhd_get_csound_prim(wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idim,csoundR)
      dmean(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundL(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2+tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*csoundR(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)+0.5d0*tmp1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp2(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*tmp3(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)**2*(wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
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
    end if

    if(stagger_grid) then
      ! calculate velocities related to different UCT schemes
      select case(type_ct)
      case('average')
      case('uct_contact')
        if(.not.allocated(vcts%vnorm)) allocate(vcts%vnorm(ixImin1:ixImax1,&
           ixImin2:ixImax2,1:ndim))
        ! get average normal velocity at cell faces
        vcts%vnorm(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)=0.5d0*(wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim))+wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idim)))
      case('uct_hll')
        if(.not.allocated(vcts%vbarC)) then
          allocate(vcts%vbarC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir,2),&
             vcts%vbarLC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir,2),&
             vcts%vbarRC(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir,2))
          allocate(vcts%cbarmin(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
             vcts%cbarmax(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim))
        end if
        ! Store magnitude of characteristics
        if(present(cmin)) then
          vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idim)=max(-cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
          vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idim)=max( cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
        else
          vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idim)=max( cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2),zero)
          vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idim)=vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)
        end if

        idimN=mod(idim,ndir)+1 ! 'Next' direction
        idimE=mod(idim+1,ndir)+1 ! Electric field direction
        ! Store velocities
        vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           1)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idimN))
        vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           1)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idimN))
        vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           1)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           1) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim))

        vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           2)=wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idimE))
        vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           2)=wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idimE))
        vcts%vbarC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           2)=(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)*vcts%vbarLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           2) +vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)*vcts%vbarRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
           1))/(vcts%cbarmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idim)+vcts%cbarmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim))
      case default
        call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
      end select
    end if

  end subroutine mhd_get_cbounds

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2),&
        AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2), b2(ixImin1:ixImax1,&
       ixImin2:ixImax2), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        gamma2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    if (mhd_boris_type == boris_reduced_force) then
      call mhd_gamma2_alfven(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, w, gamma2)
    else
      gamma2 = 1.0d0
    end if

    call mhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,csound)

    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = mhd_mag_en_all(w,ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2) * gamma2

    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = b2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cfast2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2-4.0d0*csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * mhd_mag_i_all(w,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim)**2 * inv_rho * gamma2

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(AvMinCs2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    if (.not. MHD_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
       if (mhd_boris_type == boris_simplification) then
          csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = mhd_gamma_alfven(w,&
              ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
             ixOmax2) * csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))),&
           mhd_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound

  !> Calculate fast magnetosonic wave speed
  subroutine mhd_get_csound_prim(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,csound)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: csound(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: cfast2(ixImin1:ixImax1,ixImin2:ixImax2),&
        AvMinCs2(ixImin1:ixImax1,ixImin2:ixImax2), b2(ixImin1:ixImax1,&
       ixImin2:ixImax2), kmax
    double precision :: inv_rho(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    inv_rho=1.d0/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    if (mhd_boris_type == boris_reduced_force) then
      call mhd_gamma2_alfven(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, w, gamma_A2)
    else
      gamma_A2 = 1.0d0
    end if

    if(mhd_energy) then
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,p_)*inv_rho
    else
      csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*mhd_adiab*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**gamma_1
    end if
    ! store |B|^2 in v
    b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)        = mhd_mag_en_all(w,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2) * gamma_A2
    cfast2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = b2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * inv_rho+csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = cfast2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)**2-4.0d0*csound(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) * mhd_mag_i_all(w,ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim)**2 * inv_rho * gamma_A2

    where(AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)<zero)
       AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
    end where

    AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(AvMinCs2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2))

    if (.not. MHD_Hall) then
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
       if (mhd_boris_type == boris_simplification) then
          csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = mhd_gamma_alfven(w,&
              ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
             ixOmax2) * csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    else
       ! take the Hall velocity into account:
       ! most simple estimate, high k limit:
       ! largest wavenumber supported by grid: Nyquist (in practise can reduce by some factor)
       kmax = dpi/min(dxlevel(1),dxlevel(2),bigdouble)*half
       csound(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = &
          max(sqrt(half*(cfast2(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)+AvMinCs2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))),&
           mhd_etah * sqrt(b2(ixOmin1:ixOmax1,ixOmin2:ixOmax2))*inv_rho*kmax)
    end if

  end subroutine mhd_get_csound_prim

  !> Calculate thermal pressure=(gamma-1)*(e-0.5*m**2/rho-b**2/2) within ixO^L
  subroutine mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,pth)
    use mod_global_parameters

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(out):: pth(ixImin1:ixImax1,ixImin2:ixImax2)

    if (check_small_values) then
      call mhd_find_small_values(.false., w, x, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, 'mhd_get_pthermal')
    end if

    if(mhd_energy) then
      if(mhd_internal_e) then
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gamma_1*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)
      else
        pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=gamma_1*(w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)- mhd_kin_en(w,ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2)- mhd_mag_en(w,ixImin1,ixImin2,&
           ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))
      end if
    else
      pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_adiab*w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,rho_)**mhd_gamma
    end if
  end subroutine mhd_get_pthermal

  subroutine mhd_find_small_values(primitive, w, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, subname)
    use mod_global_parameters
    use mod_small_values
    logical, intent(in)             :: primitive
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    character(len=*), intent(in)    :: subname

    double precision :: smallw(1:nw)
    integer :: idir, flag(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_check_w(primitive, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, w, flag, smallw)

    if (any(flag(ixOmin1:ixOmax1,ixOmin2:ixOmax2) /= 0)) call &
       small_values_error(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, flag, subname, smallw)

  end subroutine mhd_find_small_values

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  !> csound2=gamma*p/rho
  subroutine mhd_get_csound2(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,csound2)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: csound2(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    if(mhd_energy) then
      call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,csound2)
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*csound2(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)
    else
      csound2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=mhd_gamma*mhd_adiab*w(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**gamma_1
    end if
  end subroutine mhd_get_csound2

  !> Calculate total pressure within ixO^L including magnetic pressure
  subroutine mhd_get_p_total(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,p)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: p(ixImin1:ixImax1,ixImin2:ixImax2)

    call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,p)

    p(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = p(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2) + 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
        mag(:))**2, dim=ndim+1)

  end subroutine mhd_get_p_total

  !> Calculate fluxes within ixO^L.
  subroutine mhd_get_flux(wC,w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idim,f)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)          :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idim
    ! conservative w
    double precision, intent(in) :: wC(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    ! primitive w
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision,intent(out) :: f(ixImin1:ixImax1,ixImin2:ixImax2,nwflux)

    double precision             :: pgas(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2),tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision, allocatable:: vHall(:,:,:)
    integer                      :: idirmin, iw, idir, jdir, kdir

    if (mhd_Hall) then
      allocate(vHall(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir))
      call mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,vHall)
    end if

    if(B0field) tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum(block%B0(&
       ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,idim)*w(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,mag(:)),dim=ndim+1)

    if(mhd_energy) then
      pgas=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)
    else
      pgas=mhd_adiab*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)**mhd_gamma
    end if

    ptotal = pgas + 0.5d0*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2,&
        dim=ndim+1)

    ! Get flux of density
    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)

    ! Get flux of tracer
    do iw=1,mhd_n_tracer
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,tracer(iw))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         tracer(iw))
    end do

    ! Get flux of momentum
    ! f_i[m_k]=v_i*m_k-b_k*b_i [+ptotal if i==k]
    if (mhd_boris_type == boris_reduced_force) then
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) = pgas(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) = 0.0d0
        end if
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=f(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))
      end do
    else
      ! Normal case (no Boris approximation)
      do idir=1,ndir
        if(idim==idir) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=ptotal(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))
          if(B0field) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idir))=f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idir))+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
        else
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))= -w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mag(idir))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(idim))
        end if
        if (B0field) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=f(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mom(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(idir))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
             idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,idim)
        end if
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=f(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))
      end do
    end if

    ! Get flux of energy
    ! f_i[e]=v_i*e+v_i*ptotal-b_i*(b_k*v_k)
    if(mhd_energy) then
      if (mhd_internal_e) then
         f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,p_)
         if (mhd_Hall) then
            call mpistop("solve internal energy not implemented for Hall MHD")
         endif
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idim))*(wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           e_)+ptotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2))-w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mag(idim))*sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:)),dim=ndim+1)
        if(mhd_solve_eaux) f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           eaux_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mom(idim))*wC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)

        if (B0field) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,e_) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(idim)) * tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) - sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(:))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),&
              dim=ndim+1) * block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
              idim)
        end if

        if (mhd_Hall) then
        ! f_i[e]= f_i[e] + vHall_i*(b_k*b_k) - b_i*(vHall_k*b_k)
           if (mhd_etah>zero) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,e_) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idim) * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2,&
                 dim=ndim+1) - w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idim)) * sum(vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 :)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),dim=ndim+1)
              if (B0field) then
                 f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = f(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,e_) + vHall(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,idim) * tmp(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2) - sum(vHall(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,:)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                    mag(:)),dim=ndim+1) * block%B0(ixOmin1:ixOmax1,&
                    ixOmin2:ixOmax2,idim,idim)
              end if
           end if
        end if
      end if
    end if

    ! compute flux of magnetic field
    ! f_i[b_k]=v_i*b_k-v_k*b_i
    do idir=1,ndir
      if (idim==idir) then
        ! f_i[b_i] should be exactly 0, so we do not use the transport flux
        if (mhd_glm) then
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,psi_)
        else
           f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=zero
        end if
      else
        f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idir))-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(idim))*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))

        if (B0field) then
          f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=f(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,mag(idir))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idim))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
             idim)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mom(idir))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,idim)
        end if

        if (mhd_Hall) then
          ! f_i[b_k] = f_i[b_k] + vHall_i*b_k - vHall_k*b_i
          if (mhd_etah>zero) then
            if (B0field) then
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = f(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,mag(idir)) - vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,idir)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idim))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim,&
                 idim)) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idim)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,&
                 idim))
            else
              f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = f(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,mag(idir)) - vHall(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,idir)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(idim)) + vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))
            end if
          end if
        end if

      end if
    end do

    if (mhd_glm) then
      !f_i[psi]=Ch^2*b_{i} Eq. 24e and Eq. 38c Dedner et al 2002 JCP, 175, 645
      f(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)  = &
         cmax_global**2*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim))
    end if

  end subroutine mhd_get_flux

  !> w[iws]=w[iws]+qdt*S[iws,wCT] where S is the source based on wCT within ixO
  subroutine mhd_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    use mod_radiative_cooling, only: radiative_cooling_add_source
    use mod_viscosity, only: viscosity_add_source
    use mod_gravity, only: gravity_add_source

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)            :: active

    if (.not. qsourcesplit) then
      ! Source for solving internal energy
      if(mhd_internal_e) then
        active = .true.
        call internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
           ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,e_)
      endif

      ! Source for B0 splitting
      if (B0field) then
        active = .true.
        call add_source_B0split(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      end if

      ! Sources for resistivity in eqs. for e, B1, B2 and B3
      if (abs(mhd_eta)>smalldouble)then
        active = .true.
        call add_source_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      end if

      if (mhd_eta_hyper>0.d0)then
        active = .true.
        call add_source_hyperres(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      end if
    end if

      
    if(.not.source_split_divb .and. .not.qsourcesplit .and. istep==nstep) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
        call add_source_janhunen(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
        call add_source_powel(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
        call add_source_glm(dt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,pso(saveigrid)%w,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    else if(source_split_divb .and. qsourcesplit) then
      ! Sources related to div B
      select case (type_divb)
      case (divb_none)
        ! Do nothing
      case (divb_glm)
        active = .true.
        call add_source_glm(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_powel)
        active = .true.
        call add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_janhunen)
        active = .true.
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_linde)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_lindejanhunen)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        call add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_lindepowel)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        call add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_lindeglm)
        active = .true.
        call add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
        call add_source_glm(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
           ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
      case (divb_ct)
        continue ! Do nothing
      case (divb_multigrid)
        continue ! Do nothing
      case default
        call mpistop('Unknown divB fix')
      end select
    end if
   

    if(mhd_radiative_cooling) then
      call radiative_cooling_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
         ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    end if

    if(mhd_viscosity) then
      call viscosity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,mhd_energy,qsourcesplit,active)
    end if

    if(mhd_gravity) then
      call gravity_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,total_energy,qsourcesplit,active)
    end if

    if (mhd_boris_type == boris_reduced_force) then
      call boris_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    end if

  end subroutine mhd_add_source

  subroutine boris_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    logical, intent(in) :: qsourcesplit
    logical, intent(inout) :: active

    double precision :: JxB(ixImin1:ixImax1,ixImin2:ixImax2,3)
    double precision :: gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer          :: idir

    ! Boris source term is always unsplit
    if (qsourcesplit) return

    call get_lorentz(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,w,JxB)
    call mhd_gamma2_alfven(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, wCT, gamma_A2)

    do idir = 1, ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir)) + qdt * gamma_A2 * JxB(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, idir)
    end do

  end subroutine boris_add_source

  !> Compute the Lorentz force (JxB)
  subroutine get_lorentz(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,JxB)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: JxB(ixImin1:ixImax1,ixImin2:ixImax2,3)
    double precision                :: a(ixImin1:ixImax1,ixImin2:ixImax2,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,3)
    integer                         :: idir, idirmin
    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3)

    b=0.0d0
    do idir = 1, ndir
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2, idir) = mhd_mag_i_all(w, ixImin1,&
         ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2,idir)
    end do

    ! store J current in a
    call get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,idirmin,current)

    a=0.0d0
    do idir=7-2*ndir,3
      a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=current(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
    end do

    call cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,a,b,JxB)
  end subroutine get_lorentz

  !> Compute 1/(1+v_A^2/c^2) for Boris' approximation, where v_A is the Alfven
  !> velocity
  subroutine mhd_gamma2_alfven(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, w, gamma_A2)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(out) :: gamma_A2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (mhd_boris_c < 0.0d0) then
      ! Good for testing the non-conservative momentum treatment
      gamma_A2 = 1.0d0
    else
      ! Compute the inverse of 1 + B^2/(rho * c^2)
      gamma_A2 = 1.0d0 / (1.0d0 + mhd_mag_en_all(w, ixImin1,ixImin2,ixImax1,&
         ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2) / (w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2, rho_) * mhd_boris_c**2))
    end if
  end subroutine mhd_gamma2_alfven

  !> Compute 1/sqrt(1+v_A^2/c^2) for Boris simplification, where v_A is the
  !> Alfven velocity
  function mhd_gamma_alfven(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2) result(gamma_A)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: gamma_A(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    call mhd_gamma2_alfven(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, w, gamma_A)
    gamma_A = sqrt(gamma_A)
  end function mhd_gamma_alfven

  subroutine internal_energy_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,w,x,ie)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2,ie
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),divv(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,v)
    if(slab_uniform) then
      if(nghostcells .gt. 2) then
        call divvector(v,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,divv,sixthorder=.true.)
      else
        call divvector(v,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,divv,fourthorder=.true.)
      end if
    else
     call divvector(v,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
        ixOmax2,divv)
    end if
    call mhd_get_pthermal(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2,pth)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ie)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       ie)-qdt*pth(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*divv(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
  end subroutine internal_energy_add_source

  !> Source terms after split off time-independent magnetic field
  subroutine add_source_B0split(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: a(ixImin1:ixImax1,ixImin2:ixImax2,3),&
        b(ixImin1:ixImax1,ixImin2:ixImax2,3), axb(ixImin1:ixImax1,&
       ixImin2:ixImax2,3)
    integer :: idir

    a=0.d0
    b=0.d0
    ! for force-free field J0xB0 =0
    if(.not.B0field_forcefree) then
      ! store B0 magnetic field in b
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir)=block%B0(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:ndir,0)

      ! store J0 current in a
      do idir=7-2*ndir,3
        a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=block%J0(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,idir)
      end do
      call cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,a,b,axb)
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=axb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,:)*qdt
      ! add J0xB0 source term in momentum equations
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1:ndir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(1:ndir))+axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:ndir)
    end if

    if(total_energy) then
      a=0.d0
      ! for free-free field -(vxB0) dot J0 =0
      b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(:))
      ! store full magnetic field B0+B1 in b
      if(.not.B0field_forcefree) b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         :)=b(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)+block%B0(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,:,0)
      ! store velocity in a
      do idir=1,ndir
        a(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)=wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mom(idir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_)
      end do
      call cross_product(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,a,b,axb)
      axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)=axb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,:)*qdt
      ! add -(vxB) dot J0 source term in energy equation
      do idir=7-2*ndir,3
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,e_)-axb(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idir)*block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
      end do
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_B0')

  end subroutine add_source_B0split

  !> Add resistive source to w within ixO Uses 3 point stencil (1 neighbour) in
  !> each direction, non-conservative. If the fourthorder precompiler flag is
  !> set, uses fourth order central difference for the laplacian. Then the
  !> stencil is 5 (2 neighbours).
  subroutine add_source_res1(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idir,jdir,kdir,idirmin,idim,&
       jxOmin1,jxOmin2,jxOmax1,jxOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,ix
    integer :: lxOmin1,lxOmin2,lxOmax1,lxOmax2, kxOmin1,kxOmin2,kxOmax1,&
       kxOmax2

    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       eta(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: gradeta(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    ! Calculating resistive sources involve one extra layer
    if (mhd_4th_order) then
      ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmax1=ixOmax1+2;ixAmax2=ixOmax2+2;
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
    end if

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2) call mpistop&
       ("Error in add_source_res1: Non-conforming input limits")

    ! Calculate current density and idirmin
    call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=mhd_eta
       gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndim)=zero
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixAmin1,ixAmin2,ixAmax1,ixAmax2,idirmin,x,current,eta)
       ! assumes that eta is not function of current?
       do idim=1,ndim
          call gradient(eta,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,idim,tmp)
          gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)=tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)
       end do
    end if

    if(B0field) then
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=wCT(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndir,0)
    else
      Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=wCT(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndir))
    end if

    do idir=1,ndir
       ! Put B_idir into tmp2 and eta*Laplace B_idir into tmp
       if (mhd_4th_order) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=Bf(ixImin1:ixImax1,&
            ixImin2:ixImax2,idir)
         do idim=1,ndim
            lxOmin1=ixOmin1+2*kr(idim,1);lxOmin2=ixOmin2+2*kr(idim,2)
            lxOmax1=ixOmax1+2*kr(idim,1);lxOmax2=ixOmax2+2*kr(idim,2);
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmax1=ixOmax1+kr(idim,1);jxOmax2=ixOmax2+kr(idim,2);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
            kxOmin1=ixOmin1-2*kr(idim,1);kxOmin2=ixOmin2-2*kr(idim,2)
            kxOmax1=ixOmax1-2*kr(idim,1);kxOmax2=ixOmax2-2*kr(idim,2);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(-tmp2(lxOmin1:lxOmax1,&
               lxOmin2:lxOmax2)+16.0d0*tmp2(jxOmin1:jxOmax1,&
               jxOmin2:jxOmax2)-30.0d0*tmp2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+16.0d0*tmp2(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2)-tmp2(kxOmin1:kxOmax1,&
               kxOmin2:kxOmax2)) /(12.0d0 * dxlevel(idim)**2)
         end do
       else
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
         tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=Bf(ixImin1:ixImax1,&
            ixImin2:ixImax2,idir)
         do idim=1,ndim
            jxOmin1=ixOmin1+kr(idim,1);jxOmin2=ixOmin2+kr(idim,2)
            jxOmax1=ixOmax1+kr(idim,1);jxOmax2=ixOmax2+kr(idim,2);
            hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
            hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+(tmp2(jxOmin1:jxOmax1,&
               jxOmin2:jxOmax2)-2.0d0*tmp2(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)+tmp2(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2))/dxlevel(idim)**2
         end do
       end if

       ! Multiply by eta
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

       ! Subtract grad(eta) x J = eps_ijk d_j eta J_k if eta is non-constant
       if (mhd_eta<zero)then
          do jdir=1,ndim; do kdir=idirmin,3
             if (lvc(idir,jdir,kdir)/=0)then
                if (lvc(idir,jdir,kdir)==1)then
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2)-gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      jdir)*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
                else
                   tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                      ixOmin2:ixOmax2)+gradeta(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                      jdir)*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,kdir)
                end if
             end if
          end do; end do
       end if

       ! Add sources related to eta*laplB-grad(eta) x J to B and e
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(idir))+qdt*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if (mhd_energy) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,e_)+qdt*tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
          if(mhd_solve_eaux) then
            w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)=w(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,eaux_)+qdt*tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)*Bf(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)
          end if
       end if
    end do ! idir

    if (mhd_energy) then
       ! de/dt+=eta*J**2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
       do idir=idirmin,3
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir)**2
       end do
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          e_)+qdt*eta(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       if(mhd_solve_eaux) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,eaux_)+qdt*eta(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_res1')

  end subroutine add_source_res1

  !> Add resistive source to w within ixO
  !> Uses 5 point stencil (2 neighbours) in each direction, conservative
  subroutine add_source_res2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       eta(ixImin1:ixImax1,ixImin2:ixImax2),curlj(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:3)
    double precision :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,1:3)
    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idir,idirmin,idirmin1

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmax1=ixOmax1+2;ixAmax2=ixOmax2+2;

    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2) call mpistop&
       ("Error in add_source_res2: Non-conforming input limits")

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
    ! Calculate current density within ixL: J=curl B, thus J_i=eps_ijk*d_j B_k
    ! Determine exact value of idirmin while doing the loop.
    call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,idirmin,current)

    if (mhd_eta>zero)then
       eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=mhd_eta
    else
       call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImax1,ixImax2,&
          ixAmin1,ixAmin2,ixAmax1,ixAmax2,idirmin,x,current,eta)
    end if

    ! dB/dt= -curl(J*eta), thus B_i=B_i-eps_ijk d_j Jeta_k
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    do idir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idir)=current(ixAmin1:ixAmax1,&
          ixAmin2:ixAmax2,idir)*eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)
    end do
    curlj=0.d0
    call curlvector(tmpvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,curlj,idirmin1,1,3)
    if(stagger_grid.and.ndim==2.and.ndir==3) then
      ! if 2.5D
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(ndir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(ndir))-qdt*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         ndir)
    else
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1:ndir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(1:ndir))-qdt*curlj(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:ndir)
    end if

    if(mhd_energy) then
      ! de/dt= +div(B x Jeta) = eta J^2 - B dot curl(eta J)
      ! de1/dt= eta J^2 - B1 dot curl(eta J)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+qdt*(sum(current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)**2,&
         dim=ndim+1)*eta(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(1:ndir))*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir),&
         dim=ndim+1))
      if(mhd_solve_eaux) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,eaux_)+qdt*(sum(current(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,:)**2,dim=ndim+1)*eta(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)-sum(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mag(1:ndir))*curlj(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:ndir),&
           dim=ndim+1))
      end if
    end if

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_res2')
  end subroutine add_source_res2

  !> Add Hyper-resistive source to w within ixO
  !> Uses 9 point stencil (4 neighbours) in each direction.
  subroutine add_source_hyperres(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    !.. local ..
    double precision                :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndir:3)
    double precision                :: tmpvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3),tmpvec2(ixImin1:ixImax1,ixImin2:ixImax2,1:3),tmp(ixImin1:ixImax1,&
       ixImin2:ixImax2),ehyper(ixImin1:ixImax1,ixImin2:ixImax2,1:3)
    integer                         :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idir,&
       jdir,kdir,idirmin,idirmin1

    ixAmin1=ixOmin1-3;ixAmin2=ixOmin2-3;ixAmax1=ixOmax1+3;ixAmax2=ixOmax2+3;
    if (ixImin1>ixAmin1.or.ixImax1<ixAmax1.or.ixImin2>ixAmin2.or.&
       ixImax2<ixAmax2) call mpistop&
       ("Error in add_source_hyperres: Non-conforming input limits")

    call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,idirmin,current)
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    do jdir=idirmin,3
       tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,jdir)=current(ixAmin1:ixAmax1,&
          ixAmin2:ixAmax2,jdir)
    end do

    ixAmin1=ixOmin1-2;ixAmin2=ixOmin2-2;ixAmax1=ixOmax1+2;ixAmax2=ixOmax2+2;
    call curlvector(tmpvec,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,tmpvec2,idirmin1,1,3)

    ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
    tmpvec(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    call curlvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,tmpvec,idirmin1,1,3)
    ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir) = - tmpvec(ixAmin1:ixAmax1,&
       ixAmin2:ixAmax2,1:ndir)*mhd_eta_hyper

    ixAmin1=ixOmin1;ixAmin2=ixOmin2;ixAmax1=ixOmax1;ixAmax2=ixOmax2;
    tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
    call curlvector(ehyper,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,ixAmin2,&
       ixAmax1,ixAmax2,tmpvec2,idirmin1,1,3)

    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-tmpvec2(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*qdt
    end do

    if (mhd_energy) then
      ! de/dt= +div(B x Ehyper)
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
      tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,1:ndir)=zero
      do idir=1,ndir; do jdir=1,ndir; do kdir=idirmin,3
        tmpvec2(ixAmin1:ixAmax1,ixAmin2:ixAmax2,idir) = tmpvec(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,idir)+ lvc(idir,jdir,kdir)*wCT(ixAmin1:ixAmax1,&
           ixAmin2:ixAmax2,mag(jdir))*ehyper(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           kdir)
      end do; end do; end do
      tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
      call divvector(tmpvec2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,tmp)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*qdt
      if(mhd_solve_eaux) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,eaux_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,eaux_)+tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*qdt
      end if
    end if

    if (check_small_values)  call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_hyperres')

  end subroutine add_source_hyperres

  subroutine add_source_glm(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Dedner JCP 2002, 175, 645 _equation 24_
    ! giving the EGLM-MHD scheme
    use mod_global_parameters
    use mod_geometry

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qdt, wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision:: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    integer          :: idim,idir
    double precision :: gradPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb, mhd_divb_4thorder)

    ! dPsi/dt =  - Ch^2/Cp^2 Psi
    if (mhd_glm_alpha < zero) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = &
         abs(mhd_glm_alpha)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
    else
      ! implicit update of Psi variable
      ! equation (27) in Mignone 2010 J. Com. Phys. 229, 2117
      if(slab_uniform) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(dxlevel(:)))*w(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           psi_) = dexp(-qdt*cmax_global*mhd_glm_alpha/minval(block%ds(&
           ixOmin1:ixOmax1,ixOmin2:ixOmax2,:),dim=ndim+1))*w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,psi_)
      end if
    end if

    ! gradient of Psi
    do idim=1,ndim
       select case(typegrad)
       case("central")
          call gradient(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
             ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
             gradPsi)
       case("limited")
          call gradientS(wCT(ixImin1:ixImax1,ixImin2:ixImax2,psi_),ixImin1,&
             ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,&
             gradPsi)
       end select
       if (total_energy) then
       ! e  = e  -qdt (b . grad(Psi))
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_) = w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,e_)-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(idim))*gradPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_glm')

  end subroutine add_source_glm

  !> Add divB related sources to w within ixO corresponding to Powel
  subroutine add_source_powel(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2),&
       v(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb, mhd_divb_4thorder)

    ! calculate velocity
    call mhd_get_v(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,v)

    if (total_energy) then
      ! e = e - qdt (v . b) * div b
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)-qdt*sum(v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:)*wCT(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(:)),dim=ndim+1)*divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
    end if

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*v(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    ! m = m - qdt b div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(idir))-qdt*mhd_mag_i_all(w,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
         idir)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_powel')

  end subroutine add_source_powel

  subroutine add_source_janhunen(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add divB related sources to w within ixO
    ! corresponding to Janhunen, just the term in the induction equation.
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt,   wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision                :: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    integer                         :: idir

    ! We calculate now div B
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb, mhd_divb_4thorder)

    ! b = b - qdt v * div b
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))-qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mom(idir))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_janhunen')

  end subroutine add_source_janhunen

  subroutine add_source_linde(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    ! Add Linde's divB related sources to wnew within ixO
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, wCT(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer :: idim, idir, ixpmin1,ixpmin2,ixpmax1,ixpmax2, i1,i2, iside
    double precision :: divb(ixImin1:ixImax1,ixImin2:ixImax2),&
       graddivb(ixImin1:ixImax1,ixImin2:ixImax2)
    logical, dimension(-1:1,-1:1) :: leveljump

    ! Calculate div B
    ixpmin1=ixOmin1-1;ixpmin2=ixOmin2-1;ixpmax1=ixOmax1+1;ixpmax2=ixOmax2+1;
    call get_divb(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,ixpmax1,&
       ixpmax2,divb, mhd_divb_4thorder)

    ! for AMR stability, retreat one cell layer from the boarders of level jump
    do i2=-1,1
    do i1=-1,1
      if(i1==0.and.i2==0) cycle
      if(neighbor_type(i1,i2,saveigrid)==2 .or. neighbor_type(i1,i2,&
         saveigrid)==4) then
        leveljump(i1,i2)=.true.
      else
        leveljump(i1,i2)=.false.
      end if
    end do
    end do

    ixpmin1=ixOmin1;ixpmin2=ixOmin2;ixpmax1=ixOmax1;ixpmax2=ixOmax2;
    do idim=1,ndim
      select case(idim)
       case(1)
          do iside=1,2
            i1=kr(1,1)*(2*iside-3);i2=kr(2,1)*(2*iside-3);
            if (leveljump(i1,i2)) then
              if (iside==1) then
                ixpmin1=ixOmin1-i1
              else
                ixpmax1=ixOmax1-i1
              end if
            end if
          end do
       
       case(2)
          do iside=1,2
            i1=kr(1,2)*(2*iside-3);i2=kr(2,2)*(2*iside-3);
            if (leveljump(i1,i2)) then
              if (iside==1) then
                ixpmin2=ixOmin2-i2
              else
                ixpmax2=ixOmax2-i2
              end if
            end if
          end do
       
      end select
    end do

    ! Add Linde's diffusive terms
    do idim=1,ndim
       ! Calculate grad_idim(divb)
       select case(typegrad)
       case("central")
         call gradient(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
            ixpmax1,ixpmax2,idim,graddivb)
       case("limited")
         call gradientS(divb,ixImin1,ixImin2,ixImax1,ixImax2,ixpmin1,ixpmin2,&
            ixpmax1,ixpmax2,idim,graddivb)
       end select

       ! Multiply by Linde's eta*dt = divbdiff*(c_max*dx)*dt = divbdiff*dx**2
       if (slab_uniform) then
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2)*divbdiff/(1.0d0/dxlevel(1)**2+&
             1.0d0/dxlevel(2)**2)
       else
          graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)=graddivb(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2)*divbdiff /(1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,1)**2+1.0d0/block%ds(ixpmin1:ixpmax1,&
             ixpmin2:ixpmax2,2)**2)
       end if

       w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,mag(idim))=w(ixpmin1:ixpmax1,&
          ixpmin2:ixpmax2,mag(idim))+graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)

       if (typedivbdiff=='all' .and. total_energy) then
         ! e += B_idim*eta*grad_idim(divb)
         w(ixpmin1:ixpmax1,ixpmin2:ixpmax2,e_)=w(ixpmin1:ixpmax1,&
            ixpmin2:ixpmax2,e_)+wCT(ixpmin1:ixpmax1,ixpmin2:ixpmax2,&
            mag(idim))*graddivb(ixpmin1:ixpmax1,ixpmin2:ixpmax2)
       end if
    end do

    if (check_small_values) call mhd_handle_small_values(.false.,w,x,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       'add_source_linde')

  end subroutine add_source_linde

  !> Calculate div B within ixO
  subroutine get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,divb, fourthorder)

    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: divb(ixImin1:ixImax1,ixImin2:ixImax2)
    logical, intent(in), optional   :: fourthorder

    double precision                   :: bvec(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndir)
    double precision                   :: divb_corner(ixImin1:ixImax1,&
       ixImin2:ixImax2), sign
    double precision                   :: aux_vol(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer                            :: ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
        idir, ic1,ic2, ixmin1,ixmin2,ixmax1,ixmax2

    if(stagger_grid) then
      divb=0.d0
      do idir=1,ndim
        ixCmin1=ixOmin1-kr(idir,1);ixCmin2=ixOmin2-kr(idir,2)
        ixCmax1=ixOmax1-kr(idir,1);ixCmax2=ixOmax2-kr(idir,2);
        divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divb(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+block%ws(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idir)*block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           idir)-block%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idir)*block%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    else
      bvec(ixImin1:ixImax1,ixImin2:ixImax2,:)=w(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(:))
      select case(typediv)
      case("central")
        call divvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,divb,fourthorder)
      case("limited")
        call divvectorS(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,divb)
      end select
    end if

  end subroutine get_divb

  !> get dimensionless div B = |divB| * volume / area / |B|
  subroutine get_normalized_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,divb)

    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision                   :: divb(ixImin1:ixImax1,&
       ixImin2:ixImax2), dsurface(ixImin1:ixImax1,ixImin2:ixImax2)

    integer :: ixAmin1,ixAmin2,ixAmax1,ixAmax2,idims

    call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,divb)
    if(slab_uniform) then
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0.5d0*abs(divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/sqrt(mhd_mag_en_all(w,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2))/sum(1.d0/dxlevel(:))
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;
      ixAmax1=ixOmax1-1;ixAmax2=ixOmax2-1;
      dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)= &
         sum(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:),dim=ndim+1)
      do idims=1,ndim
        ixAmin1=ixOmin1-kr(idims,1);ixAmin2=ixOmin2-kr(idims,2)
        ixAmax1=ixOmax1-kr(idims,1);ixAmax2=ixOmax2-kr(idims,2);
        dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=dsurface(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2)+block%surfaceC(ixAmin1:ixAmax1,ixAmin2:ixAmax2,&
           idims)
      end do
      divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=abs(divb(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2))/sqrt(mhd_mag_en_all(w,ixImin1,ixImin2,ixImax1,&
         ixImax2,ixOmin1,ixOmin2,ixOmax1,&
         ixOmax2))*block%dvolume(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/dsurface(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end if

  end subroutine get_normalized_divb

  !> Calculate idirmin and the idirmin:3 components of the common current array
  !> make sure that dxlevel(^D) is set correctly.
  subroutine get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idirmin,current)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)  :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, ixImin1,ixImin2,&
       ixImax1,ixImax2
    double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    integer, intent(out) :: idirmin
    integer :: idir, idirmin0

    ! For ndir=2 only 3rd component of J can exist, ndir=1 is impossible for MHD
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3),&
       bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    idirmin0 = 7-2*ndir

    bvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,mag(1:ndir))

    call curlvector(bvec,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,current,idirmin,idirmin0,ndir)

    if(B0field) current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin0:3)=current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin0:3)+block%J0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idirmin0:3)

  end subroutine get_current

  !> If resistivity is not zero, check diffusion time limit for dt
  subroutine mhd_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    use mod_global_parameters
    use mod_usr_methods
    use mod_radiative_cooling, only: cooling_get_dt
    use mod_viscosity, only: viscosity_get_dt
    use mod_gravity, only: gravity_get_dt

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: dtnew
    double precision, intent(in)    :: dx1,dx2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    integer                       :: idirmin,idim
    double precision              :: dxarr(ndim)
    double precision              :: current(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndir:3),eta(ixImin1:ixImax1,ixImin2:ixImax2)

    dtnew = bigdouble

    dxarr(1)=dx1;dxarr(2)=dx2;
    if (mhd_eta>zero)then
       dtnew=dtdiffpar*minval(dxarr(1:ndim))**2/mhd_eta
    else if (mhd_eta<zero)then
       call get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,idirmin,current)
       call usr_special_resistivity(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,idirmin,x,current,eta)
       dtnew=bigdouble
       do idim=1,ndim
         if(slab_uniform) then
           dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)/dxarr(idim)**2)))
         else
           dtnew=min(dtnew,dtdiffpar/(smalldouble+maxval(eta(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)/block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              idim)**2)))
         end if
       end do
    end if

    if(mhd_eta_hyper>zero) then
      if(slab_uniform) then
        dtnew=min(dtdiffpar*minval(dxarr(1:ndim))**4/mhd_eta_hyper,dtnew)
      else
        dtnew=min(dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1:ndim))**4/mhd_eta_hyper,dtnew)
      end if
    end if

    if(mhd_radiative_cooling) then
      call cooling_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

    if(mhd_viscosity) then
      call viscosity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

    if(mhd_gravity) then
      call gravity_get_dt(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,dtnew,dx1,dx2,x)
    end if

  end subroutine mhd_get_dt

  ! Add geometrical source terms to w
  subroutine mhd_add_source_geom(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,wCT,w,x)
    use mod_global_parameters
    use mod_geometry

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    integer          :: iw,idir, h1xmin1,h1xmin2,h1xmax1,h1xmax2, h2xmin1,&
       h2xmin2,h2xmax1,h2xmax2
    double precision :: tmp(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp1(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    integer :: mr_,mphi_ ! Polar var. names
    integer :: br_,bphi_

    mr_=mom(1); mphi_=mom(1)-1+phi_  ! Polar var. names
    br_=mag(1); bphi_=mag(1)-1+phi_

    select case (coordinate)
    case (cylindrical)
      if (angmomfix) then
        call mpistop("angmomfix not implemented yet in MHD")
      endif
      call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,tmp)
      if(phi_>0) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mr_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)*(tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)-wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,bphi_)**2+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           mphi_)**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mphi_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)*(-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mphi_)*wCT(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mr_)/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           rho_) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,br_))
        if(.not.stagger_grid) then
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bphi_)=w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,bphi_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             1)*(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             bphi_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mr_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             br_)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mphi_)) /wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
        end if
      else
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mr_)=w(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,mr_)+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)*tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      end if
      if(mhd_glm) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,br_)=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,br_)+qdt*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         psi_)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
    case (spherical)
       h1xmin1=ixOmin1-kr(1,1);h1xmin2=ixOmin2-kr(1,2)
       h1xmax1=ixOmax1-kr(1,1);h1xmax2=ixOmax2-kr(1,2)
       h2xmin1=ixOmin1-kr(2,1);h2xmin2=ixOmin2-kr(2,2)
       h2xmax1=ixOmax1-kr(2,1);h2xmax2=ixOmax2-kr(2,2);
       call mhd_get_p_total(wCT,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
          ixOmin2,ixOmax1,ixOmax2,tmp1)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp1(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       if(B0field) then
         tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum(block%B0(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,:,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:)),&
            dim=ndim+1)
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
       ! m1
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1)-block%surfaceC(h1xmin1:h1xmax1,h1xmin2:h1xmax2,&
          1))/block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       if(ndir>1) then
         do idir=2,ndir
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(idir))**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_)-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))**2
           if(B0field) tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              idir,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))
         end do
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       ! b1
       if(mhd_glm) then
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mag(1))+qdt/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            1)*2.0d0*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
       end if

       
       ! m2
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp1(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)
       if(B0field) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       end if
       ! This will make hydrostatic p=const an exact solution
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2) *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          2)-block%surfaceC(h2xmin1:h2xmax1,h2xmin2:h2xmax2,&
          2)) /block%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-(wCT(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mom(2))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2)))
       if (B0field) then
          tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
             0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(2)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             mag(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,0)
       end if
       if(ndir==3) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(3))**2/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(3))**2)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
         if (B0field) then
            tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)-2.0d0*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               3,0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               mag(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
               2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
         end if
       end if
       w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mom(2))+qdt*tmp(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       ! b2
       if(.not.stagger_grid) then
         tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
            mag(1)))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         if(B0field) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
              0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
              0))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
         end if
         if(mhd_glm) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2) + dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)
         end if
         w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))=w(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2,mag(2))+qdt*tmp(ixOmin1:ixOmax1,&
            ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
       end if
      

       if(ndir==3) then
         ! m3
         if(.not.angmomfix) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-(wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(1))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(1)))  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3))) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
                 0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
                 0)  +(block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
                 0)*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(3)) +wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mag(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
                 0)) *dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 2))/dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(3))=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(3))+qdt*tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
         else
           call mpistop("angmomfix not implemented yet in MHD")
         end if
         ! b3
         if(.not.stagger_grid) then
           tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(wCT(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mom(1))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(1)))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_)  -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(3))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(2)) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mom(2))*wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              mag(3)))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              2)) /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              rho_)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
           if (B0field) then
              tmp(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=tmp(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2)+(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(1))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
                 0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1,&
                 0))/wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 rho_) -(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(3))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2,&
                 0) -wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 mom(2))*block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3,&
                 0))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 2)) /(wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
                 rho_)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
           end if
           w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=w(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2,mag(3))+qdt*tmp(ixOmin1:ixOmax1,&
              ixOmin2:ixOmax2)/x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
         end if
       end if
    end select
  end subroutine mhd_add_source_geom

  !> Compute 2 times total magnetic energy
  function mhd_mag_en_all(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2) result(mge)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (B0field) then
      mge = sum((w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(:))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,:,block%iw0))**2,&
          dim=ndim+1)
    else
      mge = sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2, dim=ndim+1)
    end if
  end function mhd_mag_en_all

  !> Compute full magnetic field by direction
  function mhd_mag_i_all(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idir) result(mgf)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: mgf(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    if (B0field) then
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          mag(idir))+block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir,block%iw0)
    else
      mgf = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(idir))
    end if
  end function mhd_mag_i_all

  !> Compute evolving magnetic energy
  function mhd_mag_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2) result(mge)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: mge(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    mge = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mag(:))**2,&
        dim=ndim+1)
  end function mhd_mag_en

  !> compute kinetic energy
  function mhd_kin_en(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, inv_rho) result(ke)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision              :: ke(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    double precision, intent(in), optional :: inv_rho(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    if (present(inv_rho)) then
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) * inv_rho
    else
       ke = 0.5d0 * sum(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, mom(:))**2,&
           dim=ndim+1) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
    end if
  end function mhd_kin_en

  subroutine mhd_getv_Hall(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,vHall)
    use mod_global_parameters

    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: vHall(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:3)

    integer          :: idir, idirmin
    double precision :: current(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndir:3)

    ! Calculate current density and idirmin
    call get_current(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,idirmin,current)
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:3) = zero
    vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin:3) = - mhd_etah*current(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       idirmin:3)
    do idir = idirmin, 3
       vHall(ixOmin1:ixOmax1,ixOmin2:ixOmax2,idir) = vHall(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,idir)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
    end do

  end subroutine mhd_getv_Hall

  subroutine mhd_getdt_Hall(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,dx1,dx2,dthall)
    use mod_global_parameters

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: dthall
    !.. local ..
    double precision :: dxarr(ndim)
    double precision :: bmag(ixImin1:ixImax1,ixImin2:ixImax2)

    dthall=bigdouble

    ! because we have that in cmax now:
    return

    dxarr(1)=dx1;dxarr(2)=dx2;

    if (.not. B0field) then
       bmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(sum(w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(:))**2, dim=ndim+1))
       bmag(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sqrt(sum((w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,mag(:)) + block%B0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
          1:ndir,block%iw0))**2))
    end if

    if(slab_uniform) then
      dthall=dtdiffpar*minval(dxarr(1:ndim))**&
         2.0d0/(mhd_etah*maxval(bmag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)))
    else
      dthall=dtdiffpar*minval(block%ds(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1:ndim))**2.0d0/(mhd_etah*maxval(bmag(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)))
    end if

  end subroutine mhd_getdt_Hall

  subroutine mhd_modify_wLR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,wLC,wRC,wLp,wRp,s,idir)
    use mod_global_parameters
    use mod_usr_methods
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idir
    double precision, intent(in)    :: qt
    double precision, intent(inout) :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    type(state)                     :: s
    double precision                :: dB(ixImin1:ixImax1,ixImin2:ixImax2),&
        dPsi(ixImin1:ixImax1,ixImin2:ixImax2)

    if(stagger_grid) then
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=s%ws(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=s%ws(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=s%ws(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir))=s%ws(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
    else
      ! Solve the Riemann problem for the linear 2x2 system for normal
      ! B-field and GLM_Psi according to Dedner 2002:
      ! This implements eq. (42) in Dedner et al. 2002 JcP 175
      ! Gives the Riemann solution on the interface
      ! for the normal B component and Psi in the GLM-MHD system.
      ! 23/04/2013 Oliver Porth
      dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)   = wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir)) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idir))
      dPsi(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = wRp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,psi_) - wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_)

      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idir))   = 0.5d0 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idir)) + wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         mag(idir))) - 0.5d0/cmax_global * dPsi(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         psi_)       = 0.5d0 * (wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         psi_) + wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         psi_)) - 0.5d0*cmax_global * dB(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))
      wRp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,psi_)
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))
      wRC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,psi_)
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idir)) = wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mag(idir))
      wLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,psi_) = wLp(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,psi_)
    end if

    if(associated(usr_set_wLR)) call usr_set_wLR(ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,wLC,wRC,wLp,wRp,s,idir)

  end subroutine mhd_modify_wLR

  subroutine mhd_boundary_adjust
    use mod_global_parameters
    integer :: iB, idim, iside, iigrid, igrid
    integer :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2, i1,i2

    ixGmin1=ixGlo1;ixGmin2=ixGlo2;ixGmax1=ixGhi1;ixGmax2=ixGhi2;
     do iigrid=1,igridstail; igrid=igrids(iigrid);
        if(.not.phyboundblock(igrid)) cycle
        saveigrid=igrid
        block=>ps(igrid)
        dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
        do idim=1,ndim
           ! to avoid using as yet unknown corner info in more than 1D, we
           ! fill only interior mesh ranges of the ghost cell ranges at first,
           ! and progressively enlarge the ranges to include corners later
           do iside=1,2
              i1=kr(1,idim)*(2*iside-3);i2=kr(2,idim)*(2*iside-3);
              if (neighbor_type(i1,i2,igrid)/=1) cycle
              iB=(idim-1)*2+iside
              if(.not.boundary_divbfix(iB)) cycle
              if(any(typeboundary(:,iB)=="special")) then
                ! MF nonlinear force-free B field extrapolation and data driven
                ! require normal B of the first ghost cell layer to be untouched by
                ! fixdivB=0 process, set boundary_divbfix_skip(iB)=1 in par file
                select case (idim)
                case (1)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin1=ixGmax1+1-nghostcells+&
                         boundary_divbfix_skip(2*1);ixOmin2=ixGmin2;
                      ixOmax1=ixGmax1;ixOmax2=ixGmax2;
                   else
                      ! minimal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmin2;
                      ixOmax1=ixGmin1-1+nghostcells-boundary_divbfix_skip(2*1-&
                         1);ixOmax2=ixGmax2;
                   end if 
                case (2)
                   if (iside==2) then
                      ! maximal boundary
                      ixOmin1=ixGmin1
                      ixOmin2=ixGmax2+1-nghostcells+&
                         boundary_divbfix_skip(2*2);
                      ixOmax1=ixGmax1;ixOmax2=ixGmax2;
                   else
                      ! minimal boundary
                      ixOmin1=ixGmin1;ixOmin2=ixGmin2;
                      ixOmax1=ixGmax1
                      ixOmax2=ixGmin2-1+nghostcells-boundary_divbfix_skip(2*2-&
                         1);
                   end if 
                end select
                call fixdivB_boundary(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,&
                   ixOmin2,ixOmax1,ixOmax2,ps(igrid)%w,ps(igrid)%x,iB)
              end if
           end do
        end do
     end do

  end subroutine mhd_boundary_adjust

  subroutine fixdivB_boundary(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,iB)
    use mod_global_parameters

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,iB
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)

    double precision :: dx1x2,dx1x3,dx2x1,dx2x3,dx3x1,dx3x2
    integer :: ix1,ix2,ixFmin1,ixFmin2,ixFmax1,ixFmax2

    select case(iB)
     case(1)
       ! 2nd order CD for divB=0 to set normal B component better
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1+1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=w(ix1+1,ixFmin2:ixFmax2,&
              mag(1)) +dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmax1,ixFmin1,-1
           w(ix1-1,ixFmin2:ixFmax2,mag(1))=( (w(ix1+1,ixFmin2:ixFmax2,&
              mag(1))+w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1,&
              ixFmin2:ixFmax2,1)+(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,&
              ixFmin2:ixFmax2,mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              2)-(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,&
              mag(2)))*block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,&
              2) )/block%surfaceC(ix1-1,ixFmin2:ixFmax2,1)-w(ix1,&
              ixFmin2:ixFmax2,mag(1))
         end do
       end if
      
       
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
     case(2)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1-1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx1x2=dxlevel(1)/dxlevel(2)
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=w(ix1-1,ixFmin2:ixFmax2,&
              mag(1)) -dx1x2*(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))-w(ix1,&
              ixFmin2-1:ixFmax2-1,mag(2)))
         enddo
       else
         do ix1=ixFmin1,ixFmax1
           w(ix1+1,ixFmin2:ixFmax2,mag(1))=( (w(ix1-1,ixFmin2:ixFmax2,&
              mag(1))+w(ix1,ixFmin2:ixFmax2,mag(1)))*block%surfaceC(ix1-1,&
              ixFmin2:ixFmax2,1)-(w(ix1,ixFmin2+1:ixFmax2+1,mag(2))+w(ix1,&
              ixFmin2:ixFmax2,mag(2)))*block%surfaceC(ix1,ixFmin2:ixFmax2,&
              2)+(w(ix1,ixFmin2:ixFmax2,mag(2))+w(ix1,ixFmin2-1:ixFmax2-1,&
              mag(2)))*block%surfaceC(ix1,ixFmin2-1:ixFmax2-1,&
              2) )/block%surfaceC(ix1,ixFmin2:ixFmax2,1)-w(ix1,ixFmin2:ixFmax2,&
              mag(1))
         end do
       end if
      
       
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
     case(3)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2+1
       ixFmax2=ixOmax2+1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=w(ixFmin1:ixFmax1,ix2+1,&
              mag(2)) +dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         enddo
       else
         do ix2=ixFmax2,ixFmin2,-1
           w(ixFmin1:ixFmax1,ix2-1,mag(2))=( (w(ixFmin1:ixFmax1,ix2+1,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              2)+(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              1)-(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,&
              1) )/block%surfaceC(ixFmin1:ixFmax1,ix2-1,2)-w(ixFmin1:ixFmax1,&
              ix2,mag(2))
         end do
       end if
      
       
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
     case(4)
       if(total_energy) call mhd_to_primitive(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
       
       ixFmin1=ixOmin1+1
       ixFmax1=ixOmax1-1
       ixFmin2=ixOmin2-1
       ixFmax2=ixOmax2-1
       if(slab_uniform) then
         dx2x1=dxlevel(2)/dxlevel(1)
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=w(ixFmin1:ixFmax1,ix2-1,&
              mag(2)) -dx2x1*(w(ixFmin1+1:ixFmax1+1,ix2,&
              mag(1))-w(ixFmin1-1:ixFmax1-1,ix2,mag(1)))
         end do
       else
         do ix2=ixFmin2,ixFmax2
           w(ixFmin1:ixFmax1,ix2+1,mag(2))=( (w(ixFmin1:ixFmax1,ix2-1,&
              mag(2))+w(ixFmin1:ixFmax1,ix2,&
              mag(2)))*block%surfaceC(ixFmin1:ixFmax1,ix2-1,&
              2)-(w(ixFmin1+1:ixFmax1+1,ix2,mag(1))+w(ixFmin1:ixFmax1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1:ixFmax1,ix2,&
              1)+(w(ixFmin1:ixFmax1,ix2,mag(1))+w(ixFmin1-1:ixFmax1-1,ix2,&
              mag(1)))*block%surfaceC(ixFmin1-1:ixFmax1-1,ix2,&
              1) )/block%surfaceC(ixFmin1:ixFmax1,ix2,2)-w(ixFmin1:ixFmax1,ix2,&
              mag(2))
         end do
       end if
      
       
       if(total_energy) call mhd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,&
          ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x)
     
     case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine fixdivB_boundary

  
  subroutine mhd_clean_divb_multigrid(qdt, qt, active)
    use mod_forest
    use mod_global_parameters
    use mod_multigrid_coupling
    use mod_geometry

    double precision, intent(in) :: qdt    !< Current time step
    double precision, intent(in) :: qt     !< Current time
    logical, intent(inout)       :: active !< Output if the source is active
    integer                      :: iigrid, igrid, id
    integer                      :: n, nc, lvl, ixmin1,ixmin2,ixmax1,ixmax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2, idim
    type(tree_node), pointer     :: pnode
    double precision             :: tmp(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
        grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2, ndim)
    double precision             :: res
    double precision, parameter  :: max_residual = 1d-3
    double precision, parameter  :: residual_reduction = 1d-10
    integer, parameter           :: max_its      = 50
    double precision             :: residual_it(max_its), max_divb

    mg%operator_type = mg_laplacian

    ! Set boundary conditions
    do n = 1, 2*ndim
       idim = (n+1)/2
       select case (typeboundary(mag(idim), n))
       case ('symm')
          ! d/dx B = 0, take phi = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('asymm')
          ! B = 0, so grad(phi) = 0
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('cont')
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('special')
          ! Assume Dirichlet boundary conditions, derivative zero
          mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       case ('periodic')
          ! Nothing to do here
       case default
          print *, "divb_multigrid warning: unknown b.c.: ", &
               trim(typeboundary(mag(idim), n))
          mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
          mg%bc(n, mg_iphi)%bc_value = 0.0_dp
       end select
    end do

    ixmin1=ixMlo1-1;ixmin2=ixMlo2-1;ixmax1=ixMhi1+1;ixmax2=ixMhi2+1;
    max_divb = 0.0d0

    ! Store divergence of B as right-hand side
    do iigrid = 1, igridstail
       igrid =  igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id
       lvl   =  mg%boxes(id)%lvl
       nc    =  mg%box_size_lvl(lvl)

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);

       call get_divb(ps(igrid)%w(ixGlo1:ixGhi1,ixGlo2:ixGhi2, 1:nw), ixGlo1,&
          ixGlo2,ixGhi1,ixGhi2, ixMlo1,ixMlo2,ixMhi1,ixMhi2, tmp,&
           mhd_divb_4thorder)
       mg%boxes(id)%cc(1:nc,1:nc, mg_irhs) = tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)
       max_divb = max(max_divb, maxval(abs(tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2))))
    end do

    ! Solve laplacian(phi) = divB
    if(stagger_grid) then
      call MPI_ALLREDUCE(MPI_IN_PLACE, max_divb, 1, MPI_DOUBLE_PRECISION,&
          MPI_MAX, icomm, ierrmpi)

      if (mype == 0) print *, "Performing multigrid divB cleaning"
      if (mype == 0) print *, "iteration vs residual"
      ! Solve laplacian(phi) = divB
      do n = 1, max_its
         call mg_fas_fmg(mg, n>1, max_res=residual_it(n))
         if (mype == 0) write(*, "(I4,E11.3)") n, residual_it(n)
         if (residual_it(n) < residual_reduction * max_divb) exit
      end do
      if (mype == 0 .and. n > max_its) then
         print *, "divb_multigrid warning: not fully converged"
         print *, "current amplitude of divb: ", residual_it(max_its)
         print *, "multigrid smallest grid: ", &
              mg%domain_size_lvl(:, mg%lowest_lvl)
         print *, "note: smallest grid ideally has <= 8 cells"
         print *, "multigrid dx/dy/dz ratio: ", mg%dr(:, 1)/mg%dr(1, 1)
         print *, "note: dx/dy/dz should be similar"
      end if
    else
      do n = 1, max_its
         call mg_fas_vcycle(mg, max_res=res)
         if (res < max_residual) exit
      end do
      if (res > max_residual) call mpistop("divb_multigrid: no convergence")
    end if


    ! Correct the magnetic field
    do iigrid = 1, igridstail
       igrid = igrids(iigrid);
       pnode => igrid_to_node(igrid, mype)%node
       id    =  pnode%id

       ! Geometry subroutines expect this to be set
       block => ps(igrid)
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);

       ! Compute the gradient of phi
       tmp(ixmin1:ixmax1,ixmin2:ixmax2) = mg%boxes(id)%cc(:,:, mg_iphi)

       if(stagger_grid) then
         do idim =1, ndim
           ixCmin1=ixMlo1-kr(idim,1);ixCmin2=ixMlo2-kr(idim,2);
           ixCmax1=ixMhi1;ixCmax2=ixMhi2;
           call gradientx(tmp,ps(igrid)%x,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixCmin1,&
              ixCmin2,ixCmax1,ixCmax2,idim,grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
              idim),.false.)
           ! Apply the correction B* = B - gradient(phi)
           ps(igrid)%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              idim)=ps(igrid)%ws(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              idim)-grad(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim)
         end do
         ! store cell-center magnetic energy
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = sum(ps(igrid)%w(ixMlo1:ixMhi1,&
            ixMlo2:ixMhi2, mag(1:ndim))**2, dim=ndim+1)
         ! change cell-center magnetic field
         call mhd_face_to_center(ixMlo1,ixMlo2,ixMhi1,ixMhi2,ps(igrid))
       else
         do idim = 1, ndim
            call gradient(tmp,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
               ixMhi2,idim,grad(ixGlo1:ixGhi1,ixGlo2:ixGhi2, idim))
         end do
         ! store cell-center magnetic energy
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = sum(ps(igrid)%w(ixMlo1:ixMhi1,&
            ixMlo2:ixMhi2, mag(1:ndim))**2, dim=ndim+1)
         ! Apply the correction B* = B - gradient(phi)
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             mag(1:ndim)) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             mag(1:ndim)) - grad(ixMlo1:ixMhi1,ixMlo2:ixMhi2, :)
       end if

       if(total_energy) then
         ! Determine magnetic energy difference
         tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2) = 0.5_dp * &
            (sum(ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2, mag(1:ndim))**2,&
             dim=ndim+1) - tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2))
         ! Keep thermal pressure the same
         ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             e_) = ps(igrid)%w(ixMlo1:ixMhi1,ixMlo2:ixMhi2,&
             e_) + tmp(ixMlo1:ixMhi1,ixMlo2:ixMhi2)
       end if
    end do

    active = .true.

  end subroutine mhd_clean_divb_multigrid
 

  subroutine mhd_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,s)
    use mod_global_parameters

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt,qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wprim(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)

    select case(type_ct)
    case('average')
      call update_faces_average(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,qt,qdt,fC,fE,sCT,s)
    case('uct_contact')
      call update_faces_contact(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,s)
    case('uct_hll')
      call update_faces_hll(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,qt,qdt,fE,sCT,s)
    case default
      call mpistop('choose average, uct_contact,or uct_hll for type_ct!')
    end select

  end subroutine mhd_update_faces

  !> get electric field though averaging neighors to update faces in CT
  subroutine update_faces_average(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,qdt,fC,fE,sCT,s)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt, qdt
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)

    integer                            :: hxCmin1,hxCmin2,hxCmax1,hxCmax2,&
       ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,&
       ixCmmin1,ixCmmin2,ixCmmax1,ixCmmax2
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndim:3)

    associate(bfaces=>s%ws,x=>s%x)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    ixCmax1=ixOmax1;ixCmax2=ixOmax2;
    ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,sCT,s,jce)

    fE=zero

    do idim1=1,ndim
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmax1=ixCmax1+kr(idim1,1);jxCmax2=ixCmax2+kr(idim1,2);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmax1=ixCmax1+kr(idim2,1);hxCmax2=ixCmax2+kr(idim2,2);
            ! Interpolate to edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwdim1,&
               idim2)+fC(jxCmin1:jxCmax1,jxCmin2:jxCmax2,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwdim2,&
               idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iwdim2,idim1))

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(mhd_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)+jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)

            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=qdt*s%dsC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)

            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=zero

    ! Calculate circulation on each face

    do idim1=1,ndim ! Coordinate perpendicular to face
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmax1=ixCmax1-kr(idim2,1);hxCmax2=ixCmax2-kr(idim2,2);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2);
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=zero
      end where
      ! Time update
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=bfaces(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)
    end do

    end associate

  end subroutine update_faces_average

  !> update faces using UCT contact mode by Gardiner and Stone 2005 JCP 205, 509
  subroutine update_faces_contact(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,qt,qdt,wp,fC,fE,sCT,s)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt, qdt
    ! cell-center primitive variables
    double precision, intent(in)       :: wp(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    type(state)                        :: sCT, s
    double precision, intent(in)       :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux,1:ndim)
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)

    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    ! electric field at cell centers
    double precision                   :: ECC(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)
    ! gradient of E at left and right side of a cell face
    double precision                   :: EL(ixImin1:ixImax1,ixImin2:ixImax2),&
       ER(ixImin1:ixImax1,ixImin2:ixImax2)
    ! gradient of E at left and right side of a cell corner
    double precision                   :: ELC(ixImin1:ixImax1,ixImin2:ixImax2),&
       ERC(ixImin1:ixImax1,ixImin2:ixImax2)
    ! current on cell edges
    double precision                   :: jce(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)
    ! total magnetic field at cell centers
    double precision                   :: Btot(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    integer                            :: hxCmin1,hxCmin2,hxCmax1,hxCmax2,&
       ixCmin1,ixCmin2,ixCmax1,ixCmax2,jxCmin1,jxCmin2,jxCmax1,jxCmax2,ixAmin1,&
       ixAmin2,ixAmax1,ixAmax2,ixBmin1,ixBmin2,ixBmax1,ixBmax2
    integer                            :: idim1,idim2,idir,iwdim1,iwdim2

    associate(bfaces=>s%ws,x=>s%x,w=>s%w,vnorm=>vcts%vnorm)

    if(B0field) then
      Btot(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=wp(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndim))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:ndim,0)
    else
      Btot(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=wp(ixImin1:ixImax1,&
         ixImin2:ixImax2,mag(1:ndim))
    end if
    ECC=0.d0
    ! Calculate electric field at cell centers
    do idim1=1,ndim; do idim2=1,ndim; do idir=7-2*ndim,3
      if(lvc(idim1,idim2,idir)==1)then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,idir)=ECC(ixImin1:ixImax1,&
            ixImin2:ixImax2,idir)+Btot(ixImin1:ixImax1,ixImin2:ixImax2,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,mom(idim2))
      else if(lvc(idim1,idim2,idir)==-1) then
         ECC(ixImin1:ixImax1,ixImin2:ixImax2,idir)=ECC(ixImin1:ixImax1,&
            ixImin2:ixImax2,idir)-Btot(ixImin1:ixImax1,ixImin2:ixImax2,&
            idim1)*wp(ixImin1:ixImax1,ixImin2:ixImax2,mom(idim2))
      endif
    enddo; enddo; enddo

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,sCT,s,jce)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.
    fE=zero
    ! evaluate electric field along cell edges according to equation (41)
    do idim1=1,ndim
      iwdim1 = mag(idim1)
      do idim2=1,ndim
        iwdim2 = mag(idim2)
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Allow only even permutations
          if (lvc(idim1,idim2,idir)==1) then
            ixCmax1=ixOmax1;ixCmax2=ixOmax2;
            ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1;
            ! Assemble indices
            jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
            jxCmax1=ixCmax1+kr(idim1,1);jxCmax2=ixCmax2+kr(idim1,2);
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmax1=ixCmax1+kr(idim2,1);hxCmax2=ixCmax2+kr(idim2,2);
            ! average cell-face electric field to cell edges
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)=quarter*(fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwdim1,&
               idim2)+fC(jxCmin1:jxCmax1,jxCmin2:jxCmax2,iwdim1,&
               idim2)-fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iwdim2,&
               idim1)-fC(hxCmin1:hxCmax1,hxCmin2:hxCmax2,iwdim2,idim1))

            ! add slope in idim2 direction from equation (50)
            ixAmin1=ixCmin1;ixAmin2=ixCmin2;
            ixAmax1=ixCmax1+kr(idim1,1);ixAmax2=ixCmax2+kr(idim1,2);
            EL(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=fC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,iwdim1,idim2)-ECC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,idir)
            hxCmin1=ixAmin1+kr(idim2,1);hxCmin2=ixAmin2+kr(idim2,2)
            hxCmax1=ixAmax1+kr(idim2,1);hxCmax2=ixAmax2+kr(idim2,2);
            ER(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=fC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,iwdim1,idim2)-ECC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,idir)
            where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)>0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=EL(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)
            else where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)<0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=EL(jxCmin1:jxCmax1,&
                 jxCmin2:jxCmax2)
            else where
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(EL(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)+EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
            end where
            hxCmin1=ixCmin1+kr(idim2,1);hxCmin2=ixCmin2+kr(idim2,2)
            hxCmax1=ixCmax1+kr(idim2,1);hxCmax2=ixCmax2+kr(idim2,2);
            where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idim1)>0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ER(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)
            else where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idim1)<0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ER(jxCmin1:jxCmax1,&
                 jxCmin2:jxCmax2)
            else where
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(ER(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)+ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
            end where
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)+0.25d0*(ELC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)+ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

            ! add slope in idim1 direction from equation (50)
            jxCmin1=ixCmin1+kr(idim2,1);jxCmin2=ixCmin2+kr(idim2,2)
            jxCmax1=ixCmax1+kr(idim2,1);jxCmax2=ixCmax2+kr(idim2,2);
            ixAmin1=ixCmin1;ixAmin2=ixCmin2;
            ixAmax1=ixCmax1+kr(idim2,1);ixAmax2=ixCmax2+kr(idim2,2);
            EL(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=-fC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,iwdim2,idim1)-ECC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,idir)
            hxCmin1=ixAmin1+kr(idim1,1);hxCmin2=ixAmin2+kr(idim1,2)
            hxCmax1=ixAmax1+kr(idim1,1);hxCmax2=ixAmax2+kr(idim1,2);
            ER(ixAmin1:ixAmax1,ixAmin2:ixAmax2)=-fC(ixAmin1:ixAmax1,&
               ixAmin2:ixAmax2,iwdim2,idim1)-ECC(hxCmin1:hxCmax1,&
               hxCmin2:hxCmax2,idir)
            where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim2)>0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=EL(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)
            else where(vnorm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim2)<0.d0)
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=EL(jxCmin1:jxCmax1,&
                 jxCmin2:jxCmax2)
            else where
              ELC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(EL(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)+EL(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
            end where
            hxCmin1=ixCmin1+kr(idim1,1);hxCmin2=ixCmin2+kr(idim1,2)
            hxCmax1=ixCmax1+kr(idim1,1);hxCmax2=ixCmax2+kr(idim1,2);
            where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idim2)>0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ER(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)
            else where(vnorm(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idim2)<0.d0)
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=ER(jxCmin1:jxCmax1,&
                 jxCmin2:jxCmax2)
            else where
              ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=0.5d0*(ER(ixCmin1:ixCmax1,&
                 ixCmin2:ixCmax2)+ER(jxCmin1:jxCmax1,jxCmin2:jxCmax2))
            end where
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)+0.25d0*(ELC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2)+ERC(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

            ! add current component of electric field at cell edges E=-vxB+eta J
            if(mhd_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)+jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)

            ! times time step and edge length
            fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=fE(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)*qdt*s%dsC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               idir)
            if (.not.slab) then
              where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 r_)+half*dxlevel(r_))<1.0d-9)
                fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=zero
              end where
            end if
          end if
        end do
      end do
    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=zero

    ! Calculate circulation on each face
    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmax1=ixCmax1-kr(idim2,1);hxCmax2=ixCmax2-kr(idim2,2);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idir))
        end do
      end do
      ! Divide by the area of the face to get dB/dt
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2);
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=zero
      end where
      ! Time update cell-face magnetic field component
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=bfaces(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)
    end do

    end associate

  end subroutine update_faces_contact

  !> update faces
  subroutine update_faces_hll(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,qdt,fE,sCT,s)
    use mod_global_parameters
    use mod_constrained_transport
    use mod_usr_methods

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: qt, qdt
    double precision, intent(inout)    :: fE(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)
    type(state)                        :: sCT, s

    double precision                   :: vtilL(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)
    double precision                   :: vtilR(ixImin1:ixImax1,&
       ixImin2:ixImax2,2)
    double precision                   :: bfacetot(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)
    double precision                   :: btilL(s%ixGsmin1:s%ixGsmax1,&
       s%ixGsmin2:s%ixGsmax2,ndim)
    double precision                   :: btilR(s%ixGsmin1:s%ixGsmax1,&
       s%ixGsmin2:s%ixGsmax2,ndim)
    double precision                   :: cp(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)
    double precision                   :: cm(ixImin1:ixImax1,ixImin2:ixImax2,&
       2)
    double precision                   :: circ(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    ! current on cell edges
    double precision                   :: jce(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)
    integer                            :: hxCmin1,hxCmin2,hxCmax1,hxCmax2,&
       ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixCpmin1,ixCpmin2,ixCpmax1,ixCpmax2,&
       jxCmin1,jxCmin2,jxCmax1,jxCmax2,ixCmmin1,ixCmmin2,ixCmmax1,ixCmmax2
    integer                            :: idim1,idim2,idir

    associate(bfaces=>s%ws,bfacesCT=>sCT%ws,x=>s%x,vbarC=>vcts%vbarC,&
       cbarmin=>vcts%cbarmin,cbarmax=>vcts%cbarmax)

    ! Calculate contribution to FEM of each edge,
    ! that is, estimate value of line integral of
    ! electric field in the positive idir direction.

    ! Loop over components of electric field

    ! idir: electric field component we need to calculate
    ! idim1: directions in which we already performed the reconstruction
    ! idim2: directions in which we perform the reconstruction

    ! if there is resistivity, get eta J
    if(mhd_eta/=zero) call get_resistive_electric_field(ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,sCT,s,jce)

    fE=zero

    do idir=7-2*ndim,3
      ! Indices
      ! idir: electric field component
      ! idim1: one surface
      ! idim2: the other surface
      ! cyclic permutation: idim1,idim2,idir=1,2,3
      ! Velocity components on the surface
      ! follow cyclic premutations:
      ! Sx(1),Sx(2)=y,z ; Sy(1),Sy(2)=z,x ; Sz(1),Sz(2)=x,y

      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-1+kr(idir,1);ixCmin2=ixOmin2-1+kr(idir,2);

      ! Set indices and directions
      idim1=mod(idir,3)+1
      idim2=mod(idir+1,3)+1

      jxCmin1=ixCmin1+kr(idim1,1);jxCmin2=ixCmin2+kr(idim1,2)
      jxCmax1=ixCmax1+kr(idim1,1);jxCmax2=ixCmax2+kr(idim1,2);
      ixCpmin1=ixCmin1+kr(idim2,1);ixCpmin2=ixCmin2+kr(idim2,2)
      ixCpmax1=ixCmax1+kr(idim2,1);ixCpmax2=ixCmax2+kr(idim2,2);

      ! Reconstruct transverse transport velocities
      call reconstruct(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
         ixCmax2,idim2,vbarC(ixImin1:ixImax1,ixImin2:ixImax2,idim1,1),&
         vtilL(ixImin1:ixImax1,ixImin2:ixImax2,2),vtilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,2))

      call reconstruct(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
         ixCmax2,idim1,vbarC(ixImin1:ixImax1,ixImin2:ixImax2,idim2,2),&
         vtilL(ixImin1:ixImax1,ixImin2:ixImax2,1),vtilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,1))

      ! Reconstruct magnetic fields
      ! Eventhough the arrays are larger, reconstruct works with
      ! the limits ixG.
      if(B0field) then
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,&
           idim1)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           idim1)+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,idim1,idim1)
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,&
           idim2)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,&
           idim2)+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,idim2,idim2)
      else
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,&
           idim1)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,idim1)
        bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,&
           idim2)=bfacesCT(ixImin1:ixImax1,ixImin2:ixImax2,idim2)
      end if
      call reconstruct(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
         ixCmax2,idim2,bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,idim1),&
         btilL(ixImin1:ixImax1,ixImin2:ixImax2,idim1),btilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,idim1))

      call reconstruct(ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
         ixCmax2,idim1,bfacetot(ixImin1:ixImax1,ixImin2:ixImax2,idim2),&
         btilL(ixImin1:ixImax1,ixImin2:ixImax2,idim2),btilR(ixImin1:ixImax1,&
         ixImin2:ixImax2,idim2))

      ! Take the maximum characteristic

      cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=max(cbarmin(ixCpmin1:ixCpmax1,&
         ixCpmin2:ixCpmax2,idim1),cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1))
      cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)=max(cbarmax(ixCpmin1:ixCpmax1,&
         ixCpmin2:ixCpmax2,idim1),cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1))

      cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=max(cbarmin(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,idim2),cbarmin(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim2))
      cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)=max(cbarmax(jxCmin1:jxCmax1,&
         jxCmin2:jxCmax2,idim2),cbarmax(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim2))


      ! Calculate eletric field
      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=-(cp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1)*btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim2) + cm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1)*btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim2) - cp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         1)*(btilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim2)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim2)))/(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,1)+cm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,1)) +(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)*vtilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)*btilL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1) + cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)*vtilR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)*btilR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1) - cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)*cm(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)*(btilR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1)-btilL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1)))/(cp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,2)+cm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,2))

      ! add current component of electric field at cell edges E=-vxB+eta J
      if(mhd_eta/=zero) fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idir)=fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)+jce(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idir)

      fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=qdt*s%dsC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idir)*fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)

      if (.not.slab) then
        where(abs(x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           r_)+half*dxlevel(r_)).lt.1.0d-9)
          fE(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=zero
        end where
      end if

    end do

    ! allow user to change inductive electric field, especially for boundary driven applications
    if(associated(usr_set_electric_field)) call usr_set_electric_field(ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,fE,sCT)

    circ(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)=zero

    ! Calculate circulation on each face: interal(fE dot dl)

    do idim1=1,ndim ! Coordinate perpendicular to face
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2);
      do idim2=1,ndim
        do idir=7-2*ndim,3 ! Direction of line integral
          ! Assemble indices
          hxCmin1=ixCmin1-kr(idim2,1);hxCmin2=ixCmin2-kr(idim2,2)
          hxCmax1=ixCmax1-kr(idim2,1);hxCmax2=ixCmax2-kr(idim2,2);
          ! Add line integrals in direction idir
          circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idim1)+lvc(idim1,idim2,idir)*(fE(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idir)-fE(hxCmin1:hxCmax1,hxCmin2:hxCmax2,idir))
        end do
      end do
    end do

    ! Divide by the area of the face to get dB/dt
    do idim1=1,ndim
      ixCmax1=ixOmax1;ixCmax2=ixOmax2;
      ixCmin1=ixOmin1-kr(idim1,1);ixCmin2=ixOmin2-kr(idim1,2);
      where(s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         idim1) > 1.0d-9*s%dvolume(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=circ(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idim1)/s%surfaceC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           idim1)
      elsewhere
        circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=zero
      end where
      ! Time update
      bfaces(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)=bfaces(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,idim1)-circ(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idim1)
    end do

    end associate
  end subroutine update_faces_hll

  !> calculate eta J at cell edges
  subroutine get_resistive_electric_field(ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,sCT,s,jce)
    use mod_global_parameters
    use mod_usr_methods
    use mod_geometry

    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state), intent(in)            :: sCT, s
    ! current on cell edges
    double precision :: jce(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndim:3)

    ! current on cell centers
    double precision :: jcc(ixImin1:ixImax1,ixImin2:ixImax2,7-2*ndim:3)
    ! location at cell faces
    double precision :: xs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,1:ndim)
    ! resistivity
    double precision :: eta(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: gradi(ixGslo1:ixGshi1,ixGslo2:ixGshi2)
    integer :: ix1,ix2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixAmin1,ixAmin2,ixAmax1,&
       ixAmax2,ixBmin1,ixBmin2,ixBmax1,ixBmax2,idir,idirmin,idim1,idim2

    associate(x=>s%x,dx=>s%dx,w=>s%w,wCT=>sCT%w,wCTs=>sCT%ws)
    ! calculate current density at cell edges
    jce=0.d0
    do idim1=1,ndim
      do idim2=1,ndim
        do idir=7-2*ndim,3
          if (lvc(idim1,idim2,idir)==0) cycle
          ixCmax1=ixOmax1;ixCmax2=ixOmax2;
          ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1;
          ixBmax1=ixCmax1-kr(idir,1)+1;ixBmax2=ixCmax2-kr(idir,2)+1;
          ixBmin1=ixCmin1;ixBmin2=ixCmin2;
          ! current at transverse faces
          xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,:)=x(ixBmin1:ixBmax1,&
             ixBmin2:ixBmax2,:)
          xs(ixBmin1:ixBmax1,ixBmin2:ixBmax2,idim2)=x(ixBmin1:ixBmax1,&
             ixBmin2:ixBmax2,idim2)+half*dx(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
             idim2)
          call gradientx(wCTs(ixGslo1:ixGshi1,ixGslo2:ixGshi2,idim2),xs,&
             ixGslo1,ixGslo2,ixGshi1,ixGshi2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,&
             idim1,gradi,.true.)
          if (lvc(idim1,idim2,idir)==1) then
            jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=jce(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)+gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
          else
            jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=jce(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2,idir)-gradi(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
          end if
        end do
      end do
    end do
    ! get resistivity
    if(mhd_eta>zero)then
      jce(ixImin1:ixImax1,ixImin2:ixImax2,:)=jce(ixImin1:ixImax1,&
         ixImin2:ixImax2,:)*mhd_eta
    else
      ixAmin1=ixOmin1-1;ixAmin2=ixOmin2-1;ixAmax1=ixOmax1+1;ixAmax2=ixOmax2+1;
      call get_current(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,idirmin,jcc)
      call usr_special_resistivity(wCT,ixImin1,ixImin2,ixImax1,ixImax2,ixAmin1,&
         ixAmin2,ixAmax1,ixAmax2,idirmin,x,jcc,eta)
      ! calcuate eta on cell edges
      do idir=7-2*ndim,3
        ixCmax1=ixOmax1;ixCmax2=ixOmax2;
        ixCmin1=ixOmin1+kr(idir,1)-1;ixCmin2=ixOmin2+kr(idir,2)-1;
        jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=0.d0
       do ix2=0,1
       do ix1=0,1
          if( ix1==1 .and. 1==idir  .or. ix2==1 .and. 2==idir ) cycle
          ixAmin1=ixCmin1+ix1;ixAmin2=ixCmin2+ix2;
          ixAmax1=ixCmax1+ix1;ixAmax2=ixCmax2+ix2;
          jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=jcc(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,idir)+eta(ixAmin1:ixAmax1,ixAmin2:ixAmax2)
       end do
       end do
        jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=jcc(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idir)*0.25d0
        jce(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)=jce(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,idir)*jcc(ixCmin1:ixCmax1,ixCmin2:ixCmax2,idir)
      enddo
    end if

    end associate
  end subroutine get_resistive_electric_field

  !> calculate cell-center values from face-center values
  subroutine mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,s)
    use mod_global_parameters
    ! Non-staggered interpolation range
    integer, intent(in)                :: ixOmin1,ixOmin2,ixOmax1,ixOmax2
    type(state)                        :: s

    integer                            :: fxOmin1,fxOmin2,fxOmax1,fxOmax2,&
        gxOmin1,gxOmin2,gxOmax1,gxOmax2, hxOmin1,hxOmin2,hxOmax1,hxOmax2,&
        jxOmin1,jxOmin2,jxOmax1,jxOmax2, kxOmin1,kxOmin2,kxOmax1,kxOmax2, idim

    associate(w=>s%w, ws=>s%ws)

    ! calculate cell-center values from face-center values in 2nd order
    do idim=1,ndim
      ! Displace index to the left
      ! Even if ixI^L is the full size of the w arrays, this is ok
      ! because the staggered arrays have an additional place to the left.
      hxOmin1=ixOmin1-kr(idim,1);hxOmin2=ixOmin2-kr(idim,2)
      hxOmax1=ixOmax1-kr(idim,1);hxOmax2=ixOmax2-kr(idim,2);
      ! Interpolate to cell barycentre using arithmetic average
      ! This might be done better later, to make the method less diffusive.
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(idim))=half/s%surface(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2,idim)*(ws(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idim)*s%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         idim)+ws(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
         idim)*s%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,idim))
    end do

    ! calculate cell-center values from face-center values in 4th order
    !do idim=1,ndim
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);

    !  ! Interpolate to cell barycentre using fourth order central formula
    !  w(ixO^S,mag(idim))=(0.0625d0/s%surface(ixO^S,idim))*&
    !         ( -ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !     +9.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !     +9.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !           -ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) )
    !end do

    ! calculate cell-center values from face-center values in 6th order
    !do idim=1,ndim
    !  fxO^L=ixO^L-3*kr(idim,^D);
    !  gxO^L=ixO^L-2*kr(idim,^D);
    !  hxO^L=ixO^L-kr(idim,^D);
    !  jxO^L=ixO^L+kr(idim,^D);
    !  kxO^L=ixO^L+2*kr(idim,^D);

    !  ! Interpolate to cell barycentre using sixth order central formula
    !  w(ixO^S,mag(idim))=(0.00390625d0/s%surface(ixO^S,idim))* &
    !     (  +3.0d0*ws(fxO^S,idim)*s%surfaceC(fxO^S,idim) &
    !       -25.0d0*ws(gxO^S,idim)*s%surfaceC(gxO^S,idim) &
    !      +150.0d0*ws(hxO^S,idim)*s%surfaceC(hxO^S,idim) &
    !      +150.0d0*ws(ixO^S,idim)*s%surfaceC(ixO^S,idim) &
    !       -25.0d0*ws(jxO^S,idim)*s%surfaceC(jxO^S,idim) &
    !        +3.0d0*ws(kxO^S,idim)*s%surfaceC(kxO^S,idim) )
    !end do

    end associate

  end subroutine mhd_face_to_center

  !> calculate magnetic field from vector potential
  subroutine b_from_vector_potential(ixIsmin1,ixIsmin2,ixIsmax1,ixIsmax2,&
      ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, ws, x)
    use mod_global_parameters
    use mod_constrained_transport

    integer, intent(in)                :: ixIsmin1,ixIsmin2,ixIsmax1,ixIsmax2,&
        ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout)    :: ws(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,1:nws)
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)

    double precision                   :: Adummy(ixIsmin1:ixIsmax1,&
       ixIsmin2:ixIsmax2,1:3)

    call b_from_vector_potentialA(ixIsmin1,ixIsmin2,ixIsmax1,ixIsmax2, ixImin1,&
       ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, ws, x,&
        Adummy)

  end subroutine b_from_vector_potential

end module mod_mhd_phys
