!> Reaction-diffusion module (physics routines)
!>
!> Two types of reaction-diffusion systems are included: a Gray-Scott model and
!> a model credited to Schnakenberg 1979. For the latter, defaults settings are
!> taken from section 4.4 (p. 401) of "Time-dependent advection diffusion
!> reaction systems" by Hundsdorfer & Verwer. 
!>
module mod_rd_phys

  implicit none
  private

  integer, protected, public :: u_ = 1
  integer, protected, public :: v_ = 2

  !> Whether particles module is added
  logical, public, protected              :: rd_particles = .false.

  integer            :: equation_type   = 1
  integer, parameter :: eq_gray_scott   = 1
  integer, parameter :: eq_schnakenberg = 2

  !> Diffusion coefficient for first species (u)
  double precision, public, protected :: D1 = 0.05d0
  !> Diffusion coefficient for second species (v)
  double precision, public, protected :: D2 = 1.0d0

  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_alpha = 0.1305d0
  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_beta  = 0.7695d0
  !> Parameter for Schnakenberg model
  double precision, public, protected :: sb_kappa = 100.0d0

  !> Parameter for Gray-Scott model
  double precision, public, protected :: gs_F = 0.046d0
  !> Parameter for Gray-Scott model
  double precision, public, protected :: gs_k = 0.063d0

  !> Whether to handle the explicitly handled source term in split fashion
  logical :: rd_source_split = .false.

  ! Public methods
  public :: rd_phys_init

contains

  subroutine rd_params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n
    character(len=20)            :: equation_name

    namelist /rd_list/ D1, D2, sb_alpha, sb_beta, sb_kappa, gs_F, gs_k,&
        equation_name, rd_particles, rd_source_split

    equation_name = "gray-scott"

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status='old')
       read(unitpar, rd_list, end=111)
111    close(unitpar)
    end do

    select case (equation_name)
    case ("gray-scott")
       equation_type = eq_gray_scott
    case ("schnakenberg")
       equation_type = eq_schnakenberg
    case default
       call mpistop("Unknown equation_name (not gray-scott or schnakenberg)")
    end select

  end subroutine rd_params_read

  !> Write this module's parameters to a snapsoht
  subroutine rd_write_info(fh)
    use mod_global_parameters
    integer, intent(in)                 :: fh
    integer, parameter                  :: n_par = 0
    integer, dimension(MPI_STATUS_SIZE) :: st
    integer                             :: er
    integer                             :: idim

    call MPI_FILE_WRITE(fh, n_par, 1, MPI_INTEGER, st, er)
  end subroutine rd_write_info

  subroutine rd_phys_init()
    use mod_global_parameters
    use mod_physics
    use mod_multigrid_coupling
    use mod_particles, only: particles_init

    call rd_params_read(par_files)

    physics_type = "rd"
    phys_energy  = .false.
    phys_req_diagonal = .false.
    use_particles = rd_particles

    ! Use the first variable as a density
    u_ = var_set_fluxvar("u", "u")
    v_ = var_set_fluxvar("v", "v")

    ! Disable flux conservation near AMR boundaries, since we have no fluxes
    fix_conserve_global = .false.

    phys_get_cmax     => rd_get_cmax
    phys_get_cbounds  => rd_get_cbounds
    phys_get_flux     => rd_get_flux
    phys_to_conserved => rd_to_conserved
    phys_to_primitive => rd_to_primitive
    phys_add_source   => rd_add_source
    phys_get_dt       => rd_get_dt
    phys_write_info   => rd_write_info
    phys_check_params => rd_check_params
    phys_implicit_update   => rd_implicit_update
    phys_evaluate_implicit => rd_evaluate_implicit

    ! Initialize particles module
    if (rd_particles) then
       call particles_init()
       phys_req_diagonal = .true.
    end if

  end subroutine rd_phys_init

  subroutine rd_check_params
    use mod_multigrid_coupling
    use mod_global_parameters
    integer :: n

    if (any(flux_scheme /= "source")) then
       ! there are no fluxes, only source terms in reaction-diffusion
       call mpistop("mod_rd requires flux_scheme = source")
    end if

    if (use_imex_scheme) then
       use_multigrid = .true.
       ! Set boundary conditions for the multigrid solver
       do n = 1, 2*ndim
          select case (typeboundary(u_, n))
          case ('symm')
             ! d/dx u = 0
             mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          case ('asymm')
             ! u = 0
             mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          case ('cont')
             ! d/dx u = 0
             mg%bc(n, mg_iphi)%bc_type = mg_bc_neumann
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          case ('periodic')
             ! Nothing to do here
          case default
             print *, "divb_multigrid warning: unknown b.c.: ", &
                  trim(typeboundary(u_, n))
             mg%bc(n, mg_iphi)%bc_type = mg_bc_dirichlet
             mg%bc(n, mg_iphi)%bc_value = 0.0_dp
          end select
       end do
    end if

  end subroutine rd_check_params

  subroutine rd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)

    ! Do nothing (primitive and conservative are equal for rd module)
  end subroutine rd_to_conserved

  subroutine rd_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2, nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)

    ! Do nothing (primitive and conservative are equal for rd module)
  end subroutine rd_to_primitive

  subroutine rd_get_cmax(w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, cmax)
    use mod_global_parameters
    integer, intent(in)                       :: ixImin1,ixImin2,ixImax1,&
       ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)              :: w(ixImin1:ixImax1,&
       ixImin2:ixImax2, nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(inout)           :: cmax(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
  end subroutine rd_get_cmax

  subroutine rd_get_cbounds(wLC, wRC, wLp, wRp, x, ixImin1,ixImin2,ixImax1,&
     ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim, cmax, cmin)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: wLC(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRC(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: wLp(ixImin1:ixImax1,ixImin2:ixImax2,&
        nw), wRp(ixImin1:ixImax1,ixImin2:ixImax2,nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(inout) :: cmax(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(inout), optional :: cmin(ixImin1:ixImax1,&
       ixImin2:ixImax2)

    if (present(cmin)) then
       cmin(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
    else
       cmax(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
    end if

  end subroutine rd_get_cbounds

  subroutine rd_get_dt(w, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2, dtnew, dx1,dx2, x)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2, 1:2)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(inout) :: dtnew
    double precision                :: maxrate

    ! dt < dx^2 / (2 * ndim * diffusion_coeff)
    ! use dtdiffpar < 1 for explicit and > 1 for imex/split
    dtnew = dtdiffpar * minval([ dx1,dx2 ])**2 / (2 * ndim * max(D1, D2))

    ! Estimate time step for reactions
    select case (equation_type)
    case (eq_gray_scott)
       maxrate = max(maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, v_))**2 + gs_F,&
           maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, v_) * w(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2, u_)) - gs_F - gs_k)
       dtnew = min(dtnew, 0.5d0/maxrate)
    case (eq_schnakenberg)
       maxrate = max(maxval(abs(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           v_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_) - 1)),&
           maxval(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_))**2)
       dtnew = min(dtnew, 0.5d0/(sb_kappa * maxrate))
    case default
       call mpistop("Unknown equation type")
    end select

  end subroutine rd_get_dt

  ! There is nothing to add to the transport flux in the transport equation
  subroutine rd_get_flux(wC, w, x, ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2, idim, f)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, idim
    double precision, intent(in)    :: wC(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(in)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2, 1:2)
    double precision, intent(out)   :: f(ixImin1:ixImax1,ixImin2:ixImax2,&
        nwflux)
    double precision                :: v(ixImin1:ixImax1,ixImin2:ixImax2)

    f(ixOmin1:ixOmax1,ixOmin2:ixOmax2, :) = 0.0d0
  end subroutine rd_get_flux

  ! w[iw]= w[iw]+qdt*S[wCT, qtC, x] where S is the source based on wCT within ixO
  subroutine rd_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,active)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: qdt
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2, 1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)
    double precision                :: lpl_u(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lpl_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    logical, intent(in)             :: qsourcesplit
    logical, intent(inout)          :: active

    ! here we add the reaction terms (always) and the diffusion if no imex is used
    if (qsourcesplit .eqv. rd_source_split) then
       if (.not.use_imex_scheme) then
          call rd_laplacian(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2, wCT(ixImin1:ixImax1,ixImin2:ixImax2, u_), lpl_u)
          call rd_laplacian(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2, wCT(ixImin1:ixImax1,ixImin2:ixImax2, v_), lpl_v)
       else
          ! for all IMEX scheme variants: only add the reactions
          lpl_u = 0.0d0
          lpl_v = 0.0d0
       end if

       select case (equation_type)
       case (eq_gray_scott)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, u_) + qdt * (D1 * lpl_u - wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, u_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              v_)**2 + gs_F * (1 - wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_)))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, v_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, v_) + qdt * (D2 * lpl_v + wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, u_) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              v_)**2 - (gs_F + gs_k) * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              v_))
       case (eq_schnakenberg)
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, u_) + qdt * (D1 * lpl_u + sb_kappa * (sb_alpha - &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_) + wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, u_)**2 * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              v_)))
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, v_) = w(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, v_) + qdt * (D2 * lpl_v + sb_kappa * (sb_beta - &
             wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2, u_)**2 * wCT(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2, v_)))
       case default
          call mpistop("Unknown equation type")
       end select

       ! enforce getbc call after source addition
       active = .true.
    end if

  end subroutine rd_add_source

  !> Compute the Laplacian using a standard second order scheme. For now this
  !> method only works in slab geometries. Requires one ghost cell only.
  subroutine rd_laplacian(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,var,lpl)
    use mod_global_parameters
    integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
       ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)  :: var(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision, intent(out) :: lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    integer                       :: idir, jxOmin1,jxOmin2,jxOmax1,jxOmax2,&
        hxOmin1,hxOmin2,hxOmax1,hxOmax2
    double precision              :: h_inv2

    if (slab) then
       lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = 0.0d0
       do idir = 1, ndim
          hxOmin1=ixOmin1-kr(idir,1);hxOmin2=ixOmin2-kr(idir,2)
          hxOmax1=ixOmax1-kr(idir,1);hxOmax2=ixOmax2-kr(idir,2);
          jxOmin1=ixOmin1+kr(idir,1);jxOmin2=ixOmin2+kr(idir,2)
          jxOmax1=ixOmax1+kr(idir,1);jxOmax2=ixOmax2+kr(idir,2);
          h_inv2 = 1/dxlevel(idir)**2
          lpl(ixOmin1:ixOmax1,ixOmin2:ixOmax2) = lpl(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) + h_inv2 * (var(jxOmin1:jxOmax1,&
             jxOmin2:jxOmax2) - 2 * var(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2) + var(hxOmin1:hxOmax1,hxOmin2:hxOmax2))
       end do
    else
       call mpistop("rd_laplacian not implemented in this geometry")
    end if

  end subroutine rd_laplacian

  subroutine put_laplacians_onegrid(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,w)
    use mod_global_parameters
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
        1:nw)

    double precision                :: lpl_u(ixOmin1:ixOmax1,ixOmin2:ixOmax2),&
        lpl_v(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

    call rd_laplacian(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2, w(ixImin1:ixImax1,ixImin2:ixImax2, u_), lpl_u)
    call rd_laplacian(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2, w(ixImin1:ixImax1,ixImin2:ixImax2, v_), lpl_v)

    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,u_)=D1*lpl_u
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,v_)=D2*lpl_v

  end subroutine put_laplacians_onegrid
  
  !> inplace update of psa==>F_im(psa)
  subroutine rd_evaluate_implicit(qtC,psa)
    use mod_global_parameters
    type(state), target :: psa(max_blocks)
    double precision, intent(in) :: qtC

    integer :: iigrid, igrid, level
    integer :: ixOmin1,ixOmin2,ixOmax1,ixOmax2

    ixOmin1=ixGlo1+1;ixOmin2=ixGlo2+1;ixOmax1=ixGhi1-1;ixOmax2=ixGhi2-1;
    !$OMP PARALLEL DO PRIVATE(igrid)
    do iigrid=1,igridstail; igrid=igrids(iigrid);
       dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);
       call put_laplacians_onegrid(ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,psa(igrid)%w)
    end do
    !$OMP END PARALLEL DO

  end subroutine rd_evaluate_implicit

  !> Implicit solve of psa=psb+dtfactor*dt*F_im(psa)
  subroutine rd_implicit_update(dtfactor,qdt,qtC,psa,psb)
    use mod_global_parameters
    use mod_forest
    use mod_multigrid_coupling

    type(state), target :: psa(max_blocks)
    type(state), target :: psb(max_blocks)
    double precision, intent(in) :: qdt
    double precision, intent(in) :: qtC
    double precision, intent(in) :: dtfactor

    integer                      :: n
    double precision             :: res, max_residual, lambda

    integer                        :: iw_to,iw_from
    integer                        :: iigrid, igrid, id
    integer                        :: nc, lvl
    type(tree_node), pointer       :: pnode
    real(dp)                       :: fac

    ! Avoid setting a very restrictive limit to the residual when the time step
    ! is small (as the operator is ~ 1/(D * qdt))
    if (qdt < dtmin) then
        if(mype==0)then
            print *,'skipping implicit solve: dt too small!'
            print *,'Currently at time=',global_time,' time step=',qdt,' dtmin=',dtmin
        endif
        return
    endif
    max_residual = 1d-7/qdt

    mg%operator_type = mg_helmholtz
    call mg_set_methods(mg)

    if (.not. mg%is_allocated) call mpistop(&
       "multigrid tree not allocated yet")

    ! First handle the u variable ***************************************
    lambda           = 1/(dtfactor * qdt * D1)
    call helmholtz_set_lambda(lambda)

    call mg_copy_to_tree(u_, mg_irhs, factor=-lambda, state_from=psb)
    call mg_copy_to_tree(u_, mg_iphi, state_from=psb)

    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, 10
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    call mg_copy_from_tree_gc(mg_iphi, u_, state_to=psa)
    ! Done with the u variable ***************************************

    ! Next handle the v variable ***************************************
    lambda = 1/(dtfactor * qdt * D2)
    call helmholtz_set_lambda(lambda)

    call mg_copy_to_tree(v_, mg_irhs, factor=-lambda, state_from=psb)
    call mg_copy_to_tree(v_, mg_iphi, state_from=psb)

    call mg_fas_fmg(mg, .true., max_res=res)
    do n = 1, 10
       call mg_fas_vcycle(mg, max_res=res)
       if (res < max_residual) exit
    end do

    call mg_copy_from_tree_gc(mg_iphi, v_, state_to=psa)
    ! Done with the v variable ***************************************

  end subroutine rd_implicit_update

end module mod_rd_phys
