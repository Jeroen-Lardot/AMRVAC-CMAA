!> Module for handling split source terms (split from the fluxes)
module mod_source

  implicit none
  private

  public :: add_split_source
  public :: addsource2

contains

  subroutine add_split_source(prior)
    use mod_global_parameters
    use mod_ghostcells_update
    use mod_thermal_conduction, only: phys_thermal_conduction
    use mod_physics, only: phys_req_diagonal, phys_global_source_before,&
        phys_global_source_after

    logical, intent(in) :: prior

    double precision :: qdt, qt
    integer :: iigrid, igrid, i1,i2
    logical :: src_active

    integer :: i, itdiff, Ndiff

    ! add thermal conduction
    if(associated(phys_thermal_conduction)) call phys_thermal_conduction()

    src_active = .false.

    if (prior .and. associated(phys_global_source_before)) then
       call phys_global_source_before(dt, qt, src_active)
    end if

    if ((.not.prior).and.(typesourcesplit=='sf' .or. typesourcesplit=='ssf')) &
       return

    if (prior) then
       qt=global_time
    else
       qt=global_time+dt
    end if

    !$OMP PARALLEL DO PRIVATE(igrid,qdt,i1,i2)
    do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
       qdt=dt_grid(igrid)
       block=>ps(igrid)
       call addsource1_grid(igrid,qdt,qt,ps(igrid)%w,src_active)
    end do
    !$OMP END PARALLEL DO

    if (.not. prior .and. associated(phys_global_source_after)) then
       call phys_global_source_after(dt, qt, src_active)
    end if

    if (src_active) then
       call getbc(qt,0.d0,ps,1,nwflux+nwaux,phys_req_diagonal)
    end if

  end subroutine add_split_source

  subroutine addsource1_grid(igrid,qdt,qt,w,src_active)

    use mod_global_parameters

    integer, intent(in) :: igrid
    double precision, intent(in) :: qdt, qt
    double precision, intent(inout) :: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
    logical, intent(inout) :: src_active

    double precision :: w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

    saveigrid=igrid
    typelimiter=type_limiter(node(plevel_,igrid))
    typegradlimiter=type_gradient_limiter(node(plevel_,igrid))

    dxlevel(1)=rnode(rpdx1_,igrid);dxlevel(2)=rnode(rpdx2_,igrid);

    w1(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)

    select case (typesourcesplit)
    case ('sf')
       call addsource2(qdt  ,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,1,nw,qt,w1,qt,w,ps(igrid)%x,.true.,src_active)
    case ('sfs')
       call addsource2(qdt/2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,1,nw,qt,w1,qt,w,ps(igrid)%x,.true.,src_active)
    case ('ssf')
       call addsource2(qdt/2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,ixGhi1,&
          ixGhi2,1,nw,qt,w,qt,w1,ps(igrid)%x,.true.,src_active)
       call addsource2(qdt  ,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,1,nw,qt,w1,qt,w,ps(igrid)%x,.true.,src_active)
    case ('ssfss')
       call addsource2(qdt/4,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixGlo1,ixGlo2,ixGhi1,&
          ixGhi2,1,nw,qt,w,qt,w1,ps(igrid)%x,.true.,src_active)
       call addsource2(qdt/2,ixGlo1,ixGlo2,ixGhi1,ixGhi2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,1,nw,qt,w1,qt,w,ps(igrid)%x,.true.,src_active)
    case default
       write(unitterm,*)'No such typesourcesplit=',typesourcesplit
       call mpistop("Error: Unknown typesourcesplit!")
    end select

  end subroutine addsource1_grid

  !> Add source within ixO for iws: w=w+qdt*S[wCT]
  subroutine addsource2(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x,qsourcesplit,src_active)
    use mod_global_parameters
    use mod_physics, only: phys_add_source
    use mod_usr_methods, only: usr_source
    ! differences with VAC is in iw^LIM and in declaration of ranges for wCT,w

    integer, intent(in)              :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in)     :: qdt, qtC, qt
    double precision, intent(in)     :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout)  :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    logical, intent(in)              :: qsourcesplit
    logical, intent(inout), optional :: src_active
    logical                          :: tmp_active

    tmp_active = .false.

    ! physics defined sources, typically explicitly added,
    ! along with geometrical source additions
    call phys_add_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,wCT,w,x,qsourcesplit,tmp_active)

    ! user defined sources, typically explicitly added
    if ((qsourcesplit .eqv. source_split_usr) .and. associated(usr_source)) &
       then
       tmp_active = .true.
       call usr_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
          ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    end if

    if (present(src_active)) src_active = src_active .or. tmp_active

  end subroutine addsource2

end module mod_source
