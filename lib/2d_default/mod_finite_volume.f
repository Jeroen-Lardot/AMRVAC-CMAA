!> Module with finite volume methods for fluxes
module mod_finite_volume
  implicit none
  private

  public :: finite_volume
  public :: hancock
  public :: reconstruct_LR

contains

  !> The non-conservative Hancock predictor for TVDLF
  !>
  !> on entry:
  !> input available on ixI^L=ixG^L asks for output on ixO^L=ixG^L^LSUBnghostcells
  !> one entry: (predictor): wCT -- w_n        wnew -- w_n   qdt=dt/2
  !> on exit :  (predictor): wCT -- w_n        wnew -- w_n+1/2
  subroutine hancock(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,idimsmin,idimsmax,qtC,sCT,qt,snew,dx1,dx2,x)
    use mod_physics
    use mod_global_parameters
    use mod_source, only: addsource2

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, idimsmin,idimsmax
    double precision, intent(in) :: qdt, qtC, qt, dx1,dx2, x(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:ndim)
    type(state) :: sCT, snew

    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim,&
        wLC, wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision :: fLC(ixImin1:ixImax1,ixImin2:ixImax2, nwflux),&
        fRC(ixImin1:ixImax1,ixImin2:ixImax2, nwflux)
    double precision :: dxinv(1:ndim),dxdim(1:ndim)
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2

    associate(wCT=>sCT%w,wnew=>snew%w)
    ! Expand limits in each idims direction in which fluxes are added
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-kr(idims,1);ixmin2=ixmin2-kr(idims,2)
       ixmax1=ixmax1+kr(idims,1);ixmax2=ixmax2+kr(idims,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in Hancock: Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
    dxdim(1)=dx1;dxdim(2)=dx2;
    do idims= idimsmin,idimsmax
      block%iw0=idims
      ! Calculate w_j+g_j/2 and w_j-g_j/2
      ! First copy all variables, then upwind wLC and wRC.
      ! wLC is to the left of ixO, wRC is to the right of wCT.
      hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
      hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

      wRp(hxOmin1:hxOmax1,hxOmin2:hxOmax2,1:nwflux)=wprim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:nwflux)
      wLp(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1:nwflux)=wprim(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1:nwflux)

      ! apply limited reconstruction for left and right status at cell interfaces
      call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,hxOmin1,hxOmin2,hxOmax1,hxOmax2,idims,wprim,wLC,wRC,&
         wLp,wRp,x,dxdim(idims))

      ! Calculate the fLC and fRC fluxes
      call phys_get_flux(wRC,wRp,x,ixImin1,ixImin2,ixImax1,ixImax2,hxOmin1,&
         hxOmin2,hxOmax1,hxOmax2,idims,fRC)
      call phys_get_flux(wLC,wLp,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,idims,fLC)

      ! Advect w(iw)
      do iw=1,nwflux
        if (slab_uniform) then
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,&
               ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, wnew, inv_volume)
            wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2,iw)+dxinv(idims) * inv_volume * &
               (fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2, iw)-fRC(hxOmin1:hxOmax1,&
               hxOmin2:hxOmax2, iw))
           else
             wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2,iw)+dxinv(idims)* (fLC(ixOmin1:ixOmax1,&
                ixOmin2:ixOmax2, iw)-fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2, iw))
           end if
        else
          if (associated(phys_iw_methods(iw)%inv_capacity)) then
            call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,&
               ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, wnew, inv_volume)
          else
            inv_volume = 1.0d0
          end if
          inv_volume = inv_volume/block%dvolume(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2)

          wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,iw) - qdt * inv_volume &
             *(block%surfaceC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
             idims)*fLC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
              iw) -block%surfaceC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
             idims)*fRC(hxOmin1:hxOmax1,hxOmin2:hxOmax2, iw))
        end if
      end do
    end do ! next idims
    block%iw0=0

    do iw = 1, nwflux
      if (associated(phys_iw_methods(iw)%inv_capacity)) then
        ! Copy state before adding source terms
        wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, iw) = wnew(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, iw)
      end if
    end do

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.)

    ! If there are capacity functions, now correct the added source terms
    do iw = 1, nwflux
      if (associated(phys_iw_methods(iw)%inv_capacity)) then
        call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,ixImax2,&
            ixOmin1,ixOmin2,ixOmax1,ixOmax2, wnew, inv_volume)
        wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2, iw) = wprim(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, iw) + inv_volume * (wnew(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2, iw) - wprim(ixOmin1:ixOmax1,ixOmin2:ixOmax2, iw))
      end if
    end do

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'finite_volume')
    end associate
  end subroutine hancock

  !> finite volume method
  subroutine finite_volume(method,qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,idimsmin,idimsmax, qtC,sCT,qt,snew,sold,fC,fE,dx1,&
     dx2,x)
    use mod_physics
    use mod_global_parameters
    use mod_tvd, only:tvdlimit2
    use mod_source, only: addsource2
    use mod_usr_methods

    character(len=*), intent(in)                          :: method
    double precision, intent(in)                          :: qdt, qtC, qt, dx1,&
       dx2
    integer, intent(in)                                   :: ixImin1,ixImin2,&
       ixImax1,ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, idimsmin,idimsmax
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        intent(in) ::  x
    type(state)                                           :: sCT, snew, sold
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
       1:ndim)    :: fC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       7-2*ndim:3)         :: fE

    ! primitive w at cell center
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wprim
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form, needed for better performance
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nwflux) :: fLC, fRC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cmaxC
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2)      :: cminC
    double precision, dimension(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)      :: inv_volume
    double precision, dimension(1:ndim)     :: dxinv, dxdim
    ! cell-face location coordinates
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: xi
    integer, dimension(ixImin1:ixImax1,ixImin2:ixImax2)               :: &
       patchf
    integer :: idims, iw, ixmin1,ixmin2,ixmax1,ixmax2, hxOmin1,hxOmin2,hxOmax1,&
       hxOmax2, ixCmin1,ixCmin2,ixCmax1,ixCmax2, ixCRmin1,ixCRmin2,ixCRmax1,&
       ixCRmax2, kxCmin1,kxCmin2,kxCmax1,kxCmax2, kxRmin1,kxRmin2,kxRmax1,&
       kxRmax2

    associate(wCT=>sCT%w, wnew=>snew%w, wold=>sold%w)

    fC=0.d0

    ! The flux calculation contracts by one in the idims direction it is applied.
    ! The limiter contracts the same directions by one more, so expand ixO by 2.
    ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
    do idims= idimsmin,idimsmax
       ixmin1=ixmin1-2*kr(idims,1);ixmin2=ixmin2-2*kr(idims,2)
       ixmax1=ixmax1+2*kr(idims,1);ixmax2=ixmax2+2*kr(idims,2);
    end do
    if (ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2<ixmax2) &
       call mpistop("Error in fv : Nonconforming input limits")

    wprim=wCT
    call phys_to_primitive(ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,ixImin2,&
       ixImax1,ixImax2,wprim,x)

    dxinv(1)=-qdt/dx1;dxinv(2)=-qdt/dx2;
    dxdim(1)=dx1;dxdim(2)=dx2;
    do idims= idimsmin,idimsmax
       ! use interface value of w0 at idims
       block%iw0=idims

       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       kxCmin1=ixImin1;kxCmin2=ixImin2; kxCmax1=ixImax1-kr(idims,1)
       kxCmax2=ixImax2-kr(idims,2);
       kxRmin1=kxCmin1+kr(idims,1);kxRmin2=kxCmin2+kr(idims,2)
       kxRmax1=kxCmax1+kr(idims,1);kxRmax2=kxCmax2+kr(idims,2);

       if(stagger_grid) then
          ! ct needs all transverse cells
          ixCmax1=ixOmax1+nghostcells-nghostcells*kr(idims,1)
          ixCmax2=ixOmax2+nghostcells-nghostcells*kr(idims,2)
          ixCmin1=hxOmin1-nghostcells+nghostcells*kr(idims,1)
          ixCmin2=hxOmin2-nghostcells+nghostcells*kr(idims,2);
       else
         ! ixC is centered index in the idims direction from ixOmin-1/2 to ixOmax+1/2
         ixCmax1=ixOmax1;ixCmax2=ixOmax2; ixCmin1=hxOmin1;ixCmin2=hxOmin2;
       end if

       ! wRp and wLp are defined at the same locations, and will correspond to
       ! the left and right reconstructed values at a cell face. Their indexing
       ! is similar to cell-centered values, but in direction idims they are
       ! shifted half a cell towards the 'lower' direction.
       wRp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxRmin1:kxRmax1,&
          kxRmin2:kxRmax2,1:nw)
       wLp(kxCmin1:kxCmax1,kxCmin2:kxCmax2,1:nw)=wprim(kxCmin1:kxCmax1,&
          kxCmin2:kxCmax2,1:nw)

       ! Determine stencil size
       ixCRmin1 = max(ixCmin1 - phys_wider_stencil,ixGlo1)
       ixCRmin2 = max(ixCmin2 - phys_wider_stencil,ixGlo2)
       ixCRmax1 = min(ixCmax1 + phys_wider_stencil,ixGhi1)
       ixCRmax2 = min(ixCmax2 + phys_wider_stencil,ixGhi2)

       ! get cell-face coordinates
       xi=x
       xi(ixImin1:ixImax1,ixImin2:ixImax2,idims)=xi(ixImin1:ixImax1,&
          ixImin2:ixImax2,idims)+0.5d0*sCT%dx(ixImin1:ixImax1,ixImin2:ixImax2,&
          idims)

       ! apply limited reconstruction for left and right status at cell interfaces
       call reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
          ixCRmax1,ixCRmax2,ixCRmin1,ixCRmin2,ixCRmax1,ixCRmax2,idims,wprim,&
          wLC,wRC,wLp,wRp,xi,dxdim(idims))

       ! special modification of left and right status before flux evaluation
       call phys_modify_wLR(ixImin1,ixImin2,ixImax1,ixImax2,ixCRmin1,ixCRmin2,&
          ixCRmax1,ixCRmax2,qt,wLC,wRC,wLp,wRp,sCT,idims)

       ! evaluate physical fluxes according to reconstructed status
       call phys_get_flux(wLC,wLp,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2,idims,fLC)
       call phys_get_flux(wRC,wRp,xi,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,&
          ixCmin2,ixCmax1,ixCmax2,idims,fRC)

       ! estimating bounds for the minimum and maximum signal velocities
       if(method=='tvdlf'.or.method=='tvdmu') then
         call phys_get_cbounds(wLC,wRC,wLp,wRp,xi,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC)
       else
         call phys_get_cbounds(wLC,wRC,wLp,wRp,xi,ixImin1,ixImin2,ixImax1,&
            ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,cmaxC,cminC)
       end if

       ! use approximate Riemann solver to get flux at interfaces
       select case(method)
       case('tvdmu')
         call get_Riemann_flux_tvdmu()
       case('tvdlf')
         call get_Riemann_flux_tvdlf()
       case('hll')
         call get_Riemann_flux_hll()
       case('hllc','hllcd')
         call get_Riemann_flux_hllc()
       case('hlld')
         call get_Riemann_flux_hlld()
       case default
         call mpistop('unkown Riemann flux')
       end select

       if(associated(usr_set_flux)) call usr_set_flux(ixImin1,ixImin2,ixImax1,&
          ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,qt,wLC,wRC,wLp,wRp,sCT,idims,&
          fC)

    end do ! Next idims
    block%iw0=0

    if(stagger_grid) call phys_update_faces(ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,qdt,wprim,fC,fE,sCT,snew)

    do idims= idimsmin,idimsmax
       hxOmin1=ixOmin1-kr(idims,1);hxOmin2=ixOmin2-kr(idims,2)
       hxOmax1=ixOmax1-kr(idims,1);hxOmax2=ixOmax2-kr(idims,2);

       ! Multiply the fluxes by -dt/dx since Flux fixing expects this
       if (slab_uniform) then
          fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
             idims)=dxinv(idims)*fC(ixImin1:ixImax1,ixImin2:ixImax2,1:nwflux,&
             idims)

          do iw = 1, nwflux
            if (associated(phys_iw_methods(iw)%inv_capacity)) then
              call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,&
                 ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, wnew, inv_volume)
              wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw) + inv_volume * (fC(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw,idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,&
                 iw,idims))
            else
              wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                 idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,idims))
            end if
          end do
       else
          if (.not. angmomfix) then ! default case
            if (associated(phys_iw_methods(iw)%inv_capacity)) then
              call phys_iw_methods(iw)%inv_capacity(ixImin1,ixImin2,ixImax1,&
                 ixImax2, ixOmin1,ixOmin2,ixOmax1,ixOmax2, wnew, inv_volume)
            else
              inv_volume = 1.0d0
            end if
            inv_volume = inv_volume/block%dvolume(ixOmin1:ixOmax1,&
               ixOmin2:ixOmax2)

            do iw=1,nwflux
              fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                 idims)=-qdt*fC(ixImin1:ixImax1,ixImin2:ixImax2,iw,&
                 idims)*block%surfaceC(ixImin1:ixImax1,ixImin2:ixImax2,idims)
              wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
                 ixOmin2:ixOmax2,iw) + (fC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw,&
                 idims)-fC(hxOmin1:hxOmax1,hxOmin2:hxOmax2,iw,&
                 idims)) * inv_volume
            enddo
          else
            ! If angular momentum conserving way to solve the equations,
            ! some fluxes additions need to be treated specifically
            call phys_angmomfix(fC,x,wnew,ixImin1,ixImin2,ixImax1,ixImax2,&
               ixOmin1,ixOmin2,ixOmax1,ixOmax2,idims)
          endif
       end if

       ! For the MUSCL scheme apply the characteristic based limiter
       if (method=='tvdmu') call tvdlimit2(method,qdt,ixImin1,ixImin2,ixImax1,&
          ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,ixOmin1,ixOmin2,ixOmax1,&
          ixOmax2,idims,wLC,wRC,wnew,x,fC,dx1,dx2)

    end do ! Next idims

    if (.not.slab.and.idimsmin==1) call phys_add_source_geom(qdt,ixImin1,&
       ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)

    if(stagger_grid) call phys_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
       snew)
 
    if(phys_solve_eaux) then
      call phys_energy_synchro(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
         ixOmin2,ixOmax1,ixOmax2,wCT,wnew,x)
    endif

    call addsource2(qdt*dble(idimsmax-idimsmin+1)/dble(ndim), ixImin1,ixImin2,&
       ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,1,nw,qtC,wCT,qt,wnew,x,&
       .false.)

    ! check and optionally correct unphysical values
    call phys_handle_small_values(.false.,wnew,x,ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,'finite_volume')

  end associate
  contains

    subroutine get_Riemann_flux_tvdmu()

      do iw=1,nwflux
         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=half*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))
         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)
      end do
    end subroutine get_Riemann_flux_tvdmu

    subroutine get_Riemann_flux_tvdlf()
      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      fac = -0.5d0*tvdlfeps*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux

         ! To save memory we use fLC to store (F_L+F_R)/2=half*(fLC+fRC)
         fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=0.5d0*(fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw))

         ! Add TVDLF dissipation to the flux
         if (flux_type(idims, iw) /= flux_no_dissipation) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2, iw) + fac*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
         end if

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)
      end do ! Next iw
    end subroutine get_Riemann_flux_tvdlf

    subroutine get_Riemann_flux_hll()

      double precision :: fac(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
          div(ixCmin1:ixCmax1,ixCmin2:ixCmax2)

      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      elsewhere
        patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
      endwhere

      fac = tvdlfeps*cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      div = 1/(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-cminC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))

      ! Calculate fLC=f(uL_j+1/2) and fRC=f(uR_j+1/2) for each iw
      do iw=1,nwflux
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) = half*(fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                iw) -tvdlfeps*max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                dabs(cminC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)
               ! Add hll dissipation to the flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw) = (cmaxC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw)-cminC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2) * fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                   iw) +fac*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))) * div
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)== 2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2, iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2, iw)
            endwhere
         endif

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2, iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hll

    subroutine get_Riemann_flux_hllc()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux)     :: whll, Fhll, fCD
      double precision, dimension(ixImin1:ixImax1,&
         ixImin2:ixImax2)              :: lambdaCD

      patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  1
      where(cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) >= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = -2
      elsewhere(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) <= zero)
         patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2) =  2
      endwhere
      ! Use more diffusive scheme, is actually TVDLF and selected by patchf=4
      if(method=='hllcd') call phys_diffuse_hllcd(ixImin1,ixImin2,ixImax1,&
         ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2,idims,wLC,wRC,fLC,fRC,patchf)

      !---- calculate speed lambda at CD ----!
      if(any(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==1)) call &
         phys_get_lCD(wLC,wRC,fLC,fRC,cminC,cmaxC,idims,ixImin1,ixImin2,&
         ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,ixCmax2, whll,Fhll,lambdaCD,&
         patchf)

      ! now patchf may be -1 or 1 due to phys_get_lCD
      if(any(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))== 1))then
         !======== flux at intermediate state ========!
         call phys_get_wCD(wLC,wRC,whll,fRC,fLC,Fhll,patchf,lambdaCD,cminC,&
            cmaxC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,ixCmax1,&
            ixCmax2,idims,fCD)
      endif ! Calculate the CD flux

      ! use hll flux for the auxiliary internal e
      if(mhd_energy.and.mhd_solve_eaux) then
        iw=eaux_
        fCD(ixCmin1:ixCmax1,ixCmin2:ixCmax2, iw) = (cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            iw)-cminC(ixCmin1:ixCmax1,ixCmin2:ixCmax2) * fRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2, iw) +cminC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*cmaxC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           iw)))/(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-cminC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))
      end if

      do iw=1,nwflux
         if (flux_type(idims, iw) == flux_tvdlf) then
            fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) = 0.5d0 * (fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                abs(cminC(ixCmin1:ixCmax1,&
               ixCmin2:ixCmax2))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
               iw) - wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
         else
            where(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==-2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(abs(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2))==1)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fCD(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==2)
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==3)
               ! fallback option, reducing to HLL flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=Fhll(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2,iw)
            elsewhere(patchf(ixCmin1:ixCmax1,ixCmin2:ixCmax2)==4)
               ! fallback option, reducing to TVDLF flux
               fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw) = half*((fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)+fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)) -tvdlfeps * max(cmaxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
                   dabs(cminC(ixCmin1:ixCmax1,&
                  ixCmin2:ixCmax2))) * (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                  iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))
            endwhere
         end if

         fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw,idims)=fLC(ixCmin1:ixCmax1,&
            ixCmin2:ixCmax2,iw)

      end do ! Next iw
    end subroutine get_Riemann_flux_hllc

    !> HLLD Riemann flux from Miyoshi 2005 JCP, 208, 315 and Guo 2016 JCP, 327, 543
    subroutine get_Riemann_flux_hlld()
      use mod_mhd_phys
      implicit none
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w1R,w1L,f1R,f1L,f2R,f2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
         1:nwflux) :: w2R,w2L
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: sm,s1R,&
         s1L,suR,suL,Bx
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2) :: pts,ptR,&
         ptL,signBx,r1L,r1R,tmp
      ! velocity from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: vRC,&
          vLC
      ! magnetic field from the right and the left reconstruction
      double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,ndir) :: BR,&
          BL
      integer :: ip1,ip2,ip3,idir,ix1,ix2

      associate (sR=>cmaxC,sL=>cminC)

      f1R=0.d0
      f1L=0.d0
      ip1=idims
      ip3=3
      vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLp(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(:))
      if(B0field) then
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))+block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
      else
        BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))
        BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)=wLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))
      end if
      ! HLL estimation of normal magnetic field at cell interfaces
      Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1)-fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip1))-fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip1)))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))
      ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         p_)+0.5d0*sum(BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)**2,dim=ndim+1)
      ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         p_)+0.5d0*sum(BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:)**2,dim=ndim+1)
      suR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)
      suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)
      ! Miyoshi equation (38) and Guo euqation (20)
      sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1)-ptR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+ptL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-suL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      ! Miyoshi equation (39) and Guo euqation (28)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip1))=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      ! Guo equation (22)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(B0field) then
        ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wRp(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)+0.5d0*sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(:))**2,dim=ndim+1)
        ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=wLp(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,p_)+0.5d0*sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(:))**2,dim=ndim+1)
      end if

      ! Miyoshi equation (43) and Guo equation (27)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)/(sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))

      ip2=mod(ip1+1,ndir)
      if(ip2==0) ip2=ndir
      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2
      where(r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=0.d0)
        r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      endwhere
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2
      where(r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)/=0.d0)
        r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)
      endwhere
      ! Miyoshi equation (44)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)-Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip2)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! partial solution for later usage
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(suL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(sL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         ip1))-Bx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2)*r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)
      if(ndir==3) then
        ip3=mod(ip1+2,ndir)
        if(ip3==0) ip3=ndir
        ! Miyoshi equation (46)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*BR(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)-Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*BL(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip3)*(sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip1))*r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        ! Miyoshi equation (47)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=BR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=BL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,ip3)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      end if
      ! Miyoshi equation (45)
      w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=BR(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=BL(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,ip2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))
      if(B0field) then
        ! Guo equation (26)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(:))-block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,&
           ip1)
      end if
      ! equation (48)
      if(mhd_energy) then
        ! Guo equation (25) equivalent to Miyoshi equation (41)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=suR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))+ptR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=suL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))+ptL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        if(B0field) then
          ! Guo equation (32)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             :,ip1)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))),dim=ndim+1)
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)=w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,p_)+sum(block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             :,ip1)*(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:))),dim=ndim+1)
        end if
        ! Miyoshi equation (48) and main part of Guo euqation (31)
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=((sL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1))*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)-ptL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           ip1)+w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,p_)*sm(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)+Bx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           :)*wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)))/(sL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        if(B0field) then
          ! Guo equation (31)
          w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_)+(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-sum(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*vRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1))/(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))
          w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,e_)+(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)-sum(wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             mag(:))*block%B0(ixCmin1:ixCmax1,ixCmin2:ixCmax2,:,ip1),&
             dim=ndim+1)*vLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             ip1))/(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)-sm(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2))
        end if
      end if

      ! Miyoshi equation (49) and Guo equation (35)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,rho_)=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip1))=w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip1))

      r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt(w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,rho_))
      tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=1.d0/(r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
      signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sign(1.d0,Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))
      ! Miyoshi equation (51) and Guo equation (33)
      s1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)+abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      s1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sm(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)-abs(Bx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))/r1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      ! Miyoshi equation (59) and Guo equation (41)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mom(ip2))
      ! Miyoshi equation (61) and Guo equation (43)
      w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=(r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mag(ip2))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))+r1L(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         mom(ip2)))*signBx(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip2))=w2R(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2,mag(ip2))
      if(ndir==3) then
        ! Miyoshi equation (60) and Guo equation (42)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))+(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(ip3))
        ! Miyoshi equation (62) and Guo equation (44)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=(r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mag(ip3))+r1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))+r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3))-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(ip3)))*signBx(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2))*tmp(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(ip3))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mag(ip3))
      end if
      ! Miyoshi equation (63) and Guo equation (45)
      if(mhd_energy) then
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)+r1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,e_)=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,e_)-r1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(sum(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1)-sum(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           mom(:))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mag(:)),&
           dim=ndim+1))*signBx(ixCmin1:ixCmax1,ixCmin2:ixCmax2)
      end if

      ! convert velocity to momentum
      do idir=1,ndir
        w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2R(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
        w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,mom(idir))=w2L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,mom(idir))*w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           rho_)
      end do

      ! get fluxes of intermedate states
      do iw=1,nwflux
        if (flux_type(idims, iw) == flux_tvdlf) then
          !! hll flux for normal B
          !f1L(ixC^S,iw)=(sR(ixC^S)*fLC(ixC^S, iw)-sL(ixC^S)*fRC(ixC^S, iw) &
          !          +sR(ixC^S)*sL(ixC^S)*(wRC(ixC^S,iw)-wLC(ixC^S,iw)))/(sR(ixC^S)-sL(ixC^S))
          ! tvldf flux for normal B
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)= half*(fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2, iw) + fRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
              iw) -tvdlfeps*max(sR(ixCmin1:ixCmax1,ixCmin2:ixCmax2),&
              dabs(sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))) * &
             (wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)-wLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else if(flux_type(idims, iw) == flux_special) then
          ! known flux (fLC=fRC) for normal B and psi_ in GLM method
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)
        else
          f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fLC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+sL(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=fRC(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+sR(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1L(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
          f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2,iw)+s1R(ixCmin1:ixCmax1,&
             ixCmin2:ixCmax2)*(w2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
             iw)-w1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw))
        end if
      end do

      ! use hll flux for the auxiliary internal e
      if(mhd_energy.and.mhd_solve_eaux) then
        iw=eaux_
        f1L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=(sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*fLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            iw)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2)*fRC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2, iw) +sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*sL(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)*(wRC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           iw)-wLC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)))/(sR(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2)-sL(ixCmin1:ixCmax1,ixCmin2:ixCmax2))
        f1R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
        f2L(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
        f2R(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)=f1L(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,iw)
      end if

      ! Miyoshi equation (66) and Guo equation (46)
     do ix2=ixCmin2,ixCmax2
     do ix1=ixCmin1,ixCmax1
        if(sL(ix1,ix2)>0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=fLC(ix1,ix2,1:nwflux)
        else if(s1L(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f1L(ix1,ix2,1:nwflux)
        else if(sm(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f2L(ix1,ix2,1:nwflux)
        else if(s1R(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f2R(ix1,ix2,1:nwflux)
        else if(sR(ix1,ix2)>=0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=f1R(ix1,ix2,1:nwflux)
        else if(sR(ix1,ix2)<0.d0) then
          fC(ix1,ix2,1:nwflux,ip1)=fRC(ix1,ix2,1:nwflux)
        end if
     end do
     end do

      end associate
    end subroutine get_Riemann_flux_hlld

  end subroutine finite_volume

  !> Determine the upwinded wLC(ixL) and wRC(ixR) from w.
  !> the wCT is only used when PPM is exploited.
  subroutine reconstruct_LR(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
     ixLmax1,ixLmax2,ixRmin1,ixRmin2,ixRmax1,ixRmax2,idims,w,wLC,wRC,wLp,wRp,x,&
     dxdim)
    use mod_physics
    use mod_global_parameters
    use mod_limiter

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2, ixRmin1,ixRmin2,ixRmax1,ixRmax2, idims
    double precision, intent(in) :: dxdim
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: w
    ! left and right constructed status in conservative form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLC,&
        wRC
    ! left and right constructed status in primitive form
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:nw) :: wLp,&
        wRp
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim) :: x

    integer            :: jxRmin1,jxRmin2,jxRmax1,jxRmax2, ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2, jxCmin1,jxCmin2,jxCmax1,jxCmax2, iw
    double precision   :: ldw(ixImin1:ixImax1,ixImin2:ixImax2),&
        rdw(ixImin1:ixImax1,ixImin2:ixImax2), dwC(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision   :: a2max

    select case (typelimiter)
    case (limiter_venk)
       call venklimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp)
    case (limiter_mp5)
       call MP5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,ixLmax1,&
          ixLmax2,idims,w,wLp,wRp)
    case (limiter_weno3)
       call WENO3limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,1)
    case (limiter_wenoyc3)
       call WENO3limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,2)
    case (limiter_weno5)
       call WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,1)
    case (limiter_weno5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,1)
    case (limiter_wenoz5)
       call WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,2)
    case (limiter_wenoz5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,2)
    case (limiter_wenozp5)
       call WENO5limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,3)
    case (limiter_wenozp5nm)
       call WENO5NMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,dxdim,w,wLp,wRp,3)
    case (limiter_weno5cu6)
       call WENO5CU6limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp)
    case (limiter_weno7)
       call WENO7limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp,1)
    case (limiter_mpweno7)
       call WENO7limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp,2)
    case (limiter_exeno7)
       call exENO7limiter(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
          ixLmax1,ixLmax2,idims,w,wLp,wRp)
    case (limiter_ppm)
       call PPMlimiter(ixImin1,ixImin2,ixImax1,ixImax2,ixMlo1,ixMlo2,ixMhi1,&
          ixMhi2,idims,w,w,wLp,wRp)
    case default
       jxRmin1=ixRmin1+kr(idims,1);jxRmin2=ixRmin2+kr(idims,2)
       jxRmax1=ixRmax1+kr(idims,1);jxRmax2=ixRmax2+kr(idims,2);
       ixCmax1=jxRmax1;ixCmax2=jxRmax2; ixCmin1=ixLmin1-kr(idims,1)
       ixCmin2=ixLmin2-kr(idims,2);
       jxCmin1=ixCmin1+kr(idims,1);jxCmin2=ixCmin2+kr(idims,2)
       jxCmax1=ixCmax1+kr(idims,1);jxCmax2=ixCmax2+kr(idims,2);

       do iw=1,nwflux
          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=dlog10(w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw))
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=dlog10(wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw))
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=dlog10(wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw))
          end if

          dwC(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=w(jxCmin1:jxCmax1,&
             jxCmin2:jxCmax2,iw)-w(ixCmin1:ixCmax1,ixCmin2:ixCmax2,iw)
          if(need_global_a2max) then 
            a2max=a2max_global(idims)
          else
            select case(idims)
            case(1)
              a2max=schmid_rad1
            
            case(2)
              a2max=schmid_rad2
            
            case default
              call mpistop("idims is wrong in mod_limiter")
            end select
          end if
            
          ! limit flux from left and/or right
          call dwlimiter2(dwC,ixImin1,ixImin2,ixImax1,ixImax2,ixCmin1,ixCmin2,&
             ixCmax1,ixCmax2,idims,typelimiter,ldw,rdw,a2max=a2max)
          wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)=wLp(ixLmin1:ixLmax1,&
             ixLmin2:ixLmax2,iw)+half*ldw(ixLmin1:ixLmax1,ixLmin2:ixLmax2)
          wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)=wRp(ixRmin1:ixRmax1,&
             ixRmin2:ixRmax2,iw)-half*rdw(jxRmin1:jxRmax1,jxRmin2:jxRmax2)

          if (loglimit(iw)) then
             w(ixCmin1:jxCmax1,ixCmin2:jxCmax2,iw)=10.0d0**w(ixCmin1:jxCmax1,&
                ixCmin2:jxCmax2,iw)
             wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,&
                iw)=10.0d0**wLp(ixLmin1:ixLmax1,ixLmin2:ixLmax2,iw)
             wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,&
                iw)=10.0d0**wRp(ixRmin1:ixRmax1,ixRmin2:ixRmax2,iw)
          end if
       end do
    end select

    wLC(ixLmin1:ixLmax1,ixLmin2:ixLmax2,1:nw)=wLp(ixLmin1:ixLmax1,&
       ixLmin2:ixLmax2,1:nw)
    wRC(ixRmin1:ixRmax1,ixRmin2:ixRmax2,1:nw)=wRp(ixRmin1:ixRmax1,&
       ixRmin2:ixRmax2,1:nw)
    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixLmin1,ixLmin2,&
       ixLmax1,ixLmax2,wLC,x)
    call phys_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixRmin1,ixRmin2,&
       ixRmax1,ixRmax2,wRC,x)

  end subroutine reconstruct_LR

end module mod_finite_volume
