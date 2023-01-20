!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmmat_tem(                                &
     pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpdvi,gpsp1, &
     gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
     gpadv,gpvep,gpprp,gpgrp,elaut,elapt,gpst1,gpstrm,gpstrc,h,elvel,elpre,gpgve,ielem,p1vec,gpgdv,&
     gpgde,gpdde,gpgdd,elrbu)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmmat_tem
  ! NAME 
  !    nsi_elmmat_tem
  ! DESCRIPTION
  !    Compute element matrix and RHS related to the derivative of miu w.r.t temperature
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,ntens,lnods
  use def_nastin, only       :  fvins_nsi,kfl_regim_nsi
  use def_nastin, only       :  ncomp_nsi
  use def_master, only       :  press,RhsadjTem_nsi,RhsadjNas_tem,press_forw
  use def_kermod, only       :  kfl_adj_prob
  use mod_ker_proper

  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ielem
  real(rp),    intent(in)    :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)    :: gppor(pgaus),gpgvi(ndime,pgaus),gpdvi(pnode,pgaus),gpgdv(ndime,pnode,pgaus)
  real(rp),    intent(inout) :: gpsp1(pgaus),gptt1(pgaus)
  real(rp),    intent(inout) :: gpsp2(pgaus),gptt2(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(inout) :: gpvep(ndime,pgaus)
  real(rp),    intent(in)    :: gpprp(pgaus)      
  real(rp),    intent(in)    :: gpgrp(ndime,pgaus)
  real(rp),    intent(in)    :: gplap(pnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(inout) :: elaut(pnode*ndime,pnode)
  real(rp),    intent(inout) :: elapt(pnode,pnode)
  real(rp),    intent(in)    :: gpst1(pgaus)
  real(rp),    intent(in)    :: gpstrm(ndime,pgaus)
  real(rp),    intent(in)    :: gpstrc(pgaus)
  real(rp),    intent(in)    :: h(2)
  real(rp),    intent(in)    :: elvel(ndime,pnode,ncomp_nsi),elpre(pnode)
  real(rp),    intent(in)    :: gpgde(ndime,pgaus)                    ! grad(den)
  real(rp),    intent(in)    :: gpdde(pnode,pgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp),    intent(in)    :: gpgdd(ndime,pnode,pgaus)              ! Density derivatives w.r.t nodal temperature and coordinates
  real(rp),    intent(in)    :: gpgve(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: elrbu(ndime,pnode)
  integer(ip)                :: inode,jnode,idofn
  integer(ip)                :: jdof2,idof1,idof3,idof2
  integer(ip)                :: igaus,idime,jdof1
!  integer(ip)                :: jdime
  real(rp)                   :: fact0,fact1,fact2,fact3
  real(rp)                   :: gpvno
  integer(ip)                :: ipoin
  real(rp)                   :: gpsp1_tem(pgaus, pnode),gpsp2_tem(pgaus, pnode)    ! derivatives of tau1 and tau2 w.r.t nodal temperature
  real(rp),    intent(in)    :: p1vec(pnode,pgaus)
  real(rp)                   :: gpgpr(ndime,pgaus), elpre_forw(pnode)              !!! pressure gradients and elemental pressure
  real(rp)                   :: elrhstem(pnode)

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  do idofn = 1,pevat
     do jnode = 1,pnode
        elaut(idofn,jnode) = 0.0_rp
     end do
  end do

  do inode = 1,pnode
     do jnode = 1,pnode
        elapt(inode,jnode) = 0.0_rp
     end do
  end do


  do inode = 1,pnode
     elrhstem(inode) = 0.0_rp
  end do

  do igaus = 1,pgaus
     do inode = 1,pnode
        gpsp1_tem(igaus,inode) = 0.0_rp
        gpsp2_tem(igaus,inode) = 0.0_rp
     enddo
  enddo


  !   !----------------------------------------------------------------------
  !   !
  !   ! velocity gradients
  !   !
  !   !----------------------------------------------------------------------
  !   
  !   do igaus = 1,pgaus
  !       do idime = 1,ndime
  !    do jdime = 1,ndime
  !        gpgve(jdime,idime,igaus) = 0.0_rp
  !    end do
  !       end do
  !       do inode = 1,pnode
  !    do idime = 1,ndime
  !        do jdime = 1,ndime
  !          gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
  !           + gpcar(jdime,inode,igaus) * elvel(idime,inode,1)
  !        end do
  !    end do
  !       end do
  !   end do

  !----------------------------------------------------------------------
  !
  ! pressure gradients
  !
  !----------------------------------------------------------------------
  do inode = 1,pnode
     ipoin = lnods(inode,ielem)
     if (kfl_adj_prob == 0) then
        elpre_forw(inode) = press(ipoin,1)
     else
        elpre_forw(inode) = press_forw(ipoin,1)
     endif
  end do

  do igaus = 1,pgaus
     do idime = 1,ndime
        gpgpr(idime,igaus)   = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpgpr(idime,igaus) = gpgpr(idime,igaus) + elpre_forw(inode) * gpcar(idime,inode,igaus)
        end do
     end do
  end do


  !----------------------------------------------------------------------
  !
  ! Stabilization derivatives w.r.t temperature
  !
  !----------------------------------------------------------------------

  ! derivatives related to d(miu)/dT  
  do igaus = 1,pgaus
     call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
     if (gpvno /= 0.0_rp) then
        do inode = 1,pnode
           gpsp1_tem(igaus,inode) = -4.0_rp*gpsp1(igaus)**2*gpdvi(inode,igaus)/(h(2)**2)
           gpsp2_tem(igaus,inode) =  4.0_rp*gpdvi(inode,igaus)
        enddo
     endif
  enddo

  ! derivatives related to d(rho)/dT for low mach
  if (kfl_regim_nsi == 3 ) then

     do igaus = 1,pgaus
        call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
        if (gpvno /= 0.0_rp) then
           do inode = 1,pnode
              gpsp1_tem(igaus,inode) = gpsp1_tem(igaus,inode) - 2.0_rp*gpvno*gpsp1(igaus)**2/h(1) * gpdde(inode,igaus)
              gpsp2_tem(igaus,inode) = gpsp2_tem(igaus,inode) + 2.0_rp*gpvno*h(2)**2/h(1) * gpdde(inode,igaus)
           enddo
        endif
     enddo

  endif

  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to viscosity derivatives (LAPLACIAN) => [dRm/dT] = [ d(mu)/dT dui/dxj , dvi/dxj ]
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then

     do igaus = 1,pgaus

        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1

           do jnode = 1,pnode              
              fact0 = gpvol(igaus) * gpdvi(jnode,igaus)
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact0*(gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus))
              elaut(idof2,jnode) = elaut(idof2,jnode) + fact0*(gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus))
           end do

        end do
     end do

  else

     do igaus = 1,pgaus

        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           idof3 = idof2+1

           do jnode = 1,pnode              
              fact0 = gpvol(igaus) * gpdvi(jnode,igaus)
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact0*&
                   (gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus) + gpgve(3,1,igaus)*gpcar(3,inode,igaus))
              elaut(idof2,jnode) = elaut(idof2,jnode) + fact0*&
                   (gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus) + gpgve(3,2,igaus)*gpcar(3,inode,igaus))
              elaut(idof3,jnode) = elaut(idof3,jnode) + fact0*&
                   (gpgve(1,3,igaus)*gpcar(1,inode,igaus) + gpgve(2,3,igaus)*gpcar(2,inode,igaus) + gpgve(3,3,igaus)*gpcar(3,inode,igaus))
           end do

        end do
     end do


  end if

  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to gradient viscosity derivatives (LAPLACIAN) => - [ d(d(mu)/dxj)/dT dui/dxj , (p1vec-v)]
  !
  !----------------------------------------------------------------------
  if( ndime == 2 ) then
     do igaus = 1,pgaus
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           do jnode = 1,pnode 
              fact1 = ( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,1 ,igaus)* gpgdv(2,jnode,igaus) )
              fact2 = ( gpgve(1,2 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) )
              elaut(idof1,jnode) = elaut(idof1,jnode) - fact1 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))          ! Aut_x
              elaut(idof2,jnode) = elaut(idof2,jnode) - fact2 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))          ! Aut_y         
           end do
        end do
     end do

  else

     do igaus = 1,pgaus
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           idof3 = idof2+1
           do jnode = 1,pnode 
              fact1 = ( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,1 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(3,1 ,igaus)* gpgdv(3,jnode,igaus))
              fact2 = ( gpgve(1,2 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(3,2 ,igaus)* gpgdv(3,jnode,igaus))
              fact3 = ( gpgve(1,3 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,3 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(3,3 ,igaus)* gpgdv(3,jnode,igaus))

              elaut(idof1,jnode) = elaut(idof1,jnode) - fact1 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))          ! Aut_x
              elaut(idof2,jnode) = elaut(idof2,jnode) - fact2 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))          ! Aut_y         
              elaut(idof3,jnode) = elaut(idof3,jnode) - fact3 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))          ! Aut_y 
           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to viscosity derivatives (DIVERGENCE) => [dRm/dT] = [ d(mu)/dT duj/dxi , dvi/dxj ]
  !
  !----------------------------------------------------------------------

  if( fvins_nsi > 0.9_rp ) then

     if( ndime == 2 ) then

        do igaus = 1,pgaus

           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1

              do jnode = 1,pnode              
                 fact0 = gpvol(igaus) * gpdvi(jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact0*(gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(1,2,igaus)*gpcar(2,inode,igaus))
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact0*(gpgve(2,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus))
              end do

           end do
        end do

     else

        do igaus = 1,pgaus

           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1

              do jnode = 1,pnode              
                 fact0 = gpvol(igaus) * gpdvi(jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact0* &
                      (gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(1,2,igaus)*gpcar(2,inode,igaus) + gpgve(1,3,igaus)*gpcar(3,inode,igaus))
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact0* &
                      (gpgve(2,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus) + gpgve(2,3,igaus)*gpcar(3,inode,igaus))
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact0* &
                      (gpgve(3,1,igaus)*gpcar(1,inode,igaus) + gpgve(3,2,igaus)*gpcar(2,inode,igaus) + gpgve(3,3,igaus)*gpcar(3,inode,igaus))
              end do

           end do
        end do

     end if

  end if

  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to gradient viscosity derivatives (DIVERGENCE) => -A [ d(dmu/dxj)dT duj/dxi , (p1vec-v) ]  
  ! 
  !----------------------------------------------------------------------
  if( fvins_nsi > 0.9_rp ) then

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode 
                 fact1 = ( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(1,2 ,igaus)* gpgdv(2,jnode,igaus) )
                 fact2 = ( gpgve(2,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) )
                 elaut(idof1,jnode) = elaut(idof1,jnode) - fact1 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))     ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) - fact2 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))     ! Aut_y         
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode 
                 fact1 = ( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(1,2 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(1,3 ,igaus)* gpgdv(3,jnode,igaus))
                 fact2 = ( gpgve(2,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(2,3 ,igaus)* gpgdv(3,jnode,igaus))
                 fact3 = ( gpgve(3,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(3,2 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(3,3 ,igaus)* gpgdv(3,jnode,igaus))

                 elaut(idof1,jnode) = elaut(idof1,jnode) - fact1 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))    ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) - fact2 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))    ! Aut_y         
                 elaut(idof3,jnode) = elaut(idof3,jnode) - fact3 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))    ! Aut_y 
              end do
           end do
        end do
     endif !ndime

  endif
  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to gradient viscosity derivatives (COMPLETED) => + 2/3 B [ d(dmu/dxi (div u))dT , (p1vec-v) ]
  ! 
  !----------------------------------------------------------------------

  if( fvins_nsi > 1.9_rp ) then

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode 
                 fact1 = 2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus)) * gpgdv(1,jnode,igaus)
                 fact2 = 2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus)) * gpgdv(2,jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))       ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))       ! Aut_y         
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode 
                 fact1 = 2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3 ,igaus) ) * gpgdv(1,jnode,igaus)
                 fact2 = 2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3 ,igaus) ) * gpgdv(2,jnode,igaus)
                 fact3 = 2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3 ,igaus) ) * gpgdv(3,jnode,igaus)

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))     ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))     ! Aut_y         
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact3 * (p1vec(inode,igaus)-gpsha(inode,igaus)*gpvol(igaus))     ! Aut_y 
              end do
           end do
        end do
     endif !ndime

  endif

  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to gradient viscosity derivatives (COMPLETED) => - 2/3 * B * ( d(mu)/dT (div u) , dvi/dxi )
  ! 
  !----------------------------------------------------------------------

  if( fvins_nsi > 1.9_rp ) then

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode 
                 fact1 = -2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus)) * gpdvi(jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * gpcar(1,inode,igaus) *gpvol(igaus)       ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact1 * gpcar(2,inode,igaus) *gpvol(igaus)       ! Aut_y         
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode 
                 fact1 = -2.0_rp/3.0_rp * ( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3,igaus)) * gpdvi(jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * gpcar(1,inode,igaus) *gpvol(igaus)       ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact1 * gpcar(2,inode,igaus) *gpvol(igaus)       ! Aut_y         
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact1 * gpcar(3,inode,igaus) *gpvol(igaus)       ! Aut_y         
              end do
           end do
        end do
     endif !ndime

  endif


  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to rho*d(tau1)/dTj*adv*dvi/dx*( Res_m )
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then

     do igaus = 1,pgaus
        fact0 = gpvol(igaus)*gpden(igaus)
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           fact1 = fact0 *( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) )
           do jnode = 1,pnode              
              fact2 = fact1 * gpsp1_tem(igaus,jnode)
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact2 * ( gpstrm(1,igaus))          ! Aut_x
              elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * ( gpstrm(2,igaus))          ! Aut_x
           end do
        end do
     end do
  else
     do igaus = 1,pgaus
        fact0 = gpvol(igaus)*gpden(igaus)
        do inode = 1,pnode
           idof1 = 3*inode-2
           idof2 = idof1+1
           idof3 = idof2+1
           fact1 = fact0 *( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) & 
                + gpadv(3,igaus) * gpcar(3,inode,igaus) )
           do jnode = 1,pnode              
              fact2 = fact1 * gpsp1_tem(igaus,jnode)
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact2 * ( gpstrm(1,igaus))          ! Aut_x
              elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * ( gpstrm(2,igaus))          ! Aut_y
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact2 * ( gpstrm(3,igaus))          ! Aut_z
           end do
        end do
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to d(tau2)/dTj*dvi/dx*Res_c
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then

     do igaus = 1,pgaus

        fact0 = gpvol(igaus)*gpstrc(igaus)
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           fact1 = fact0 * gpcar(1,inode,igaus)
           fact2 = fact0 * gpcar(2,inode,igaus)
           do jnode = 1,pnode              
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * gpsp2_tem(igaus,jnode)          ! Aut_x
              elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * gpsp2_tem(igaus,jnode)          ! Aut_y
           end do
        end do
     end do

  else

     do igaus = 1,pgaus

        fact0 = gpvol(igaus)*gpstrc(igaus)
        do inode = 1,pnode
           idof1 = 3*inode-2
           idof2 = idof1+1
           idof3 = idof2+1
           fact1 = fact0 * gpcar(1,inode,igaus)
           fact2 = fact0 * gpcar(2,inode,igaus)
           fact3 = fact0 * gpcar(3,inode,igaus)
           do jnode = 1,pnode              
              elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * gpsp2_tem(igaus,jnode)          ! Aut_x
              elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * gpsp2_tem(igaus,jnode)          ! Aut_y
              elaut(idof3,jnode) = elaut(idof3,jnode) + fact3 * gpsp2_tem(igaus,jnode)          ! Aut_z
           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! elapt: jacobian matrix related to d(tau1)/dTj*dqi/dxk*Res_mk
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then
     do igaus = 1,pgaus 
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           do jnode = 1,pnode
              fact1 = (gpstrm(1,igaus)*gpcar(1,jnode,igaus)+gpstrm(2,igaus)*gpcar(2,jnode,igaus))*fact0

              elapt(jnode,inode) = elapt(jnode,inode) + fact1 * gpsp1_tem(igaus,inode)   ! Apt

           end do
        end do
     end do
  else
     do igaus = 1,pgaus   
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           do jnode = 1,pnode
              fact1 = (gpstrm(1,igaus)*gpcar(1,jnode,igaus)+gpstrm(2,igaus)*gpcar(2,jnode,igaus) &
                   + gpstrm(3,igaus)*gpcar(3,jnode,igaus) )*fact0

              elapt(jnode,inode) = elapt(jnode,inode) + fact1 * gpsp1_tem(igaus,inode)   ! Apt              

           end do
        end do
     end do
  endif

  !----------------------------------------------------------------------
  !
  ! elapt: jacobian matrix related to viscosity derivatives (LAPLACIAN) p2vec * d(dk/dxj)/dT dui/dxj
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then
     do igaus = 1,pgaus
        do inode = 1,pnode
           do jnode = 1,pnode 
              fact1 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,1 ,igaus)* gpgdv(2,jnode,igaus) )
              fact2 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,2 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) )
              elapt(inode,jnode) = elapt(inode,jnode) - fact1 * gpcar(1,inode,igaus) - fact2 * gpcar(2,inode,igaus)         ! Apt
           end do
        end do
     end do

  else

     do igaus = 1,pgaus
        do inode = 1,pnode
           do jnode = 1,pnode 
              fact1 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,1 ,igaus)* gpgdv(2,jnode,igaus) )
              fact2 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,2 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) )
              fact3 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,3 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,3 ,igaus)* gpgdv(2,jnode,igaus) )
              elapt(inode,jnode) = elapt(inode,jnode) - fact1 * gpcar(1,inode,igaus) - fact2 * gpcar(2,inode,igaus) &
                   - fact3 * gpcar(3,inode,igaus)                                       ! Apt
           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! elapt: jacobian matrix related to viscosity derivatives (DIVERGENCE) -A [ d(dmu/dxj)/dT duj/dxi , p2vec ]
  !
  !---------------------------------------------------------------------

  if( fvins_nsi > 0.9_rp ) then

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode 
                 fact1 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(1,2 ,igaus)* gpgdv(2,jnode,igaus) )
                 fact2 = gpvol(igaus)*gpsp1(igaus)*( gpgve(2,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) )
                 elapt(inode,jnode) = elapt(inode,jnode) - fact1 * gpcar(1,inode,igaus) - fact2 * gpcar(2,inode,igaus)         ! Apt
              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode 
                 fact1 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(1,2 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(1,3 ,igaus)* gpgdv(3,jnode,igaus))
                 fact2 = gpvol(igaus)*gpsp1(igaus)*( gpgve(2,1 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,2 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(2,3 ,igaus)* gpgdv(3,jnode,igaus))
                 fact3 = gpvol(igaus)*gpsp1(igaus)*( gpgve(1,3 ,igaus)* gpgdv(1,jnode,igaus) + gpgve(2,3 ,igaus)* gpgdv(2,jnode,igaus) + gpgve(3,3 ,igaus)* gpgdv(3,jnode,igaus))
                 elapt(inode,jnode) = elapt(inode,jnode) - fact1 * gpcar(1,inode,igaus) - fact2 * gpcar(2,inode,igaus) &
                      - fact3 * gpcar(3,inode,igaus)                                       ! Apt
              end do
           end do
        end do

     end if !ndime

  endif


  !----------------------------------------------------------------------
  !
  ! elaut: jacobian matrix related to gradient viscosity derivatives (COMPLETED) => +2/3 B [ d(dmu/dxi (div u))dT , p2vec ]
  ! 
  !----------------------------------------------------------------------
  if( fvins_nsi > 1.9_rp ) then

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode 
                 fact1 = 2.0_rp/3.0_rp*gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus))* gpgdv(1,jnode,igaus)
                 fact2 = 2.0_rp/3.0_rp*gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus))* gpgdv(2,jnode,igaus)
                 elapt(inode,jnode) = elapt(inode,jnode) + fact1 * gpcar(1,inode,igaus) + fact2 * gpcar(2,inode,igaus)         ! Apt
              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode 
                 fact1 = 2.0_rp/3.0_rp*gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3 ,igaus))* gpgdv(1,jnode,igaus)
                 fact2 = 2.0_rp/3.0_rp*gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3 ,igaus))* gpgdv(2,jnode,igaus)
                 fact3 = 2.0_rp/3.0_rp*gpvol(igaus)*gpsp1(igaus)*( gpgve(1,1 ,igaus) + gpgve(2,2 ,igaus) + gpgve(3,3 ,igaus))* gpgdv(3,jnode,igaus)
                 elapt(inode,jnode) = elapt(inode,jnode) + fact1 * gpcar(1,inode,igaus) + fact2 * gpcar(2,inode,igaus) &
                      + fact3 * gpcar(3,inode,igaus)                                       ! Apt
              end do
           end do
        end do

     end if !ndime

  endif

  !----------------------------------------------------------------------
  !
  ! modification of elaut and elapt due to LOW MACH (variable density)
  !
  !----------------------------------------------------------------------  

  if (kfl_regim_nsi == 3 ) then

     !----------------------------------------------------------------------
     !
     ! elaut: tau2 dN/dx * d( u*grad(rho) / rho)/dT = tau2 * dN/dx * ()
     !
     !----------------------------------------------------------------------  

     if( ndime == 2 ) then

        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode 
                 fact1 = gpvol(igaus)/gpden(igaus) *      ( gpadv(1,igaus)*gpgdd(1,jnode,igaus) + gpadv(2,igaus)*gpgdd(2,jnode,igaus)) - &
                      gpvol(igaus)/(gpden(igaus)**2) * ( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) ) * gpdde(jnode,igaus)

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * gpsp2(igaus) * gpcar(1,inode,igaus)          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact1 * gpsp2(igaus) * gpcar(2,inode,igaus)          ! Aut_y         
              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode 
                 fact1 = gpvol(igaus)/gpden(igaus) *      ( gpadv(1,igaus)*gpgdd(1,jnode,igaus) + gpadv(2,igaus)*gpgdd(2,jnode,igaus) + &
                      gpadv(3,igaus)*gpgdd(3,jnode,igaus)) - &
                      gpvol(igaus)/(gpden(igaus)**2) * ( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) + &
                      gpadv(3,igaus)*gpgde(3,igaus) ) * gpdde(jnode,igaus)

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * gpsp2(igaus) * gpcar(1,inode,igaus)          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact1 * gpsp2(igaus) * gpcar(2,inode,igaus)          ! Aut_y         
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact1 * gpsp2(igaus) * gpcar(3,inode,igaus)          ! Aut_z         

              end do
           end do
        end do

     end if

     !-------------------------------------------------------------
     !    
     !     elaut: d(tau2)/dT dN/dx * u*grad(rho) / rho = 
     !
     !----------------------------------------------------------------------  

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           fact0 = gpvol(igaus)*( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) )/gpden(igaus)

           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact0 * gpsp2_tem(igaus,jnode) * gpcar(1,inode,igaus)          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact0 * gpsp2_tem(igaus,jnode) * gpcar(2,inode,igaus)          ! Aut_y         

              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           fact0 = gpvol(igaus)*( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) + gpadv(3,igaus)*gpgde(3,igaus))/gpden(igaus)

           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode 

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact0 * gpsp2_tem(igaus,jnode) * gpcar(1,inode,igaus)          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact0 * gpsp2_tem(igaus,jnode) * gpcar(2,inode,igaus)          ! Aut_y         
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact0 * gpsp2_tem(igaus,jnode) * gpcar(3,inode,igaus)          ! Aut_z         

              end do
           end do
        end do

     end if

     !----------------------------------------------------------------------
     !
     ! elapt: jacobian matrix related to q * d( u*grad(rho) / rho )/dT
     !
     !----------------------------------------------------------------------

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode 
                 fact1 = gpvol(igaus)/gpden(igaus) *      ( gpadv(1,igaus)*gpgdd(1,jnode,igaus) + gpadv(2,igaus)*gpgdd(2,jnode,igaus)) - &
                      gpvol(igaus)/(gpden(igaus)**2) * ( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) ) * gpdde(jnode,igaus)

                 elapt(inode,jnode) = elapt(inode,jnode) + fact1 * gpsha(inode,igaus)         ! Apt

              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           do inode = 1,pnode
              do jnode = 1,pnode 
                 fact1 = gpvol(igaus)/gpden(igaus) *      ( gpadv(1,igaus)*gpgdd(1,jnode,igaus) + gpadv(2,igaus)*gpgdd(2,jnode,igaus) + &
                      gpadv(3,igaus)*gpgdd(3,jnode,igaus)) - &
                      gpvol(igaus)/(gpden(igaus)**2) * ( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) + &
                      gpadv(3,igaus)*gpgde(3,igaus) ) * gpdde(jnode,igaus)

                 elapt(inode,jnode) = elapt(inode,jnode) + fact1 * gpsha(inode,igaus)         ! Apt

              end do
           end do
        end do

     end if

     !----------------------------------------------------------------------
     !
     ! elaut: jacobian matrix related to tau1*d(rho)/dTj*adv*dvi/dx*( Res_m )
     !
     !----------------------------------------------------------------------

     if( ndime == 2 ) then

        do igaus = 1,pgaus
           fact0 = gpvol(igaus)*gpsp1(igaus)
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              fact1 = fact0 *( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) )
              do jnode = 1,pnode
                 fact2 = fact1 * gpdde(jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact2 * ( gpstrm(1,igaus))          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * ( gpstrm(2,igaus))          ! Aut_x
              end do
           end do
        end do
     else
        do igaus = 1,pgaus
           fact0 = gpvol(igaus)*gpsp1(igaus)
           do inode = 1,pnode
              idof1 = 3*inode-2
              idof2 = idof1+1
              idof3 = idof2+1
              fact1 = fact0 *( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) & 
                   + gpadv(3,igaus) * gpcar(3,inode,igaus) )
              do jnode = 1,pnode
                 fact2 = fact1 * gpdde(jnode,igaus)
                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact2 * ( gpstrm(1,igaus))          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * ( gpstrm(2,igaus))          ! Aut_y
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact2 * ( gpstrm(3,igaus))          ! Aut_z
              end do
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! elaut: jacobian matrix related to p1vec * d(rho)/dT * uj dui/dxj
     !
     !----------------------------------------------------------------------

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode 
                 fact1 = gpdde(jnode,igaus) * ( gpgve(1,1 ,igaus)* gpadv(1,igaus) + gpgve(2,1 ,igaus)* gpadv(2,igaus) )
                 fact2 = gpdde(jnode,igaus) * ( gpgve(1,2 ,igaus)* gpadv(1,igaus) + gpgve(2,2 ,igaus)* gpadv(2,igaus) )

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * p1vec(inode,igaus)          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * p1vec(inode,igaus)          ! Aut_y         
              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode 
                 fact1 = gpdde(jnode,igaus) * ( gpgve(1,1 ,igaus)* gpadv(1,igaus) + gpgve(2,1 ,igaus)* gpadv(2,igaus) + gpgve(3,1 ,igaus)* gpadv(3,igaus))
                 fact2 = gpdde(jnode,igaus) * ( gpgve(1,2 ,igaus)* gpadv(1,igaus) + gpgve(2,2 ,igaus)* gpadv(2,igaus) + gpgve(3,2 ,igaus)* gpadv(3,igaus))
                 fact3 = gpdde(jnode,igaus) * ( gpgve(1,3 ,igaus)* gpadv(1,igaus) + gpgve(2,3 ,igaus)* gpadv(2,igaus) + gpgve(3,3 ,igaus)* gpadv(3,igaus))

                 elaut(idof1,jnode) = elaut(idof1,jnode) + fact1 * p1vec(inode,igaus)          ! Aut_x
                 elaut(idof2,jnode) = elaut(idof2,jnode) + fact2 * p1vec(inode,igaus)          ! Aut_y         
                 elaut(idof3,jnode) = elaut(idof3,jnode) + fact3 * p1vec(inode,igaus)          ! Aut_y 
              end do
           end do
        end do

     end if

     !----------------------------------------------------------------------
     !
     ! elapt: jacobian matrix related to tau1 dq/dxi * d(rho)/dT * uj dui/dxj
     !
     !----------------------------------------------------------------------

     if( ndime == 2 ) then
        do igaus = 1,pgaus
           fact1 = ( gpgve(1,1 ,igaus)* gpadv(1,igaus) + gpgve(2,1 ,igaus)* gpadv(2,igaus) )*gpvol(igaus)
           fact2 = ( gpgve(1,2 ,igaus)* gpadv(1,igaus) + gpgve(2,2 ,igaus)* gpadv(2,igaus) )*gpvol(igaus)
           do inode = 1,pnode
              do jnode = 1,pnode 
                 elapt(inode,jnode) = elapt(inode,jnode) + gpdde(jnode,igaus)*gpsp1(igaus)* &
                      (gpcar(1,inode,igaus)*fact1 + gpcar(2,inode,igaus)*fact2) ! Aut_x
              end do
           end do
        end do

     else

        do igaus = 1,pgaus
           fact1 = ( gpgve(1,1 ,igaus)* gpadv(1,igaus) + gpgve(2,1 ,igaus)* gpadv(2,igaus) + gpgve(3,1 ,igaus)* gpadv(3,igaus))*gpvol(igaus)
           fact2 = ( gpgve(1,2 ,igaus)* gpadv(1,igaus) + gpgve(2,2 ,igaus)* gpadv(2,igaus) + gpgve(3,2 ,igaus)* gpadv(3,igaus))*gpvol(igaus)
           fact3 = ( gpgve(1,3 ,igaus)* gpadv(1,igaus) + gpgve(2,3 ,igaus)* gpadv(2,igaus) + gpgve(3,3 ,igaus)* gpadv(3,igaus))*gpvol(igaus)
           do inode = 1,pnode
              do jnode = 1,pnode 
                 elapt(inode,jnode) = elapt(inode,jnode) + gpdde(jnode,igaus)*gpsp1(igaus)* &
                      (gpcar(1,inode,igaus)*fact1+gpcar(2,inode,igaus)*fact2+gpcar(3,inode,igaus)*fact3) ! Aut_x
              end do
           end do
        end do

     end if


  endif !low mach

  !----------------------------------------------------------------------
  !
  ! rhsid vector related to viscosity derivative matrix   Transpose[dRm/dT]*[Lambda_u] + Transpose[dRc/dT]*[Lambda_p]
  !
  !----------------------------------------------------------------------

  do inode = 1, pnode
     do jnode = 1, pnode
        jdof1 = 2*jnode-1
        jdof2 = jdof1+1
        elrhstem(inode) = elrhstem(inode) + elaut(jdof1,inode)*elvel(1,jnode,2) + elaut(jdof2,inode)*elvel(2,jnode,2)
     enddo
  enddo

  do inode = 1, pnode
     do jnode = 1, pnode
        elrhstem(inode) = elrhstem(inode) + elapt(jnode,inode)*elpre(jnode)
     enddo
  enddo

  !----------------------------------------------------------------------
  !
  !                 RhsadjTem_nsi   <------   elrhstem
  !
  !----------------------------------------------------------------------

  do inode = 1,pnode
     RhsadjTem_nsi(ielem)%a(inode) = elrhstem(inode)
  enddo

  !----------------------------------------------------------------------
  !
  !                 !  ELRHS = ELRHS - [dRt/du] [Lambda_t]  sent by temper 
  !
  !----------------------------------------------------------------------

  do idime = 1,ndime
     do inode = 1,pnode
        elrbu(idime,inode) = elrbu(idime,inode) - RhsadjNas_tem(ielem)%a(idime, inode) 
     enddo
  enddo


end subroutine nsi_elmmat_tem
