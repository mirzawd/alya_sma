!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmmat_st(                                &
     pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
     gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
     gpadv,gpvep,gpprp,gpgrp,elauu,elaup,elapp,elapu,gpst1,&
     gpstrm,gpstrc,h,elrbu,elrbp,elvel,ielem,gptem,gpgde,dgpmut_dvel,gpgve)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmmat_st
  ! NAME 
  !    nsi_elmmat_st
  ! DESCRIPTION
  !    Compute element matrix and RHS respect to the derivative of stabilization term
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,ntens
  use def_nastin, only       :  nbdfp_nsi
  use def_nastin, only       :  kfl_regim_nsi
  use def_nastin, only       :  ncomp_nsi
  use def_master, only       :  kfl_coupl,ID_NASTIN,ID_TURBUL
  use def_kermod, only       :  kfl_adj_prob
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ielem
  real(rp),    intent(in)    :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)    :: gppor(pgaus),gpgvi(ndime,pgaus)
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
  !   real(rp),    intent(in)    :: gprhs(ndime,pgaus)
  !   real(rp),    intent(inout) :: gprhs_sgs(ndime,pgaus)
  !   real(rp),    intent(in)    :: gprhc(pgaus)
  !   real(rp),    intent(in)    :: gprh2(pgaus)
  !   real(rp),    intent(in)    :: rmomu(pnode,pgaus)
  !   real(rp),    intent(in)    :: rcont(ndime,pnode,pgaus)
  !   real(rp),    intent(out)   :: p1vec(pnode,pgaus)
  !   real(rp),    intent(out)   :: p2vec(ndime,pnode,pgaus)
  !   real(rp),    intent(out)   :: p2sca(pnode,pgaus)
  !   real(rp),    intent(out)   :: wgrgr(pnode,pnode,pgaus)
  !   real(rp),    intent(out)   :: wgrvi(pnode,pgaus)
  real(rp),    intent(inout)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(inout)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(inout)   :: elapp(pnode,pnode)
  real(rp),    intent(inout)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(inout)   :: elrbu(ndime,pnode)
  real(rp),    intent(inout)   :: elrbp(pnode)
  real(rp),    intent(in)    :: gpst1(pgaus)
  real(rp),    intent(in)    :: gpstrm(ndime,pgaus)
  real(rp),    intent(in)    :: gpstrc(pgaus)
  real(rp),    intent(in)    :: h(2)
  real(rp),    intent(in)    :: elvel(ndime,pnode,ncomp_nsi)
  real(rp),    intent(in)    :: gptem(pgaus, nbdfp_nsi)
  real(rp),    intent(in)    :: gpgde(ndime,pgaus)
  real(rp),    intent(in)    :: dgpmut_dvel(ndime,pnode,pgaus)
  real(rp),    intent(in)    :: gpgve(ndime,ndime,pgaus)
  integer(ip)                :: inode,jnode,idofn,jdofn
  integer(ip)                :: jdof2,jdof3,idof1,idof3,idof2
  integer(ip)                :: igaus,idime,jdof1
  integer(ip)                :: ind
  real(rp)                   :: fact2,fact3,fact4,fact5,fact6
  real(rp)                   :: fact0
  real(rp)                   :: fact1
  real(rp)                   :: gpsp1der(pgaus,pnode,ndime),gpsp2der(pgaus,pnode,ndime),gpvno
  real(rp)                   :: eljacuu(pnode*ndime,pnode*ndime)
  real(rp)                   :: eljacpu(pnode,pnode*ndime)
  real(rp)                   :: elvelvec(pnode*ndime)
  real(rp)                   :: jacuuprodu(ndime*pnode),jacpuprodu(pnode)

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------
  do igaus = 1,pgaus
     do inode = 1,pnode
        do idime = 1,ndime
           gpsp1der(igaus,inode,idime) = 0.0_rp
           gpsp2der(igaus,inode,idime) = 0.0_rp
        enddo
     enddo
  enddo

  do idofn = 1,pevat
     do jdofn = 1,pevat
        eljacuu(jdofn,idofn) = 0.0_rp
     end do
     do jnode = 1,pnode
        eljacpu(jnode,idofn) = 0.0_rp
     end do
  end do

  do inode = 1,pnode
     do idime = 1,ndime
        ind = (inode-1)*ndime + idime
        if (kfl_adj_prob == 0) then
           elvelvec(ind) = elvel(idime,inode,1)
        endif
     enddo
  enddo

  do inode = 1,pnode
     jacpuprodu(inode) = 0.0_rp
     do idime = 1,ndime
        ind = (inode-1)*ndime + idime
        jacuuprodu(ind) = 0.0_rp
     enddo
  enddo

  !----------------------------------------------------------------------
  !
  ! Stabilization derivatives
  !
  !----------------------------------------------------------------------


  do igaus = 1,pgaus
     call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
     if (gpvno /= 0.0_rp) then
        do inode = 1,pnode
           do idime = 1,ndime
              gpsp1der(igaus,inode,idime) = -2.0_rp*gpden(igaus)*gpsp1(igaus)*gpsp1(igaus)*gpsha(inode,igaus)*gpadv(idime,igaus)/(gpvno*h(1))
              gpsp2der(igaus,inode,idime) =  2.0_rp*gpden(igaus)*h(2)*h(2)*gpsha(inode,igaus)*gpadv(idime,igaus)/(gpvno*h(1))
           enddo
        enddo
     endif
  enddo

  if (kfl_coupl(ID_NASTIN,ID_TURBUL) >= 1 ) then

     do igaus = 1,pgaus
        call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
        if (gpvno /= 0.0_rp) then
           do inode = 1,pnode
              do idime = 1,ndime
                 gpsp1der(igaus,inode,idime) = gpsp1der(igaus,inode,idime) &
                      -4.0_rp*dgpmut_dvel(idime,inode,igaus)*gpsp1(igaus)*gpsp1(igaus)/(h(2)*h(2))
                 gpsp2der(igaus,inode,idime) = gpsp2der(igaus,inode,idime) + 4.0_rp*dgpmut_dvel(idime,inode,igaus)
              enddo
           enddo
        endif
     enddo


  endif

  !----------------------------------------------------------------------
  !
  ! Jacuu: tau1*dv/dx*Res_m
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then

     do igaus = 1,pgaus
        fact0 = gpvol(igaus) * gpsp1(igaus)
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           fact1 = fact0 * gpcar(1,inode,igaus) 
           fact2 = fact0 * gpcar(2,inode,igaus)
           do jnode = 1,pnode              
              jdof1              =   2*jnode-1
              jdof2              =   jdof1+1
              fact3 = fact1 * gpsha(jnode,igaus)
              fact4 = fact2 * gpsha(jnode,igaus) 
              eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact3 * gpstrm(1,igaus)          ! Auu_xx
              eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact3 * gpstrm(2,igaus)          ! Auu_yx
              eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact4 * gpstrm(1,igaus)          ! Auu_xy
              eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact4 * gpstrm(2,igaus)          ! Auu_yy
           end do
        end do
     end do

  else

     do igaus = 1,pgaus
        fact0 = gpvol(igaus) * gpsp1(igaus)
        do inode = 1,pnode
           idof1 = 3*inode-2
           idof2 = idof1+1
           idof3 = idof2+1
           fact1 = fact0 * gpcar(1,inode,igaus) 
           fact2 = fact0 * gpcar(2,inode,igaus)
           fact3 = fact0 * gpcar(3,inode,igaus)
           do jnode = 1,pnode              
              jdof1              =   3*jnode-2
              jdof2              =   jdof1+1
              jdof3              =   jdof2+1
              fact4 = fact1 * gpsha(jnode,igaus)
              fact5 = fact2 * gpsha(jnode,igaus)
              fact6 = fact3 * gpsha(jnode,igaus) 

              eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact4 * gpstrm(1,igaus)          ! Auu_xx
              eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact4 * gpstrm(2,igaus)          ! Auu_yx
              eljacuu(idof3,jdof1) = eljacuu(idof3,jdof1) + fact4 * gpstrm(3,igaus)          ! Auu_zx              
              eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact5 * gpstrm(1,igaus)          ! Auu_xy
              eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact5 * gpstrm(2,igaus)          ! Auu_yy
              eljacuu(idof3,jdof2) = eljacuu(idof3,jdof2) + fact5 * gpstrm(3,igaus)          ! Auu_zy
              eljacuu(idof1,jdof3) = eljacuu(idof1,jdof3) + fact6 * gpstrm(1,igaus)          ! Auu_xz
              eljacuu(idof2,jdof3) = eljacuu(idof2,jdof3) + fact6 * gpstrm(2,igaus)          ! Auu_yz
              eljacuu(idof3,jdof3) = eljacuu(idof3,jdof3) + fact6 * gpstrm(3,igaus)          ! Auu_zz   

           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Jacuu: d(tau1)/du*adv*dv/dx*Res_m
  !
  !----------------------------------------------------------------------
  if( ndime == 2 ) then

     do igaus = 1,pgaus
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           fact1 = fact0 *( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) )

           do jnode = 1,pnode              
              jdof1              =   2*jnode-1
              jdof2              =   jdof1+1
              fact2 = fact1 * gpsp1der(igaus,jnode,1)
              fact3 = fact1 * gpsp1der(igaus,jnode,2)

              eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact2 * gpstrm(1,igaus)          ! Auu_xx
              eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact2 * gpstrm(2,igaus)          ! Auu_yx
              eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact3 * gpstrm(1,igaus)          ! Auu_xy
              eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact3 * gpstrm(2,igaus)          ! Auu_yy

           end do
        end do
     end do

  else

     do igaus = 1,pgaus
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           idof1 = 3*inode-2
           idof2 = idof1+1
           idof3 = idof2+1
           fact1 = fact0 *( gpadv(1,igaus) * gpcar(1,inode,igaus) + gpadv(2,igaus) * gpcar(2,inode,igaus) & 
                + gpadv(3,igaus) * gpcar(3,inode,igaus) )

           do jnode = 1,pnode              
              jdof1              =   3*jnode-2
              jdof2              =   jdof1+1
              jdof3              =   jdof2+1
              fact2 = fact1 * gpsp1der(igaus,jnode,1)
              fact3 = fact1 * gpsp1der(igaus,jnode,2)
              fact4 = fact1 * gpsp1der(igaus,jnode,3)

              eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact2 * gpstrm(1,igaus)          ! Auu_xx
              eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact2 * gpstrm(2,igaus)          ! Auu_yx
              eljacuu(idof3,jdof1) = eljacuu(idof3,jdof1) + fact2 * gpstrm(3,igaus)          ! Auu_zx
              eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact3 * gpstrm(1,igaus)          ! Auu_xy
              eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact3 * gpstrm(2,igaus)          ! Auu_yy
              eljacuu(idof3,jdof2) = eljacuu(idof3,jdof2) + fact3 * gpstrm(3,igaus)          ! Auu_zy
              eljacuu(idof1,jdof3) = eljacuu(idof1,jdof3) + fact4 * gpstrm(1,igaus)          ! Auu_xz
              eljacuu(idof2,jdof3) = eljacuu(idof2,jdof3) + fact4 * gpstrm(2,igaus)          ! Auu_yz
              eljacuu(idof3,jdof3) = eljacuu(idof3,jdof3) + fact4 * gpstrm(3,igaus)          ! Auu_zz

           end do
        end do
     end do

  end if
  ! 
  !   !----------------------------------------------------------------------
  !   !
  !   ! Jacuu: d(tau2)/du*dv/dxi*Res_c
  !   !
  !   !----------------------------------------------------------------------
  ! 
  if( ndime == 2 ) then

     do igaus = 1,pgaus

        fact0 = gpvol(igaus)*gpstrc(igaus)

        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           fact1 = fact0 * gpcar(1,inode,igaus)
           fact2 = fact0 * gpcar(2,inode,igaus)

           do jnode = 1,pnode              
              jdof1              =   2*jnode-1
              jdof2              =   jdof1+1

              eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact1 * gpsp2der(igaus,jnode,1)          ! Auu_xx
              eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact2 * gpsp2der(igaus,jnode,1)          ! Auu_yx
              eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact1 * gpsp2der(igaus,jnode,2)          ! Auu_xy
              eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact2 * gpsp2der(igaus,jnode,2)          ! Auu_yy

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
              jdof1              =   3*jnode-2
              jdof2              =   jdof1+1
              jdof3              =   jdof2+1

              eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact1 * gpsp2der(igaus,jnode,1)          ! Auu_xx
              eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact2 * gpsp2der(igaus,jnode,1)          ! Auu_yx
              eljacuu(idof3,jdof1) = eljacuu(idof3,jdof1) + fact3 * gpsp2der(igaus,jnode,1)          ! Auu_zx
              eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact1 * gpsp2der(igaus,jnode,2)          ! Auu_xy
              eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact2 * gpsp2der(igaus,jnode,2)          ! Auu_yy
              eljacuu(idof3,jdof2) = eljacuu(idof3,jdof2) + fact3 * gpsp2der(igaus,jnode,2)          ! Auu_zy
              eljacuu(idof1,jdof3) = eljacuu(idof1,jdof3) + fact1 * gpsp2der(igaus,jnode,3)          ! Auu_xz
              eljacuu(idof2,jdof3) = eljacuu(idof2,jdof3) + fact2 * gpsp2der(igaus,jnode,3)          ! Auu_yz
              eljacuu(idof3,jdof3) = eljacuu(idof3,jdof3) + fact3 * gpsp2der(igaus,jnode,3)          ! Auu_zz

           end do

        end do
     end do

  end if
  !  
  !   !----------------------------------------------------------------------
  !   !
  !   ! Jacpu: d(tau1)/du*dq/dxi*Res_mi
  !   !
  !   !----------------------------------------------------------------------
  !  
  if( ndime == 2 ) then
     do igaus = 1,pgaus 
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1
           do jnode = 1,pnode
              fact1 = (gpstrm(1,igaus)*gpcar(1,jnode,igaus)+gpstrm(2,igaus)*gpcar(2,jnode,igaus))*fact0

              eljacpu(jnode,idof1) = eljacpu(jnode,idof1) + fact1 * gpsp1der(igaus,inode,1)   ! Apu_x
              eljacpu(jnode,idof2) = eljacpu(jnode,idof2) + fact1 * gpsp1der(igaus,inode,2)   ! Apu_y

           end do
        end do
     end do
  else
     do igaus = 1,pgaus   
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           idof1 = 3*inode-2
           idof2 = idof1+1
           idof3 = idof2+1
           do jnode = 1,pnode
              fact1 = (gpstrm(1,igaus)*gpcar(1,jnode,igaus)+gpstrm(2,igaus)*gpcar(2,jnode,igaus) &
                   + gpstrm(3,igaus)*gpcar(3,jnode,igaus) )*fact0

              eljacpu(jnode,idof1) = eljacpu(jnode,idof1) + fact1 * gpsp1der(igaus,inode,1)   ! Apu_x
              eljacpu(jnode,idof2) = eljacpu(jnode,idof2) + fact1 * gpsp1der(igaus,inode,2)   ! Apu_y
              eljacpu(jnode,idof3) = eljacpu(jnode,idof3) + fact1 * gpsp1der(igaus,inode,3)   ! Apu_z

           end do
        end do
     end do
  endif

  !----------------------------------------------------------------------
  !
  ! modification of elauu and elapu due to turbulence viscousity
  !
  !----------------------------------------------------------------------  

  if (kfl_coupl(ID_NASTIN,ID_TURBUL) >= 1 ) then

     !----------------------------------------------------------------------
     !
     ! eljacuu: jacobian matrix related to viscosity derivatives (LAPLACIAN) => [dRm/du] = [ d(mut)/du dui/dxj , dvi/dxj ]
     !
     !----------------------------------------------------------------------

     if( ndime == 2 ) then

        do igaus = 1,pgaus
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1

              do jnode = 1,pnode              
                 jdof1 = 2*jnode-1
                 jdof2 = jdof1+1

                 fact0 = gpvol(igaus) * dgpmut_dvel(1,jnode,igaus)
                 fact1 = gpvol(igaus) * dgpmut_dvel(2,jnode,igaus)

                 eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact0*(gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus))
                 eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact0*(gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus))
                 eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact1*(gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus))
                 eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact1*(gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus))

              end do

           end do
        end do

     else

        do igaus = 1,pgaus

           do inode = 1,pnode
              idof1 = 3*inode-2
              idof2 = idof1+1
              idof3 = idof2+1

              do jnode = 1,pnode              
                 jdof1 = 3*jnode-2
                 jdof2 = jdof1+1
                 jdof3 = jdof2+1

                 fact0 = gpvol(igaus) * dgpmut_dvel(1,jnode,igaus)
                 fact1 = gpvol(igaus) * dgpmut_dvel(2,jnode,igaus)
                 fact2 = gpvol(igaus) * dgpmut_dvel(3,jnode,igaus)

                 eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact0*&
                      (gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus) + gpgve(3,1,igaus)*gpcar(3,inode,igaus))
                 eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact0*&
                      (gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus) + gpgve(3,2,igaus)*gpcar(3,inode,igaus))
                 eljacuu(idof3,jdof1) = eljacuu(idof3,jdof1) + fact0*&
                      (gpgve(1,3,igaus)*gpcar(1,inode,igaus) + gpgve(2,3,igaus)*gpcar(2,inode,igaus) + gpgve(3,3,igaus)*gpcar(3,inode,igaus))

                 eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact1*&
                      (gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus) + gpgve(3,1,igaus)*gpcar(3,inode,igaus))
                 eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact1*&
                      (gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus) + gpgve(3,2,igaus)*gpcar(3,inode,igaus))
                 eljacuu(idof3,jdof2) = eljacuu(idof3,jdof2) + fact1*&
                      (gpgve(1,3,igaus)*gpcar(1,inode,igaus) + gpgve(2,3,igaus)*gpcar(2,inode,igaus) + gpgve(3,3,igaus)*gpcar(3,inode,igaus))

                 eljacuu(idof1,jdof3) = eljacuu(idof1,jdof3) + fact2*&
                      (gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus) + gpgve(3,1,igaus)*gpcar(3,inode,igaus))
                 eljacuu(idof2,jdof3) = eljacuu(idof2,jdof3) + fact2*&
                      (gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus) + gpgve(3,2,igaus)*gpcar(3,inode,igaus))
                 eljacuu(idof3,jdof3) = eljacuu(idof3,jdof3) + fact2*&
                      (gpgve(1,3,igaus)*gpcar(1,inode,igaus) + gpgve(2,3,igaus)*gpcar(2,inode,igaus) + gpgve(3,3,igaus)*gpcar(3,inode,igaus))

              end do

           end do
        end do


     end if


  endif

  !----------------------------------------------------------------------
  !
  ! modification of elauu and elapu due to LOW MACH
  !
  !----------------------------------------------------------------------  

  if (kfl_regim_nsi == 3 ) then
     !----------------------------------------------------------------------
     !
     ! Jacuu: tau2/rho * d( u*d(rho)/dx )/du
     !
     !----------------------------------------------------------------------  
     if( ndime == 2 ) then

        do igaus = 1,pgaus
           fact0 = gpsp2(igaus) * gpvol(igaus)       
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1   
              do jnode = 1,pnode
                 jdof1              =   2*jnode-1
                 jdof2              =   jdof1+1
                 fact1  = fact0 * gpsha(jnode,igaus) * gpgde(1,igaus) / gpden(igaus)
                 fact2  = fact0 * gpsha(jnode,igaus) * gpgde(2,igaus) / gpden(igaus)

                 eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact1 * gpcar(1,inode,igaus)          ! Auu_xx
                 eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact1 * gpcar(2,inode,igaus)          ! Auu_yx
                 eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact2 * gpcar(1,inode,igaus)          ! Auu_xy
                 eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact2 * gpcar(2,inode,igaus)          ! Auu_yy      

              enddo !jnode
           enddo !inode
        end do !igaus

     else

        do igaus = 1,pgaus
           fact0 = gpsp2(igaus) * gpvol(igaus)
           do inode = 1,pnode
              idof1 = 3*inode-2
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode
                 jdof1              =   3*jnode-2
                 jdof2              =   jdof1+1
                 jdof3              =   jdof2+1
                 fact1  = fact0 * gpsha(jnode,igaus) * gpgde(1,igaus) / gpden(igaus)
                 fact2  = fact0 * gpsha(jnode,igaus) * gpgde(2,igaus) / gpden(igaus)
                 fact3  = fact0 * gpsha(jnode,igaus) * gpgde(3,igaus) / gpden(igaus)

                 eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact1 * gpcar(1,inode,igaus)          ! Auu_xx
                 eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact1 * gpcar(2,inode,igaus)          ! Auu_yx
                 eljacuu(idof3,jdof1) = eljacuu(idof3,jdof1) + fact1 * gpcar(3,inode,igaus)          ! Auu_zx
                 eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact2 * gpcar(1,inode,igaus)          ! Auu_xy
                 eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact2 * gpcar(2,inode,igaus)          ! Auu_yy      
                 eljacuu(idof3,jdof2) = eljacuu(idof3,jdof2) + fact2 * gpcar(3,inode,igaus)          ! Auu_zy      
                 eljacuu(idof1,jdof3) = eljacuu(idof1,jdof3) + fact3 * gpcar(1,inode,igaus)          ! Auu_xz
                 eljacuu(idof2,jdof3) = eljacuu(idof2,jdof3) + fact3 * gpcar(2,inode,igaus)          ! Auu_yz      
                 eljacuu(idof3,jdof3) = eljacuu(idof3,jdof3) + fact3 * gpcar(3,inode,igaus)          ! Auu_zz      

              enddo !jnode
           enddo !inode      
        end do !igaus

     end if
     !----------------------------------------------------------------------
     !
     ! Jacuu: d(tau2)/du * ui*d(rho)/dxi /rho
     !
     !----------------------------------------------------------------------  
     if( ndime == 2 ) then

        do igaus = 1,pgaus
           fact0 = gpvol(igaus)*( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) )/gpden(igaus)
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1    
              do jnode = 1,pnode
                 jdof1              =   2*jnode-1
                 jdof2              =   jdof1+1

                 eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact0 * gpsp2der(igaus,jnode,1) * gpcar(1,inode,igaus)          ! Auu_xx
                 eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact0 * gpsp2der(igaus,jnode,1) * gpcar(2,inode,igaus)          ! Auu_yx
                 eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact0 * gpsp2der(igaus,jnode,2) * gpcar(1,inode,igaus)          ! Auu_xy
                 eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact0 * gpsp2der(igaus,jnode,2) * gpcar(2,inode,igaus)          ! Auu_yy      

              enddo !jnode
           enddo !inode
        end do !igaus

     else

        do igaus = 1,pgaus
           fact0 = gpvol(igaus)*( gpadv(1,igaus)*gpgde(1,igaus) + gpadv(2,igaus)*gpgde(2,igaus) + gpadv(3,igaus)*gpgde(3,igaus))/gpden(igaus)
           do inode = 1,pnode
              idof1 = 3*inode-2
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode
                 jdof1              =   3*jnode-2
                 jdof2              =   jdof1+1
                 jdof3              =   jdof2+1

                 eljacuu(idof1,jdof1) = eljacuu(idof1,jdof1) + fact0 * gpsp2der(igaus,jnode,1) * gpcar(1,inode,igaus)          ! Auu_xx
                 eljacuu(idof2,jdof1) = eljacuu(idof2,jdof1) + fact0 * gpsp2der(igaus,jnode,1) * gpcar(2,inode,igaus)          ! Auu_yx
                 eljacuu(idof3,jdof1) = eljacuu(idof3,jdof1) + fact0 * gpsp2der(igaus,jnode,1) * gpcar(3,inode,igaus)          ! Auu_zx
                 eljacuu(idof1,jdof2) = eljacuu(idof1,jdof2) + fact0 * gpsp2der(igaus,jnode,2) * gpcar(1,inode,igaus)          ! Auu_xy
                 eljacuu(idof2,jdof2) = eljacuu(idof2,jdof2) + fact0 * gpsp2der(igaus,jnode,2) * gpcar(2,inode,igaus)          ! Auu_yy
                 eljacuu(idof3,jdof2) = eljacuu(idof3,jdof2) + fact0 * gpsp2der(igaus,jnode,2) * gpcar(3,inode,igaus)          ! Auu_zy      
                 eljacuu(idof1,jdof3) = eljacuu(idof1,jdof3) + fact0 * gpsp2der(igaus,jnode,3) * gpcar(1,inode,igaus)          ! Auu_xz
                 eljacuu(idof2,jdof3) = eljacuu(idof2,jdof3) + fact0 * gpsp2der(igaus,jnode,3) * gpcar(2,inode,igaus)          ! Auu_yz
                 eljacuu(idof3,jdof3) = eljacuu(idof3,jdof3) + fact0 * gpsp2der(igaus,jnode,3) * gpcar(3,inode,igaus)          ! Auu_zz      

              enddo !jnode
           enddo !inode
        end do !igaus

     end if



     !----------------------------------------------------------------------
     !
     ! Jacpu: q * d[ u*d(rho)/dx ]/du /rho
     !
     !----------------------------------------------------------------------

     if( ndime == 2 ) then
        do igaus = 1,pgaus 
           fact0 = gpvol(igaus)
           do inode = 1,pnode
              idof1 = 2*inode-1
              idof2 = idof1+1
              do jnode = 1,pnode
                 fact1  = fact0 * gpsha(jnode,igaus) * gpgde(1,igaus) / gpden(igaus)
                 fact2  = fact0 * gpsha(jnode,igaus) * gpgde(2,igaus) / gpden(igaus)

                 eljacpu(jnode,idof1) = eljacpu(jnode,idof1) + fact1 * gpsha(inode,igaus)   ! Apu_x
                 eljacpu(jnode,idof2) = eljacpu(jnode,idof2) + fact2 * gpsha(inode,igaus)   ! Apu_y

              end do
           end do
        end do
     else
        do igaus = 1,pgaus   
           fact0 = gpvol(igaus)
           do inode = 1,pnode
              idof1 = 3*inode-2
              idof2 = idof1+1
              idof3 = idof2+1
              do jnode = 1,pnode
                 fact1  = fact0 * gpsha(jnode,igaus) * gpgde(1,igaus) / gpden(igaus)
                 fact2  = fact0 * gpsha(jnode,igaus) * gpgde(2,igaus) / gpden(igaus)
                 fact3  = fact0 * gpsha(jnode,igaus) * gpgde(3,igaus) / gpden(igaus)

                 eljacpu(jnode,idof1) = eljacpu(jnode,idof1) + fact1 * gpsha(inode,igaus)   ! Apu_x
                 eljacpu(jnode,idof2) = eljacpu(jnode,idof2) + fact2 * gpsha(inode,igaus)   ! Apu_y
                 eljacpu(jnode,idof3) = eljacpu(jnode,idof3) + fact3 * gpsha(inode,igaus)   ! Apu_z

              end do
           end do
        end do
     endif

  endif ! end low mach


  !----------------------------------------------------------------------
  !
  ! Update Auu and Apu for forward and adjoint regarding the exact linearization:
  !
  ! Auu: elauu = elauu + eljacuu
  ! Apu: elapu = elapu + eljacpu
  !
  !----------------------------------------------------------------------

  do idofn = 1,pevat
     do jdofn = 1,pevat
        elauu(jdofn,idofn) = elauu(jdofn,idofn) + eljacuu(jdofn,idofn)
     end do
     do jnode = 1,pnode
        elapu(jnode,idofn) = elapu(jnode,idofn) + eljacpu(jnode,idofn)
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Update bu and bp for adjoint only regarding the time terms
  !
  ! bu = 0.0
  ! bp = 0.0
  !
  !----------------------------------------------------------------------

  !   if(kfl_adj_prob == 1) then  
  !     ! put elrbu and elrbp to zero and then add times terms 
  !     elrbu = 0.0_rp
  !     elrbp = 0.0_rp
  !   endif

  !----------------------------------------------------------------------
  !
  ! Update bu and bp for forward only regarding the exact linearization:
  !
  ! bu: elrbu = elrbu + eljacuu * elvel
  ! bp: elrbp = elrbp + eljacpu * elvel
  !
  !----------------------------------------------------------------------

  !   if(kfl_adj_prob == 0) then  
  !     ! eljacuu * elvel = jacuuprodu
  !     call mbvab0(jacuuprodu,eljacuu,elvelvec,pevat,pevat)
  !     ! eljacpu * elvel = jacpuprodu
  !     call mbvab0(jacpuprodu,eljacpu,elvelvec,pnode,pevat)
  !     do inode = 1,pnode
  !       elrbp(inode) = elrbp(inode) + jacpuprodu(inode)
  !       do idime = 1,ndime
  !    ind = (inode-1)*ndime + idime
  !    elrbu(idime,inode) = elrbu(idime,inode) + jacuuprodu(ind)
  !       enddo
  !     end do    
  !   endif



end subroutine nsi_elmmat_st
