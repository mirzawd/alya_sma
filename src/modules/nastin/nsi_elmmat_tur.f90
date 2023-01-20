!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_elmmat_tur(                                &
     pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpdvi,gpsp1, &
     gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
     gpadv,gpvep,gpprp,gpgrp,gpst1,gpstrm,gpstrc,h,elvel,elpre,gpgve,ielem,p1vec,gpgdv,&
     gpgde,gpdde,gpgdd,elrbu,dgpmut_dtur)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmmat_tur
  ! NAME 
  !    nsi_elmmat_tur
  ! DESCRIPTION
  !    Compute element matrix and RHS related to the derivative of miu w.r.t turbulence
  ! USES
  ! USED BY
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode,ntens,lnods
  use def_nastin, only       :  ncomp_nsi
!  use def_nastin, only       :  kfl_fixno_nsi
  use def_master, only       :  press,RhsadjTur_nsi,RhsadjNas_tur,press_forw
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
  real(rp),    intent(in)    :: gpst1(pgaus)
  real(rp),    intent(in)    :: gpstrm(ndime,pgaus)
  real(rp),    intent(in)    :: gpstrc(pgaus)
  real(rp),    intent(in)    :: h(2)
  real(rp),    intent(in)    :: elvel(ndime,pnode,ncomp_nsi),elpre(pnode)
  real(rp),    intent(in)    :: gpgde(ndime,pgaus)                    ! grad(den)
  real(rp),    intent(in)    :: gpdde(pnode,pgaus)                    ! Density derivatives w.r.t nodal temperature
  real(rp),    intent(in)    :: gpgdd(ndime,pnode,pgaus)              ! Density derivatives w.r.t nodal temperature and coordinates
  real(rp),    intent(in)    :: dgpmut_dtur(1,pnode,pgaus)            ! Turbulence viscosity derivatives w.r.t turbulence unk
  real(rp),    intent(in)    :: gpgve(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: elrbu(ndime,pnode)
  integer(ip)                :: inode,jnode,idofn
  integer(ip)                :: jdof2,jdof3,idof1,idof3,idof2
  integer(ip)                :: igaus,idime,jdof1
  integer(ip)                :: iturb
  real(rp)                   :: fact0,fact1,fact2,fact3
  real(rp)                   :: gpvno
  integer(ip)                :: ipoin,nturb
  real(rp)                   :: dgpsp1_dtur(1,pgaus, pnode)    ! derivatives of tau1 w.r.t nodal k and w
  real(rp)                   :: dgpsp2_dtur(1,pgaus, pnode)    ! derivatives of tau2 w.r.t nodal k and w
  real(rp),    intent(in)    :: p1vec(pnode,pgaus)
  real(rp)                   :: gpgpr(ndime,pgaus), elpre_forw(pnode)              !!! pressure gradients and elemental pressure
  real(rp)                   :: elrhstur(1,pnode)
  real(rp)                   :: elaut(1,pnode*ndime,pnode)
  real(rp)                   :: elapt(1,pnode,pnode)

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  ! has to be defined in kernel
  nturb = 1_ip    

  do idofn = 1,pevat
     do jnode = 1,pnode
        do iturb = 1, nturb 
           elaut(iturb,idofn,jnode) = 0.0_rp
        enddo
     end do
  end do
  do inode = 1,pnode
     do jnode = 1,pnode
        do iturb = 1, nturb
           elapt(iturb,inode,jnode) = 0.0_rp
        enddo
     end do
  end do
  do inode = 1,pnode
     do iturb = 1, nturb
        elrhstur(iturb,inode) = 0.0_rp
     enddo
  end do
  do igaus = 1,pgaus
     do inode = 1,pnode
        do iturb = 1, nturb
           dgpsp1_dtur(iturb,igaus,inode) = 0.0_rp
           dgpsp2_dtur(iturb,igaus,inode) = 0.0_rp
        enddo
     enddo
  enddo

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
  ! Stabilization derivatives w.r.t turbulence
  !
  !----------------------------------------------------------------------

  ! derivatives related to d(miu)/dk and d(miu)/dw
  do igaus = 1,pgaus
     call vecnor(gpadv(1,igaus),ndime,gpvno,2_ip)
     if (gpvno /= 0.0_rp ) then
        do inode = 1,pnode
           do iturb = 1, nturb
              dgpsp1_dtur(iturb,igaus,inode) = -4.0_rp*gpsp1(igaus)**2*dgpmut_dtur(iturb,inode,igaus)/(h(2)**2)
              dgpsp2_dtur(iturb,igaus,inode) =  4.0_rp*dgpmut_dtur(iturb,inode,igaus)       
           enddo
        enddo
     endif
  enddo

  !----------------------------------------------------------------------
  !
  ! Jacobian matrix related to viscosity derivatives (LAPLACIAN) => elauk: [dRm/dk] = [ d(mut)/dk dui/dxj , dvi/dxj ]
  !                                                                 elauw: [dRm/dw] = [ d(mut)/dw dui/dxj , dvi/dxj ]
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then

     do igaus = 1,pgaus

        do inode = 1,pnode
           idof1 = 2*inode-1
           idof2 = idof1+1

           do jnode = 1,pnode              
              do iturb = 1, nturb
                 fact0 = gpvol(igaus) * dgpmut_dtur(iturb,jnode,igaus)
                 elaut(iturb,idof1,jnode) = elaut(iturb,idof1,jnode) + fact0*(gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus))
                 elaut(iturb,idof2,jnode) = elaut(iturb,idof2,jnode) + fact0*(gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus))
              enddo
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
              do iturb = 1, nturb
                 fact0 = gpvol(igaus) * dgpmut_dtur(iturb,jnode,igaus)
                 elaut(iturb,idof1,jnode) = elaut(iturb,idof1,jnode) + fact0*&
                      (gpgve(1,1,igaus)*gpcar(1,inode,igaus) + gpgve(2,1,igaus)*gpcar(2,inode,igaus) + gpgve(3,1,igaus)*gpcar(3,inode,igaus))
                 elaut(iturb,idof2,jnode) = elaut(iturb,idof2,jnode) + fact0*&
                      (gpgve(1,2,igaus)*gpcar(1,inode,igaus) + gpgve(2,2,igaus)*gpcar(2,inode,igaus) + gpgve(3,2,igaus)*gpcar(3,inode,igaus))
                 elaut(iturb,idof3,jnode) = elaut(iturb,idof3,jnode) + fact0*&
                      (gpgve(1,3,igaus)*gpcar(1,inode,igaus) + gpgve(2,3,igaus)*gpcar(2,inode,igaus) + gpgve(3,3,igaus)*gpcar(3,inode,igaus))
              enddo
           end do

        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! elauk: jacobian matrix related to rho*d(tau1)/dk*adv*dvi/dx*( Res_m )
  ! elauw: jacobian matrix related to rho*d(tau1)/dw*adv*dvi/dx*( Res_m )
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
              do iturb = 1, nturb
                 fact2 = fact1 * dgpsp1_dtur(iturb,igaus,jnode)
                 elaut(iturb,idof1,jnode) = elaut(iturb,idof1,jnode) + fact2 * ( gpstrm(1,igaus))          ! Auk_x
                 elaut(iturb,idof2,jnode) = elaut(iturb,idof2,jnode) + fact2 * ( gpstrm(2,igaus))          ! Auk_x
              enddo
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
              do iturb = 1, nturb
                 fact2 = fact1 * dgpsp1_dtur(iturb,igaus,jnode)
                 elaut(iturb,idof1,jnode) = elaut(iturb,idof1,jnode) + fact2 * ( gpstrm(1,igaus))          ! Auk_x
                 elaut(iturb,idof2,jnode) = elaut(iturb,idof2,jnode) + fact2 * ( gpstrm(2,igaus))          ! Auk_y
                 elaut(iturb,idof3,jnode) = elaut(iturb,idof3,jnode) + fact2 * ( gpstrm(3,igaus))          ! Auk_z
              enddo
           end do
        end do
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! elauk: jacobian matrix related to d(tau2)/dk*dvi/dx*Res_c
  ! elauw: jacobian matrix related to d(tau2)/dw*dvi/dx*Res_c
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
              do iturb = 1, nturb
                 elaut(iturb,idof1,jnode) = elaut(iturb,idof1,jnode) + fact1 * dgpsp2_dtur(iturb,igaus,jnode)          ! Auk_x
                 elaut(iturb,idof2,jnode) = elaut(iturb,idof2,jnode) + fact2 * dgpsp2_dtur(iturb,igaus,jnode)          ! Auk_y
              enddo
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
              do iturb = 1, nturb
                 elaut(iturb,idof1,jnode) = elaut(iturb,idof1,jnode) + fact1 * dgpsp2_dtur(iturb,igaus,jnode)          ! Auk_x
                 elaut(iturb,idof2,jnode) = elaut(iturb,idof2,jnode) + fact2 * dgpsp2_dtur(iturb,igaus,jnode)          ! Auk_y
                 elaut(iturb,idof3,jnode) = elaut(iturb,idof3,jnode) + fact3 * dgpsp2_dtur(iturb,igaus,jnode)          ! Auk_z
              enddo
           end do
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! elapk: jacobian matrix related to d(tau1)/dk*dqi/dxk*Res_mk
  ! elapw: jacobian matrix related to d(tau1)/dw*dqi/dxk*Res_mk
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then
     do igaus = 1,pgaus 
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           do jnode = 1,pnode
              fact1 = (gpstrm(1,igaus)*gpcar(1,jnode,igaus)+gpstrm(2,igaus)*gpcar(2,jnode,igaus))*fact0
              do iturb = 1, nturb
                 elapt(iturb,jnode,inode) = elapt(iturb,jnode,inode) + fact1 * dgpsp1_dtur(iturb,igaus,inode)   ! Apt
              enddo
           end do
        end do
     end do
  else
     do igaus = 1,pgaus   
        fact0 = gpvol(igaus)
        do inode = 1,pnode
           do jnode = 1,pnode
              fact1 = (gpstrm(1,igaus)*gpcar(1,jnode,igaus)+gpstrm(2,igaus)*gpcar(2,jnode,igaus)+gpstrm(3,igaus)*gpcar(3,jnode,igaus) )*fact0
              do iturb = 1, nturb
                 elapt(iturb,jnode,inode) = elapt(iturb,jnode,inode) + fact1 * dgpsp1_dtur(iturb,igaus,inode)   ! Apt
              enddo
           end do
        end do
     end do
  endif


  !----------------------------------------------------------------------
  !
  ! Apply B.C. on elaut
  !
  !----------------------------------------------------------------------

  !   do inode = 1,pnode
  !      idof1 = (inode-1) * ndime
  !      ipoin = lnods(inode,ielem)
  !      do idime = 1,ndime
  !         idof1 = idof1+1
  !         if(      kfl_fixno_nsi(idime,ipoin) ==  1 &
  !              .or.kfl_fixno_nsi(idime,ipoin) ==  8 &
  !              .or.kfl_fixno_nsi(idime,ipoin) ==  9 &
  !              .or.kfl_fixno_nsi(idime,ipoin) ==  5 &
  !              .or.kfl_fixno_nsi(idime,ipoin) ==  6 &
  !              .or.kfl_fixno_nsi(idime,ipoin) ==  7 &
  !              .or.kfl_fixno_nsi(idime,ipoin) == 11 ) elaut(iturb,idof1,:) = 0.0_rp
  !      enddo
  !   enddo

  !----------------------------------------------------------------------
  !
  ! rhsid vector related to viscosity derivative matrix   Transpose[dRm/dT]*[Lambda_u] + Transpose[dRc/dT]*[Lambda_p]
  !
  !----------------------------------------------------------------------

  if( ndime == 2 ) then
     do inode = 1, pnode
        do jnode = 1, pnode
           jdof1 = 2*jnode-1
           jdof2 = jdof1+1
           do iturb = 1, nturb
              elrhstur(iturb,inode) = elrhstur(iturb,inode) + elaut(iturb,jdof1,inode)*elvel(1,jnode,2) + elaut(iturb,jdof2,inode)*elvel(2,jnode,2)
           enddo
        enddo
     enddo
  else
     do inode = 1, pnode
        do jnode = 1, pnode
           jdof1 = 3*jnode-2
           jdof2 = jdof1+1
           jdof3 = jdof2+1
           do iturb = 1, nturb
              elrhstur(iturb,inode) = elrhstur(iturb,inode) + elaut(iturb,jdof1,inode)*elvel(1,jnode,2) + elaut(iturb,jdof2,inode)*elvel(2,jnode,2) + &
                   elaut(iturb,jdof3,inode)*elvel(3,jnode,2)
           enddo
        enddo
     enddo
  endif

  do inode = 1, pnode
     do jnode = 1, pnode
        do iturb = 1, nturb
           elrhstur(iturb,inode) = elrhstur(iturb,inode) + elapt(iturb,jnode,inode)*elpre(jnode)
        enddo
     enddo
  enddo

  !----------------------------------------------------------------------
  !
  !                 RhsadjTur_nsi(1)   <------   elrhskin
  !                 RhsadjTur_nsi(2)   <------   elrhsome
  !
  !----------------------------------------------------------------------  
  do inode = 1,pnode
     do iturb = 1, nturb
        RhsadjTur_nsi(ielem)%a(iturb,inode) = elrhstur(iturb,inode)
     enddo
  enddo

  !----------------------------------------------------------------------
  !
  !                 !  ELRHS = ELRHS - [dRt/du] [Lambda_t]  sent by turbul
  !
  !----------------------------------------------------------------------
  do idime = 1,ndime
     do inode = 1,pnode
        elrbu(idime,inode) = elrbu(idime,inode) - RhsadjNas_tur(ielem)%a(idime, inode) 
     enddo
  enddo


end subroutine nsi_elmmat_tur
