!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_elmpre(&
     ielem,pnode,pgaus,pmate,eldis,eldip,elvel,elacc,eldix,elvex,elcod,&
     elmof,elrst,gpsha,gpcar,gphea,gpgi0,gpdet,gpvol,&
     elepo,elepo_new_inverse,elepo_new_determinant,&
     gpvel,gpacc,&
     gprat,gpdis,gpcau,gpigd,gpcod,gpgdi,&
     gpigd_eps,gpgdi_eps,gpdet_eps,&
     gpmof,gprestre,gpdet_min,&
     ielem_min,trpn_sld_temp)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_elmpre
  ! NAME
  !    sld_elmpre
  ! DESCRIPTION
  !    Compute some Gauss values
  ! INPUT
  !    GPGI0 ... Previous deformation gradient tensor ............ F(n)
  ! OUTPUT
  !    GPDIS ... Deformation ..................................... phi (phidot is elvel)
  !    GPVEL ... Deformation velocity............................. phidot
  !    GPGDI ... Updated deformation gradient .....................F(n+1) = grad(phi)
  !    GPIGD ... Inverse of updated deformation gradient tensor .. F^{-1}(n+1)
  !    GPRAT ... Rate of Deformation   ........................... Fdot = grad(phidot)
  !    GPCAU ... Right Cauchy-Green deformation tensor ........... C = F^t x F
  !    GPCOD ... Physical coordinates of the gauss points
  !    GPRESTRE ... Residual stresses
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp,lg
  use def_elmtyp, only       :  QUA04, HEX08
  use def_domain, only       :  mnode,ndime,mgaus,lnods,ltype
  use def_solidz, only       :  volum_sld
  use def_solidz, only       :  kfl_serei_sld
  use def_solidz, only       :  kfl_xfeme_sld, kfl_ninex_sld
  use def_solidz, only       :  kfl_gdepo
  use def_solidz, only       :  kfl_moduf_sld,nvgij_inv_sld,kfl_restr_sld
  use def_master, only       :  ITER_K, TIME_N

  implicit none

  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: pgaus
  integer(ip), intent(in)    :: pmate
  real(rp),    intent(in)    :: eldis(ndime,pnode,*)
  real(rp),    intent(in)    :: eldip(ndime,pnode  )
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: elacc(ndime,pnode,3)
  real(rp),    intent(in)    :: elrst(    6,pnode)
  real(rp),    intent(in)    :: eldix(ndime,pnode,*)
  real(rp),    intent(in)    :: elvex(ndime,pnode)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elmof(pnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gphea(pnode,pgaus)
  real(rp),    intent(in)    :: gpgi0(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpdet(pgaus)
  real(rp),    intent(in)    :: gpdet_eps(pgaus)
  real(rp),    intent(inout) :: gpvol(pgaus)
  real(rp),    intent(out)   :: elepo(ndime,ndime,pnode)
  real(rp),    intent(inout) :: trpn_sld_temp(pnode)
  real(rp),    intent(out)   :: elepo_new_inverse(ndime,ndime,pnode)
  real(rp),    intent(out)   :: elepo_new_determinant(pnode)
  real(rp),    intent(out)   :: gpvel(ndime,pgaus,*)
  real(rp),    intent(out)   :: gpacc(ndime,pgaus,3)
  real(rp),    intent(out)   :: gprat(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpdis(ndime,pgaus,*)       !< Deformation: phi
  real(rp),    intent(out)   :: gpcau(ndime,ndime,pgaus)   !< Right Cauchy-Green deformation tensor: C = F^t x F
  real(rp),    intent(out)   :: gpcod(ndime,mgaus)
  real(rp),    intent(out)   :: gpigd(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpigd_eps(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpgdi_eps(ndime,ndime,pgaus)
  real(rp),    intent(out)   :: gpmof(pgaus)
  real(rp),    intent(out)   :: gprestre(ndime,ndime,pgaus)
  real(rp),    intent(inout) :: gpdet_min
  integer(ip), intent(inout) :: ielem_min
  real(rp)                   :: detme,volel,xmean
  integer(ip)                :: inode,igaus,idime,jdime
  integer(ip)                :: kdime,ipoin,kauxi,ivoig
  real(rp)                   :: detf0,hglas,gpgdi_aux(ndime,ndime,pgaus)
  !
  ! GPDIS and GPVEL: displacement and velocity
  !
  kauxi= 0
  gpacc = 0.0_rp
  elepo = 0.0_rp
  elepo_new_determinant = 0.0_rp
  elepo_new_inverse = 0.0_rp

  do idime = 1,ndime
     do igaus = 1,pgaus
        gpcod(idime,igaus)   = 0.0_rp
        gpdis(idime,igaus,ITER_K) = 0.0_rp
        gpdis(idime,igaus,TIME_N) = 0.0_rp
        gpvel(idime,igaus,ITER_K) = 0.0_rp
        do inode = 1,pnode
           gpdis(idime,igaus,ITER_K) = gpdis(idime,igaus,ITER_K) + eldis(idime,inode,ITER_K)*gpsha(inode,igaus)
           gpdis(idime,igaus,TIME_N) = gpdis(idime,igaus,TIME_N) + eldis(idime,inode,TIME_N)*gpsha(inode,igaus)
           gpvel(idime,igaus,ITER_K) = gpvel(idime,igaus,ITER_K) + elvel(idime,inode)*gpsha(inode,igaus)
           gpacc(idime,igaus,ITER_K) = gpacc(idime,igaus,ITER_K) + elacc(idime,inode,ITER_K)*gpsha(inode,igaus)
           gpacc(idime,igaus,TIME_N) = gpacc(idime,igaus,TIME_N) + elacc(idime,inode,TIME_N)*gpsha(inode,igaus)
           gpcod(idime,igaus)   = gpcod(idime,igaus)   + elcod(idime,inode)*gpsha(inode,igaus)
        end do
        do jdime = 1,ndime
           gprestre(idime,jdime,igaus) = 0.0_rp
           gpgdi(idime,jdime,igaus) = 0.0_rp
           gpgdi_eps(idime,jdime,igaus) = 0.0_rp
           gprat(idime,jdime,igaus) = 0.0_rp
        end do
     end do
  end do
  !
  ! GPDIS and GPVEL: displacement and velocity
  ! to include the contribution of the enrichment dofs
  !
  if (kfl_xfeme_sld == 1) then
     do idime=1,ndime
        do igaus=1,pgaus
           do inode=1,pnode
              gpdis(idime,igaus,ITER_K)=gpdis(idime,igaus,ITER_K)&
                   +eldix(idime,inode,ITER_K)*gpsha(inode,igaus)*gphea(inode,igaus)
              gpdis(idime,igaus,TIME_N)=gpdis(idime,igaus,TIME_N)&
                   +eldix(idime,inode,TIME_N)*gpsha(inode,igaus)*gphea(inode,igaus)
              gpvel(idime,igaus,ITER_K)=gpvel(idime,igaus,ITER_K)&
                   +elvex(idime,inode)*gpsha(inode,igaus)*gphea(inode,igaus)
           end do
        end do
     end do
  end if

  !
  ! GPGDI: Deformation gradients
  !
  do igaus=1,pgaus
     do idime=1,ndime
        do jdime=1,ndime
           do inode=1,pnode
              gpgdi(idime,jdime,igaus)=gpgdi(idime,jdime,igaus)&
                   +eldis(idime,inode,ITER_K)*gpcar(jdime,inode,igaus)
              gprat(idime,jdime,igaus)=gprat(idime,jdime,igaus)&
                   +elvel(idime,inode)*gpcar(jdime,inode,igaus)
              gpgdi_eps(idime,jdime,igaus)=gpgdi_eps(idime,jdime,igaus)&
                   +(eldis(idime,inode,ITER_K)+eldip(idime,inode))*gpcar(jdime,inode,igaus)
           end do
        end do
     end do
  end do

  !
  ! GPRESTRE: Residual stresses
  !
  if (kfl_restr_sld < 0) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              do inode=1,pnode
                 ivoig= nvgij_inv_sld(idime,jdime)
                 gprestre(idime,jdime,igaus)=gprestre(idime,jdime,igaus)&
                      +elrst(ivoig,inode)*gpsha(inode,igaus)
              end do
           end do
        end do
     end do
  end if

  !
  ! GPGDI: Deformation gradients for xfem case
  ! to include the contribution of the enrichment dofs  !
  ! MARIANO: REVISAR ESTO PARA EL NINEX
  if (kfl_xfeme_sld == 1) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              do inode=1,pnode
                 gpgdi(idime,jdime,igaus)=gpgdi(idime,jdime,igaus)&
                      +eldix(idime,inode,1)*gpcar(jdime,inode,igaus)*gphea(inode,igaus)
                 gprat(idime,jdime,igaus)=gprat(idime,jdime,igaus)&
                      +elvex(idime,inode)*gpcar(jdime,inode,igaus)*gphea(inode,igaus)
              end do
           end do
        end do
     end do
  end if

  !
  ! GPMOF: Modulating fields, elmof is set to 1 when there are no fields (then, gpmof will be 1)
  !
  if( kfl_moduf_sld(1) > 0 ) then
     do igaus = 1,pgaus
        gpmof(igaus) = 0.0_rp
        do inode = 1,pnode
           gpmof(igaus) = gpmof(igaus) + elmof(inode) * gpsha(inode,igaus)
        end do
     end do
  end if

  detme= 0.0_rp
  volel= 0.0_rp
  gpigd= 0.0_rp
  gpigd_eps= 0.0_rp

  do igaus=1,pgaus

     gpgdi(1,1,igaus)= gpgdi(1,1,igaus) + 1.0_rp
     gpgdi(2,2,igaus)= gpgdi(2,2,igaus) + 1.0_rp
     if (ndime==3) gpgdi(3,3,igaus)= gpgdi(3,3,igaus) + 1.0_rp
     gpgdi_eps(1,1,igaus)= gpgdi_eps(1,1,igaus) + 1.0_rp
     gpgdi_eps(2,2,igaus)= gpgdi_eps(2,2,igaus) + 1.0_rp
     if (ndime==3) gpgdi_eps(3,3,igaus)= gpgdi_eps(3,3,igaus) + 1.0_rp
     !Compute J
     call invmtx(gpgdi(1,1,igaus),gpigd(1,1,igaus),gpdet(igaus),ndime)
     if (kfl_ninex_sld == 1) call invmtx(gpgdi_eps(1,1,igaus),gpigd_eps(1,1,igaus),gpdet_eps(igaus),ndime)

     if( gpdet(igaus) <= 0.0_rp .and. gpdet(igaus) < gpdet_min ) then
        !$OMP CRITICAL (compute_gpdet_min)
        if( gpdet(igaus) < gpdet_min ) then
           !
           ! GPDET < GPDET_MIN must be checked just in case another thread has just modified it!
           !
           gpdet_min = gpdet(igaus)
           ielem_min = ielem
        end if
        !$OMP END CRITICAL (compute_gpdet_min)

     end if

     volel = volel + gpvol(igaus)
     detme = detme + gpvol(igaus) * gpdet(igaus)

     !
     ! GDEPO: Push forward
     !
     if( kfl_gdepo /= 0 ) then
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           xmean = gpsha(inode,igaus) * gpvol(igaus)
           elepo_new_determinant(inode) = elepo_new_determinant(inode) + xmean * gpdet(igaus)
           do idime= 1,ndime
              do jdime= 1,ndime
                 elepo(idime,jdime,inode) = elepo(idime,jdime,inode) + xmean * gpgdi(idime,jdime,igaus)
                 elepo_new_inverse(idime,jdime,inode) = elepo_new_inverse(idime,jdime,inode) + xmean * gpigd(idime,jdime,igaus)
              end do
           end do
        end do
     end if

  end do
  volum_sld(1) = volum_sld(1) + volel         ! reference
  volum_sld(2) = volum_sld(2) + detme         ! deformed
  !
  ! Selective reduced integration
  !
  if( kfl_serei_sld == 1_ip .and. &
       (ltype(ielem) == QUA04 .or. ltype(ielem) == HEX08) ) then
     detf0 = 0.0_rp
     do igaus = 1,pgaus
        detf0 = detf0 + 0.125_rp*gpdet(igaus)
     end do
     ! Selective reduced integration: (see de Souza Neto et al. (1996))
     !  replace F = F * ( mean(J) / J )^{1/ndime}
     !  where mean(J) = J_0 (evaluate at the centroid)
     do igaus = 1,pgaus
        hglas = (detf0/gpdet(igaus))**(1.0_rp/real(ndime,rp))
        gpgdi_aux(:,:,igaus) = hglas*gpgdi(:,:,igaus)
        gpgdi(:,:,igaus)     = gpgdi_aux(:,:,igaus)
        ! Compute inverse and detF
        call invmtx(gpgdi(1,1,igaus),gpigd(1,1,igaus),gpdet(igaus),ndime)
     end do
  end if
  !
  ! GPCAU: Cauchy tensor
  !
  do igaus=1,pgaus
     do kdime=1,ndime
        do jdime=1,ndime
           gpcau(jdime,kdime,igaus)=0.0_rp
           do idime=1,ndime
              gpcau(jdime,kdime,igaus)= gpcau(jdime,kdime,igaus) + gpgdi(idime,jdime,igaus)*gpgdi(idime,kdime,igaus)
           end do
        end do
     end do
  end do
  
end subroutine sld_elmpre

