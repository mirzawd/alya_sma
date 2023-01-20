!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_crapro(&
     ielem,pgaus,pnode,gpcod,elcod,gpgdi,gppio,gpdet)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_crapro
  ! NAME
  !    sld_crapro_VM
  ! DESCRIPTION
  !    Projection to nodes
  ! INPUT
  !    GPGDI ....... Deformation tensor ...................... F = grad(phi)
  !    GPPIO ....... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPDET ....... Deformation jacobian determinant ........ J = |F|
  ! OUTPUT
  !    CRAPX_SLD ... Position of crack plane in element
  !    CRANX_SLD ... Normal of crack plane in element
  !    LEENR_SLD ... List of enriched/cracked elements
  ! USES
  !    sld_xfemic(1)
  ! USED BY
  !    sld_elmope(6)
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_elmtyp, only       :  ELEXT,ELCUT
  use def_domain, only       :  ndime,lelch,nmate,lnods,kfl_elcoh
  use def_solidz, only       :  kfl_xfeme_sld,leenr_sld,parch_sld,celen_sld
  use def_solidz, only       :  lecoh_sld,cranx_sld,crapx_sld,lmate_sld,kfl_cohes_sld
  use def_solidz, only       :  nopio_sld,sgmax_sld,sgult_sld
  implicit none

  integer(ip), intent(in)    :: ielem,pgaus,pnode
  real(rp),    intent(in)    :: gpcod(ndime,pgaus),elcod(ndime,pnode)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpdet(pgaus)
  real(rp),    intent(in)    :: gppio(ndime,ndime,pgaus)
  real(rp)                   :: eval(ndime),evec(ndime,ndime)
  integer(ip)                :: idime,jdime,ipoin,inode,i,k
  integer(ip)                :: nrank,pmate,eltip
  real(rp)                   :: sigma(ndime,ndime),sigma3(3,3)
  real(rp)                   :: eval3(3),evec3(3,3),dummr,frori(3,3),finor(3)
  real(rp)                   :: sgmax,sgdir(ndime),tcrit,eleng
  real(rp)                   :: weigh,dista,fibor,fiber,towei,elref
  real(rp)                   :: tcomp,sgsur,sgin1,sgdj2  !compressive strenght, VMsurface, I1, J2
  
  logical                    :: debugging
  

  fiber = 400.0_rp ! * (2.0_rp * rand() - 1.0_rp) ! fiber orientation in degree
  elref = 4.0_rp * atan(1.0_rp)
  fibor = fiber * atan(1.0_rp)/45.0_rp
  finor(1) = -sin(fibor)
  finor(2) =  cos(fibor)
  finor(3) =  0.0_rp
  
  debugging =  .false.
  !
  ! Only if XFEM enrichment is used
  !
  if( kfl_xfeme_sld == 1 .and. leenr_sld(ielem) == 0) then
  
     pmate = 1
     if( nmate > 1 ) then
       pmate = lmate_sld(ielem)
     end if
     tcrit = parch_sld(1,pmate)
     eleng = celen_sld(ielem)
     !
     ! Obtain the crack position: CRAPX_SLD(ielem)
     !
     eltip = 0
     call sld_xfemic(1_ip,pnode,pgaus,dummr,gpcod,elcod,ielem,eltip)
     !
     ! Calculate sigma (weighted averaged Cauchy stress tensor)
     !
     sgmax = 0.0_rp
     sigma = 0.0_rp
     towei = 0.0_rp
     do inode = 1,pnode
        dista = 0.0_rp
        do jdime = 1,ndime
           dista = dista + (elcod(jdime,inode) - crapx_sld(jdime,ielem))* &
                           (elcod(jdime,inode) - crapx_sld(jdime,ielem))
        end do
        dista = sqrt(dista)
        weigh = 1.0_rp / real(pnode,rp) ! exp(-elref*dista/eleng) ! 
        towei = towei + weigh
        ipoin = lnods(inode,ielem)           
        do idime = 1,ndime
           do jdime = 1,ndime
              sigma(idime,jdime) = sigma(idime,jdime) +  weigh * nopio_sld(jdime+(idime-1)*ndime,ipoin)
           end do
        end do
     end do
     
!     do idime = 1,ndime
!        do jdime = 1,ndime
!           sigma(idime,jdime) = sigma(idime,jdime)/towei
!        end do
!     end do
     
     if (fiber > 360.0_rp) then 
        !
        ! Compute the principal tensile direction of averaged sigma
        !
        if( ndime <= 2_ip ) then
           sigma3 = 0.0_rp
           do i = 1,3
              sigma3(i,i) = 1.0_rp
           end do
           do idime = 1,ndime
              do jdime = 1,ndime
                 sigma3(idime,jdime) = sigma(idime,jdime)
              end do
           end do
           call spcdec(sigma3,eval3,evec3,nrank,1_ip,'SLD_CRAPRO, NDIME=2, COMPUTING AVG SIGMA PPAL DIRECTION')
           do k = ndime+1,3
              idime = 0
              do i = 1,3
                 if(abs(evec3(k,i)-1.0_rp) > 1.0e-12_rp) then
                    idime = idime + 1
                    eval(idime) = eval3(i)
                    do jdime = 1,ndime
                       evec(jdime,idime) = evec3(jdime,i)
                    end do
                 endif
              end do
           end do
        else if( ndime == 3_ip ) then
           call spcdec(sigma,eval,evec,nrank,1_ip,'SLD_CRAPRO, NDIME=3, COMPUTING AVG SIGMA PPAL DIRECTION')
        endif
     
        sgmax = 0.0_rp
        do idime = 1,ndime
           if (eval(idime)>=sgmax) then
              sgmax = eval(idime)
              do jdime = 1,ndime
                 sgdir(jdime) = evec(idime,jdime)
              end do
           end if
        end do
        
     else
        !
        ! Compute the stress projection on the prescribed fiber direction
        !
        do idime = 1,ndime
           sgdir(idime) = finor(idime)
        end do
        frori = 0.0_rp
        sgmax = 0.0_rp
        do idime = 1,ndime
           do jdime = 1,ndime
              frori(idime,jdime) = sgdir(idime)*sgdir(jdime)
              sgmax = sgmax + sigma(idime,jdime)*frori(idime,jdime)
           end do
        end do
        
     end if

     cranx_sld(:,ielem) = sgdir
     sgmax_sld(ielem)   = sgmax

#ifdef NO_COLORING
     !$OMP ATOMIC
#endif
     sgult_sld          = max(sgult_sld,sgmax) 
      
     ! Von Misses surface
     tcomp = 10.0_rp*tcrit
     sgin1 = 0.0_rp
     sgdj2 = 0.0_rp
     do idime = 1,ndime
        sgin1 = sgin1 + sigma(idime,idime)
        sgdj2 = sgdj2 + (sigma(idime,idime))*(sigma(idime,idime))
     end do
     sgdj2 = 3*sgdj2 - sgin1*sgin1
 
     !sgsur = 0.5_rp*sgdj2 + (tcomp-tcrit)*sgin1 - tcomp*tcrit
     if (ndime==2) then
        sgsur = sqrt(0.5_rp*((sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))))
     else
        sgsur = sqrt(0.5_rp*((sigma(1,1)-sigma(2,2))*(sigma(1,1)-sigma(2,2))+(sigma(2,2)-sigma(3,3))*(sigma(2,2)-sigma(3,3))+(sigma(3,3)-sigma(1,1))*(sigma(3,3)-sigma(1,1))))
     end if     

     if (eltip == 1 .and. sgmax_sld(ielem) > tcrit .and. sgsur > 0) then
        !
        !  Introduce crack on element when stress exceeds critical level (threshold)
        !
        leenr_sld(ielem) = -1
        if( kfl_cohes_sld /= 0 .and. kfl_elcoh > 0 ) lecoh_sld(ielem) = 1
        lelch(ielem) = ELCUT

        if (debugging) write(*,'(a30)')             '--| ALYA     PROPAGATION CRACK'
        if (debugging) write(*,'(a20,1x,i8,1x,a20)')'--| ALYA     ELEMENT',ielem,'EXCEEDS CRIT. STRESS'
        if (debugging) write(*,'(a27,1x,e12.5)')    '--| ALYA     ELEMENT STRESS',sgmax
        if (debugging) write(*,'(a30,3(1x,f12.7))') '--| ALYA     CRACK ORIENTATION',(cranx_sld(idime,ielem),idime=1,ndime)
        if (debugging) write(*,'(a36,1x,f12.5)')    '--| ALYA     CRACK ORIENTATION ANGLE',&
                       atan(-cranx_sld(1,ielem)/cranx_sld(2,ielem))*45.0_rp/atan(1.0_rp) 
        if (debugging) write(*,'(a30,3(1x,f12.7))') '--| ALYA     CRACK POSITION   ',(crapx_sld(idime,ielem),idime=1,ndime)
         
     end if
        
  end if

end subroutine sld_crapro
