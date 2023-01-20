!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_updunk
  ! NAME
  !    nsi_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates for the velocity and
  !    pressure.
  ! USED BY
  !    nsi_begste (itask=1)
  !    nsi_begite (itask=2)
  !    nsi_endite (itask=3, inner loop)
  !    nsi_endite (itask=4, outer loop)
  !    nsi_endste (itask=5)
  !    nsi_endste (itask=6)               Updates Pressure
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,           only : ip,rp
  use def_master
  use def_master,           only : INOTMASTER,veloc,press,densi,INOTEMPTY
  use def_master,           only : tempe,unkno,vesgs,solve,rbbou
  use def_master,           only : coupling,ID_NASTIN,lzone,npoi1,npoi2,npoi3
  use def_domain,           only : npoin,nelem,ndime
  use def_domain,           only : ngaus,ltype,nrbod
  use def_nastin,           only : NSI_MONOLITHIC
  use def_nastin,           only : relax_nsi,kfl_tiaor_nsi
  use def_nastin,           only : kfl_tiacc_nsi,kfl_tisch_nsi
  use def_nastin,           only : kfl_sgsac_nsi,kfl_resid_nsi
  use def_nastin,           only : veold_nsi
  use def_nastin,           only : ivari_nsi,kfl_dttyp_nsi
  use def_nastin,           only : ivari_nsi_mom,ivari_nsi_cont
  use def_nastin,           only : kfl_relax_nsi
  use def_nastin,           only : nprev_nsi,relap_nsi
  use def_nastin,           only : kfl_sgsti_nsi,kfl_relap_nsi
  use def_nastin,           only : kfl_updpr_nsi,ndbgs_nsi
  use def_nastin,           only : kfl_confi_nsi
  use def_nastin,           only : kfl_ini_ts_guess_order_nsi
  use def_nastin,           only : kfl_bubbl_nsi
  use def_nastin,           only : kfl_fscon_nsi
  use def_nastin,           only : kfl_stabi_nsi, NSI_FRACTIONAL_STEP, NSI_ASGS
  use def_nastin,           only : vefix
  use mod_nsi_bubble,       only : nsi_bubble_update
  use mod_array_operations, only : array_operations_axpby
  use mod_array_operations, only : array_operations_copy
  use mod_array_operations, only : array_operations_axoz
  use mod_communications,   only : PAR_SUM
  use mod_nsi_efect,        only : nsi_set_fringe_node_bcs
  use def_coupli,           only : kfl_efect 
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,             only : arrays_number
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,itotv,idime,itime,ielem
  integer(ip)             :: pelty,pgaus,iimbo,icomp
  real(rp)                :: rela1_nsi,paver,preno
  real(rp),    pointer    :: vpfor(:,:)
  real(rp),    pointer    :: vptor(:,:)
  
  if( INOTEMPTY ) then

     select case (itask)

     case ( ITASK_INIUNK )
        !
        ! (:,3) <= (:,1)
        !
        do icomp = 1,size(veloc,DIM=3,KIND=ip)
           do ipoin = 1,npoin
              veloc(:,ipoin,icomp) = veloc(:,ipoin,1)
           end do
        end do
        do icomp = 1,size(press,DIM=2,KIND=ip)
           do ipoin = 1,npoin
              press(ipoin,icomp) = press(ipoin,1) 
           end do
        end do
        
     case( ITASK_BEGSTE )
        !
        ! Extrapolate (:,1)
        ! (:,2) <= (:,1)
        !
        if( kfl_ini_ts_guess_order_nsi == 2 ) then
           do ipoin = 1,npoin
              veloc(:,ipoin,1) = 2.0_rp * veloc(:,ipoin,1) - veloc(:,ipoin,nprev_nsi+1)
           end do
        end if
        do ipoin = 1,npoin        
           veloc(:,ipoin,2) = veloc(:,ipoin,1)
           press(ipoin,2)   = press(ipoin,1)
        end do

     case( ITASK_BEGITE )
        !
        ! UNKNO <= (:,1)
        !
        veloc(1:ndime,1:npoin,1) = veloc(1:ndime,1:npoin,2)
        press(1:npoin,1)         = press(1:npoin,2)
        !
        ! Dirichlet boundary conditions on velocity for Embedded Finite Element
        ! Coupling Technique (EFECT) 
        !
        if( kfl_efect ) call nsi_set_fringe_node_bcs(2_ip, vefix)

        if( NSI_MONOLITHIC ) then
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = (ipoin-1) * solve(1) % ndofn + idime
                 unkno(itotv) = veloc(idime,ipoin,1)
              end do
              unkno(itotv+1) = press(ipoin,1)
           end do

        else
           !
           ! Schur complement: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           itotv = 0
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = (ipoin-1) * ndime + idime
                 unkno(itotv) = veloc(idime,ipoin,1)
              end do
              unkno(ndbgs_nsi+ipoin) = press(ipoin,1)
           end do
        end if

     case( ITASK_INNITE )
        !
        ! (:,1) <= UNKNO
        !
        do ipoin = 1,npoin
           do idime = 1,ndime
              itotv= (ipoin-1)*ndime + idime
              veloc(idime,ipoin,1) = unkno(itotv)
           end do
           press(ipoin,1) = unkno(ndbgs_nsi+ipoin)
        end do
       
     case( ITASK_ENDINN )
        !
        ! (:,1) <= UNKNO
        !
        if( NSI_MONOLITHIC ) then
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = (ipoin-1) * solve(1) % ndofn + idime
                 veloc(idime,ipoin,1) = unkno(itotv)
              end do
              press(ipoin,1) = unkno(itotv+1)
           end do
        else
           !
           ! Schur complement: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv= (ipoin-1)*ndime + idime
                 veloc(idime,ipoin,1) = unkno(itotv)
              end do
              press(ipoin,1) = unkno(ndbgs_nsi+ipoin)
           end do
        end if
        !
        ! Average pressure to zero
        !
        if( kfl_confi_nsi == -2 ) then
           paver = 0.0_rp
           preno = 0.0_rp
           do ipoin = 1,npoin
              if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then
                 paver = paver + press(ipoin,1)
                 preno = preno + 1.0_rp
              end if
           end do
           call PAR_SUM(paver,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
           call PAR_SUM(preno,'IN MY CODE WITHOUT MASTER',INCLUDE_ROOT=.true.)
           paver = paver / preno
           press(1:npoin,1) = press(1:npoin,1) - paver
        end if
        !
        ! Pressure bubble
        !
        if( kfl_bubbl_nsi /= 0 ) call nsi_bubble_update()

     case( ITASK_ENDITE )
        !
        ! Assign u(n,i-1,*) <-- u(n,i,*)
        !
        ! If velocity residual is required, save previous interation in veold_nsi
        ! 
        if( kfl_resid_nsi == 1 ) then
           veold_nsi(1:ndime,1:npoin) = veloc(1:ndime,1:npoin,2)
        end if
        veloc(1:ndime,1:npoin,2) = veloc(1:ndime,1:npoin,1)
        press(1:npoin,2)         = press(1:npoin,1)

     case( ITASK_ENDSTE )
        !
        ! u(n-1,*,*) <-- u(n,*,*)
        !
        if( kfl_tisch_nsi == 1 .and. kfl_tiacc_nsi == 2 ) then
           !
           ! Crank-Nicolson method
           !
           veloc(1:ndime,1:npoin,1) = 2.0_rp*veloc(1:ndime,1:npoin,1)-veloc(1:ndime,1:npoin,3)
           if( kfl_updpr_nsi== 1 ) then
              press(1:npoin,1) = 2.0_rp*press(1:npoin,1)-press(1:npoin,3)
           end if

        else if( kfl_tisch_nsi == 2 ) then
           !
           ! BDF scheme
           !
           do itime = 2+kfl_tiaor_nsi,4,-1
              veloc(1:ndime,1:npoin,itime) = veloc(1:ndime,1:npoin,itime-1)
           end do

        else if( kfl_tisch_nsi == 4 ) then
           !
           ! Multistep FS scheme
           !
           if(kfl_fscon_nsi == 0) then   ! NON-CONSISTENT - Capuano approximation this is where we need the n-1 pressure
              press(1:npoin,4) = press(1:npoin,3)
           end if
           if (NSI_FRACTIONAL_STEP.and.kfl_stabi_nsi == NSI_ASGS)  &
                veloc(1:ndime,1:npoin,4) = veloc(1:ndime,1:npoin,3)
        end if
        press(1:npoin,3) = press(1:npoin,1)
        !
        ! Adaptive time step
        !
        if( kfl_dttyp_nsi == 2 &  ! ADAPTIVE TIME STEP
             .or. output_postprocess_check_variable_postprocess(arrays_number('DUDT ')) .or. &
             &    output_postprocess_check_variable_postprocess(arrays_number('D2UDT')) ) then
           veloc(1:ndime,1:npoin,5) = veloc(1:ndime,1:npoin,4)
           veloc(1:ndime,1:npoin,4) = veloc(1:ndime,1:npoin,3)
           veloc(1:ndime,1:npoin,3) = veloc(1:ndime,1:npoin,1)
        else
           veloc(1:ndime,1:npoin,3) = veloc(1:ndime,1:npoin,1)
        end if

        if( kfl_sgsti_nsi == 1 ) then
           !
           ! Time tracking of the subscales
           !
           if( kfl_sgsac_nsi == 2 .and. kfl_tiacc_nsi == 2 ) then
              !$OMP PARALLEL DO SCHEDULE (STATIC)          &
              !$OMP DEFAULT  ( NONE )                      &
              !$OMP PRIVATE  ( ielem, pelty, pgaus )       &
              !$OMP SHARED   ( ltype, ngaus, nelem,        &
#ifndef NDIMEPAR
              !$OMP            ndime,                      &
#endif
              !$OMP            vesgs )
              do ielem = 1,nelem
                 pelty = ltype(ielem)
                 pgaus = ngaus(pelty)
                 vesgs(ielem) % a(1:ndime,1:pgaus,1) =             &
                      & 2.0_rp*vesgs(ielem) % a(1:ndime,1:pgaus,1) &
                      &       -vesgs(ielem) % a(1:ndime,1:pgaus,2)
              end do
              !$OMP END PARALLEL DO
           end if
           !$OMP PARALLEL DO SCHEDULE (STATIC)             &
           !$OMP DEFAULT  ( NONE )                         &
           !$OMP PRIVATE  ( ielem, pelty, pgaus )          &
           !$OMP SHARED   ( nelem, ltype, ngaus,           &
#ifndef NDIMEPAR
           !$OMP            ndime,                         &
#endif
           !$OMP            vesgs )
           do ielem = 1,nelem
              pelty = abs(ltype(ielem))
              pgaus = ngaus(pelty)
              vesgs(ielem) % a(1:ndime,1:pgaus,2) = vesgs(ielem) % a(1:ndime,1:pgaus,1)
           end do
           !$OMP END PARALLEL DO
        end if
        if ( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod /= 0 ) ) then   ! rigid body
           do iimbo = 1,nrbod
              vpfor         => rbbou(iimbo) % vpfor
              vptor         => rbbou(iimbo) % vptor
              vpfor(1:3,4)  =  vpfor(1:3,3)
              vpfor(1:3,3)  =  vpfor(1:3,1)
              vptor(1:3,4)  =  vptor(1:3,3)
              vptor(1:3,3)  =  vptor(1:3,1)
           end do
        end if

     case ( 1100_ip )
        !
        ! Equal all variables
        !
        do ipoin = 1,npoin
           veloc(:,ipoin,nprev_nsi) = veloc(:,ipoin,1)
           veloc(:,ipoin,2)         = veloc(:,ipoin,1)
           press(ipoin,nprev_nsi)   = press(ipoin,1)
           press(ipoin,2)           = press(ipoin,1)
        end do

        if( NSI_FRACTIONAL_STEP.and.kfl_stabi_nsi == NSI_ASGS ) then
           do ipoin = 1,npoin
              veloc(:,ipoin,4) = veloc(:,ipoin,1)
           end do
        end if
        
        if( kfl_sgsti_nsi == 1 ) then
           !$OMP PARALLEL DO SCHEDULE (STATIC)              &
           !$OMP DEFAULT  ( NONE )                          &
           !$OMP PRIVATE  ( ielem, pelty, pgaus )           &
           !$OMP SHARED   ( nelem, ltype, ngaus,            &
#ifndef NDIMEPAR
           !$OMP            ndime,                          &
#endif
           !$OMP            vesgs )
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              vesgs(ielem) % a(1:ndime,1:pgaus,2) = vesgs(ielem) % a(1:ndime,1:pgaus,1)
           end do
           !$OMP END PARALLEL DO
        end if

     case ( 1400_ip )
        !
        ! Assign xx(n,i-1,*) <-- xx(n,i,*) for rigid body variables
        !
        if ( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod /= 0 ) ) then   ! rigid body
           do iimbo = 1,nrbod
              vpfor         => rbbou(iimbo) % vpfor
              vptor         => rbbou(iimbo) % vptor
              vpfor(1:3,2)  =  vpfor(1:3,1)
              vptor(1:3,2)  =  vptor(1:3,1)
           end do
        end if

     case ( 1500_ip )
        !
        ! Relax UNKNO
        !
        if( .not. NSI_MONOLITHIC ) then
           !
           ! Block Gauss-Seidel: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           if( ivari_nsi == ivari_nsi_mom .and. kfl_relax_nsi /= 0 ) then
              rela1_nsi = 1.0_rp - relax_nsi
              do ipoin = 1,npoin
                 do idime = 1,ndime
                    itotv = (ipoin-1)*ndime+idime
                    unkno(itotv) = relax_nsi * unkno(itotv)&
                         &       + rela1_nsi * veloc(idime,ipoin,1)
                 end do
              end do
           end if
           if( ivari_nsi == ivari_nsi_cont .and. kfl_relap_nsi /= 0 ) then
              rela1_nsi = 1.0_rp-relap_nsi
              do ipoin = 1,npoin
                 itotv = ndbgs_nsi+ipoin
                 unkno(itotv) = relap_nsi * unkno(itotv)&
                      &       + rela1_nsi * press(ipoin,1)
              end do
           end if

        else
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           if( kfl_relax_nsi /= 0 ) then
              rela1_nsi = 1.0_rp - relax_nsi
              do ipoin = 1,npoin
                 do idime = 1,ndime
                    itotv = (ipoin-1) * solve(1) % ndofn + idime
                    unkno(itotv) = relax_nsi * unkno(itotv)&
                         &       + rela1_nsi * veloc(idime,ipoin,1)
                 end do
              end do
           end if
           if( kfl_relap_nsi /= 0 ) then
              rela1_nsi = 1.0_rp - relap_nsi
              do ipoin = 1,npoin
                 itotv = ipoin * solve(1) % ndofn
                 unkno(itotv) = relax_nsi * unkno(itotv)&
                      &       + rela1_nsi * press(ipoin,1)
              end do
           end if

        end if

     case ( 2000_ip )
        !
        ! UNKNO <= PRESS(:,:,1)
        ! Hydrostatic pressure
        !
        unkno(1:npoin) = press(1:npoin,1)

     case ( 2200_ip )
        !
        ! UNKNO <= VELOC(:,:,1)
        !
        if( NSI_MONOLITHIC ) then
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = (ipoin-1) * solve(1) % ndofn + idime
                 unkno(itotv) = veloc(idime,ipoin,1)
              end do
              unkno(itotv+1) = press(ipoin,1)
           end do

        else
           !
           ! Schur complement: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           itotv     = 0
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = itotv + 1
                 unkno(itotv) = veloc(idime,ipoin,nprev_nsi)
              end do
              unkno(ndbgs_nsi+ipoin) = press(ipoin,nprev_nsi)
           end do
        end if

     case( 2300_ip )
        !
        ! VELOC(:,:,1) <= UNKNO
        !
        if( .not. NSI_MONOLITHIC ) then
           !
           ! Schur complement: UNKNO=[ u1 v1 u2 v2 ... p1 p2 ... ]
           !
           itotv = 0
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = itotv+1
                 veloc(idime,ipoin,1) = unkno(itotv)
              end do
              press(ipoin,1) = unkno(ndbgs_nsi+ipoin)
           end do

        else
           !
           ! Monolithic: UNKNO=[ u1 v1 p1 u2 v2 p2 ... ]
           !
           do ipoin = 1,npoin
              do idime = 1,ndime
                 itotv = (ipoin-1) * solve(1) % ndofn + idime
                 veloc(idime,ipoin,1) = unkno(itotv)
              end do
              press(ipoin,1) = unkno(itotv+1)
           end do

        end if

     end select

  else  ! the RB variables must also be updated by the master

     select case (itask)

     case( ITASK_ENDSTE )

        if ( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod /= 0 ) ) then   ! rigid body
           do iimbo = 1,nrbod
              vpfor         => rbbou(iimbo) % vpfor
              vptor         => rbbou(iimbo) % vptor
              vpfor(1:3,4)  =  vpfor(1:3,3)
              vpfor(1:3,3)  =  vpfor(1:3,1)
              vptor(1:3,4)  =  vptor(1:3,3)
              vptor(1:3,3)  =  vptor(1:3,1)
           end do
        end if

     case(1400_ip)
        !
        ! Assign xx(n,i-1,*) <-- xx(n,i,*) for rigid body variables
        !
        if ( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod /= 0 ) ) then   ! rigid body
           do iimbo = 1,nrbod
              vpfor         => rbbou(iimbo) % vpfor
              vptor         => rbbou(iimbo) % vptor
              vpfor(1:3,2)  =  vpfor(1:3,1)
              vptor(1:3,2)  =  vptor(1:3,1)
           end do
        end if

     end select

  end if

end subroutine nsi_updunk

