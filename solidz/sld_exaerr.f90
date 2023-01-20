!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine sld_exaerr(itask)
  !-----------------------------------------------------------------------
  !****f* Nstinc/sld_exaerr
  ! NAME
  !    sld_exacso
  ! DESCRIPTION
  !    This routine computes the FEM errors (referred to an analytical
  !    solution defined in exacso.f90). The errors are normalized by the
  !    appropriate norm of the exact solution, except when this norm is zero.
  ! USES
  ! USED BY
  !    sld_output
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_solidz
  use def_elmtyp
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_MAX
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)    :: dummr(3)
  integer(ip) :: ielem,ipoin,idime                    ! Indices and dimensions
  integer(ip) :: ibopo
  real(rp)    :: u(ndime)                             ! u(i) = ui
  real(rp)    :: dudx(ndime,ndime)                    ! dudx(i,j) = duj/dxi
  real(rp)    :: d2udx2(ndime,ndime,ndime)            ! d2udx2(i,j,k) = d^2uk/dxidxj

  integer(ip) :: ielem_min
  real(rp)    :: gpdet_min

  if( kfl_exacs_sld /= 0 ) then

     if( itask == 1 .and. INOTMASTER ) then
        !
        ! Impose exact Dirichlet boundary and initial condition
        !
        do ipoin = 1,npoin
           ibopo = lpoty(ipoin)
           call sld_exacso(1_ip,coord(1,ipoin),u,dudx,d2udx2,dummr,dummr)
           do idime = 1,ndime
              if( kfl_fixno_sld(idime,ipoin) > 0 ) then
                 bvess_sld(idime,ipoin,1) = u(idime)
                 if( ibopo > 0 ) kfl_fixno_sld(idime,ipoin) = 1
              end if
           end do
        end do

     else
        !
        ! Compute errors
        !
        err01_sld = 0.0_rp
        err02_sld = 0.0_rp
        err0i_sld = 0.0_rp
        err11_sld = 0.0_rp
        err12_sld = 0.0_rp
        err1i_sld = 0.0_rp
        if( INOTMASTER ) then
           do ielem = 1,nelem
              call sld_element_operations(&
                   8_ip,ielem,ielem_min,gpdet_min)
           end do
        end if

        call PAR_SUM(2_ip,err01_sld)
        call PAR_SUM(2_ip,err02_sld)
        call PAR_MAX(2_ip,err0i_sld)
        call PAR_SUM(2_ip,err11_sld)
        call PAR_SUM(2_ip,err12_sld)
        call PAR_MAX(2_ip,err1i_sld)

        if( INOTSLAVE ) then

           err02_sld(1) = sqrt(err02_sld(1))
           err12_sld(1) = sqrt(err12_sld(1))
           err02_sld(2) = sqrt(err02_sld(2))
           err12_sld(2) = sqrt(err12_sld(2))

           if( err01_sld(2) > 0.0_rp ) err01_sld(1) = err01_sld(1) / err01_sld(2)
           if( err02_sld(2) > 0.0_rp ) err02_sld(1) = err02_sld(1) / err02_sld(2)
           if( err0i_sld(2) > 0.0_rp ) err0i_sld(1) = err0i_sld(1) / err0i_sld(2)
           if( err11_sld(2) > 0.0_rp ) err11_sld(1) = err11_sld(1) / err11_sld(2)
           if( err12_sld(2) > 0.0_rp ) err12_sld(1) = err12_sld(1) / err12_sld(2)
           if( err1i_sld(2) > 0.0_rp ) err1i_sld(1) = err1i_sld(1) / err1i_sld(2)

           write(momod(modul)%lun_outpu,100)                    &
                &       ittim,itinn(modul),                     &
                &       err01_sld(1),err02_sld(1),err0i_sld(1), &
                &       err11_sld(1),err12_sld(1),err1i_sld(1)
        end if

     end if
  end if

100 format(///,10X,'FINITE ELEMENT ERRORS',                                    &
       &              /,10X,'=====================',//,                        &
       &  '          TIME STEP NO.',I5,',  ITERATION NO. ',I5,/,10X,23('-'),/, &
       &  '          NORM    DISPLACEMENT               ',/,10X,23('-'),/,     &
       &  '          W(0,1) ',es16.8e3,/, &
       &  '          W(0,2) ',es16.8e3,/, &
       &  '          W(0,i) ',es16.8e3,/, &
       &  '          W(1,1) ',es16.8e3,/, &
       &  '          W(1,2) ',es16.8e3,/, &
       &  '          W(1,i) ',es16.8e3,/,10X,23('-'))

end subroutine sld_exaerr
