!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tur_endite(itask)
!-----------------------------------------------------------------------
!****f* Turbul/tur_endite
! NAME 
!    tur_endite
! DESCRIPTION
!    This routine checks convergence and performs updates of the
!    turbulence variables at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    tur_cvgunk
!    tur_updunk
! USED BY
!    tur_doiter
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_turbul
  use def_kermod
  use mod_messages,   only : livinf
  use mod_ker_updpro, only : ker_updpro
  implicit none
  integer(ip) :: itask

  select case ( itask )

  case ( ITASK_ENDINN )
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || f(n,i,j) - f(n,i,j-1)|| / ||f(n,i,j)||) and update unknowns:
     !  f(n,i,j-1) <-- f(n,i,j) 
     !
     call tur_updrel()     ! Compute relaxation factor
     ! do not relax k-equation of k-omega model because is linear
     if (.not.tur_k_omega.or.(iunkn_tur.ne.1)) call tur_updunk(888_ip) ! Relax UNKNO 
     call tur_cvgunk(ITASK_ENDINN) ! Compute residual=UNTUR-UNKN
     call tur_updunk(ITASK_ENDINN) ! Actualize UNTUR=UNKNO and TURMU
     !
     ! Residual projection and subgrid scale equation
     !
     call tur_solsgs()
    
     call ker_updpro(ITASK_ENDITE)
     
  case ( ITASK_ENDITE )
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || f(n,i,*) - f(n,i-1,*)|| / ||f(n,i,*)||) and update unknowns:
     !  f(n,i-1,*) <-- f(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call tur_cvgunk(ITASK_ENDITE)
     call tur_updunk(ITASK_ENDITE)
     
  end select

end subroutine tur_endite
