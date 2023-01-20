!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine ale_smodef()
  !-----------------------------------------------------------------------
  !****f* domain/ale_smodef
  ! NAME
  !    domain
  ! DESCRIPTION
  !    This routines 
  !    KFL_FIXNO_ALE(:,:) = -1 ... FMALE free
  !                       =  0 ... ALE free
  !                       =  1 ... ALE fixed
  !                       =  3 ... FMALE fixed
  ! USED BY
  !    Turnon 
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_parame
  use def_elmtyp
  use def_domain
  use def_alefor
  use mod_memchk
  use def_kermod,         only : kfl_adj_prob
  use mod_messages, only : livinf

  implicit none
  integer(ip)    :: idime,ipoin
  character(150) :: messa

  if( kfl_defor_ale /= 0 .or. kfl_smoot_ale /= 0 .or. kfl_smobo_ale /= 0 ) then

     !-------------------------------------------------------------------
     !
     ! Copy new coordinates
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        do ipoin = 1,npoin
           do idime = 1,ndime
              coord_ale(idime,ipoin,1) = coord(idime,ipoin) 
           end do
        end do
     end if

     !-------------------------------------------------------------------
     !
     ! Boundary smoothing: Sharp edges have KFL_FIXNO_ALE(1,IPOIN) = 2
     !
     !-------------------------------------------------------------------

     if( kfl_smobo_ale /= 0 ) then
        call runend('CHECK ASSEMBLY, NO LONGER SUPPORTED')
        messa = intost(nsmob_ale) 
        messa = 'BOUNDARY SMOOTHING USING '//trim(messa)//' STEPS'
        call livinf(59_ip,trim(messa),modul)
        call ale_smobou()
     end if

     !-------------------------------------------------------------------
     !                   
     ! Mesh deformation
     !
     ! Update displacement and mesh velocity
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !  
     !-------------------------------------------------------------------

     if( kfl_defor_ale /= 0 ) then
        messa = intost(ndefo_ale) 
        messa = 'MESH DEFORMATION USING '//trim(messa)//' LOADING STEPS'
        call livinf(59_ip,trim(messa),modul)
        call ale_deform()
        call ale_updunk(ITASK_ENDINN)
     end if

     !-------------------------------------------------------------------
     !                                                    
     ! Mesh smoothing
     !
     ! Update displacement and mesh velocity
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !                                                                   
     !-------------------------------------------------------------------

     if( kfl_smoot_ale /= 0 ) then
        messa = intost(nsmoo_ale) 
        messa = 'MESH SMOOTHING USING '//trim(messa)//' STEP'
        call livinf(59_ip,trim(messa),modul)
        call ale_smooth()
        call ale_updunk(ITASK_ENDINN)
     end if

     !-------------------------------------------------------------------
     ! 
     ! Update displacement and mesh velocity
     ! FMALE: 1. does not update coordinate
     !        2. Invert mesh velocity
     !
     !-------------------------------------------------------------------
     !call ale_updunk(3_ip)

  end if
  
  if (kfl_adj_prob == 1) then
    
     !-------------------------------------------------------------------
     ! 
     ! Update mesh sensitivities 
     ! 
     !-------------------------------------------------------------------

     call ale_mshsen
     return
  endif

end subroutine ale_smodef
