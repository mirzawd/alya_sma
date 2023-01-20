!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_begite
!-----------------------------------------------------------------------
!****f* Exmedi/exm_begite
! NAME 
!    exm_begite
! DESCRIPTION
!    This routine starts the internal iteration 
! USES
!    exm_inisol
!    exm_updunk
! USED BY
!    exm_doiter
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain    
  use      def_exmedi
  use      mod_eccoupling
  use mod_messages, only : livinf
  use mod_biofibers, only : kfl_biofibers, biofibers

  implicit none
  integer (ip) :: imate,kmodel_maxvalue
!
! Initializations
!
  kfl_goite_exm = 1
  itinn(modul)  = 0
  if(miinn_exm==0) kfl_goite_exm=0
!  if(itcou==1) call exm_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Obtain the initial guess for inner iterations: y(2) <- y(1)
  !

  !
  ! check if any TT-like model is present
  !
  kmodel_maxvalue= 0
  do imate= 1,nmate
     kmodel_maxvalue = max(kfl_cellmod(imate),kmodel_maxvalue)
  end do
  
  call exm_updunk(ITASK_BEGITE)    ! u(,ITER_K) <-- u(,ITER_AUX)
 
  if( kfl_biofibers )then
     if(kfl_gcoup_exm == 1_ip )then 
        call biofibers % get_current_fibers_at_nodes() 
        fiber => biofibers % nod_lng(:,:,1)
     else
        call biofibers % get_reference_fibers_at_nodes() 
        fiber => biofibers % nod_lng(:,:,2)
     endif
     
  endif

  
end subroutine exm_begite
