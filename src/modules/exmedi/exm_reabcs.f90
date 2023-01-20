!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine exm_reabcs
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_reabcs
  ! NAME 
  !    exm_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions... 
  !
  ! USES
  !    exm_bcntoe
  !    ecoute
  !    memchk
  !    runend
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      mod_opebcs
  use      def_exmedi
    use mod_ecoute, only      : ecoute
  implicit none
  !
  ! Allocate memory
  !
  if( kfl_icodn > 0 ) then
     call opnbcs(0_ip,1_ip,ndofn_exm,0_ip,tncod_exm) 
  end if
  if( kfl_icodb > 0 ) then
     call opbbcs(0_ip,1_ip,1_ip,tbcod_exm)      
  end if

  if( INOTSLAVE ) then

     !
     ! Reach the nodal-wise section.
     !

     rewind(lisda)
     call ecoute('exm_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('exm_reabcs')
     end do
     !
     ! Read data.
     !
     call ecoute('exm_reabcs')
!     if(kfl_rstar==1) then
!        if(words(1)/='THIS ') call runend('THIS IS NOT A RE-START FILE')
!        call ecoute('exm_reabcs')
!     else
!        if(words(1)=='THIS ') call runend('THIS IS A RE-START FILE')
!     end if
     do while(words(1)/='ENDBO')

        if (exists('RETURN')) then
           kfl_exboc_exm = 0        
        else
           kfl_exboc_exm = 1
        end if

        if(words(1)=='CODES'.and.exists('NODES')) then 
           !
           ! User-defined codes on nodes
           !
           tncod => tncod_exm
           call reacod(1_ip)

        else if(words(1)=='CODES'.and.exists('BOUND')) then
           !
           ! User-defined codes on boundaries
           !          
           tbcod     => tbcod_exm(1:)
           call reacod(2_ip)
           
        else if( words(1) == 'CODES' ) then
           call runend('CODES section without NODES or BOUNDARIES')

        end if

        call ecoute('exm_reabcs')

     end do

  end if

end subroutine exm_reabcs
