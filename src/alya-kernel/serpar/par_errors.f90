!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine par_errors(itask)
!------------------------------------------------------------------------
!****f* Parall/par_errors
! NAME
!    par_errors
! DESCRIPTION
!    Check errors
! OUTPUT
! USED BY
!    Parall
!***
!------------------------------------------------------------------------
  use      def_kintyp
  use      def_domain
  use      def_master
  use      def_parall
  implicit none
  integer(ip), intent(in) :: itask

  if(kfl_paral==0) then

     select case(itask)

     case(1)
        if(rp/=8) then
           call runend('PARALL: ONLY WORKS WITH DOUBLE PRECISION')
        end if
        
        !if(ip/=4) then
        !   call runend('PARALL: ONLY WORKS WITH 4-BYTE INTEGERS')
        !end if
        
        if(kfl_autbo==1) then
           call runend('PARALL: CANNOT USE AUTOMATIC BOUNDARY GENERATION')
        end if
        
        !if(kfl_ptask==2.and.kfl_postp_par==1) then
        !   call runend('PARALL: WHEN DOING RESTART, MASTER CANNOT BE IN CHARGE OF POSTPROCESS')
        !end if
        
     case(2)

     end select

  end if

end subroutine par_errors
