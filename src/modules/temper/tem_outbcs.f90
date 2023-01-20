!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine tem_outbcs
!-----------------------------------------------------------------------
!****f* Nastin/tem_outbcs
! NAME 
!    tem_outbcs
! DESCRIPTION
!    Postprocess boundary conditions. This could be useful to include
!    this new file when running the same problem over again.
! USES
! USED BY
!    tem_output
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_temper
  use      mod_iofile
  implicit none
  integer(ip) :: ipoin

  do ipoin=1,npoin
     write(lun_bound_tem,1,advance='no') ipoin
     write(lun_bound_tem,2,advance='no') kfl_fixno_tem(1,ipoin)        
     write(lun_bound_tem,5,advance='no') therm(ipoin,1)
     if(kfl_conbc_tem==0) then
        write(lun_bound_tem,4,advance='no') kfl_funno_tem(ipoin)
     end if
  end do
 
  close(lun_bound_tem)

1 format(1x,i7)
2 format(1x,i1)
4 format(1x,i2)
5 format(1x,e16.8E3)

end subroutine tem_outbcs
