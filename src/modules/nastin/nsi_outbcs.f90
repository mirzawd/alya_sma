!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_outbcs
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_outbcs
  ! NAME 
  !    nsi_outbcs
  ! DESCRIPTION
  !    Postprocess boundary conditions on nodes in a format directly
  !    readable by Alya. This could be useful to include this new file
  !    when running the same problem over again, and to avoid:
  !    1. Recomputing the automatic intersections when using slip and
  !       wall law boundary conditions.
  !    2. Reinterpolating from a background mesh.
  ! USES
  ! USED BY
  !    nsi_output
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastin
  use      mod_iofile
  implicit none
  integer(ip) :: ipoin,idime,ibopo

  if( ISEQUEN ) then
     write(lun_bound_nsi,6) '  ON_NODES'
     do ipoin=1,npoin
        write(lun_bound_nsi,1,advance='no') ipoin
        if(ndime==2) then
           write(lun_bound_nsi,2,advance='no') (kfl_fixno_nsi(idime,ipoin),idime=1,ndime)
        else
           write(lun_bound_nsi,3,advance='no') (kfl_fixno_nsi(idime,ipoin),idime=1,ndime)
        end if
        do idime=1,ndime
           write(lun_bound_nsi,5,advance='no') veloc(idime,ipoin,1)
        end do
        if(kfl_conbc_nsi==0) then
           write(lun_bound_nsi,4,advance='no') kfl_funno_nsi(ipoin)
        end if
        ibopo=lpoty(ipoin)
        if(ibopo>0) then
           write(lun_bound_nsi,1,advance='no') kfl_fixrs_nsi(ipoin)
        else
           write(lun_bound_nsi,1,advance='no') zero              
        end if
        write(lun_bound_nsi,*)
     end do
     write(lun_bound_nsi,6) '  END_ON_NODES'

     close(lun_bound_nsi)
  end if

1 format(1x,i7)
2 format(1x,2(i1))
3 format(1x,3(i1))
4 format(1x,i2)
5 format(1x,e16.8E3)
6 format(a)

end subroutine nsi_outbcs
