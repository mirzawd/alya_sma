!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine nsi_out_res_eocoe(xx,nsize,ipass_aux)
  !
  ! outputs results for eocoe project
  !
  use def_kintyp
  implicit none

  integer(ip), intent(in)    :: nsize
  real(rp),    intent(in)    :: xx(nsize)
  integer(ip), intent(in)    :: ipass_aux

  integer(ip)        :: i
  character(50)      :: char_aux


  write(char_aux,'(i1)')ipass_aux
  char_aux='result'//trim(char_aux)//'.txt'
  open(777,file=char_aux)
  ! xx
  write(777,'(i15)')nsize
  do i=1,nsize
     write(777,'(e17.10)')xx(i)
  end do
  close(777)

 end subroutine nsi_out_res_eocoe
