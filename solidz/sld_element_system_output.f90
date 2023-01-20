!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @{file   sld_element_system_output.f90
!> @author  Guillaume Houzeaux version of nastin (adapted by ECR)
!> @date    April, 2018
!>          - Subroutine written
!> @brief   Output elemental matrix-vector system
!> @details Output elemental matrix-vector system
!>
!>          \verbatim
!>
!>          Without enrichement
!>          -------------------
!>            +-        +  +- -+     +-  -+
!>            | Kuu Kuv |  | u |     | fu |
!>            |         |  |   |  =  |    |
!>            | Kvu Kvv |  | v |     | fv |
!>            +-       -+  +- -+     +-  -+
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_element_system_output(pnode,elmat,elrhs)

  use def_kintyp, only :  ip,rp
  use def_master, only :  intost, modul, itinn, ittim, cutim
  use def_domain, only :  ndime

  implicit none

  integer(ip), intent(in)           :: pnode
  ! Element matrices
  real(rp),    intent(in)           :: elmat(ndime*pnode,ndime*pnode)
  real(rp),    intent(in)           :: elrhs(ndime*pnode)

  integer(ip),   parameter          :: lreal=9
  character(13)                     :: FMT1
  integer(ip)                       :: inode,idime,idofn,jnode,jdime,jdofn
  character(3)                      :: vanam

  !
  ! Format
  !
  FMT1 = '(1x,e'//trim(intost(lreal))//'.'//trim(intost(lreal-7))//',)'
  !
  ! Time info
  !
  write(99,100) cutim
  write(99,101) ittim
  write(99,102) itinn(modul)
  !
  ! First line
  !
  write(99,'(a)',advance='no') '+-'
  do idofn = 1,pnode*ndime
     do idime = 1,lreal+1
        write(99,'(a)',advance='no') ' '
     end do
  end do

  write(99,'(a)',advance='no') '-+'
  write(99,'(a)',advance='no') '  +-   -+'

  write(99,'(a)',advance='no') '     +-'
  do idime = 1,lreal
     write(99,'(a)',advance='no') ' '
  end do
  write(99,'(a)',advance='no') '-+'

  write(99,*)
  !
  ! Matrix - vector
  !
  do inode = 1,pnode
     do idime = 1,ndime
        !
        ! K
        !
        write(99,'(a)',advance='no') '|'
        idofn = (inode-1)*ndime+idime
        do jnode = 1,pnode
           do jdime = 1,ndime
              jdofn = (jnode-1)*ndime+jdime
              write(99,FMT1,advance='no') elmat(idofn,jdofn)
           end do
        end do
        !
        ! u
        !
        write(99,'(a)',advance='no') ' '
        write(99,'(a)',advance='no') ' |'
        if( idime == 1 ) then
           vanam = 'u'
        else if( idime == 2 ) then
           vanam = 'v'
        else if( idime == 3 ) then
           vanam = 'w'
        end if
        vanam = trim(vanam) // trim(intost(inode))
        write(99,'(a)',advance='no') '  | '//vanam//' |'
        !
        ! f
        !
        write(99,'(a)',advance='no') '     |'
        write(99,FMT1,advance='no') elrhs(idofn)
        write(99,'(a)',advance='no') ' |'

        write(99,*)
     end do
  end do
  !
  ! Last line
  !
  write(99,'(a)',advance='no') '+-'
  do idofn = 1,pnode*ndime
     do idime = 1,lreal+1
        write(99,'(a)',advance='no') ' '
     end do
  end do

  write(99,'(a)',advance='no') '-+'
  write(99,'(a)',advance='no') '  +-   -+'

  write(99,'(a)',advance='no') '     +-'
  do idime = 1,lreal
     write(99,'(a)',advance='no') ' '
  end do
  write(99,'(a)',advance='no') '-+'

  write(99,*)
  flush(99)

100 format ('Time =', e16.8e3)
101 format ('Time step =', i9)
102 format ('N-R iteration =', i9)

end subroutine sld_element_system_output

