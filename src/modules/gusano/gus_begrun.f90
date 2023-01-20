!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Gusano
!> @{
!> @file    gus_begrun.f90
!> @author  houzeaux
!> @date    2020-10-23
!> @brief   Begin run
!> @details Begin the run
!> @} 
!-----------------------------------------------------------------------

subroutine gus_begrun()

  use def_master
  use def_domain
  use def_gusano
  use mod_memory_basic
  use mod_maths_basic, only : maths_normalize_vector
  implicit none
  integer(ip)          :: ipoin,ipoin1,ipoin2,ielem
  integer(ip), pointer :: num_angles(:)
  real(rp),    pointer :: angles(:,:,:)
  real(rp)             :: vec1(3)

  nullify(num_angles)
  nullify(angles)
  !
  ! Angles
  !
  call memory_alloca(mem_modul(1:2,modul),'NUM_ANGLES','gus_begrun',num_angles,npoin)
  call memory_alloca(mem_modul(1:2,modul),'ANGLES'    ,'gus_begrun',angles    ,2_ip,ndime,npoin)
  do ielem = 1,nelem
     if( ltype(ielem) > 0 ) then
        ipoin1                              = lnods(1,ielem)
        ipoin2                              = lnods(2,ielem)
        num_angles(ipoin1)                  = num_angles(ipoin1) + 1
        num_angles(ipoin2)                  = num_angles(ipoin2) + 1
        vec1(1:ndime)                       = coord(1:ndime,ipoin2)-coord(1:ndime,ipoin1)
        call maths_normalize_vector(ndime,vec1)
        angles(num_angles(ipoin1),:,ipoin1) = vec1(1:ndime)
        angles(num_angles(ipoin2),:,ipoin2) = vec1(1:ndime)
     end if
  end do
  do ipoin = 1,npoin
     if( num_angles(ipoin1) == 2 ) then
        angle_gus(ipoin) = dot_product(angles(1,1:ndime,ipoin),angles(2,1:ndime,ipoin))
     end if
  end do
  call memory_deallo(mem_modul(1:2,modul),'NUM_ANGLES','gus_begrun',num_angles)
  call memory_deallo(mem_modul(1:2,modul),'ANGLES'    ,'gus_begrun',angles)

end subroutine gus_begrun
