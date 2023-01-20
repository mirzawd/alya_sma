!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



subroutine renelm()
  !-----------------------------------------------------------------------
  !****f* domain/renelm
  ! NAME
  !    domain
  ! DESCRIPTION
  !    Renumber the elements and the element arrays
  ! OUTPUT
  ! USED BY
  !    Turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp,      only : ip
!  use def_domain,      only : ompss_domains
  use def_domain,      only : meshe
  use def_master,      only : ISLAVE
  use def_kermod,      only : ndivi
  use def_kermod,      only : kfl_renumbering_nelem
  use mod_renumbering, only : renumbering_elements
  use mod_renumbering, only : renumbering_element_arrays
  use mod_messages,    only : livinf
  implicit none
!  integer(ip)          :: isubd,offset_element
  integer(ip), pointer :: permr(:)

  if( ISLAVE .and. kfl_renumbering_nelem /= 0 .and. meshe(ndivi) % nelem > 0 ) then

     call livinf(0_ip,'RENUMBER ELEMENTS',0_ip)
     nullify(permr)
     !
     ! Compute permutation arrays: uncomment following lines to renumber
     ! each OMPSS subdomain independently
     ! 
     !#ifdef ALYA_OMPSS
     !offset_element = 0
     !do isubd = 1,size(ompss_domains)
     !   call renumbering_elements(1_ip,meshe(ndivi),permr,ompss_domains(isubd) % elements,offset_element) 
     !   offset_element = offset_element + size(ompss_domains(isubd) % elements)
     !end do
     !#else
     call renumbering_elements(kfl_renumbering_nelem,meshe(ndivi),permr,NAME='PERM2') 
     !#endif
     !
     ! Update element arrays using permutation
     !
     call renumbering_element_arrays(permr,NAME='PERM2')
     
  end if

end subroutine renelm
