!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_output_miniapp.f90
!> @author  guillaume
!> @date    2021-07-14
!> @brief   miniapp output
!> @details Output for the Nastin miniapp
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_output_miniapp()

#include "def_vector_size.inc"
  use def_kintyp_basic, only : ip,rp
  use mod_memory_basic, only : memory_size
  use mod_elmgeo,       only : element_type
  use mod_strings,      only : integer_to_string
  use def_master,       only : kfl_paral 
  use def_master,       only : INOTMASTER
  use def_master,       only : veloc
  use def_master,       only : rhsid
  use def_master,       only : namda
  use mod_iofile,       only : iofile_available_unit
  use def_domain
  use mod_parall
  use def_nastin
  implicit none
  integer(ip)   :: ielty,isubd,pnode,pgaus,pdime
  integer(ip)   :: n_rhsid,n_veloc
  integer(ip)   :: n_num_pack_par
  integer(ip)   :: n_packs
  integer(ip)   :: n_neighbours
  integer(ip)   :: n_elements,ii,ivect
  integer(4)    :: unit4  
  character(11) :: wstrat
  
  if( INOTMASTER .and. kfl_miniapp_nsi /= 0 ) then
 
     unit4 = iofile_available_unit(90_ip)
     open(unit4,file=trim(namda)//'-miniapp-'//integer_to_string(kfl_paral)//'.bin',status='unknown')
     !
     ! Some keywords
     !
     ivect = int(VECTOR_SIZE,ip)
     write(unit4,*) ivect
#ifdef ALYA_OMPSS         
     wstrat = 'ALYA_OMPSS '
#elif NO_COLORING 
     wstrat = 'NO_COLORING'
#else                  
     wstrat = 'NULL       '
#endif
     write(unit4,*) wstrat
     !
     ! Dimensions
     !
     n_rhsid         = memory_size(rhsid)
     n_veloc         = memory_size(veloc,3_ip)
     n_num_pack_par  = memory_size(num_pack_par)
     write(unit4,*) n_veloc,n_rhsid
     write(unit4,*) n_num_pack_par     
     write(unit4,*) par_omp_nelem_chunk 
     !
     ! elmar
     !
     write(unit4,*) ndime,mnode,npoin,nelem     
     write(unit4,*) lexis(1:nelty)
     do ielty = 1,nelty
        if( lexis(ielty) /= 0 ) then
           pdime = ldime(ielty)
           pnode = element_type(ielty) % number_nodes
           pgaus = elmar(ielty) % pgaus
           write(unit4,*) ielty,pgaus,pnode,pdime
           if( associated(elmar(ielty) % shape) ) write(unit4,*) elmar(ielty) % shape(1:pnode,1:pgaus)
           if( associated(elmar(ielty) % deriv) ) write(unit4,*) elmar(ielty) % deriv(1:pdime,1:pnode,1:pgaus)
           if( associated(elmar(ielty) % weigp) ) write(unit4,*) elmar(ielty) % weigp(1:pgaus)
        end if
     end do
     !
     ! Mesh arrays
     !
     write(unit4,*) lnods(1:mnode,1:nelem)
     write(unit4,*) ltype(1:nelem)
     write(unit4,*) lmate(1:nelem)
     write(unit4,*) lnnod(1:nelem)
     write(unit4,*) lgaus(1:nelem)
     write(unit4,*) coord(1:ndime,1:npoin)
     !
     ! def_master
     !
     write(unit4,*) veloc(1:ndime,1:npoin,1:n_veloc)
     write(unit4,*) rhsid(1:n_rhsid)
     !
     ! mod_parall
     !
     write(unit4,*) num_subd_par
     if( associated(num_pack_par) ) write(unit4,*) num_pack_par(1:num_subd_par)
     if( associated(list_elements_par) ) then 
        do isubd = 1,num_subd_par
           write(unit4,*) memory_size(list_elements_par(isubd) % packs)           
           if( associated(list_elements_par(isubd) % packs) ) then
              do ii = 1,num_pack_par(isubd)
                 n_packs = memory_size(list_elements_par(isubd) % packs(ii) % l)
                 write(unit4,*) n_packs 
                 if( n_packs > 0 ) write(unit4,*) list_elements_par(isubd) % packs(ii) % l(1:n_packs)
              end do
           end if
        end do
     end if
     !
     ! ompss domains
     !
     if( associated(ompss_domains) ) then
        write(unit4,*) 1_ip
     else
        write(unit4,*) 0_ip
     end if

     if( associated(ompss_domains) ) then
        do isubd = 1,num_subd_par
           n_neighbours = memory_size(ompss_domains(isubd) % neighbours)
           n_elements   = memory_size(ompss_domains(isubd) % elements)
           write(unit4,*) ompss_domains(isubd) % neighIdx,ompss_domains(isubd) % elemIdx,n_neighbours,n_elements
           if( n_neighbours > 0 ) write(unit4,*) ompss_domains(isubd) % neighbours(1:n_neighbours)
           if( n_elements   > 0 ) write(unit4,*) ompss_domains(isubd) % elements  (1:n_elements)
        end do
     end if
     !
     ! Close unit 
     !
     flush(unit4)
     close(unit4)
  end if
  
end subroutine nsi_output_miniapp
