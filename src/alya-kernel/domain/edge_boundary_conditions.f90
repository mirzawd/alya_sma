!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    edge_boundary_conditions.f90
!> @date    04/02/2016
!> @author  Guillaume Houzeaux
!> @brief   Compute edge codes
!> @details Compute edge codes KFL_CODED by extrapolating boundary codes
!>
!> @} 
!-----------------------------------------------------------------------

subroutine edge_boundary_conditions()
  use def_kintyp, only : ip
  use def_kermod, only : kfl_edge_elements
  use def_kermod, only : ndivi 
  use def_domain, only : meshe
  use def_domain, only : mcodb
  use def_domain, only : kfl_codbo
  use def_domain, only : kfl_coded
  use def_domain, only : kfl_icodb
  use def_master,   only : INOTMASTER
  use mod_messages, only : livinf
  use mod_domain,   only : domain_memory_allocate
  implicit none
  integer(ip) :: iboun,pblty,pedge,icodb,iedgg,iedge,code1,code2
  integer(ip) :: mcodb1

  if( kfl_edge_elements == 1 .and. kfl_icodb > 0 ) then

     call livinf(0_ip,'COMPUTE BOUNDARY CONDITIONS ON EDGES',0_ip)

     if( INOTMASTER ) then
        !
        ! Allocate memory for KFL_CODED(:,:)
        !
        call domain_memory_allocate('KFL_CODED')
        mcodb1    = mcodb+1
        kfl_coded = mcodb1
        !
        ! Extrapolate from bondary to edge
        !
        do iboun = 1,meshe(ndivi) % nboun
           pblty = meshe(ndivi) % ltypb(iboun) 
           pedge = meshe(ndivi) % lnneb(iboun)
           icodb = kfl_codbo(iboun)
           do iedge = 1,pedge 
              iedgg = meshe(ndivi) % ledgb(iedge,iboun)
              code1 = kfl_coded(1,iedgg)
              code2 = kfl_coded(2,iedgg)
              if( icodb /= code1 .and. icodb /= code2 ) then
                 if( kfl_coded(1,iedgg) == mcodb1 ) then
                    kfl_coded(1,iedgg) = icodb
                 else
                    kfl_coded(2,iedgg) = icodb
                 end if
              end if
              if( kfl_coded(1,iedgg) > kfl_coded(2,iedgg) ) then
                 code1              = kfl_coded(1,iedgg)
                 kfl_coded(1,iedgg) = kfl_coded(2,iedgg)
                 kfl_coded(2,iedgg) = code1               
              end if
           end do
        end do

     end if

  end if

end subroutine edge_boundary_conditions
