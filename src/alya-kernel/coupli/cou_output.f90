!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!----------------------------------------------------------------------
!> @addtogroup Coupling
!> @{
!> @file    cou_output.f90
!> @author  Guillaume Houzeaux
!> @date    18/06/2014
!> @brief   Output coupling 
!> @details Output some information about the couplinga
!> @} 
!----------------------------------------------------------------------

subroutine cou_output()
  use def_kintyp
  use def_parame
  use mod_parall
  use def_domain
  use def_master
  use def_coupli
  use mod_communications
  use mod_iofile, only : iofile_flush_unit
  use def_mpi
#include "def_mpi.inc"
  
  implicit none  
  integer(ip)                  :: ielem,ipoin
  integer(ip)                  :: ipart,dummi,ineig
  integer(ip)                  :: ipart_world,icoup
  integer(ip)                  :: kelem
  integer(ip)                  :: jneig
  integer(ip)                  :: number_wet_points,mneig
  integer(ip)                  :: npoin_wet         
  integer(ip)                  :: nboun_wet         
  integer(ip)                  :: nneig_source     
  integer(ip)                  :: istep,nstep,part1,part2     
  type(comm_data_par), pointer :: commu
  MY_MPI_COMM                  :: PAR_COMM_TO_USE
  real(rp)                     :: y

  integer(ip)                  :: nneig
  integer(ip), pointer         :: lneig(:)
  real(rp),    pointer         :: xcoor_gat(:,:)
  integer(ip), pointer         :: nneig_gat(:)
  integer(4),  pointer         :: nneig4_gat(:)
  integer(ip), pointer         :: lneig_gat(:)
  integer(ip), pointer         :: lresu_gat(:)
  integer(ip), pointer         :: lcode_gat(:)
  type(i1p),   pointer         :: ledge(:)
  integer(ip), pointer         :: comm_last(:)
  integer(ip), pointer         :: comm_matrix(:,:)

  integer(ip) :: lun_coupli_msh
  integer(ip) :: lun_coupli_res

  return

  lun_coupli_msh = 101
  lun_coupli_res = 102

  if( ISEQUEN ) return

  call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
  
  nneig         =  0
  nullify( lneig      )
  nullify( xcoor_gat  )
  nullify( nneig_gat  )
  nullify( nneig4_gat )
  nullify( lneig_gat  )
  nullify( lresu_gat  )
  nullify( lcode_gat  )
  nullify( ledge      )
  nullify( comm_last  )
  nullify( comm_matrix)

  if( PAR_MY_WORLD_RANK == 0 ) then
     allocate( xcoor_gat (3,0:PAR_WORLD_SIZE-1) )
     allocate( nneig_gat (  0:PAR_WORLD_SIZE-1) )
     allocate( nneig4_gat(  0:PAR_WORLD_SIZE-1) )
     allocate( lresu_gat (  0:PAR_WORLD_SIZE-1) )
     allocate( lcode_gat (  0:PAR_WORLD_SIZE-1) )
     allocate( ledge     (  0:PAR_WORLD_SIZE-1) )
     do ipart = 0,PAR_WORLD_SIZE-1 
        nullify(ledge(ipart) % l)
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! Couplings
  !
  !----------------------------------------------------------------------

  do icoup = 1,mcoup
     color_target      = coupling_type(icoup) % color_target
     color_source      = coupling_type(icoup) % color_source
     number_wet_points = coupling_type(icoup) % wet % number_wet_points
     npoin_wet         = coupling_type(icoup) % wet % npoin_wet
     nboun_wet         = coupling_type(icoup) % wet % nboun_wet
     nneig_source      = coupling_type(icoup) % commd % nneig
     !
     ! Geometry
     !
     call PAR_GATHER(nneig_source,nneig_gat,'IN THE WORLD') 

     if( PAR_MY_WORLD_RANK == 0 ) then        
        dummi = 0
        do ipart = 0,PAR_WORLD_SIZE-1 
           dummi = dummi + nneig_gat(ipart)
           nneig4_gat(ipart) = int(nneig_gat(ipart),4)
        end do
        if( associated(lneig_gat) ) deallocate( lneig_gat )
        allocate( lneig_gat(dummi) )
     else
        if( INOTMASTER ) then
           if( associated(lneig) ) deallocate( lneig )
           allocate( lneig(coupling_type(icoup) % commd % nneig) )
           do ineig = 1,coupling_type(icoup) % commd % nneig   
              ipart        = coupling_type(icoup) % commd % neights(ineig)
              ipart_world  = PAR_COMM_COLOR_PERM(color_target,color_source,ipart)
              lneig(ineig) = ipart_world
           end do
        end if
     end if

     call PAR_GATHERV(lneig,lneig_gat,nneig4_gat,'IN THE WORLD')

     if( PAR_MY_WORLD_RANK == 0 ) then        
        write(lun_coupli_msh,'(a,i1,a)') 'MESH COUPLING'//trim(intost(icoup))//' dimension ',ndime,' Elemtype Linear Nnode 2'
        write(lun_coupli_msh,'(a)') 'coordinates'
        do ipart = 0,PAR_WORLD_SIZE-1
           mneig = mneig + nneig_gat(ipart)
        end do
        allocate( comm_matrix(mneig,0:PAR_WORLD_SIZE-1) ) 
        allocate( comm_last(0:PAR_WORLD_SIZE-1) )
        do ipart = 0,PAR_WORLD_SIZE-1
           comm_last(ipart) = 1
           do ineig = 1,mneig
              comm_matrix(ineig,ipart) = 0
           end do
        end do

        jneig = 0
        nstep = 0
        do while( jneig < mneig )
           nstep = nstep + 1
           kelem = 0
           do ipart = 0,PAR_WORLD_SIZE-1
              ineig = comm_last(ipart)
              if( ineig <= nneig_gat(ipart) ) then
                 !
                 ! ipart+1 ========= lneig_gat(kelem)+1
                 !
                 kelem = kelem + ineig
                 part1 = ipart 
                 part2 = lneig_gat(kelem) 
                 if( comm_matrix(nstep,part1) == 0 .and. comm_matrix(nstep,part2) == 0 ) then
                    jneig                    = jneig + 2
                    comm_matrix(nstep,part1) = part2
                    comm_matrix(nstep,part2) = part1             
                    comm_last(ipart)         = comm_last(ipart) + 1
                 end if
              end if
           end do
        end do

        write(lun_coupli_msh,'(a)') 'coordinates'
        ipoin = 0
        y     = 0.0_rp
        do istep = 1,nstep
           do ipart = 0,PAR_WORLD_SIZE-1
              if( comm_matrix(nstep,part1) /= 0 ) then
                 part2 = comm_matrix(istep,part1)
                 comm_matrix(istep,part2) = 0
                 y     = y - 1.0_rp
                 ipoin = ipoin + 1 ; write(lun_coupli_msh,'(i5,3(1x,e12.6))') ipoin,real(part1,rp),y
                 ipoin = ipoin + 1 ; write(lun_coupli_msh,'(i5,3(1x,e12.6))') ipoin,real(part2,rp),y
              end if
           end do
           y = y - 2.0_rp
        end do
        write(lun_coupli_msh,'(a)') 'end coordinates'

        write(lun_coupli_msh,'(a,i1,a)') 'MESH COUPLING'//trim(intost(icoup))//' dimension ',ndime,' Elemtype Linear Nnode 2'
        write(lun_coupli_msh,'(a)') 'elements'
        ielem = 0
        ipoin = 0
        do istep = 1,nstep
           do ipart = 0,PAR_WORLD_SIZE-1
              if( comm_matrix(istep,part1) /= 0 ) then
                 part2 = comm_matrix(istep,part1)
                 ielem = ielem + 1
                 write(lun_coupli_msh,'(3(1x,i7))') ielem,ipoin+1,ipoin+2
                 ipoin = ipoin + 2
              end if
           end do
        end do

        write(lun_coupli_msh,'(a)') 'end elements' 

     end if
  end do
  !
  ! Deallocate
  !
  if( associated(lneig)      ) deallocate( lneig      )
  if( associated(xcoor_gat)  ) deallocate( xcoor_gat  )
  if( associated(nneig_gat)  ) deallocate( nneig_gat  )
  if( associated(nneig4_gat) ) deallocate( nneig4_gat )
  if( associated(lneig_gat)  ) deallocate( lneig_gat  )
  if( associated(lresu_gat)  ) deallocate( lresu_gat  )
  if( associated(lcode_gat)  ) deallocate( lcode_gat  )
  if( associated(ledge)      ) deallocate( ledge      )

  if( PAR_MY_WORLD_RANK == 0 ) then
     close(lun_coupli_msh)
     close(lun_coupli_res)     
  end if
 
end subroutine cou_output
