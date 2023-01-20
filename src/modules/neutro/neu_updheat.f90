!-----------------------------------------------------------------!
! Copyright 2005 - 2022 Barcelona Supercomputing Center.          !
! Distributed under the ALYA AVAILABLE SOURCE ("ALYA AS") LICENSE !
! for nonprofit scientific purposes only.                         !
! See companion file LICENSE.txt.                                 !
!-----------------------------------------------------------------!



!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_updheat.f90
!> @date    04/05/2022
!> @author  Ezequiel Goldberg (EG)
!> @brief   Calculate heat after neutron flux
!> @details Calculate heat after neutron flux (currently, stationary). 
!>          Heat can be included in the output and can be used as heat source in TEMPER
!> @}
!------------------------------------------------------------------------

subroutine neu_updheat()
    use def_master
    use def_domain
    use def_neutro
    implicit none
    integer(ip)              :: ipoin,iener, ielem,idire!,idime
    integer(ip), allocatable :: elem_x_nodes(:)

    if( INOTMASTER ) then
        ! Calculate kerma factor per energy group for each node (originally per material)
        ! For nodes shared between different materials, compute average kerma (for each energy group)
        if( num_materials_neu > 1 ) then
           allocate(elem_x_nodes(npoin))
           elem_x_nodes = 0_ip
           kerma_poin_neu = 0.0_rp
           do ielem=1,nelem
              do ipoin=1,lnnod(ielem)
                 do iener=1,num_energies_neu
                    kerma_poin_neu(lnods(ipoin,ielem), iener) = &
                        kerma_poin_neu(lnods(ipoin,ielem), iener) + kerma_neu(lmate(ielem), iener)
                 end do
                 elem_x_nodes(lnods(ipoin,ielem)) = elem_x_nodes(lnods(ipoin,ielem)) + 1
              end do
           end do
           do ipoin=1,npoin
              do iener=1,num_energies_neu
                 kerma_poin_neu(ipoin, iener) = kerma_poin_neu(ipoin, iener)/real(elem_x_nodes(ipoin),rp)
              end do
           end do
           DEALLOCATE(elem_x_nodes)
        else
           do iener=1,num_energies_neu
              do ipoin=1,npoin
                 kerma_poin_neu(ipoin, iener) = kerma_neu(1, iener)
              end do
           end do
        end if

        !CALCULATE HEAT
        do ipoin = 1,npoin  
            heat_sink(ipoin) = 0.0_rp 
            do iener = 1,num_energies_neu !> Iterating over the energies in each point
               do idire = 1,num_directions_neu !> Iterating over the directions in each point
                  ! if (neutr(iener,idire,ipoin,1) > 0.0_rp) then ! EG: why was this here? 
                     heat_sink(ipoin) = heat_sink(ipoin) + &
                        neutr(iener,idire,ipoin,1) * weigd_neu(idire) * kerma_poin_neu(ipoin, iener)
                  ! end if
               end do
            enddo 
         end do
    end if
end subroutine neu_updheat